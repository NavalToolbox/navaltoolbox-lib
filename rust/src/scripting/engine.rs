// Copyright (C) 2026 Antoine ANCEAU
//
// This file is part of navaltoolbox.
//
// navaltoolbox is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.

//! Rhai script execution engine.

use std::fs;

use rhai::{Engine, Scope, AST};
use thiserror::Error;

use super::context::CriteriaContext;
use super::result::{CriteriaResult, CriteriaStatus, CriterionResult, PlotData, PlotElement};

/// Errors that can occur during script execution.
#[derive(Debug, Error)]
pub enum ScriptError {
    #[error("Failed to read script file: {0}")]
    FileRead(#[from] std::io::Error),

    #[error("Script compilation error: {0}")]
    Compile(String),

    #[error("Script execution error: {0}")]
    Runtime(String),

    #[error("Script did not return a valid result")]
    InvalidResult,
}

/// Rhai script execution engine for criteria verification.
pub struct ScriptEngine {
    engine: Engine,
}

impl Default for ScriptEngine {
    fn default() -> Self {
        Self::new()
    }
}

impl ScriptEngine {
    /// Create a new script engine with navaltoolbox functions registered.
    pub fn new() -> Self {
        let mut engine = Engine::new();

        // Register CriteriaContext type and methods
        engine
            .register_type_with_name::<CriteriaContext>("CriteriaContext")
            // GZ curve methods
            .register_fn("get_heels", CriteriaContext::get_heels)
            .register_fn("get_gz_values", CriteriaContext::get_gz_values)
            .register_fn("area_under_curve", CriteriaContext::area_under_curve)
            .register_fn("gz_at_angle", CriteriaContext::gz_at_angle)
            .register_fn("find_max_gz", CriteriaContext::find_max_gz)
            .register_fn(
                "find_angle_of_vanishing_stability",
                CriteriaContext::find_angle_of_vanishing_stability,
            )
            .register_fn(
                "get_first_flooding_angle",
                CriteriaContext::get_first_flooding_angle,
            )
            .register_fn(
                "find_equilibrium_angle",
                CriteriaContext::find_equilibrium_angle,
            )
            .register_fn(
                "find_second_intercept",
                CriteriaContext::find_second_intercept,
            )
            .register_fn("get_limiting_angle", CriteriaContext::get_limiting_angle)
            // Hydrostatic properties
            .register_fn("get_gm0", CriteriaContext::get_gm0)
            .register_fn("get_gm0_dry", CriteriaContext::get_gm0_dry)
            .register_fn("get_draft", CriteriaContext::get_draft)
            .register_fn("get_trim", CriteriaContext::get_trim)
            .register_fn("get_displacement", CriteriaContext::get_displacement)
            .register_fn("get_cog", CriteriaContext::get_cog)
            // Form coefficients
            .register_fn("get_cb", CriteriaContext::get_cb)
            .register_fn("get_cm", CriteriaContext::get_cm)
            .register_fn("get_cp", CriteriaContext::get_cp)
            .register_fn("get_lwl", CriteriaContext::get_lwl)
            .register_fn("get_bwl", CriteriaContext::get_bwl)
            .register_fn("get_vcb", CriteriaContext::get_vcb)
            // Wind data
            .register_fn("has_wind_data", CriteriaContext::has_wind_data)
            .register_fn("get_emerged_area", CriteriaContext::get_emerged_area)
            .register_fn("get_wind_lever_arm", CriteriaContext::get_wind_lever_arm)
            .register_fn(
                "calculate_wind_heeling_lever",
                CriteriaContext::calculate_wind_heeling_lever,
            )
            // External parameters
            .register_fn("get_param", CriteriaContext::get_param)
            .register_fn("has_param", CriteriaContext::has_param)
            // Metadata
            .register_fn("get_vessel_name", CriteriaContext::get_vessel_name)
            .register_fn(
                "get_loading_condition",
                CriteriaContext::get_loading_condition,
            );

        // Register helper functions
        engine.register_fn("criterion", Self::create_criterion);

        Self { engine }
    }

    /// Helper to create a criterion result map from Rust.
    ///
    /// Usage in Rhai: `criterion("Name", "Desc", required, actual, "unit")`
    fn create_criterion(
        name: &str,
        description: &str,
        required: f64,
        actual: f64,
        unit: &str,
    ) -> rhai::Map {
        let margin = actual - required;
        let status = if actual >= required { "PASS" } else { "FAIL" };

        let mut map = rhai::Map::new();
        map.insert("name".into(), rhai::Dynamic::from(name.to_string()));
        map.insert(
            "description".into(),
            rhai::Dynamic::from(description.to_string()),
        );
        map.insert("required".into(), rhai::Dynamic::from(required));
        map.insert("actual".into(), rhai::Dynamic::from(actual));
        map.insert("unit".into(), rhai::Dynamic::from(unit.to_string()));
        map.insert("status".into(), rhai::Dynamic::from(status.to_string()));
        map.insert("margin".into(), rhai::Dynamic::from(margin));
        map
    }

    /// Compile a script from a string.
    pub fn compile(&self, script: &str) -> Result<AST, ScriptError> {
        self.engine
            .compile(script)
            .map_err(|e| ScriptError::Compile(e.to_string()))
    }

    /// Compile a script from a file.
    pub fn compile_file(&self, path: &str) -> Result<AST, ScriptError> {
        let script = fs::read_to_string(path)?;
        self.compile(&script)
    }

    /// Run a script from a file path.
    pub fn run_script_file(
        &self,
        path: &str,
        context: CriteriaContext,
    ) -> Result<CriteriaResult, ScriptError> {
        let script = fs::read_to_string(path)?;
        self.run_script(&script, context)
    }

    /// Run a script from a string.
    pub fn run_script(
        &self,
        script: &str,
        context: CriteriaContext,
    ) -> Result<CriteriaResult, ScriptError> {
        let ast = self.compile(script)?;
        self.run_ast(&ast, context)
    }

    /// Run a pre-compiled AST.
    pub fn run_ast(
        &self,
        ast: &AST,
        context: CriteriaContext,
    ) -> Result<CriteriaResult, ScriptError> {
        let mut scope = Scope::new();

        // Clone context for the function call argument
        let ctx_for_call = context.clone();
        scope.push("ctx", context);

        // Call the check function with the cloned context
        let result: rhai::Map = self
            .engine
            .call_fn(&mut scope, ast, "check", (ctx_for_call,))
            .map_err(|e| ScriptError::Runtime(e.to_string()))?;

        // Convert the map to CriteriaResult
        self.map_to_criteria_result(&result)
    }

    /// Convert a Rhai map to CriteriaResult.
    fn map_to_criteria_result(&self, map: &rhai::Map) -> Result<CriteriaResult, ScriptError> {
        let mut result = CriteriaResult::default();

        // Extract fields from the map
        if let Some(v) = map.get("regulation_name") {
            result.regulation_name = v.clone().into_string().unwrap_or_default();
        }
        if let Some(v) = map.get("regulation_reference") {
            result.regulation_reference = v.clone().into_string().unwrap_or_default();
        }
        if let Some(v) = map.get("vessel_name") {
            result.vessel_name = v.clone().into_string().unwrap_or_default();
        }
        if let Some(v) = map.get("loading_condition") {
            result.loading_condition = v.clone().into_string().unwrap_or_default();
        }
        if let Some(v) = map.get("displacement") {
            result.displacement = v.as_float().unwrap_or(0.0);
        }
        if let Some(v) = map.get("cog") {
            if let Some(arr) = v.clone().try_cast::<rhai::Array>() {
                if arr.len() >= 3 {
                    result.cog = [
                        arr[0].as_float().unwrap_or(0.0),
                        arr[1].as_float().unwrap_or(0.0),
                        arr[2].as_float().unwrap_or(0.0),
                    ];
                }
            }
        }
        if let Some(v) = map.get("notes") {
            result.notes = v.clone().into_string().unwrap_or_default();
        }
        if let Some(v) = map.get("overall_pass") {
            result.overall_pass = v.as_bool().unwrap_or(false);
        }

        // Extract criteria
        if let Some(criteria_val) = map.get("criteria") {
            if let Some(criteria_arr) = criteria_val.clone().try_cast::<rhai::Array>() {
                for crit_dyn in criteria_arr {
                    if let Some(crit_map) = crit_dyn.try_cast::<rhai::Map>() {
                        if let Some(crit) = self.map_to_criterion_result(&crit_map) {
                            if crit.status == CriteriaStatus::Pass {
                                result.pass_count += 1;
                            } else if crit.status == CriteriaStatus::Fail {
                                result.fail_count += 1;
                            }
                            result.criteria.push(crit);
                        }
                    }
                }
            }
        }

        // Extract plots
        if let Some(plots_val) = map.get("plots") {
            if let Some(plots_arr) = plots_val.clone().try_cast::<rhai::Array>() {
                for plot_dyn in plots_arr {
                    if let Some(plot_map) = plot_dyn.try_cast::<rhai::Map>() {
                        if let Some(plot) = self.map_to_plot_data(&plot_map) {
                            result.plots.push(plot);
                        }
                    }
                }
            }
        }

        Ok(result)
    }

    /// Convert a Rhai map to CriterionResult.
    fn map_to_criterion_result(&self, map: &rhai::Map) -> Option<CriterionResult> {
        let name = map.get("name")?.clone().into_string().unwrap_or_default();
        let description = map
            .get("description")
            .and_then(|v| v.clone().into_string().ok())
            .unwrap_or_default();
        let required_value = map
            .get("required")
            .or_else(|| map.get("required_value"))
            .and_then(|v| v.as_float().ok())
            .unwrap_or(0.0);
        let actual_value = map
            .get("actual")
            .or_else(|| map.get("actual_value"))
            .and_then(|v| v.as_float().ok())
            .unwrap_or(0.0);
        let unit = map
            .get("unit")
            .and_then(|v| v.clone().into_string().ok())
            .unwrap_or_default();

        // Parse status
        let status_str = map
            .get("status")
            .and_then(|v| v.clone().into_string().ok())
            .unwrap_or_else(|| "N/A".to_string());

        let status = match status_str.to_uppercase().as_str() {
            "PASS" => CriteriaStatus::Pass,
            "FAIL" => CriteriaStatus::Fail,
            "WARNING" | "WARN" => CriteriaStatus::Warning,
            _ => CriteriaStatus::NotApplicable,
        };

        let margin = map
            .get("margin")
            .and_then(|v| v.as_float().ok())
            .unwrap_or(actual_value - required_value);

        let notes = map.get("notes").and_then(|v| v.clone().into_string().ok());

        let plot_id = map
            .get("plot_id")
            .and_then(|v| v.clone().into_string().ok());

        Some(CriterionResult {
            name,
            description,
            required_value,
            actual_value,
            unit,
            status,
            margin,
            notes,
            plot_id,
        })
    }

    /// Convert a Rhai map to PlotData.
    fn map_to_plot_data(&self, map: &rhai::Map) -> Option<PlotData> {
        let id = map
            .get("id")
            .and_then(|v| v.clone().into_string().ok())
            .unwrap_or_else(|| "main".to_string());
        let title = map
            .get("title")
            .and_then(|v| v.clone().into_string().ok())
            .unwrap_or_default();
        let x_label = map
            .get("x_label")
            .and_then(|v| v.clone().into_string().ok())
            .unwrap_or_else(|| "Heel (°)".to_string());
        let y_label = map
            .get("y_label")
            .and_then(|v| v.clone().into_string().ok())
            .unwrap_or_else(|| "GZ (m)".to_string());

        // Parse elements
        let mut elements = Vec::new();
        if let Some(elements_val) = map.get("elements") {
            if let Some(elements_arr) = elements_val.clone().try_cast::<rhai::Array>() {
                for elem_dyn in elements_arr {
                    if let Some(elem_map) = elem_dyn.try_cast::<rhai::Map>() {
                        if let Some(elem) = self.map_to_plot_element(&elem_map) {
                            elements.push(elem);
                        }
                    }
                }
            }
        }

        Some(PlotData {
            id,
            title,
            x_label,
            y_label,
            elements,
        })
    }

    /// Convert a Rhai map to PlotElement.
    fn map_to_plot_element(&self, map: &rhai::Map) -> Option<PlotElement> {
        let elem_type = map
            .get("type")
            .and_then(|v| v.clone().into_string().ok())
            .unwrap_or_default();

        match elem_type.as_str() {
            "Curve" => {
                let name = map
                    .get("name")
                    .and_then(|v| v.clone().into_string().ok())
                    .unwrap_or_default();
                let x = self.extract_f64_array(map.get("x")?)?;
                let y = self.extract_f64_array(map.get("y")?)?;
                let color = map.get("color").and_then(|v| v.clone().into_string().ok());
                let style = map.get("style").and_then(|v| v.clone().into_string().ok());
                Some(PlotElement::Curve {
                    name,
                    x,
                    y,
                    color,
                    style,
                })
            }
            "HorizontalLine" => {
                let name = map
                    .get("name")
                    .and_then(|v| v.clone().into_string().ok())
                    .unwrap_or_default();
                let y = map.get("y")?.as_float().ok()?;
                Some(PlotElement::HorizontalLine {
                    name,
                    y,
                    x_min: map.get("x_min").and_then(|v| v.as_float().ok()),
                    x_max: map.get("x_max").and_then(|v| v.as_float().ok()),
                    color: map.get("color").and_then(|v| v.clone().into_string().ok()),
                    style: map.get("style").and_then(|v| v.clone().into_string().ok()),
                })
            }
            "VerticalLine" => {
                let name = map
                    .get("name")
                    .and_then(|v| v.clone().into_string().ok())
                    .unwrap_or_default();
                let x = map.get("x")?.as_float().ok()?;
                Some(PlotElement::VerticalLine {
                    name,
                    x,
                    y_min: map.get("y_min").and_then(|v| v.as_float().ok()),
                    y_max: map.get("y_max").and_then(|v| v.as_float().ok()),
                    color: map.get("color").and_then(|v| v.clone().into_string().ok()),
                    style: map.get("style").and_then(|v| v.clone().into_string().ok()),
                })
            }
            "Point" => {
                let name = map
                    .get("name")
                    .and_then(|v| v.clone().into_string().ok())
                    .unwrap_or_default();
                let x = map.get("x")?.as_float().ok()?;
                let y = map.get("y")?.as_float().ok()?;
                Some(PlotElement::Point {
                    name,
                    x,
                    y,
                    marker: map.get("marker").and_then(|v| v.clone().into_string().ok()),
                    color: map.get("color").and_then(|v| v.clone().into_string().ok()),
                })
            }
            _ => None,
        }
    }

    /// Extract f64 array from a Dynamic.
    fn extract_f64_array(&self, dyn_val: &rhai::Dynamic) -> Option<Vec<f64>> {
        let arr = dyn_val.clone().try_cast::<rhai::Array>()?;
        Some(arr.iter().filter_map(|v| v.as_float().ok()).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_engine_creation() {
        let _engine = ScriptEngine::new();
    }
    #[test]
    fn test_criterion_helper() {
        let engine = ScriptEngine::new();
        let script = r#"
            let res = criterion("Test", "Desc", 1.0, 1.5, "m");
            res
        "#;

        let result: rhai::Map = engine.engine.eval(script).unwrap();
        assert_eq!(result["name"].clone().into_string().unwrap(), "Test");
        assert_eq!(result["status"].clone().into_string().unwrap(), "PASS");
        assert_eq!(result["margin"].as_float().unwrap(), 0.5);
    }
}
