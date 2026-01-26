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

//! Criteria context for Rhai scripts.
//!
//! Wraps `CompleteStabilityResult` and provides calculation methods
//! accessible from Rhai scripts.

use std::collections::HashMap;

use rhai::{Array, Dynamic, Map};

use crate::stability::CompleteStabilityResult;

/// Context for Rhai scripts, wrapping CompleteStabilityResult.
///
/// Provides access to GZ curve data, hydrostatic properties, form coefficients,
/// and user-defined external parameters.
#[derive(Debug, Clone)]
pub struct CriteriaContext {
    /// The complete stability result
    result: CompleteStabilityResult,
    /// Vessel name
    vessel_name: String,
    /// Loading condition description
    loading_condition: String,
    /// External parameters provided by the user
    params: HashMap<String, Dynamic>,
}

impl CriteriaContext {
    /// Create a new criteria context from a CompleteStabilityResult.
    pub fn new(
        result: CompleteStabilityResult,
        vessel_name: String,
        loading_condition: String,
    ) -> Self {
        Self {
            result,
            vessel_name,
            loading_condition,
            params: HashMap::new(),
        }
    }

    /// Set an external parameter accessible to scripts.
    pub fn set_param(&mut self, key: &str, value: Dynamic) {
        self.params.insert(key.to_string(), value);
    }

    /// Set a float parameter.
    pub fn set_param_f64(&mut self, key: &str, value: f64) {
        self.params.insert(key.to_string(), Dynamic::from(value));
    }

    /// Set a string parameter.
    pub fn set_param_str(&mut self, key: &str, value: &str) {
        self.params
            .insert(key.to_string(), Dynamic::from(value.to_string()));
    }

    /// Set a bool parameter.
    pub fn set_param_bool(&mut self, key: &str, value: bool) {
        self.params.insert(key.to_string(), Dynamic::from(value));
    }

    // =========================================================================
    // GZ Curve Analysis Methods
    // =========================================================================

    /// Get heel angles as array.
    pub fn get_heels(&self) -> Array {
        self.result
            .gz_curve
            .heels()
            .into_iter()
            .map(Dynamic::from)
            .collect()
    }

    /// Get GZ values as array.
    pub fn get_gz_values(&self) -> Array {
        self.result
            .gz_curve
            .values()
            .into_iter()
            .map(Dynamic::from)
            .collect()
    }

    /// Calculate area under GZ curve between two angles (m·rad).
    ///
    /// Uses trapezoidal integration.
    pub fn area_under_curve(&self, from_angle: f64, to_angle: f64) -> f64 {
        let heels = self.result.gz_curve.heels();

        if heels.is_empty() {
            return 0.0;
        }

        let mut area = 0.0;

        for i in 0..heels.len() - 1 {
            let h1 = heels[i];
            let h2 = heels[i + 1];

            // Skip if entirely outside range
            if h2 <= from_angle || h1 >= to_angle {
                continue;
            }

            // Clip to range
            let x1 = h1.max(from_angle);
            let x2 = h2.min(to_angle);

            // Interpolate values at clipped points
            let y1 = self.gz_at_angle(x1);
            let y2 = self.gz_at_angle(x2);

            // Trapezoid area (convert degrees to radians)
            let dx = (x2 - x1).to_radians();
            area += 0.5 * (y1 + y2) * dx;
        }

        area
    }

    /// Get GZ value at a specific angle by linear interpolation.
    pub fn gz_at_angle(&self, angle: f64) -> f64 {
        self.result.gz_curve.interpolate(angle).unwrap_or(0.0)
    }

    /// Find maximum GZ and its angle.
    ///
    /// Returns a map with keys "angle" and "value".
    pub fn find_max_gz(&self) -> Map {
        let mut map = Map::new();

        if let Some(point) = self.result.gz_curve.max_value() {
            map.insert("angle".into(), Dynamic::from(point.heel));
            map.insert("value".into(), Dynamic::from(point.value));
        } else {
            map.insert("angle".into(), Dynamic::UNIT);
            map.insert("value".into(), Dynamic::UNIT);
        }

        map
    }

    /// Find angle where GZ returns to 0 after maximum (vanishing stability).
    ///
    /// Returns the angle or () if not found.
    pub fn find_angle_of_vanishing_stability(&self) -> Dynamic {
        let heels = self.result.gz_curve.heels();
        let values = self.result.gz_curve.values();

        if heels.is_empty() {
            return Dynamic::UNIT;
        }

        // Find max GZ index
        let max_idx = values
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0);

        // Look for zero crossing after max
        for i in max_idx..values.len() - 1 {
            if values[i] > 0.0 && values[i + 1] <= 0.0 {
                // Linear interpolation to find exact crossing
                let t = values[i] / (values[i] - values[i + 1]);
                let angle = heels[i] + t * (heels[i + 1] - heels[i]);
                return Dynamic::from(angle);
            }
        }

        Dynamic::UNIT
    }

    /// Get first flooding angle from GZ curve.
    ///
    /// Returns the first heel angle where is_flooding is true, or () if none.
    pub fn get_first_flooding_angle(&self) -> Dynamic {
        for point in &self.result.gz_curve.points {
            if point.is_flooding {
                return Dynamic::from(point.heel);
            }
        }
        Dynamic::UNIT
    }

    /// Find equilibrium angle (first positive intercept) for a given heeling arm.
    pub fn find_equilibrium_angle(&self, heeling_arm: f64) -> Dynamic {
        let heels = self.result.gz_curve.heels();
        let values = self.result.gz_curve.values();

        if heels.is_empty() {
            return Dynamic::UNIT;
        }

        // Search for first crossing where GZ increases above lever or comes from below
        // Usually equilibrium is where Residual Righting Amount = GZ - HA crosses zero from pos to neg?
        // No, GZ must equal HA.
        // Stable equilibrium: GZ curve intersects HA curve from below (slope of GZ > slope of HA).
        // Since HA is constant, we look for GZ(phi) crossing HA from below.

        for i in 0..values.len() - 1 {
            let v1 = values[i];
            let v2 = values[i + 1];

            // Check if we cross the heeling arm line
            if (v1 <= heeling_arm && v2 >= heeling_arm) || (v1 >= heeling_arm && v2 <= heeling_arm)
            {
                // Linear interpolation
                let t = (heeling_arm - v1) / (v2 - v1);
                let angle = heels[i] + t * (heels[i + 1] - heels[i]);

                // Only return "stable" equilibrium? Or first one?
                // Typically the first positive intercept is the equilibrium.
                if angle >= 0.0 {
                    return Dynamic::from(angle);
                }
            }
        }
        Dynamic::UNIT
    }

    /// Find second intercept angle (where GZ becomes less than heeling arm).
    ///
    /// This is typically the unstable equilibrium or capsize limit.
    pub fn find_second_intercept(&self, heeling_arm: f64) -> Dynamic {
        let heels = self.result.gz_curve.heels();
        let values = self.result.gz_curve.values();

        if heels.is_empty() {
            return Dynamic::UNIT;
        }

        // We want the SECOND positive intercept.
        // Or specifically where GZ CROSSES HA from ABOVE. (Unstable)

        let mut intercept_count = 0;

        for i in 0..values.len() - 1 {
            let v1 = values[i];
            let v2 = values[i + 1];

            if (v1 <= heeling_arm && v2 >= heeling_arm) || (v1 >= heeling_arm && v2 <= heeling_arm)
            {
                // Linear interpolation
                let t = (heeling_arm - v1) / (v2 - v1);
                let angle = heels[i] + t * (heels[i + 1] - heels[i]);

                if angle >= 0.0 {
                    intercept_count += 1;
                    // If this is the crossing from ABOVE (unstable), it's likely the second one
                    // provided we started from upright stable.
                    // Or we just return the 2nd one we find.
                    if intercept_count == 2 {
                        return Dynamic::from(angle);
                    }

                    // Also check if we passed peak GZ and go down?
                    // The standard definition usually implies the geometric 2nd intercept.
                }
            }
        }

        Dynamic::UNIT
    }

    /// Get limiting angle (minimum of: default, flooding angle, vanishing stability).
    pub fn get_limiting_angle(&self, default: f64) -> f64 {
        let mut limit = default;

        if let Some(flooding) = self.get_first_flooding_angle().try_cast::<f64>() {
            limit = limit.min(flooding);
        }

        if let Some(vanishing) = self.find_angle_of_vanishing_stability().try_cast::<f64>() {
            limit = limit.min(vanishing);
        }

        limit
    }

    // =========================================================================
    // Hydrostatic Properties
    // =========================================================================

    /// Get initial metacentric height GM0 (with free surface correction).
    pub fn get_gm0(&self) -> Dynamic {
        match self.result.gm0() {
            Some(v) => Dynamic::from(v),
            None => Dynamic::UNIT,
        }
    }

    /// Get initial metacentric height without free surface correction.
    pub fn get_gm0_dry(&self) -> Dynamic {
        match self.result.gm0_dry() {
            Some(v) => Dynamic::from(v),
            None => Dynamic::UNIT,
        }
    }

    /// Get draft in meters.
    pub fn get_draft(&self) -> f64 {
        self.result.hydrostatics.draft
    }

    /// Get trim in degrees.
    pub fn get_trim(&self) -> f64 {
        self.result.hydrostatics.trim
    }

    /// Get displacement in kg.
    pub fn get_displacement(&self) -> f64 {
        self.result.displacement
    }

    /// Get center of gravity as array [LCG, TCG, VCG].
    pub fn get_cog(&self) -> Array {
        self.result.cog.iter().map(|&v| Dynamic::from(v)).collect()
    }

    // =========================================================================
    // Form Coefficients
    // =========================================================================

    /// Get block coefficient.
    pub fn get_cb(&self) -> f64 {
        self.result.hydrostatics.cb
    }

    /// Get midship section coefficient.
    pub fn get_cm(&self) -> f64 {
        self.result.hydrostatics.cm
    }

    /// Get prismatic coefficient.
    pub fn get_cp(&self) -> f64 {
        self.result.hydrostatics.cp
    }

    /// Get length at waterline in meters.
    pub fn get_lwl(&self) -> f64 {
        self.result.hydrostatics.lwl
    }

    /// Get beam at waterline in meters.
    pub fn get_bwl(&self) -> f64 {
        self.result.hydrostatics.bwl
    }

    /// Get vertical center of buoyancy in meters.
    pub fn get_vcb(&self) -> f64 {
        self.result.hydrostatics.vcb()
    }

    // =========================================================================
    // Wind Data
    // =========================================================================

    /// Check if wind heeling data is available.
    pub fn has_wind_data(&self) -> bool {
        self.result.has_wind_data()
    }

    /// Get emerged lateral area in m².
    pub fn get_emerged_area(&self) -> Dynamic {
        match &self.result.wind_data {
            Some(wind) => Dynamic::from(wind.emerged_area),
            None => Dynamic::UNIT,
        }
    }

    /// Get wind lever arm in meters.
    pub fn get_wind_lever_arm(&self) -> Dynamic {
        match &self.result.wind_data {
            Some(wind) => Dynamic::from(wind.wind_lever_arm),
            None => Dynamic::UNIT,
        }
    }

    /// Calculate wind heeling lever lw1 for a given wind pressure.
    ///
    /// lw1 = (P × A × Z) / (Δ × g)
    ///
    /// Returns the lever in meters, or () if wind data not available.
    pub fn calculate_wind_heeling_lever(&self, wind_pressure: f64) -> Dynamic {
        match &self.result.wind_data {
            Some(wind) => {
                let g = 9.81;
                let lw1 = (wind_pressure * wind.emerged_area * wind.wind_lever_arm)
                    / (self.result.displacement * g);
                Dynamic::from(lw1)
            }
            None => Dynamic::UNIT,
        }
    }

    // =========================================================================
    // External Parameters
    // =========================================================================

    /// Get an external parameter set by the user.
    ///
    /// Returns the value or () if not set.
    pub fn get_param(&self, key: &str) -> Dynamic {
        self.params.get(key).cloned().unwrap_or(Dynamic::UNIT)
    }

    /// Check if a parameter exists.
    pub fn has_param(&self, key: &str) -> bool {
        self.params.contains_key(key)
    }

    // =========================================================================
    // Metadata
    // =========================================================================

    /// Get vessel name.
    pub fn get_vessel_name(&self) -> String {
        self.vessel_name.clone()
    }

    /// Get loading condition.
    pub fn get_loading_condition(&self) -> String {
        self.loading_condition.clone()
    }

    /// Get the underlying CompleteStabilityResult.
    pub fn result(&self) -> &CompleteStabilityResult {
        &self.result
    }
}
