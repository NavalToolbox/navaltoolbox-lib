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

//! Result structures for criteria scripts.

use serde::{Deserialize, Serialize};

/// Result of a criteria verification script.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CriteriaResult {
    /// Name of the regulation (e.g., "IMO A.749(18)")
    pub regulation_name: String,
    /// Document reference
    pub regulation_reference: String,
    /// Vessel name
    pub vessel_name: String,
    /// Loading condition description
    pub loading_condition: String,
    /// Displacement in kg
    pub displacement: f64,
    /// Center of gravity [LCG, TCG, VCG]
    pub cog: [f64; 3],
    /// Individual criterion results
    pub criteria: Vec<CriterionResult>,
    /// Overall pass status (true if no FAIL criteria)
    pub overall_pass: bool,
    /// Number of passing criteria
    pub pass_count: usize,
    /// Number of failing criteria
    pub fail_count: usize,
    /// General notes
    pub notes: String,
    /// Plot data for GUI visualization
    pub plots: Vec<PlotData>,
}

impl Default for CriteriaResult {
    fn default() -> Self {
        Self {
            regulation_name: String::new(),
            regulation_reference: String::new(),
            vessel_name: String::new(),
            loading_condition: String::new(),
            displacement: 0.0,
            cog: [0.0, 0.0, 0.0],
            criteria: Vec::new(),
            overall_pass: false,
            pass_count: 0,
            fail_count: 0,
            notes: String::new(),
            plots: Vec::new(),
        }
    }
}

impl std::fmt::Display for CriteriaResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "======================================================================"
        )?;
        writeln!(f, "  {}", self.regulation_name)?;
        writeln!(f, "  Reference: {}", self.regulation_reference)?;
        writeln!(
            f,
            "======================================================================"
        )?;
        writeln!(f, "Vessel: {}", self.vessel_name)?;
        writeln!(f, "Loading Condition: {}", self.loading_condition)?;
        writeln!(f, "Displacement: {:.0} kg", self.displacement)?;
        writeln!(
            f,
            "COG: ({:.2}, {:.2}, {:.2}) m",
            self.cog[0], self.cog[1], self.cog[2]
        )?;
        writeln!(f)?;
        writeln!(f, "CRITERIA RESULTS")?;
        writeln!(
            f,
            "----------------------------------------------------------------------"
        )?;

        for crit in &self.criteria {
            let icon = match crit.status {
                CriteriaStatus::Pass => "✓",
                CriteriaStatus::Fail => "✗",
                CriteriaStatus::NotApplicable => "-",
                CriteriaStatus::Warning => "!",
            };
            writeln!(f, "[{}] {}", icon, crit.name)?;
            writeln!(f, "    Required: {:.4} {}", crit.required_value, crit.unit)?;
            writeln!(f, "    Actual:   {:.4} {}", crit.actual_value, crit.unit)?;
            if crit.required_value != 0.0 {
                let margin_percent = (crit.margin / crit.required_value.abs()) * 100.0;
                writeln!(
                    f,
                    "    Margin:   {:+.4} {} ({:+.1}%)",
                    crit.margin, crit.unit, margin_percent
                )?;
            }
            if let Some(ref notes) = crit.notes {
                if !notes.is_empty() {
                    writeln!(f, "    Note: {}", notes)?;
                }
            }
            writeln!(f)?;
        }

        writeln!(
            f,
            "----------------------------------------------------------------------"
        )?;
        let status = if self.overall_pass { "PASS" } else { "FAIL" };
        writeln!(
            f,
            "OVERALL: {} ({}/{} criteria passed)",
            status,
            self.pass_count,
            self.criteria.len()
        )?;
        writeln!(
            f,
            "======================================================================"
        )?;
        Ok(())
    }
}

/// Result of a single criterion check.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CriterionResult {
    /// Criterion name (e.g., "Area 0-30°")
    pub name: String,
    /// Description
    pub description: String,
    /// Required/threshold value
    pub required_value: f64,
    /// Actual calculated value
    pub actual_value: f64,
    /// Unit (e.g., "m·rad", "°", "m")
    pub unit: String,
    /// Pass/Fail status
    pub status: CriteriaStatus,
    /// Margin (actual - required)
    pub margin: f64,
    /// Optional notes
    pub notes: Option<String>,
    /// Associated plot ID (for GUI grouping)
    pub plot_id: Option<String>,
}

impl Default for CriterionResult {
    fn default() -> Self {
        Self {
            name: String::new(),
            description: String::new(),
            required_value: 0.0,
            actual_value: 0.0,
            unit: String::new(),
            status: CriteriaStatus::NotApplicable,
            margin: 0.0,
            notes: None,
            plot_id: None,
        }
    }
}

/// Status of a criterion check.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CriteriaStatus {
    Pass,
    Fail,
    NotApplicable,
    Warning,
}

impl std::fmt::Display for CriteriaStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CriteriaStatus::Pass => write!(f, "PASS"),
            CriteriaStatus::Fail => write!(f, "FAIL"),
            CriteriaStatus::NotApplicable => write!(f, "N/A"),
            CriteriaStatus::Warning => write!(f, "WARN"),
        }
    }
}

// ============================================================================
// Plot data for GUI visualization
// ============================================================================

/// Container for plot elements related to a criteria check.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PlotData {
    /// Unique identifier for the plot
    pub id: String,
    /// Plot title
    pub title: String,
    /// X-axis label
    pub x_label: String,
    /// Y-axis label
    pub y_label: String,
    /// Graphic elements (curves, lines, areas, points)
    pub elements: Vec<PlotElement>,
}

impl Default for PlotData {
    fn default() -> Self {
        Self {
            id: "main".to_string(),
            title: String::new(),
            x_label: "Heel (°)".to_string(),
            y_label: "GZ (m)".to_string(),
            elements: Vec::new(),
        }
    }
}

/// Individual plot element.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum PlotElement {
    /// A curve/line series
    Curve {
        name: String,
        x: Vec<f64>,
        y: Vec<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        color: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        style: Option<String>,
    },

    /// A horizontal line
    HorizontalLine {
        name: String,
        y: f64,
        #[serde(skip_serializing_if = "Option::is_none")]
        x_min: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        x_max: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        color: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        style: Option<String>,
    },

    /// A vertical line
    VerticalLine {
        name: String,
        x: f64,
        #[serde(skip_serializing_if = "Option::is_none")]
        y_min: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        y_max: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        color: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        style: Option<String>,
    },

    /// A filled area between curves
    FilledArea {
        name: String,
        x: Vec<f64>,
        y_lower: Vec<f64>,
        y_upper: Vec<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        color: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        alpha: Option<f64>,
    },

    /// A point marker
    Point {
        name: String,
        x: f64,
        y: f64,
        #[serde(skip_serializing_if = "Option::is_none")]
        marker: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        color: Option<String>,
    },

    /// An annotation/label
    Annotation {
        text: String,
        x: f64,
        y: f64,
        #[serde(skip_serializing_if = "Option::is_none")]
        color: Option<String>,
    },
}
