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

//! Stability dataclasses.

/// A point on a stability curve (KN or GZ).
#[derive(Debug, Clone)]
pub struct StabilityPoint {
    /// Heel angle in degrees
    pub heel: f64,
    /// Draft at equilibrium in meters
    pub draft: f64,
    /// Trim angle at equilibrium in degrees
    pub trim: f64,
    /// KN or GZ value in meters
    pub value: f64,
    /// Whether any downflooding opening is submerged at this heel
    pub is_flooding: bool,
    /// Names of submerged openings (empty if none)
    pub flooded_openings: Vec<String>,
}

/// A complete stability curve.
#[derive(Debug, Clone)]
pub struct StabilityCurve {
    /// Type of curve ("KN" or "GZ")
    pub curve_type: String,
    /// Displacement in kg
    pub displacement: f64,
    /// Center of gravity (LCG, TCG, VCG) - for GZ curves
    pub cog: Option<[f64; 3]>,
    /// Points on the curve
    pub points: Vec<StabilityPoint>,
}

impl StabilityCurve {
    /// Creates a new KN curve.
    pub fn new_kn(displacement: f64, points: Vec<StabilityPoint>) -> Self {
        Self {
            curve_type: "KN".to_string(),
            displacement,
            cog: None,
            points,
        }
    }

    /// Creates a new GZ curve.
    pub fn new_gz(displacement: f64, cog: [f64; 3], points: Vec<StabilityPoint>) -> Self {
        Self {
            curve_type: "GZ".to_string(),
            displacement,
            cog: Some(cog),
            points,
        }
    }

    /// Returns the heel angles.
    pub fn heels(&self) -> Vec<f64> {
        self.points.iter().map(|p| p.heel).collect()
    }

    /// Returns the values (KN or GZ).
    pub fn values(&self) -> Vec<f64> {
        self.points.iter().map(|p| p.value).collect()
    }

    /// Finds the maximum value on the curve.
    pub fn max_value(&self) -> Option<&StabilityPoint> {
        self.points
            .iter()
            .max_by(|a, b| a.value.partial_cmp(&b.value).unwrap())
    }

    /// Interpolates a value at a given heel angle.
    pub fn interpolate(&self, heel: f64) -> Option<f64> {
        if self.points.is_empty() {
            return None;
        }

        // Find bracketing points
        for i in 0..self.points.len() - 1 {
            if self.points[i].heel <= heel && heel <= self.points[i + 1].heel {
                let t =
                    (heel - self.points[i].heel) / (self.points[i + 1].heel - self.points[i].heel);
                return Some(
                    self.points[i].value + t * (self.points[i + 1].value - self.points[i].value),
                );
            }
        }

        None
    }
}

/// A stability curve with wind heeling moment data.
///
/// Per IMO 2008 IS Code (MSC.267), wind heeling levers are constant at all angles.
#[derive(Debug, Clone)]
pub struct StabilityCurveWithWind {
    /// Type of curve (always "GZ")
    pub curve_type: String,
    /// Displacement in kg
    pub displacement: f64,
    /// Center of gravity (LCG, TCG, VCG)
    pub cog: [f64; 3],
    /// Points on the curve
    pub points: Vec<StabilityPoint>,
    /// Wind pressure used in N/m²
    pub wind_pressure: f64,
    /// Steady wind heeling lever (lw1) in meters - constant per IMO
    pub wind_lever_lw1: f64,
    /// Gust wind heeling lever (lw2 = 1.5 × lw1) in meters
    pub wind_lever_lw2: f64,
    /// Emerged lateral area above waterline in m²
    pub wind_area: f64,
    /// Z coordinate of emerged area center in meters
    pub wind_center_z: f64,
}

impl StabilityCurveWithWind {
    /// Creates a new GZ curve with wind data.
    pub fn new(
        displacement: f64,
        cog: [f64; 3],
        points: Vec<StabilityPoint>,
        wind_pressure: f64,
        wind_lever_lw1: f64,
        wind_area: f64,
        wind_center_z: f64,
    ) -> Self {
        Self {
            curve_type: "GZ".to_string(),
            displacement,
            cog,
            points,
            wind_pressure,
            wind_lever_lw1,
            wind_lever_lw2: wind_lever_lw1 * 1.5,
            wind_area,
            wind_center_z,
        }
    }

    /// Returns the heel angles.
    pub fn heels(&self) -> Vec<f64> {
        self.points.iter().map(|p| p.heel).collect()
    }

    /// Returns the GZ values.
    pub fn gz_values(&self) -> Vec<f64> {
        self.points.iter().map(|p| p.value).collect()
    }

    /// Returns the GZ values corrected for wind (GZ - lw1).
    pub fn gz_corrected(&self) -> Vec<f64> {
        self.points
            .iter()
            .map(|p| p.value - self.wind_lever_lw1)
            .collect()
    }

    /// Returns the minimum heel angle.
    pub fn min_heel(&self) -> Option<f64> {
        self.points
            .iter()
            .map(|p| p.heel)
            .min_by(|a, b| a.partial_cmp(b).unwrap())
    }

    /// Returns the maximum heel angle.
    pub fn max_heel(&self) -> Option<f64> {
        self.points
            .iter()
            .map(|p| p.heel)
            .max_by(|a, b| a.partial_cmp(b).unwrap())
    }

    /// Finds the maximum GZ value on the curve.
    pub fn max_gz(&self) -> Option<&StabilityPoint> {
        self.points
            .iter()
            .max_by(|a, b| a.value.partial_cmp(&b.value).unwrap())
    }
}
