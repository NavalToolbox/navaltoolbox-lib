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

//! Complete stability calculation combining hydrostatics, GZ curve, and wind data.

use super::StabilityCurve;
use crate::hydrostatics::HydrostaticState;

/// Wind heeling data from silhouette calculations.
///
/// Contains emerged area and centroid for wind moment calculations
/// per IMO 2008 IS Code (MSC.267).
#[derive(Debug, Clone)]
pub struct WindHeelingData {
    /// Total emerged lateral area above waterline (m²)
    pub emerged_area: f64,
    /// Centroid of emerged area [x, z] in meters
    pub emerged_centroid: [f64; 2],
    /// Lever arm from waterline to centroid z coordinate (m)
    pub wind_lever_arm: f64,
    /// Waterline Z at which calculations were performed
    pub waterline_z: f64,
}

impl WindHeelingData {
    /// Create new wind heeling data.
    pub fn new(emerged_area: f64, emerged_centroid: [f64; 2], waterline_z: f64) -> Self {
        let wind_lever_arm = emerged_centroid[1] - waterline_z;
        Self {
            emerged_area,
            emerged_centroid,
            wind_lever_arm,
            waterline_z,
        }
    }
}

/// Complete stability calculation result.
///
/// Combines hydrostatic properties, GZ curve, and wind heeling data
/// for a single loading condition using the same vessel.
#[derive(Debug, Clone)]
pub struct CompleteStabilityResult {
    /// Hydrostatic state at equilibrium (draft, trim, GM0, etc.)
    pub hydrostatics: HydrostaticState,
    /// GZ stability curve for the loading condition
    pub gz_curve: StabilityCurve,
    /// Wind heeling data (if silhouettes are defined)
    pub wind_data: Option<WindHeelingData>,
    /// Displacement mass in kg
    pub displacement: f64,
    /// Center of gravity [LCG, TCG, VCG] in meters
    pub cog: [f64; 3],
}

impl CompleteStabilityResult {
    /// Create a new complete stability result.
    pub fn new(
        hydrostatics: HydrostaticState,
        gz_curve: StabilityCurve,
        wind_data: Option<WindHeelingData>,
        displacement: f64,
        cog: [f64; 3],
    ) -> Self {
        Self {
            hydrostatics,
            gz_curve,
            wind_data,
            displacement,
            cog,
        }
    }

    /// Returns the initial transverse metacentric height (GM0).
    /// This is the GM at 0° heel with free surface correction.
    pub fn gm0(&self) -> Option<f64> {
        self.hydrostatics.gmt
    }

    /// Returns the initial transverse metacentric height without FSC (dry).
    pub fn gm0_dry(&self) -> Option<f64> {
        self.hydrostatics.gmt_dry
    }

    /// Returns the maximum GZ value.
    pub fn max_gz(&self) -> Option<f64> {
        self.gz_curve.max_value().map(|p| p.value)
    }

    /// Returns the heel angle at maximum GZ.
    pub fn heel_at_max_gz(&self) -> Option<f64> {
        self.gz_curve.max_value().map(|p| p.heel)
    }

    /// Returns true if wind heeling data is available.
    pub fn has_wind_data(&self) -> bool {
        self.wind_data.is_some()
    }
}
