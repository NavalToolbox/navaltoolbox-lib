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
    /// Centroid of submerged lateral area [x, z] in meters.
    ///
    /// Note: If the silhouettes represent only the emerged windage area
    /// (i.e. the submerged area is negligible, < 1% of emerged area),
    /// this falls back to the IMO approximation: z = T/2 (half the draft).
    pub submerged_centroid: [f64; 2],
    /// Lever arm Z per IS Code 2008 §2.3.2:
    /// vertical distance from centre of emerged area to centre of underwater lateral area (m)
    pub wind_lever_arm: f64,
    /// Waterline Z at which calculations were performed
    pub waterline_z: f64,
}

impl WindHeelingData {
    /// Create new wind heeling data.
    ///
    /// Computes `wind_lever_arm` as the exact vertical distance between the
    /// centroid of the emerged lateral area and the centroid of the submerged
    /// lateral area, per IMO 2008 IS Code §2.3.2:
    ///
    /// > Z = vertical distance from the centre of A to the centre of the
    /// > underwater lateral area.
    pub fn new(
        emerged_area: f64,
        emerged_centroid: [f64; 2],
        submerged_centroid: [f64; 2],
        waterline_z: f64,
    ) -> Self {
        let wind_lever_arm = emerged_centroid[1] - submerged_centroid[1];
        Self {
            emerged_area,
            emerged_centroid,
            submerged_centroid,
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
    /// Moulded breadth at amidships in meters (IS Code §2.13).
    ///
    /// Maximum breadth of the ship measured at the midship section (x = (AP+FP)/2)
    /// to the moulded line of the frame. `None` when hull geometry is unavailable
    /// (e.g. in unit-test contexts that construct this struct directly).
    pub moulded_breadth: Option<f64>,
}

impl CompleteStabilityResult {
    /// Create a new complete stability result.
    pub fn new(
        hydrostatics: HydrostaticState,
        gz_curve: StabilityCurve,
        wind_data: Option<WindHeelingData>,
        displacement: f64,
        cog: [f64; 3],
        moulded_breadth: Option<f64>,
    ) -> Self {
        Self {
            hydrostatics,
            gz_curve,
            wind_data,
            displacement,
            cog,
            moulded_breadth,
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
