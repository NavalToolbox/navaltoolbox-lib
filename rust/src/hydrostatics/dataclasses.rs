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

//! Hydrostatic state dataclass.

/// Result of hydrostatic calculations at a given draft/trim/heel.
#[derive(Debug, Clone)]
pub struct HydrostaticState {
    /// Draft at reference point in meters
    pub draft: f64,
    /// Trim angle in degrees
    pub trim: f64,
    /// Heel angle in degrees
    pub heel: f64,
    /// Submerged volume in m³
    pub volume: f64,
    /// Displacement mass in kg
    pub displacement: f64,
    /// Center of buoyancy X (LCB) in meters
    pub lcb: f64,
    /// Center of buoyancy Y (TCB) in meters
    pub tcb: f64,
    /// Center of buoyancy Z (VCB) in meters
    pub vcb: f64,
    /// Waterplane area in m²
    pub waterplane_area: f64,
    /// Waterplane centroid X (LCF) in meters
    pub lcf: f64,
    /// Transverse metacentric radius BM_t in meters
    pub bmt: f64,
    /// Longitudinal metacentric radius BM_l in meters
    pub bml: f64,
    /// Transverse metacentric height GM_t in meters (requires VCG)
    pub gmt: f64,
    /// Longitudinal metacentric height GM_l in meters (requires VCG)
    pub gml: f64,
    /// Length at waterline in meters
    pub lwl: f64,
    /// Beam at waterline in meters
    pub bwl: f64,
}

impl Default for HydrostaticState {
    fn default() -> Self {
        Self {
            draft: 0.0,
            trim: 0.0,
            heel: 0.0,
            volume: 0.0,
            displacement: 0.0,
            lcb: 0.0,
            tcb: 0.0,
            vcb: 0.0,
            waterplane_area: 0.0,
            lcf: 0.0,
            bmt: 0.0,
            bml: 0.0,
            gmt: 0.0,
            gml: 0.0,
            lwl: 0.0,
            bwl: 0.0,
        }
    }
}
