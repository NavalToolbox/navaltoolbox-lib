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

/// Options for tank handling in hydrostatic calculations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TankOptions {
    /// Include tank fluid mass in displacement calculation
    pub include_mass: bool,
    /// Include Free Surface Moment correction in GM calculations
    pub include_fsm: bool,
}

impl Default for TankOptions {
    fn default() -> Self {
        Self {
            include_mass: false,
            include_fsm: true,
        }
    }
}

impl TankOptions {
    /// Create options with all tank effects disabled
    pub fn none() -> Self {
        Self {
            include_mass: false,
            include_fsm: false,
        }
    }

    /// Create options with all tank effects enabled
    pub fn all() -> Self {
        Self {
            include_mass: true,
            include_fsm: true,
        }
    }

    /// Create options with only mass included (no FSM)
    pub fn mass_only() -> Self {
        Self {
            include_mass: true,
            include_fsm: false,
        }
    }

    /// Create options with only FSM correction (no mass)
    pub fn fsm_only() -> Self {
        Self {
            include_mass: false,
            include_fsm: true,
        }
    }
}

/// Result of hydrostatic calculations at a given draft/trim/heel.
#[derive(Debug, Clone)]
pub struct HydrostaticState {
    /// Draft at reference point in meters
    pub draft: f64,
    /// Trim angle in degrees
    pub trim: f64,
    /// Heel angle in degrees
    pub heel: f64,
    /// Draft at Aft Perpendicular in meters
    pub draft_ap: f64,
    /// Draft at Forward Perpendicular in meters
    pub draft_fp: f64,
    /// Draft at Midship (MP) in meters
    pub draft_mp: f64,

    /// Submerged volume in m³
    pub volume: f64,
    /// Total displacement mass in kg (Volume * water_density).
    /// This represents the total weight of the floating system (Vessel + Tanks).
    pub displacement: f64,
    /// Vessel displacement mass in kg (Total - Tank Contents).
    /// Corresponds to the input displacement (Lightship + Deadweight excluding tanks).
    pub vessel_displacement: f64,
    /// Tank fluid mass in kg (sum of all tank fluid masses).
    /// If TankOptions.include_mass is false, this will be 0.0.
    pub tank_displacement: f64,

    /// Center of buoyancy [LCB, TCB, VCB] in meters
    pub cob: [f64; 3],

    /// Center of gravity [LCG, TCG, VCG] in meters (if specified)
    pub cog: Option<[f64; 3]>,

    /// Waterplane area in m²
    pub waterplane_area: f64,
    /// Waterplane centroid X (LCF) in meters
    pub lcf: f64,

    /// Transverse metacentric radius BM_t in meters
    pub bmt: f64,
    /// Longitudinal metacentric radius BM_l in meters
    pub bml: f64,

    /// Transverse metacentric height GM_t in meters (requires VCG)
    /// Includes free surface correction from tanks (wet - conservative)
    pub gmt: Option<f64>,
    /// Longitudinal metacentric height GM_l in meters (requires VCG)
    /// Includes free surface correction from tanks (wet - conservative)
    pub gml: Option<f64>,

    /// Transverse metacentric height without free surface correction (dry)
    pub gmt_dry: Option<f64>,
    /// Longitudinal metacentric height without free surface correction (dry)
    pub gml_dry: Option<f64>,

    /// Length at waterline in meters
    pub lwl: f64,
    /// Beam at waterline in meters
    pub bwl: f64,

    /// Wetted surface area in m²
    pub wetted_surface_area: f64,
    /// Midship section area in m²
    pub midship_area: f64,
    /// Midship section coefficient
    pub cm: f64,
    /// Block coefficient
    pub cb: f64,
    /// Prismatic coefficient
    pub cp: f64,

    /// Transverse free surface correction (meters)
    pub free_surface_correction_t: f64,
    /// Longitudinal free surface correction (meters)
    pub free_surface_correction_l: f64,

    /// 6x6 Hydrostatic stiffness matrix (flattened row-major)
    pub stiffness_matrix: [f64; 36],

    /// Length overall submerged in meters
    pub los: f64,

    /// Sectional area curve: (x_position, area) pairs from AP to FP
    pub sectional_areas: Vec<(f64, f64)>,

    /// Freeboard at deck edge (min distance from waterline to deck) in meters
    pub freeboard: Option<f64>,
}

impl HydrostaticState {
    /// Returns the longitudinal center of buoyancy (LCB) in meters
    pub fn lcb(&self) -> f64 {
        self.cob[0]
    }

    /// Returns the transverse center of buoyancy (TCB) in meters
    pub fn tcb(&self) -> f64 {
        self.cob[1]
    }

    /// Returns the vertical center of buoyancy (VCB) in meters
    pub fn vcb(&self) -> f64 {
        self.cob[2]
    }

    /// Returns the longitudinal center of gravity (LCG) if specified
    pub fn lcg(&self) -> Option<f64> {
        self.cog.map(|c| c[0])
    }

    /// Returns the transverse center of gravity (TCG) if specified
    pub fn tcg(&self) -> Option<f64> {
        self.cog.map(|c| c[1])
    }

    /// Returns the vertical center of gravity (VCG) if specified
    pub fn vcg(&self) -> Option<f64> {
        self.cog.map(|c| c[2])
    }
}

impl Default for HydrostaticState {
    fn default() -> Self {
        Self {
            draft: 0.0,
            trim: 0.0,
            heel: 0.0,
            volume: 0.0,
            displacement: 0.0,
            vessel_displacement: 0.0,
            tank_displacement: 0.0,
            cob: [0.0, 0.0, 0.0],
            cog: None,
            waterplane_area: 0.0,
            lcf: 0.0,
            bmt: 0.0,
            bml: 0.0,
            gmt: None,
            gml: None,
            gmt_dry: None,
            gml_dry: None,
            lwl: 0.0,
            bwl: 0.0,
            wetted_surface_area: 0.0,
            midship_area: 0.0,
            cm: 0.0,
            cb: 0.0,
            cp: 0.0,
            free_surface_correction_t: 0.0,
            free_surface_correction_l: 0.0,
            stiffness_matrix: [0.0; 36],
            los: 0.0,
            draft_ap: 0.0,
            draft_fp: 0.0,
            draft_mp: 0.0,
            sectional_areas: Vec::new(),
            freeboard: None,
        }
    }
}
