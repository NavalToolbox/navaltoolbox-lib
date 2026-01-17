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

//! Tank state dataclass.

/// State of a tank at a given fill level.
#[derive(Debug, Clone)]
pub struct TankState {
    /// Fill level as fraction (0.0 to 1.0)
    pub fill_level: f64,
    /// Filled volume in m³
    pub fill_volume: f64,
    /// Fluid mass in kg
    pub fluid_mass: f64,
    /// Center of gravity [x, y, z] in meters
    pub center_of_gravity: [f64; 3],
    /// Transverse free surface moment in m⁴
    pub free_surface_moment_t: f64,
    /// Longitudinal free surface moment in m⁴
    pub free_surface_moment_l: f64,
}

impl Default for TankState {
    fn default() -> Self {
        Self {
            fill_level: 0.0,
            fill_volume: 0.0,
            fluid_mass: 0.0,
            center_of_gravity: [0.0, 0.0, 0.0],
            free_surface_moment_t: 0.0,
            free_surface_moment_l: 0.0,
        }
    }
}
