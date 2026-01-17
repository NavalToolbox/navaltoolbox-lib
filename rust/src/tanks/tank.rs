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

//! Tank implementation.
//!
//! Represents a tank with fluid management capabilities.

use nalgebra::Point3;
use parry3d_f64::shape::{Shape, TriMesh};

use crate::hull::Hull;
use crate::mesh::clip_at_waterline;

/// Represents a tank with fluid management capabilities.
#[derive(Clone)]
pub struct Tank {
    /// Tank name
    name: String,
    /// Tank geometry (mesh)
    mesh: TriMesh,
    /// Total volume in m³
    total_volume: f64,
    /// Fluid density in kg/m³
    fluid_density: f64,
    /// Current fill level (0.0 to 1.0)
    fill_level: f64,
    /// Seawater density for FSC calculation
    water_density: f64,
    /// Bounds (xmin, xmax, ymin, ymax, zmin, zmax)
    bounds: (f64, f64, f64, f64, f64, f64),
}

impl Tank {
    /// Creates a new Tank from a mesh.
    pub fn new(name: &str, mesh: TriMesh, fluid_density: f64) -> Self {
        let mass_props = mesh.mass_properties(1.0);
        let total_volume = mass_props.mass().abs();
        let aabb = mesh.local_aabb();
        let bounds = (
            aabb.mins.x, aabb.maxs.x,
            aabb.mins.y, aabb.maxs.y,
            aabb.mins.z, aabb.maxs.z,
        );

        Self {
            name: name.to_string(),
            mesh,
            total_volume,
            fluid_density,
            fill_level: 0.0,
            water_density: 1025.0,
            bounds,
        }
    }

    /// Creates a box-shaped tank from min/max coordinates.
    pub fn from_box(
        name: &str,
        x_min: f64, x_max: f64,
        y_min: f64, y_max: f64,
        z_min: f64, z_max: f64,
        fluid_density: f64,
    ) -> Self {
        let vertices = vec![
            Point3::new(x_min, y_min, z_min),
            Point3::new(x_max, y_min, z_min),
            Point3::new(x_max, y_max, z_min),
            Point3::new(x_min, y_max, z_min),
            Point3::new(x_min, y_min, z_max),
            Point3::new(x_max, y_min, z_max),
            Point3::new(x_max, y_max, z_max),
            Point3::new(x_min, y_max, z_max),
        ];

        let indices = vec![
            [0, 2, 1], [0, 3, 2],
            [4, 5, 6], [4, 6, 7],
            [0, 1, 5], [0, 5, 4],
            [2, 3, 7], [2, 7, 6],
            [0, 4, 7], [0, 7, 3],
            [1, 2, 6], [1, 6, 5],
        ];

        let mesh = TriMesh::new(vertices, indices).expect("Failed to create tank mesh");
        Self::new(name, mesh, fluid_density)
    }

    /// Returns the tank name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns the total volume in m³.
    pub fn total_volume(&self) -> f64 {
        self.total_volume
    }

    /// Returns the fill level as a fraction (0.0 to 1.0).
    pub fn fill_level(&self) -> f64 {
        self.fill_level
    }

    /// Sets the fill level as a fraction (0.0 to 1.0).
    pub fn set_fill_level(&mut self, level: f64) {
        self.fill_level = level.clamp(0.0, 1.0);
    }

    /// Returns the fill level as a percentage (0 to 100).
    pub fn fill_percent(&self) -> f64 {
        self.fill_level * 100.0
    }

    /// Sets the fill level as a percentage (0 to 100).
    pub fn set_fill_percent(&mut self, percent: f64) {
        self.fill_level = (percent / 100.0).clamp(0.0, 1.0);
    }

    /// Returns the filled volume in m³.
    pub fn fill_volume(&self) -> f64 {
        self.total_volume * self.fill_level
    }

    /// Returns the fluid mass in kg.
    pub fn fluid_mass(&self) -> f64 {
        self.fill_volume() * self.fluid_density
    }

    /// Returns the center of gravity of the fluid [x, y, z].
    pub fn center_of_gravity(&self) -> [f64; 3] {
        if self.fill_level <= 0.0 {
            return [0.0, 0.0, 0.0];
        }

        // Find Z for this fill level
        let z = self.find_z_for_fill_level(self.fill_level);

        // Clip mesh at this Z
        if let Some(clipped) = clip_at_waterline(&self.mesh, z) {
            let mass_props = clipped.mass_properties(1.0);
            let com = mass_props.local_com;
            [com.x, com.y, com.z]
        } else {
            // Full tank
            let mass_props = self.mesh.mass_properties(1.0);
            let com = mass_props.local_com;
            [com.x, com.y, com.z]
        }
    }

    /// Returns the transverse free surface moment (I_t) in m⁴.
    pub fn free_surface_moment_t(&self) -> f64 {
        if self.fill_level <= 0.0 || self.fill_level >= 1.0 {
            return 0.0;
        }

        // Simplified: use box approximation
        let length = self.bounds.1 - self.bounds.0;
        let breadth = self.bounds.3 - self.bounds.2;
        
        // I_t = L * B³ / 12
        length * breadth.powi(3) / 12.0
    }

    /// Returns the longitudinal free surface moment (I_l) in m⁴.
    pub fn free_surface_moment_l(&self) -> f64 {
        if self.fill_level <= 0.0 || self.fill_level >= 1.0 {
            return 0.0;
        }

        let length = self.bounds.1 - self.bounds.0;
        let breadth = self.bounds.3 - self.bounds.2;
        
        // I_l = B * L³ / 12
        breadth * length.powi(3) / 12.0
    }

    /// Returns the transverse free surface correction.
    pub fn free_surface_correction_t(&self) -> f64 {
        self.free_surface_moment_t() * (self.fluid_density / self.water_density)
    }

    /// Returns the longitudinal free surface correction.
    pub fn free_surface_correction_l(&self) -> f64 {
        self.free_surface_moment_l() * (self.fluid_density / self.water_density)
    }

    /// Find Z coordinate for a given fill level using bisection.
    fn find_z_for_fill_level(&self, target_fraction: f64) -> f64 {
        let z_min = self.bounds.4;
        let z_max = self.bounds.5;
        let target_volume = self.total_volume * target_fraction;

        let tolerance = target_volume * 1e-4;
        let max_iter = 50;

        let mut low = z_min;
        let mut high = z_max;

        for _ in 0..max_iter {
            let mid = (low + high) / 2.0;

            let volume = if let Some(clipped) = clip_at_waterline(&self.mesh, mid) {
                clipped.mass_properties(1.0).mass().abs()
            } else {
                0.0
            };

            let diff = volume - target_volume;

            if diff.abs() < tolerance {
                return mid;
            }

            if diff > 0.0 {
                high = mid;
            } else {
                low = mid;
            }
        }

        (low + high) / 2.0
    }
}

impl std::fmt::Debug for Tank {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Tank")
            .field("name", &self.name)
            .field("total_volume", &self.total_volume)
            .field("fill_level", &self.fill_level)
            .field("fluid_density", &self.fluid_density)
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_box_tank_volume() {
        let tank = Tank::from_box("Test", 0.0, 10.0, 0.0, 5.0, 0.0, 2.0, 1000.0);
        
        // Volume should be 10 * 5 * 2 = 100 m³
        assert!((tank.total_volume() - 100.0).abs() < 0.1, "Volume was {}", tank.total_volume());
    }

    #[test]
    fn test_tank_fill() {
        let mut tank = Tank::from_box("Test", 0.0, 10.0, 0.0, 5.0, 0.0, 2.0, 1000.0);
        
        tank.set_fill_percent(50.0);
        
        assert!((tank.fill_level() - 0.5).abs() < 1e-6);
        assert!((tank.fill_volume() - 50.0).abs() < 0.1);
        assert!((tank.fluid_mass() - 50000.0).abs() < 100.0);
    }
}
