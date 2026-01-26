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
use std::path::Path;

use crate::hull::Hull;
use crate::mesh::{clip_at_waterline, clip_by_axis_aligned_plane, load_stl, load_vtk, Axis};

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
            aabb.mins.x,
            aabb.maxs.x,
            aabb.mins.y,
            aabb.maxs.y,
            aabb.mins.z,
            aabb.maxs.z,
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

    /// Creates a Tank from a file (STL or VTK).
    pub fn from_file<P: AsRef<Path>>(path: P, fluid_density: f64) -> Result<Self, String> {
        let path = path.as_ref();
        if !path.exists() {
            return Err(format!("File not found: {}", path.display()));
        }

        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.to_lowercase())
            .unwrap_or_default();

        let mesh = match ext.as_str() {
            "stl" => load_stl(path).map_err(|e| format!("Failed to load STL: {}", e))?,
            "vtk" | "vtp" | "vtu" => {
                load_vtk(path).map_err(|e| format!("Failed to load VTK: {}", e))?
            }
            _ => return Err(format!("Unsupported file format: {}", ext)),
        };

        if mesh.vertices().is_empty() {
            return Err("Loaded mesh has no vertices".to_string());
        }

        let name = path.file_stem().and_then(|s| s.to_str()).unwrap_or("Tank");

        Ok(Self::new(name, mesh, fluid_density))
    }

    /// Creates a Tank as the intersection of a box with a hull geometry.
    ///
    /// # Arguments
    /// * `name` - Tank name
    /// * `hull` - Hull to intersect with
    /// * `x_min`, `x_max` - Longitudinal bounds
    /// * `y_min`, `y_max` - Transverse bounds
    /// * `z_min`, `z_max` - Vertical bounds
    /// * `fluid_density` - Fluid density (kg/m³)
    #[allow(clippy::too_many_arguments)]
    pub fn from_box_hull_intersection(
        name: &str,
        hull: &Hull,
        x_min: f64,
        x_max: f64,
        y_min: f64,
        y_max: f64,
        z_min: f64,
        z_max: f64,
        fluid_density: f64,
    ) -> Result<Self, String> {
        let mut mesh = hull.mesh().clone();

        // Check for invalid bounds
        if x_min >= x_max || y_min >= y_max || z_min >= z_max {
            return Err("Invalid box dimensions: min must be less than max".to_string());
        }

        // Apply clips sequentially
        // Note: The order should not strictly matter for the final result,
        // but checking for emptiness earlier is better.

        mesh = clip_by_axis_aligned_plane(&mesh, Axis::X, x_min, false).0
            .ok_or_else(|| format!("Tank '{}': No geometry remaining after X_min ({:.2}) clip. Box is fully aft of hull?", name, x_min))?;

        mesh = clip_by_axis_aligned_plane(&mesh, Axis::X, x_max, true).0
            .ok_or_else(|| format!("Tank '{}': No geometry remaining after X_max ({:.2}) clip. Box is fully fwd of hull?", name, x_max))?;

        mesh = clip_by_axis_aligned_plane(&mesh, Axis::Y, y_min, false).0
            .ok_or_else(|| format!("Tank '{}': No geometry remaining after Y_min ({:.2}) clip. Box is fully stbd of hull?", name, y_min))?;

        mesh = clip_by_axis_aligned_plane(&mesh, Axis::Y, y_max, true).0
            .ok_or_else(|| format!("Tank '{}': No geometry remaining after Y_max ({:.2}) clip. Box is fully port of hull?", name, y_max))?;

        mesh = clip_by_axis_aligned_plane(&mesh, Axis::Z, z_min, false).0
            .ok_or_else(|| format!("Tank '{}': No geometry remaining after Z_min ({:.2}) clip. Box is fully below hull?", name, z_min))?;

        mesh = clip_by_axis_aligned_plane(&mesh, Axis::Z, z_max, true).0
            .ok_or_else(|| format!("Tank '{}': No geometry remaining after Z_max ({:.2}) clip. Box is fully above hull?", name, z_max))?;

        let mass_props = mesh.mass_properties(1.0);
        let volume = mass_props.mass().abs();

        if volume < 1e-6 {
            return Err(format!("Tank '{}': Intersection resulted in near-zero volume ({:.2e} m³). Check that box overlaps hull interior.", name, volume));
        }

        Ok(Self::new(name, mesh, fluid_density))
    }

    /// Creates a box-shaped tank from min/max coordinates.
    #[allow(clippy::too_many_arguments)]
    pub fn from_box(
        name: &str,
        x_min: f64,
        x_max: f64,
        y_min: f64,
        y_max: f64,
        z_min: f64,
        z_max: f64,
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
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [2, 3, 7],
            [2, 7, 6],
            [0, 4, 7],
            [0, 7, 3],
            [1, 2, 6],
            [1, 6, 5],
        ];

        let mesh = TriMesh::new(vertices, indices).expect("Failed to create tank mesh");
        Self::new(name, mesh, fluid_density)
    }

    /// Returns the tank name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Sets the tank name.
    pub fn set_name(&mut self, name: &str) {
        self.name = name.to_string();
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

    /// Returns the fluid density in kg/m³.
    pub fn fluid_density(&self) -> f64 {
        self.fluid_density
    }

    /// Returns the filled volume in m³.
    pub fn fill_volume(&self) -> f64 {
        self.total_volume * self.fill_level
    }

    /// Returns the fluid mass in kg.
    pub fn fluid_mass(&self) -> f64 {
        self.fill_volume() * self.fluid_density
    }

    /// Returns the mesh of the fluid at the current fill level.
    pub fn get_fluid_mesh(&self) -> Option<TriMesh> {
        self.get_fluid_mesh_at(0.0, 0.0)
    }

    /// Returns the mesh of the fluid at a specific heel and trim.
    pub fn get_fluid_mesh_at(&self, heel: f64, trim: f64) -> Option<TriMesh> {
        if self.fill_level <= 1e-6 {
            return None;
        }
        if self.fill_level >= 1.0 - 1e-6 {
            return Some(self.mesh.clone());
        }

        // 1. Transform mesh to align water parallel to XY plane
        let pivot = nalgebra::Point3::new(0.0, 0.0, 0.0);
        let transformed_mesh = crate::mesh::transform_mesh(&self.mesh, heel, trim, pivot);

        // 2. Find Z level for this volume in transformed orientation
        let z = self.find_z_for_mesh(&transformed_mesh, self.total_volume * self.fill_level);

        // 3. Clip
        if let Some(clipped) = clip_at_waterline(&transformed_mesh, z).0 {
            // 4. Inverse transform back to ship frame
            use nalgebra::Rotation3;
            let roll = heel.to_radians();
            let pitch = trim.to_radians();
            let rotation = Rotation3::from_euler_angles(roll, pitch, 0.0);
            let pivot_pt = nalgebra::Point3::new(0.0, 0.0, 0.0);

            // We need to inverse transform the result mesh.
            // Since mesh crate doesn't have `inverse_transform_mesh`, we can reuse transform_mesh
            // but with inverse rotation.
            // Note: transform_mesh(mesh, roll, pitch, pivot) applies R * (v - p) + p
            // We want R^-1 * (v - p) + p.
            // Since R(roll, pitch) isn't just single axis rotations, we can't just pass -roll, -pitch easily
            // if order matters. However, for small angles or standard Euler:
            // R = Ry(pitch) * Rx(roll).
            // R^-1 = Rx(-roll) * Ry(-pitch).
            // Our transform_mesh might assume a specific order.
            // A better way is to implement inverse transform or use the fact that we can just
            // manually transform vertices here.

            let inverse_rotation = rotation.inverse();
            let new_vertices: Vec<Point3<f64>> = clipped
                .vertices()
                .iter()
                .map(|v| {
                    let p = Point3::new(v.x, v.y, v.z);
                    pivot_pt + inverse_rotation * (p - pivot_pt)
                })
                .collect();

            Some(TriMesh::new(new_vertices, clipped.indices().to_vec()).unwrap())
        } else {
            None
        }
    }

    /// Returns the underlying TriMesh of the tank container.
    pub fn mesh(&self) -> &TriMesh {
        &self.mesh
    }

    /// Returns the center of gravity of the fluid [x, y, z].
    pub fn center_of_gravity(&self) -> [f64; 3] {
        self.center_of_gravity_at(0.0, 0.0)
    }

    /// Returns the center of gravity of the fluid at a specific heel and trim.
    pub fn center_of_gravity_at(&self, heel: f64, trim: f64) -> [f64; 3] {
        if self.fill_level <= 0.0 {
            return [0.0, 0.0, 0.0];
        }

        // 1. Transform mesh to align water parallel to XY plane
        let pivot = nalgebra::Point3::new(0.0, 0.0, 0.0);
        let transformed_mesh = crate::mesh::transform_mesh(&self.mesh, heel, trim, pivot);

        // 2. Find Z level for this volume in transformed orientation
        let z = self.find_z_for_mesh(&transformed_mesh, self.total_volume * self.fill_level);

        // 3. Calculate centroid in transformed frame
        let transformed_cog = if let Some(clipped) = clip_at_waterline(&transformed_mesh, z).0 {
            let mass_props = clipped.mass_properties(1.0);
            let com = mass_props.local_com;
            nalgebra::Point3::new(com.x, com.y, com.z)
        } else {
            // Full tank
            let mass_props = transformed_mesh.mass_properties(1.0);
            let com = mass_props.local_com;
            nalgebra::Point3::new(com.x, com.y, com.z)
        };

        // 4. Inverse transform back to ship frame
        use nalgebra::Rotation3;
        let roll = heel.to_radians();
        let pitch = trim.to_radians();

        // Rotation used in transform_mesh matches:
        let rotation = Rotation3::from_euler_angles(roll, pitch, 0.0);
        let inverse_rotation = rotation.inverse();

        let original_cog = inverse_rotation * (transformed_cog - pivot) + pivot.coords;

        [original_cog.x, original_cog.y, original_cog.z]
    }

    /// Helper to find Z level for a specific mesh
    fn find_z_for_mesh(&self, mesh: &TriMesh, target_volume: f64) -> f64 {
        let aabb = mesh.local_aabb();
        let z_min = aabb.mins.z;
        let z_max = aabb.maxs.z;

        let tolerance = target_volume.max(0.001) * 1e-4;
        let max_iter = 20;

        let mut low = z_min;
        let mut high = z_max;

        for _ in 0..max_iter {
            let mid = (low + high) / 2.0;

            let volume = if let Some(clipped) = clip_at_waterline(mesh, mid).0 {
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

    /// Returns the transverse free surface moment (I_t) in m⁴.
    pub fn free_surface_moment_t(&self) -> f64 {
        if self.fill_level <= 0.0 || self.fill_level >= 1.0 {
            return 0.0;
        }

        // Simplified: use box approximation of bounding box at fill level
        // TODO: Implement exact waterplane inertia calculation
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
        assert!(
            (tank.total_volume() - 100.0).abs() < 0.1,
            "Volume was {}",
            tank.total_volume()
        );
    }

    #[test]
    fn test_tank_fill() {
        let mut tank = Tank::from_box("Test", 0.0, 10.0, 0.0, 5.0, 0.0, 2.0, 1000.0);

        tank.set_fill_percent(50.0);

        assert!((tank.fill_level() - 0.5).abs() < 1e-6);
        assert!((tank.fill_volume() - 50.0).abs() < 0.1);
        assert!((tank.fluid_mass() - 50000.0).abs() < 100.0);
    }

    #[test]
    fn test_from_box_hull_intersection() {
        use crate::hull::Hull;
        // Create a large box hull: L=20, B=10, D=10
        // Bounds: X[0,20], Y[-5,5], Z[0,10]
        let hull = Hull::from_box(20.0, 10.0, 10.0);

        // Intersect with smaller box inside: X[5,15], Y[-2,2], Z[0,5]
        // Expected Volume: 10 * 4 * 5 = 200
        let tank = Tank::from_box_hull_intersection(
            "InnerTank",
            &hull,
            5.0,
            15.0,
            -2.0,
            2.0,
            0.0,
            5.0,
            1025.0,
        )
        .expect("Intersection failed");

        assert!(
            (tank.total_volume() - 200.0).abs() < 1e-3,
            "Volume was {}",
            tank.total_volume()
        );

        // Test non-intersecting (Outside Hull)
        let res = Tank::from_box_hull_intersection(
            "Outside", &hull, 30.0, 40.0, 0.0, 1.0, 0.0, 1.0, 1025.0,
        );
        assert!(res.is_err());
    }
}
