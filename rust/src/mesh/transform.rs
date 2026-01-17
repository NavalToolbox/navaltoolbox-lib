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

//! Mesh transformation utilities.

use nalgebra::{Point3, Rotation3, Vector3};

use parry3d_f64::shape::TriMesh;

/// Returns the axis-aligned bounding box of the mesh.
/// Returns (min_x, max_x, min_y, max_y, min_z, max_z)
pub fn get_bounds(mesh: &TriMesh) -> (f64, f64, f64, f64, f64, f64) {
    let aabb = mesh.local_aabb();
    (
        aabb.mins.x,
        aabb.maxs.x,
        aabb.mins.y,
        aabb.maxs.y,
        aabb.mins.z,
        aabb.maxs.z,
    )
}

/// Transforms a mesh by applying rotation around a pivot point.
///
/// # Arguments
/// * `mesh` - The mesh to transform
/// * `heel` - Heel angle in degrees (rotation around X axis)
/// * `trim` - Trim angle in degrees (rotation around Y axis)
/// * `pivot` - Pivot point for rotation
pub fn transform_mesh(mesh: &TriMesh, heel: f64, trim: f64, pivot: Point3<f64>) -> TriMesh {
    let heel_rad = heel.to_radians();
    let trim_rad = trim.to_radians();

    // Create rotation: first X (heel), then Y (trim)
    let rot_x = Rotation3::from_axis_angle(&Vector3::x_axis(), heel_rad);
    let rot_y = Rotation3::from_axis_angle(&Vector3::y_axis(), trim_rad);
    let rotation = rot_y * rot_x;

    // Transform all vertices
    let new_vertices: Vec<Point3<f64>> = mesh
        .vertices()
        .iter()
        .map(|v| {
            let relative = v - pivot;
            let rotated = rotation * relative;
            pivot + rotated
        })
        .collect();

    // Keep same indices
    let indices: Vec<[u32; 3]> = mesh
        .indices()
        .iter()
        .map(|idx| [idx[0], idx[1], idx[2]])
        .collect();

    TriMesh::new(new_vertices, indices).expect("Failed to create transformed mesh")
}

/// Transforms a point by rotation around a pivot.
pub fn transform_point(
    point: Point3<f64>,
    heel: f64,
    trim: f64,
    pivot: Point3<f64>,
) -> Point3<f64> {
    let heel_rad = heel.to_radians();
    let trim_rad = trim.to_radians();

    let rot_x = Rotation3::from_axis_angle(&Vector3::x_axis(), heel_rad);
    let rot_y = Rotation3::from_axis_angle(&Vector3::y_axis(), trim_rad);
    let rotation = rot_y * rot_x;

    let relative = point - pivot;
    let rotated = rotation * relative;
    pivot + rotated
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_transform_point_identity() {
        let point = Point3::new(1.0, 2.0, 3.0);
        let pivot = Point3::new(0.0, 0.0, 0.0);
        let result = transform_point(point, 0.0, 0.0, pivot);
        assert!((result - point).norm() < 1e-10);
    }

    #[test]
    fn test_transform_point_heel_90() {
        let point = Point3::new(0.0, 1.0, 0.0);
        let pivot = Point3::new(0.0, 0.0, 0.0);
        let result = transform_point(point, 90.0, 0.0, pivot);
        assert!((result.x).abs() < 1e-6);
        assert!((result.y).abs() < 1e-6);
        assert!((result.z - 1.0).abs() < 1e-6);
    }
}
