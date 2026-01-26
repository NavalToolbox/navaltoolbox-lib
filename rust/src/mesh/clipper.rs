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

//! Axis-aligned plane mesh clipping with watertight cap generation.
//!
//! Custom Sutherland-Hodgman style clipper that:
//! 1. Clips triangles at a specified plane (X, Y, or Z)
//! 2. Reconstructs cut edge loops
//! 3. Triangulates loops with earcutr for watertight caps

use nalgebra::{Point2, Point3};
use ordered_float::OrderedFloat;
use parry3d_f64::shape::TriMesh;
use std::collections::{HashMap, HashSet};

const EPSILON: f64 = 1e-6;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Axis {
    X,
    Y,
    Z,
}

/// Type alias for vertex map with OrderedFloat coordinates
type VertexKey = (OrderedFloat<f64>, OrderedFloat<f64>, OrderedFloat<f64>);
type VertexMap = HashMap<VertexKey, u32>;

/// Helper to convert f64 Point3 to OrderedFloat point for hashing
fn to_ordered(p: &Point3<f64>) -> Point3<OrderedFloat<f64>> {
    Point3::new(OrderedFloat(p.x), OrderedFloat(p.y), OrderedFloat(p.z))
}

fn get_component(p: &Point3<f64>, axis: Axis) -> f64 {
    match axis {
        Axis::X => p.x,
        Axis::Y => p.y,
        Axis::Z => p.z,
    }
}

/// Canonical plane intersection for bit-exact shared edge handling.
fn intersect_segment_plane_canonical(
    p1: &Point3<f64>,
    p2: &Point3<f64>,
    axis: Axis,
    plane_val: f64,
) -> Point3<f64> {
    // Sort to ensure consistency regardless of edge direction
    let (a, b) = if p1.x < p2.x || (p1.x == p2.x && (p1.y < p2.y || (p1.y == p2.y && p1.z < p2.z)))
    {
        (p1, p2)
    } else {
        (p2, p1)
    };

    let val_a = get_component(a, axis);
    let val_b = get_component(b, axis);

    // Avoid division by zero for parallel edges
    if (val_b - val_a).abs() < 1e-12 {
        return *a;
    }

    let t = (plane_val - val_a) / (val_b - val_a);

    Point3::new(
        a.x + t * (b.x - a.x),
        a.y + t * (b.y - a.y),
        a.z + t * (b.z - a.z),
    )
}

/// Clips a mesh at the given Z (draft) plane, returning a closed (watertight) mesh and the waterplane area.
/// Keeps volume below `draft`.
pub fn clip_at_waterline(mesh: &TriMesh, draft: f64) -> (Option<TriMesh>, f64) {
    clip_by_axis_aligned_plane(mesh, Axis::Z, draft, true)
}

/// Clips a mesh by an axis-aligned plane.
/// Returns (Clipped Mesh, Cap Area).
pub fn clip_by_axis_aligned_plane(
    mesh: &TriMesh,
    axis: Axis,
    value: f64,
    keep_lower: bool,
) -> (Option<TriMesh>, f64) {
    let vertices = mesh.vertices();
    let indices = mesh.indices();

    let mut new_vertices = Vec::new();
    let mut new_indices = Vec::new();
    let mut cut_segments: Vec<(Point3<f64>, Point3<f64>)> = Vec::new();

    let mut vertex_map: VertexMap = HashMap::new();

    let mut get_or_add_vertex = |p: Point3<f64>| -> u32 {
        let key = (OrderedFloat(p.x), OrderedFloat(p.y), OrderedFloat(p.z));
        if let Some(&idx) = vertex_map.get(&key) {
            idx
        } else {
            let idx = new_vertices.len() as u32;
            new_vertices.push(p);
            vertex_map.insert(key, idx);
            idx
        }
    };

    let is_inside = |p: &Point3<f64>| -> bool {
        let v = get_component(p, axis);
        if keep_lower {
            v <= value + EPSILON
        } else {
            v >= value - EPSILON
        }
    };

    for tri in indices {
        let v0 = vertices[tri[0] as usize];
        let v1 = vertices[tri[1] as usize];
        let v2 = vertices[tri[2] as usize];

        let b0 = is_inside(&v0);
        let b1 = is_inside(&v1);
        let b2 = is_inside(&v2);

        match (b0, b1, b2) {
            (true, true, true) => {
                let i0 = get_or_add_vertex(v0);
                let i1 = get_or_add_vertex(v1);
                let i2 = get_or_add_vertex(v2);
                new_indices.push([i0, i1, i2]);
            }
            (false, false, false) => {}
            (true, false, false) => {
                let int1 = intersect_segment_plane_canonical(&v0, &v1, axis, value);
                let int2 = intersect_segment_plane_canonical(&v0, &v2, axis, value);
                let i0 = get_or_add_vertex(v0);
                let i1 = get_or_add_vertex(int1);
                let i2 = get_or_add_vertex(int2);
                new_indices.push([i0, i1, i2]);
                cut_segments.push((int1, int2));
            }
            (false, true, false) => {
                let int0 = intersect_segment_plane_canonical(&v1, &v0, axis, value);
                let int2 = intersect_segment_plane_canonical(&v1, &v2, axis, value);
                let i0 = get_or_add_vertex(int0);
                let i1 = get_or_add_vertex(v1);
                let i2 = get_or_add_vertex(int2);
                new_indices.push([i0, i1, i2]);
                cut_segments.push((int2, int0));
            }
            (false, false, true) => {
                let int0 = intersect_segment_plane_canonical(&v2, &v0, axis, value);
                let int1 = intersect_segment_plane_canonical(&v2, &v1, axis, value);
                let i0 = get_or_add_vertex(int0);
                let i1 = get_or_add_vertex(int1);
                let i2 = get_or_add_vertex(v2);
                new_indices.push([i0, i1, i2]);
                cut_segments.push((int0, int1));
            }
            (false, true, true) => {
                let int1 = intersect_segment_plane_canonical(&v0, &v1, axis, value);
                let int2 = intersect_segment_plane_canonical(&v0, &v2, axis, value);
                let iv1 = get_or_add_vertex(v1);
                let iv2 = get_or_add_vertex(v2);
                let i_int1 = get_or_add_vertex(int1);
                let i_int2 = get_or_add_vertex(int2);
                new_indices.push([i_int1, iv1, iv2]);
                new_indices.push([i_int1, iv2, i_int2]);
                cut_segments.push((int2, int1));
            }
            (true, false, true) => {
                let int0 = intersect_segment_plane_canonical(&v1, &v0, axis, value);
                let int2 = intersect_segment_plane_canonical(&v1, &v2, axis, value);
                let iv0 = get_or_add_vertex(v0);
                let iv2 = get_or_add_vertex(v2);
                let i_int0 = get_or_add_vertex(int0);
                let i_int2 = get_or_add_vertex(int2);
                new_indices.push([iv0, i_int0, iv2]);
                new_indices.push([i_int0, i_int2, iv2]);
                cut_segments.push((int0, int2));
            }
            (true, true, false) => {
                let int0 = intersect_segment_plane_canonical(&v2, &v0, axis, value);
                let int1 = intersect_segment_plane_canonical(&v2, &v1, axis, value);
                let iv0 = get_or_add_vertex(v0);
                let iv1 = get_or_add_vertex(v1);
                let i_int0 = get_or_add_vertex(int0);
                let i_int1 = get_or_add_vertex(int1);
                new_indices.push([iv0, iv1, i_int1]);
                new_indices.push([iv0, i_int1, i_int0]);
                cut_segments.push((int1, int0));
            }
        }
    }

    if new_indices.is_empty() {
        return (None, 0.0);
    }

    // --- CAP GENERATION ---
    let mut total_cap_area = 0.0;
    let mut adjacency: HashMap<Point3<OrderedFloat<f64>>, Point3<OrderedFloat<f64>>> =
        HashMap::new();
    for (start, end) in &cut_segments {
        let s = to_ordered(start);
        let e = to_ordered(end);
        adjacency.insert(s, e);
    }

    let mut visited: HashSet<Point3<OrderedFloat<f64>>> = HashSet::new();
    let keys: Vec<_> = adjacency.keys().cloned().collect();

    for start_node in keys {
        if visited.contains(&start_node) {
            continue;
        }

        let mut loop_pts_2d = Vec::new();
        let mut loop_pts_3d = Vec::new();
        let mut curr = start_node;
        let mut closed = false;
        // Limit iterations to prevent infinite loops on malformed meshes
        let max_iter = adjacency.len() + 1;

        for _ in 0..max_iter {
            if visited.contains(&curr) {
                if curr == start_node {
                    closed = true;
                }
                break;
            }
            visited.insert(curr);

            let p3d = Point3::new(curr.x.0, curr.y.0, curr.z.0);
            loop_pts_3d.push(p3d);

            // Project to 2D
            let p2d = match axis {
                Axis::Z => Point2::new(curr.x.0, curr.y.0), // XY
                Axis::X => Point2::new(curr.y.0, curr.z.0), // YZ
                Axis::Y => Point2::new(curr.z.0, curr.x.0), // ZX
            };
            loop_pts_2d.push(p2d);

            if let Some(&next) = adjacency.get(&curr) {
                curr = next;
            } else {
                break;
            }
        }

        if closed && loop_pts_2d.len() >= 3 {
            let mut flat_verts = Vec::with_capacity(loop_pts_2d.len() * 2);
            // Calculate polygon area using shoelace formula
            let mut area = 0.0;
            let n = loop_pts_2d.len();
            for i in 0..n {
                let p1 = loop_pts_2d[i];
                let p2 = loop_pts_2d[(i + 1) % n];
                area += (p1.x * p2.y) - (p2.x * p1.y);

                flat_verts.push(p1.x);
                flat_verts.push(p1.y);
            }
            total_cap_area += (area / 2.0).abs();

            let hole_indices: Vec<usize> = vec![];
            if let Ok(indices) = earcutr::earcut(&flat_verts, &hole_indices, 2) {
                for i in (0..indices.len()).step_by(3) {
                    let idx0 = indices[i];
                    let idx1 = indices[i + 1];
                    let idx2 = indices[i + 2];

                    let p0 = loop_pts_3d[idx0];
                    let p1 = loop_pts_3d[idx1];
                    let p2 = loop_pts_3d[idx2];

                    // Note: Depending on keep_lower/keep_upper and axis, the winding might need reversal.
                    // Earcutr usually outputs CCW.
                    // If keep_lower=true, Cap Normal is +Axis.
                    // Loop projection:
                    // Z (XY): U=x, V=y. Normal +Z. CCW in XY is correct.
                    // X (YZ): U=y, V=z. Normal +X. CCW in YZ (viewed from +X) is correct.
                    // Y (ZX): U=z, V=x. Normal +Y. CCW in ZX (viewed from +Y) is correct?
                    //   (z, x) -> (0,0) (1,0) (0,1).
                    //   3D: (0,0,0) (0,1,0)? No:
                    //   (0,0,0) -> (z=0, x=0).
                    //   (1,0,0) -> (z=0, x=1).
                    //   (0,0,1) -> (z=1, x=0).
                    //   CCW in (z,x): (0,0)->(1,0)->(0,1).
                    //   Points: (0, ?, 0), (1, ?, 0), (0, ?, 1).
                    //   Cross product?
                    //   v1-v0 = (1, 0, 0). v2-v0 = (0, 0, 1).
                    //   (1,0,0) x (0,0,1) = (0, -1, 0).
                    //   Normal points -Y.
                    //   But we want +Y (if keep_lower=true).
                    //   So for Y axis using (z,x) projection, we might need to swap?

                    // Actually, let's just create the triangle and rely on `TriMesh` fix or consistency.
                    // Or check `earcutr` output.
                    // For Y axis specifically:
                    // If we use (z, x): Normal (-Y) comes out.
                    // If we use (x, z): Normal (+Y) comes out.
                    // X x Z = -Y.
                    // Z x X = Y.
                    // So (z, x) corresponds to Y normal.
                    // Wait: Base vectors i, j, k.
                    // k x i = j.
                    // So (z, x) is indeed the pair giving +Y normal if CCW.
                    // Example: origin, pt on Z axis (0,0,1), pt on X axis (1,0,0).
                    // 2D: (0,0) -> (1,0) -> (0,1) [CCW].
                    // 3D: (0,0,0) -> (0,0,1) -> (1,0,0).
                    // v1=(0,0,1). v2=(1,0,0).
                    // v1 x v2 = (0, 1, 0). +Y.
                    // Correct!

                    // So projection (z, x) is correct for +Y normal.

                    // However, we need to respect the `keep_lower` parameter.
                    // If keep_lower=true: Cap Normal +Axis. (Standard winding).
                    // If keep_lower=false: Cap Normal -Axis. (Reverse winding).

                    let (p0, p1, p2) = if keep_lower {
                        (p0, p1, p2)
                    } else {
                        (p0, p2, p1)
                    };

                    let i0 = get_or_add_vertex(p0);
                    let i1 = get_or_add_vertex(p1);
                    let i2 = get_or_add_vertex(p2);
                    new_indices.push([i0, i1, i2]);
                }
            }
        }
    }

    (TriMesh::new(new_vertices, new_indices).ok(), total_cap_area)
}

#[cfg(test)]
mod tests {
    use super::*;
    use parry3d_f64::shape::Shape;

    fn create_unit_cube() -> TriMesh {
        let vertices = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(0.0, 1.0, 1.0),
        ];

        let indices = vec![
            [0, 2, 1],
            [0, 3, 2], // Bottom
            [4, 5, 6],
            [4, 6, 7], // Top
            [0, 1, 5],
            [0, 5, 4], // Front
            [2, 3, 7],
            [2, 7, 6], // Back
            [0, 4, 7],
            [0, 7, 3], // Left
            [1, 2, 6],
            [1, 6, 5], // Right
        ];

        TriMesh::new(vertices, indices).expect("Failed to create test cube")
    }

    #[test]
    fn test_clip_cube_at_half_z() {
        let cube = create_unit_cube();
        let (output, area) = clip_by_axis_aligned_plane(&cube, Axis::Z, 0.5, true);
        if let Some(clipped) = output {
            let mass_props = clipped.mass_properties(1.0);
            let volume = mass_props.mass();
            assert!((volume - 0.5).abs() < 1e-6, "Volume was {}", volume);
            assert!((mass_props.local_com.z - 0.25).abs() < 1e-6);
            assert!((area - 1.0).abs() < 1e-6, "Waterplane Area was {}", area);
        } else {
            panic!("Clipping failed");
        }
    }

    #[test]
    fn test_clip_cube_at_half_x() {
        let cube = create_unit_cube();
        let (output, _area) = clip_by_axis_aligned_plane(&cube, Axis::X, 0.5, true);
        if let Some(clipped) = output {
            let mass_props = clipped.mass_properties(1.0);
            let volume = mass_props.mass();
            assert!((volume - 0.5).abs() < 1e-6, "Volume was {}", volume);
        } else {
            panic!("Clipping failed");
        }
    }
}
