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

//! Z-plane mesh clipping with watertight cap generation.
//!
//! Custom Sutherland-Hodgman style clipper that:
//! 1. Clips triangles at Z = draft plane
//! 2. Reconstructs waterline edge loops
//! 3. Triangulates loops with earcutr for watertight caps

use nalgebra::{Point2, Point3};
use ordered_float::OrderedFloat;
use parry3d_f64::shape::TriMesh;
use std::collections::{HashMap, HashSet};

/// Type alias for vertex map with OrderedFloat coordinates
type VertexKey = (OrderedFloat<f64>, OrderedFloat<f64>, OrderedFloat<f64>);
type VertexMap = HashMap<VertexKey, u32>;

/// Helper to convert f64 Point3 to OrderedFloat point for hashing
fn to_ordered(p: &Point3<f64>) -> Point3<OrderedFloat<f64>> {
    Point3::new(OrderedFloat(p.x), OrderedFloat(p.y), OrderedFloat(p.z))
}

/// Canonical Z-plane intersection for bit-exact shared edge handling.
fn intersect_segment_z_canonical(p1: &Point3<f64>, p2: &Point3<f64>, z_plane: f64) -> Point3<f64> {
    let (a, b) = if p1.x < p2.x || (p1.x == p2.x && (p1.y < p2.y || (p1.y == p2.y && p1.z < p2.z)))
    {
        (p1, p2)
    } else {
        (p2, p1)
    };

    let t = (z_plane - a.z) / (b.z - a.z);
    Point3::new(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y), z_plane)
}

/// Clips a mesh at the given Z (draft) plane, returning a closed (watertight) mesh.
///
/// Returns `None` if no geometry remains below the waterline.
pub fn clip_at_waterline(mesh: &TriMesh, draft: f64) -> Option<TriMesh> {
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

    for tri in indices {
        let v0 = vertices[tri[0] as usize];
        let v1 = vertices[tri[1] as usize];
        let v2 = vertices[tri[2] as usize];

        let b0 = v0.z <= draft;
        let b1 = v1.z <= draft;
        let b2 = v2.z <= draft;

        match (b0, b1, b2) {
            // All below
            (true, true, true) => {
                let i0 = get_or_add_vertex(v0);
                let i1 = get_or_add_vertex(v1);
                let i2 = get_or_add_vertex(v2);
                new_indices.push([i0, i1, i2]);
            }
            // All above
            (false, false, false) => {}

            // 1 vertex below (keep tip)
            (true, false, false) => {
                let int1 = intersect_segment_z_canonical(&v0, &v1, draft);
                let int2 = intersect_segment_z_canonical(&v0, &v2, draft);
                let i0 = get_or_add_vertex(v0);
                let i1 = get_or_add_vertex(int1);
                let i2 = get_or_add_vertex(int2);
                new_indices.push([i0, i1, i2]);
                cut_segments.push((int1, int2));
            }
            (false, true, false) => {
                let int0 = intersect_segment_z_canonical(&v1, &v0, draft);
                let int2 = intersect_segment_z_canonical(&v1, &v2, draft);
                let i0 = get_or_add_vertex(int0);
                let i1 = get_or_add_vertex(v1);
                let i2 = get_or_add_vertex(int2);
                new_indices.push([i0, i1, i2]);
                cut_segments.push((int2, int0));
            }
            (false, false, true) => {
                let int0 = intersect_segment_z_canonical(&v2, &v0, draft);
                let int1 = intersect_segment_z_canonical(&v2, &v1, draft);
                let i0 = get_or_add_vertex(int0);
                let i1 = get_or_add_vertex(int1);
                let i2 = get_or_add_vertex(v2);
                new_indices.push([i0, i1, i2]);
                cut_segments.push((int0, int1));
            }

            // 2 vertices below (quad -> 2 triangles)
            (false, true, true) => {
                let int1 = intersect_segment_z_canonical(&v0, &v1, draft);
                let int2 = intersect_segment_z_canonical(&v0, &v2, draft);
                let iv1 = get_or_add_vertex(v1);
                let iv2 = get_or_add_vertex(v2);
                let i_int1 = get_or_add_vertex(int1);
                let i_int2 = get_or_add_vertex(int2);
                new_indices.push([i_int1, iv1, iv2]);
                new_indices.push([i_int1, iv2, i_int2]);
                cut_segments.push((int2, int1));
            }
            (true, false, true) => {
                let int0 = intersect_segment_z_canonical(&v1, &v0, draft);
                let int2 = intersect_segment_z_canonical(&v1, &v2, draft);
                let iv0 = get_or_add_vertex(v0);
                let iv2 = get_or_add_vertex(v2);
                let i_int0 = get_or_add_vertex(int0);
                let i_int2 = get_or_add_vertex(int2);
                new_indices.push([iv0, i_int0, iv2]);
                new_indices.push([i_int0, i_int2, iv2]);
                cut_segments.push((int0, int2));
            }
            (true, true, false) => {
                let int0 = intersect_segment_z_canonical(&v2, &v0, draft);
                let int1 = intersect_segment_z_canonical(&v2, &v1, draft);
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
        return None;
    }

    // --- CAP GENERATION ---
    let mut adjacency: HashMap<Point3<OrderedFloat<f64>>, Point3<OrderedFloat<f64>>> =
        HashMap::new();
    for (start, end) in &cut_segments {
        let s = to_ordered(start);
        let e = to_ordered(end);
        adjacency.insert(s, e);
    }

    let mut visited: HashSet<Point3<OrderedFloat<f64>>> = HashSet::new();
    let mut loops: Vec<Vec<Point2<f64>>> = Vec::new();
    let keys: Vec<_> = adjacency.keys().cloned().collect();

    for start_node in keys {
        if visited.contains(&start_node) {
            continue;
        }

        let mut loop_pts = Vec::new();
        let mut curr = start_node;
        let mut closed = false;
        let max_iter = adjacency.len() + 1;

        for _ in 0..max_iter {
            if visited.contains(&curr) {
                if curr == start_node {
                    closed = true;
                }
                break;
            }
            visited.insert(curr);
            loop_pts.push(Point2::new(curr.x.0, curr.y.0));

            if let Some(&next) = adjacency.get(&curr) {
                curr = next;
            } else {
                break;
            }
        }

        if closed && loop_pts.len() >= 3 {
            loops.push(loop_pts);
        }
    }

    for poly in loops {
        let mut flat_verts = Vec::with_capacity(poly.len() * 2);
        for p in &poly {
            flat_verts.push(p.x);
            flat_verts.push(p.y);
        }

        let hole_indices: Vec<usize> = vec![];
        if let Ok(indices) = earcutr::earcut(&flat_verts, &hole_indices, 2) {
            for i in (0..indices.len()).step_by(3) {
                let idx0 = indices[i];
                let idx1 = indices[i + 1];
                let idx2 = indices[i + 2];

                let p0 = Point3::new(poly[idx0].x, poly[idx0].y, draft);
                let p1 = Point3::new(poly[idx1].x, poly[idx1].y, draft);
                let p2 = Point3::new(poly[idx2].x, poly[idx2].y, draft);

                let i0 = get_or_add_vertex(p0);
                let i1 = get_or_add_vertex(p1);
                let i2 = get_or_add_vertex(p2);

                new_indices.push([i0, i1, i2]);
            }
        }
    }

    TriMesh::new(new_vertices, new_indices).ok()
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

        TriMesh::new(vertices, indices).expect("Failed to create test cube")
    }

    #[test]
    fn test_clip_cube_at_half_volume() {
        let cube = create_unit_cube();
        if let Some(clipped) = clip_at_waterline(&cube, 0.5) {
            let mass_props = clipped.mass_properties(1.0);
            let volume = mass_props.mass();
            assert!((volume - 0.5).abs() < 1e-6, "Volume was {}", volume);
            assert!((mass_props.local_com.z - 0.25).abs() < 1e-6);
        } else {
            panic!("Clipping failed");
        }
    }
}
