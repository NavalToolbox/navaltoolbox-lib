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

//! Waterplane property calculations.
//!
//! This module provides functions to calculate waterplane properties such as:
//! - Waterplane area
//! - Centroid (LCF, TCF)
//! - Second moments of area (for metacentric calculations)

use nalgebra::Point3;
use parry3d_f64::shape::TriMesh;

/// Properties of a waterplane at a given draft/trim/heel
#[derive(Debug, Clone)]
pub struct WaterplaneProperties {
    /// Waterplane area in m²
    pub area: f64,

    /// Centroid of waterplane [LCF, TCF] in meters
    pub centroid: [f64; 2],

    /// Second moment of area about transverse axis (for rolling)
    /// I_t = ∫∫ y² dA
    pub i_transverse: f64,

    /// Second moment of area about longitudinal axis (for pitching)
    /// I_l = ∫∫ x² dA
    pub i_longitudinal: f64,

    /// Minimum X coordinate of the waterplane (aft-most point)
    pub min_x: f64,
    /// Maximum X coordinate of the waterplane (fore-most point)
    pub max_x: f64,
    /// Minimum Y coordinate of the waterplane
    pub min_y: f64,
    /// Maximum Y coordinate of the waterplane
    pub max_y: f64,
}

/// Calculate waterplane properties from a mesh at a given draft
pub fn calculate_waterplane_properties(mesh: &TriMesh, draft: f64) -> Option<WaterplaneProperties> {
    // Extract waterline contour(s) - may have multiple for catamarans
    let contours = extract_waterline_contours(mesh, draft)?;

    if contours.is_empty() {
        return None;
    }

    // Calculate properties for all contours combined
    let mut total_area = 0.0;
    let mut moment_x = 0.0; // For LCF: ∫ x dA
    let mut moment_y = 0.0; // For TCF: ∫ y dA
    let mut i_xx = 0.0; // ∫∫ y² dA (about x-axis, for transverse)
    let mut i_yy = 0.0; // ∫∫ x² dA (about y-axis, for longitudinal)

    let mut min_x = f64::MAX;
    let mut max_x = f64::MIN;
    let mut min_y = f64::MAX;
    let mut max_y = f64::MIN;

    for contour in &contours {
        for p in contour {
            min_x = min_x.min(p.x);
            max_x = max_x.max(p.x);
            min_y = min_y.min(p.y);
            max_y = max_y.max(p.y);
        }

        // Triangulate this contour using fan triangulation from first point
        let triangles = triangulate_polygon(contour);

        for tri in triangles {
            let [p0, p1, p2] = tri;

            // Triangle area (2D in xy plane)
            let area = triangle_area_2d(&p0, &p1, &p2);
            if area <= 0.0 {
                continue;
            }

            // Triangle centroid
            let cx = (p0.x + p1.x + p2.x) / 3.0;
            let cy = (p0.y + p1.y + p2.y) / 3.0;

            // Accumulate first moments
            total_area += area;
            moment_x += area * cx;
            moment_y += area * cy;

            // Second moments about origin
            // For a triangle: I = (1/12) * A * (p0² + p1² + p2² + p0*p1 + p1*p2 + p2*p0)
            let i_xx_tri = area
                * (p0.y * p0.y
                    + p1.y * p1.y
                    + p2.y * p2.y
                    + p0.y * p1.y
                    + p1.y * p2.y
                    + p2.y * p0.y)
                / 6.0;
            let i_yy_tri = area
                * (p0.x * p0.x
                    + p1.x * p1.x
                    + p2.x * p2.x
                    + p0.x * p1.x
                    + p1.x * p2.x
                    + p2.x * p0.x)
                / 6.0;

            i_xx += i_xx_tri;
            i_yy += i_yy_tri;
        }
    }

    if total_area <= 1e-9 {
        return None;
    }

    // Centroid
    let lcf = moment_x / total_area;
    let tcf = moment_y / total_area;

    // Transfer second moments to centroidal axes using parallel axis theorem
    // I_c = I_o - A * d²
    let i_transverse = i_xx - total_area * tcf * tcf;
    let i_longitudinal = i_yy - total_area * lcf * lcf;

    Some(WaterplaneProperties {
        area: total_area,
        centroid: [lcf, tcf],
        i_transverse: i_transverse.max(0.0),
        i_longitudinal: i_longitudinal.max(0.0),
        min_x,
        max_x,
        min_y,
        max_y,
    })
}

/// Extract waterline contour(s) from mesh at given draft
fn extract_waterline_contours(mesh: &TriMesh, draft: f64) -> Option<Vec<Vec<Point3<f64>>>> {
    let vertices = mesh.vertices();
    let indices = mesh.indices();

    if vertices.is_empty() || indices.is_empty() {
        return None;
    }

    let tolerance = 1e-6;

    // Build edges that cross the waterline with their endpoints
    let mut crossing_segments: Vec<(Point3<f64>, Point3<f64>)> = Vec::new();

    for tri_indices in indices {
        let v0 = vertices[tri_indices[0] as usize];
        let v1 = vertices[tri_indices[1] as usize];
        let v2 = vertices[tri_indices[2] as usize];

        // Check each edge of the triangle
        let triangle_edges = [(v0, v1, 0), (v1, v2, 1), (v2, v0, 2)];

        for (p0, p1, _edge_id) in triangle_edges {
            let z0 = p0.z;
            let z1 = p1.z;

            // Check if edge crosses the waterline (one above, one below)
            if (z0 - draft) * (z1 - draft) < -tolerance {
                // Edge crosses waterline - compute intersection point
                let t = (draft - z0) / (z1 - z0);
                let x = p0.x + t * (p1.x - p0.x);
                let y = p0.y + t * (p1.y - p0.y);

                let intersection = Point3::new(x, y, draft);

                // Store as a directed segment
                // We need to determine which direction based on triangle winding
                if z0 < draft {
                    // p0 is below, p1 is above
                    crossing_segments.push((intersection, intersection));
                } else {
                    // p0 is above, p1 is below
                    crossing_segments.push((intersection, intersection));
                }
            }
        }
    }

    if crossing_segments.is_empty() {
        return None;
    }

    // Collect unique intersection points
    let mut points: Vec<Point3<f64>> = crossing_segments.iter().map(|(p, _)| *p).collect();

    // Remove duplicates
    points.sort_by(|a, b| {
        a.x.partial_cmp(&b.x)
            .unwrap()
            .then(a.y.partial_cmp(&b.y).unwrap())
    });
    points.dedup_by(|a, b| (a.x - b.x).abs() < tolerance && (a.y - b.y).abs() < tolerance);

    if points.is_empty() {
        return None;
    }

    // For now, form a single contour by sorting points by angle from centroid
    // This works for convex and simple concave shapes
    let cx: f64 = points.iter().map(|p| p.x).sum::<f64>() / points.len() as f64;
    let cy: f64 = points.iter().map(|p| p.y).sum::<f64>() / points.len() as f64;

    points.sort_by(|a, b| {
        let angle_a = (a.y - cy).atan2(a.x - cx);
        let angle_b = (b.y - cy).atan2(b.x - cx);
        angle_a.partial_cmp(&angle_b).unwrap()
    });

    Some(vec![points])
}

/// Separate points into connected contours
#[allow(dead_code)]
fn separate_into_contours(mut points: Vec<Point3<f64>>, _tolerance: f64) -> Vec<Vec<Point3<f64>>> {
    if points.is_empty() {
        return vec![];
    }

    // Simple approach: sort by angle from centroid
    // This works for convex shapes and simple concave cases
    let mut contours = vec![];

    while !points.is_empty() {
        // Take first point as seed
        let seed = points.remove(0);
        let mut contour = vec![seed];

        // Find centroid of remaining points to sort by angle
        if points.is_empty() {
            contours.push(contour);
            break;
        }

        let cx: f64 = points.iter().map(|p| p.x).sum::<f64>() / points.len() as f64;
        let cy: f64 = points.iter().map(|p| p.y).sum::<f64>() / points.len() as f64;

        // Sort by angle from centroid
        points.sort_by(|a, b| {
            let angle_a = (a.y - cy).atan2(a.x - cx);
            let angle_b = (b.y - cy).atan2(b.x - cx);
            angle_a.partial_cmp(&angle_b).unwrap()
        });

        contour.append(&mut points);
        contours.push(contour);
    }

    contours
}

/// Triangulate a 2D polygon using fan triangulation
fn triangulate_polygon(contour: &[Point3<f64>]) -> Vec<[Point3<f64>; 3]> {
    if contour.len() < 3 {
        return vec![];
    }

    let mut triangles = Vec::new();
    let first = contour[0];

    // Fan triangulation from first vertex
    for i in 1..contour.len() - 1 {
        triangles.push([first, contour[i], contour[i + 1]]);
    }

    triangles
}

/// Calculate 2D triangle area (ignoring z coordinate)
fn triangle_area_2d(p0: &Point3<f64>, p1: &Point3<f64>, p2: &Point3<f64>) -> f64 {
    // 2D cross product: |AB × AC| / 2
    let ab_x = p1.x - p0.x;
    let ab_y = p1.y - p0.y;
    let ac_x = p2.x - p0.x;
    let ac_y = p2.y - p0.y;

    0.5 * (ab_x * ac_y - ab_y * ac_x).abs()
}

#[cfg(test)]
mod tests {
    use super::*;
    use parry3d_f64::math::Point;
    use parry3d_f64::shape::TriMesh;

    fn create_box_mesh(loa: f64, boa: f64, depth: f64) -> TriMesh {
        let hb = boa / 2.0;
        let vertices = vec![
            Point::new(0.0, -hb, 0.0),
            Point::new(loa, -hb, 0.0),
            Point::new(loa, hb, 0.0),
            Point::new(0.0, hb, 0.0),
            Point::new(0.0, -hb, depth),
            Point::new(loa, -hb, depth),
            Point::new(loa, hb, depth),
            Point::new(0.0, hb, depth),
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
        TriMesh::new(vertices, indices).unwrap()
    }

    #[test]
    fn test_box_waterplane_area() {
        // 10m × 10m × 10m box
        let mesh = create_box_mesh(10.0, 10.0, 10.0);
        let draft = 5.0;
        let wp = calculate_waterplane_properties(&mesh, draft).unwrap();

        // Expected area = L × B = 100 m²
        let expected_area = 100.0;
        assert!(
            (wp.area - expected_area).abs() < 1.0,
            "Area: got {:.2}, expected {:.2}",
            wp.area,
            expected_area
        );
    }

    #[test]
    fn test_box_waterplane_centroid() {
        let mesh = create_box_mesh(10.0, 10.0, 10.0);
        let draft = 5.0;
        let wp = calculate_waterplane_properties(&mesh, draft).unwrap();

        // Expected LCF = 5.0 (center at x=5)
        // Expected TCF = 0.0 (symmetric)
        assert!(
            (wp.centroid[0] - 5.0).abs() < 0.1,
            "LCF: got {:.2}, expected 5.0",
            wp.centroid[0]
        );
        assert!(
            wp.centroid[1].abs() < 0.1,
            "TCF: got {:.2}, expected 0.0",
            wp.centroid[1]
        );
    }
}
