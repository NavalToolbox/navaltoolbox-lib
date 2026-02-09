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
    // Extract waterline contour(s)
    let contours = extract_waterline_contours(mesh, draft)?;

    if contours.is_empty() {
        return None;
    }

    // Accumulators for Green's theorem
    let mut total_area = 0.0;
    let mut moment_x = 0.0; // ∫ x dA
    let mut moment_y = 0.0; // ∫ y dA
    let mut i_xx = 0.0; // ∫ y² dA (about x-axis)
    let mut i_yy = 0.0; // ∫ x² dA (about y-axis)

    let mut min_x = f64::MAX;
    let mut max_x = f64::MIN;
    let mut min_y = f64::MAX;
    let mut max_y = f64::MIN;

    for contour in &contours {
        if contour.len() < 3 {
            continue;
        }

        // Green's theorem integration
        for i in 0..contour.len() {
            let p0 = contour[i];
            let p1 = contour[(i + 1) % contour.len()];

            min_x = min_x.min(p0.x);
            max_x = max_x.max(p0.x);
            min_y = min_y.min(p0.y);
            max_y = max_y.max(p0.y);

            // Cross product term for area (x_i * y_{i+1} - x_{i+1} * y_i)
            let cp = p0.x * p1.y - p1.x * p0.y;

            // Area accumulator
            total_area += cp;

            // First moments
            moment_y += (p0.y + p1.y) * cp; // Cy * 6A
            moment_x += (p0.x + p1.x) * cp; // Cx * 6A

            // Second moments
            // I_x (about X axis) = ∫ y² dA = 1/12 * Σ (y_i² + y_i*y_{i+1} + y_{i+1}²) * cp
            i_xx += (p0.y * p0.y + p0.y * p1.y + p1.y * p1.y) * cp;

            // I_y (about Y axis) = ∫ x² dA = 1/12 * Σ (x_i² + x_i*x_{i+1} + x_{i+1}²) * cp
            i_yy += (p0.x * p0.x + p0.x * p1.x + p1.x * p1.x) * cp;
        }
    }

    // Finalize values

    // Green's theorem usually gives 2 * Area
    total_area *= 0.5;

    // Handle negative area (clockwise winding)
    // We take absolute value, assuming standard properties
    // However, for holes, signed area is important.
    // If the mesh is properly oriented (CCW outside), outside contours are positive, holes result in negative area.
    // So we should just keep the sign until the end?
    // If we assume a set of disjoint solid bodies, total area should be positive.

    // For mass properties, we usually want positive area.
    // But moments should be divided by the SIGNED area to get correct centroid.
    // Let's assume consistent winding from the mesh.

    // NOTE: TriMesh triangles are usually CCW from outside.
    // The intersection segments are generated (p0->p1) such that water is on the left?
    // In extract_waterline_contours, we set direction:
    // if z0 < draft (p0 below) -> intersection -> (p0 is below, p1 is above).
    // The "water" is the solid part of the hull? No, "waterplane" is the slice of the hull.
    // "Waterplane" usually means the cross section of the hull.
    // So "solid" is "waterplane".
    // If triangle normal points OUT, and we go CCW around the solid.
    // z0 < draft < z1: Edge goes UP through waterplane. Inside is LEFT?
    // Standard rule: Looking from above (Z+), traverse boundary CCW.
    // Interior is on Left.
    // If edge p0->p1 crosses UP (z0 < z_plane < z1):
    // Intersection happens. If we define the segment vector according to Right Hand Rule around normal?
    // Let's verify winding later. For now, we normalize at the end.

    if total_area.abs() <= 1e-9 {
        return None;
    }

    // Centroid formula: Cx = (1/6A) * Σ (xi + xi+1)(xi yi+1 - xi+1 yi)
    // My accumulator `moment_x` stores Σ (xi + xi+1)(cp).
    // So Actual Moment Mx = moment_x / 6.0? NO.
    // Cx = moment_x / (3 * (2*Area)) ??
    // Cx = (moment_x / 6.0) / (total_area) ? No wait.

    // Formula check:
    // A = 0.5 * Σ cp
    // Cx = (1 / (6 * A)) * Σ (x_i + x_{i+1}) * cp
    // Mx (First Moment about Y) = A * Cx = (1/6) * Σ ...

    // So:
    // Actual First Moment_Y (Area * Cx) = moment_x / 6.0. ? No.
    // Let's re-verify variable names.
    // `moment_x` accumulator usually accumulates `My` (moment about Y axis, which gives x-coordinate).
    // My code: `moment_x += (p0.x + p1.x) * cp`. That is related to X-coordinate.
    // So this is `Moment_Y` (integral of x dA).

    // Let A_signed = total_area (which includes 0.5 factor).
    // Cx = (1/6 * sum_x_cp) / A_signed.
    // =>  A_signed * Cx = (1/6 * sum_x_cp) = (sum_x_cp / 6.0).
    // So `moment_x_val` = `moment_x_acc / 6.0`.

    // So:
    let my_val = moment_x / 6.0; // Moment about Y (gives LCF)
    let mx_val = moment_y / 6.0; // Moment about X (gives TCF)

    // If area is negative (CW winding), moments are also negative (consistent).
    // Dividing gives correct positive Centroid.
    let lcf = my_val / total_area;
    let tcf = mx_val / total_area;

    // Second moments
    // I_xx (about origin X axis) = (1/12) * sum_y_sq_cp
    // My accumulator `i_xx` is exactly `sum_y_sq_cp`.
    let i_xx_origin = i_xx / 12.0;
    let i_yy_origin = i_yy / 12.0;

    // We need positive Area and Inertia for results
    let _sign = total_area.signum();
    let area_final = total_area.abs();

    // Correct moments to be about Centroid
    // parallel axis theorem: I_c = I_o - A * d²
    // I_transverse (about LCF, TCF axis parallel to X) = I_xx_origin - Area * TCF²
    // Wait. TCF is y-coordinate. Distance from X-axis is y. Correct.

    // Note: If area was negative, `i_xx_origin` is also negative.
    // `total_area * tcf * tcf`: `total_area` is negative. `tcf` is correct coord involved.
    // So `(-I_o) - (-A) * tcf²` = `- (I_o - A*tcf²)`.
    // So we should compute `I_c_signed` then take abs.

    let i_transverse_signed = i_xx_origin - total_area * tcf * tcf;
    let i_longitudinal_signed = i_yy_origin - total_area * lcf * lcf;

    Some(WaterplaneProperties {
        area: area_final,
        centroid: [lcf, tcf],
        i_transverse: i_transverse_signed.abs(),
        i_longitudinal: i_longitudinal_signed.abs(),
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

    // Collect cut segments
    // A segment is p1->p2
    let mut segments: Vec<(Point3<f64>, Point3<f64>)> = Vec::new();

    for tri_indices in indices {
        let v0 = vertices[tri_indices[0] as usize];
        let v1 = vertices[tri_indices[1] as usize];
        let v2 = vertices[tri_indices[2] as usize];
        let tri_vs = [v0, v1, v2];

        // Check intersections for each edge
        // We want to form a directed segment inside the triangle
        // The tri is v0->v1->v2.
        // If we find intersection points P_a and P_b on edges.
        // We need to order them consistent with CCW winding of the "underwater" (or solid) part?
        // Actually, we just need consistent winding relative to the cut.
        // If we keep the "part below draft" on the left, we get a consistent contour.

        // Find edges crossing the plane
        let mut intersections = Vec::with_capacity(2);

        for i in 0..3 {
            let p_start = tri_vs[i];
            let p_end = tri_vs[(i + 1) % 3];

            let d1 = p_start.z - draft;
            let d2 = p_end.z - draft;

            // Check sign change
            if d1 * d2 < 0.0 {
                // Intersects
                let t = d1 / (d1 - d2); // 0..1
                let pt = Point3::new(
                    p_start.x + (p_end.x - p_start.x) * t,
                    p_start.y + (p_end.y - p_start.y) * t,
                    draft,
                );
                intersections.push((i, pt));
            }
        }

        if intersections.len() == 2 {
            let (idx1, p1) = intersections[0];
            let (idx2, p2) = intersections[1];

            // Determine direction.
            // We want the segment to have "submerged part" on the left?
            // Or just consistent with triangle winding.
            // Triangle normal points OUT.
            // Winding v0->v1->v2 is CCW around normal.
            // If we enter at edge 1 and exit at edge 2?
            // Let's look at the vertex being "isolated".
            // One vertex is on one side, two on the other. Or vice versa.

            // Case 1: 1 vertex below (submerged), 2 above.
            // The submerged region is a triangle corner.
            // Boundary goes: Enter edge -> Corner -> Exit edge -> ...
            // The cut segment closes this corner.
            // To keep "submerged" on Left, we must go from Exit Edge intersection to Enter Edge intersection?
            // Wait.
            // CCW around submerged corner: Corner -> Exit Edge Pt -> Enter Edge Pt -> Corner.
            // So the cut segment is Exit->Enter.

            // Which one is "Enter" vs "Exit"?
            // We follow triangle edges.
            // If edge v0->v1 goes Above->Below (enters water), intersection is "Enter".
            // If edge v1->v2 goes Below->Above (exits water), intersection is "Exit".

            // We want segment Exit->Enter.

            // Let's identify Enter/Exit.
            let get_status = |idx: usize| {
                // Edge starts at idx.
                // d1 = z[idx] - draft.
                // d2 = z[idx+1] - draft.
                // If d1 > 0 and d2 < 0: Above->Below (Enter water).
                // If d1 < 0 and d2 > 0: Below->Above (Exit water).
                let start = tri_vs[idx].z;
                let end = tri_vs[(idx + 1) % 3].z;
                if start > draft && end < draft {
                    1
                }
                // Enter
                else if start < draft && end > draft {
                    -1
                }
                // Exit
                else {
                    0
                }
            };

            let s1 = get_status(idx1);
            let s2 = get_status(idx2);

            if s1 == -1 && s2 == 1 {
                // idx1 is Exit, idx2 is Enter.
                // Segment: p1 -> p2.
                segments.push((p1, p2));
            } else if s1 == 1 && s2 == -1 {
                // idx1 is Enter, idx2 is Exit.
                // Segment: p2 -> p1.
                segments.push((p2, p1));
            }
        }
    }

    // Chain segments
    chain_segments(segments, tolerance)
}

/// Helper to chain directed segments into closed loops
fn chain_segments(
    segments: Vec<(Point3<f64>, Point3<f64>)>,
    tolerance: f64,
) -> Option<Vec<Vec<Point3<f64>>>> {
    if segments.is_empty() {
        return Some(vec![]);
    }

    // We need to match end of one segment to start of another.
    // Use a spatial index or sort.
    // Given the small number of points in a waterplane, simple sort is fast enough.

    #[derive(Clone, Copy)]
    struct Segment {
        p_start: Point3<f64>,
        p_end: Point3<f64>,
        original_idx: usize,
    }

    let mut indexed_segments: Vec<Segment> = segments
        .iter()
        .enumerate()
        .map(|(i, &(s, e))| Segment {
            p_start: s,
            p_end: e,
            original_idx: i,
        })
        .collect();

    // Sort by start point for fast lookup
    // Lexicographical sort on (x, y)
    indexed_segments.sort_by(|a, b| {
        a.p_start
            .x
            .partial_cmp(&b.p_start.x)
            .unwrap()
            .then(a.p_start.y.partial_cmp(&b.p_start.y).unwrap())
    });

    // Helper to find segment starting near p
    let find_next = |p: Point3<f64>, used: &[bool]| -> Option<usize> {
        // Binary search for X
        let found = indexed_segments.binary_search_by(|seg| {
            if seg.p_start.x < p.x - tolerance {
                std::cmp::Ordering::Less
            } else if seg.p_start.x > p.x + tolerance {
                std::cmp::Ordering::Greater
            } else {
                std::cmp::Ordering::Equal
            }
        });

        // Search neighborhood
        let start_idx = match found {
            Ok(i) => i,
            Err(i) => i,
        };

        // Scan backwards
        let mut i = start_idx;
        while i > 0 && indexed_segments[i - 1].p_start.x >= p.x - tolerance {
            i -= 1;
        }

        // Scan fwd
        while i < indexed_segments.len() && indexed_segments[i].p_start.x <= p.x + tolerance {
            let seg = &indexed_segments[i];
            if !used[seg.original_idx] {
                // Check Y and full distance
                if (seg.p_start - p).norm_squared() < tolerance * tolerance {
                    return Some(i);
                }
            }
            i += 1;
        }
        None
    };

    let count = segments.len();
    let mut used = vec![false; count];
    let mut contours = Vec::new();

    // Iterate efficiently. We need to visit all segments.
    // We can iterate through indexed_segments, start a chain if generic segment !used.

    for i in 0..count {
        let root_seg_idx_in_sorted = i;
        // Note: the loop iterates 0..count, but we should check `indexed_segments[i]`

        let root_seg = &indexed_segments[root_seg_idx_in_sorted];
        if used[root_seg.original_idx] {
            continue;
        }

        // Start new contour
        let mut contour = Vec::new();
        let mut current_seg_idx_sorted = root_seg_idx_in_sorted;

        loop {
            let seg = &indexed_segments[current_seg_idx_sorted];
            used[seg.original_idx] = true;
            contour.push(seg.p_start);

            // Find next
            let next_pt = seg.p_end;

            // Try to find a segment starting at next_pt
            if let Some(next_idx) = find_next(next_pt, &used) {
                current_seg_idx_sorted = next_idx;
            } else {
                // Chain broken or closed?
                // If closed, next_pt should be close to contour[0]
                if (next_pt - contour[0]).norm_squared() < tolerance * tolerance {
                    // Closed loop
                    break;
                } else {
                    // Open loop? Waterplane contours should be closed for valid meshes.
                    // But if mesh has holes or is open, we might stop.
                    // For now, assume it ends.
                    contour.push(next_pt);
                    break;
                }
            }
        }

        if contour.len() >= 3 {
            contours.push(contour);
        }
    }

    Some(contours)
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
