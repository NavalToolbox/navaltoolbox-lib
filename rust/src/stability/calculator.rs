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

//! Stability calculator.
//!
//! Calculates KN and GZ curves.

use crate::vessel::Vessel;
use crate::hydrostatics::HydrostaticsCalculator;
use crate::mesh::{clip_at_waterline, transform_mesh, transform_point, get_bounds};
use super::{StabilityPoint, StabilityCurve};
use nalgebra::Point3;
use parry3d_f64::shape::Shape;

/// Calculator for stability analysis (KN and GZ curves).
pub struct StabilityCalculator<'a> {
    vessel: &'a Vessel,
    water_density: f64,
}

impl<'a> StabilityCalculator<'a> {
    /// Creates a new stability calculator.
    pub fn new(vessel: &'a Vessel, water_density: f64) -> Self {
        Self {
            vessel,
            water_density,
        }
    }

    /// Calculates the GZ curve for a specific loading condition.
    ///
    /// # Arguments
    /// * `displacement_mass` - Target displacement in kg
    /// * `cog` - Center of gravity (LCG, TCG, VCG)
    /// * `heels` - List of heel angles in degrees
    pub fn calculate_gz_curve(
        &self,
        displacement_mass: f64,
        cog: [f64; 3],
        heels: &[f64],
    ) -> StabilityCurve {
        let target_volume = displacement_mass / self.water_density;
        let bounds = self.vessel.get_bounds();
        let center_x = (bounds.0 + bounds.1) / 2.0;
        let center_y = (bounds.2 + bounds.3) / 2.0;
        let z_min = bounds.4;
        let z_max = bounds.5;

        // Find upright equilibrium draft
        let upright_draft = self.find_draft_for_volume(target_volume, 0.0, 0.0, z_min, z_max);
        let _ = upright_draft; // Used for reference

        let mut points = Vec::with_capacity(heels.len());
        let mut prev_trim = 0.0;

        for &heel in heels {
            // Find equilibrium at this heel
            let (draft, trim, gz) = self.find_equilibrium_at_heel(
                target_volume,
                cog,
                heel,
                prev_trim,
                center_x,
                center_y,
                z_min,
                z_max,
            );

            prev_trim = trim;

            // Check downflooding openings
            let pivot = [center_x, center_y, draft];
            let flooded_openings = crate::downflooding::check_openings_submerged(
                self.vessel.downflooding_openings(),
                heel,
                trim,
                pivot,
                draft,
            );
            let is_flooding = !flooded_openings.is_empty();

            points.push(StabilityPoint {
                heel,
                draft,
                trim,
                value: gz,
                is_flooding,
                flooded_openings,
            });
        }

        StabilityCurve::new_gz(displacement_mass, cog, points)
    }

    /// Find draft for target volume at given heel and trim.
    fn find_draft_for_volume(
        &self,
        target_volume: f64,
        heel: f64,
        trim: f64,
        z_min: f64,
        z_max: f64,
    ) -> f64 {
        let bounds = self.vessel.get_bounds();
        let center_x = (bounds.0 + bounds.1) / 2.0;
        let center_y = (bounds.2 + bounds.3) / 2.0;

        let tolerance = target_volume * 1e-4;
        let max_iter = 50;

        let mut low = z_min;
        let mut high = z_max;

        for _ in 0..max_iter {
            let mid = (low + high) / 2.0;
            let pivot = Point3::new(center_x, center_y, mid);

            let mut total_volume = 0.0;
            for hull in self.vessel.hulls() {
                let transformed = transform_mesh(hull.mesh(), heel, trim, pivot);
                if let Some(clipped) = clip_at_waterline(&transformed, mid) {
                    total_volume += clipped.mass_properties(1.0).mass();
                }
            }

            let diff = total_volume - target_volume;

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

    /// Find equilibrium state at a specific heel angle.
    fn find_equilibrium_at_heel(
        &self,
        target_volume: f64,
        cog: [f64; 3],
        heel: f64,
        initial_trim: f64,
        center_x: f64,
        center_y: f64,
        z_min: f64,
        z_max: f64,
    ) -> (f64, f64, f64) {
        let lcg = cog[0];
        let vcg = cog[2];
        let lcb_tolerance = 0.5;
        let max_iter = 15;

        let mut trim = initial_trim;
        let mut best_draft = (z_min + z_max) / 2.0;
        let mut best_trim = trim;
        let mut best_gz = 0.0;
        let mut best_error = f64::INFINITY;

        for _ in 0..max_iter {
            let draft = self.find_draft_for_volume(target_volume, heel, trim, z_min, z_max);
            let pivot = Point3::new(center_x, center_y, draft);

            let mut total_volume = 0.0;
            let mut total_moment = [0.0, 0.0, 0.0];

            for hull in self.vessel.hulls() {
                let transformed = transform_mesh(hull.mesh(), heel, trim, pivot);
                if let Some(clipped) = clip_at_waterline(&transformed, draft) {
                    let mass_props = clipped.mass_properties(1.0);
                    let vol = mass_props.mass();
                    let cob = mass_props.local_com;

                    total_volume += vol;
                    total_moment[0] += vol * cob.x;
                    total_moment[1] += vol * cob.y;
                    total_moment[2] += vol * cob.z;
                }
            }

            if total_volume <= 0.0 {
                continue;
            }

            let lcb = total_moment[0] / total_volume;
            let tcb = total_moment[1] / total_volume;

            // Transform CoG
            let g_ship = Point3::new(cog[0], cog[1], cog[2]);
            let g_transformed = transform_point(g_ship, heel, trim, pivot);

            // GZ = -(B_y - G_y)
            let gz = -(tcb - g_transformed.y);

            // LCB error
            let lcb_error = (lcb - g_transformed.x).abs();

            if lcb_error < best_error {
                best_error = lcb_error;
                best_draft = draft;
                best_trim = trim;
                best_gz = gz;
            }

            if lcb_error < lcb_tolerance {
                return (draft, trim, gz);
            }

            // Adjust trim
            let trim_gain = 0.05;
            trim += (lcb - g_transformed.x) * trim_gain;
            trim = trim.clamp(-10.0, 10.0);
        }

        (best_draft, best_trim, best_gz)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hull::Hull;
    use nalgebra::Point3;
    use parry3d_f64::shape::TriMesh;

    fn create_box_hull(loa: f64, boa: f64, depth: f64) -> Hull {
        let hb = boa / 2.0;
        let vertices = vec![
            Point3::new(0.0, -hb, 0.0),
            Point3::new(loa, -hb, 0.0),
            Point3::new(loa, hb, 0.0),
            Point3::new(0.0, hb, 0.0),
            Point3::new(0.0, -hb, depth),
            Point3::new(loa, -hb, depth),
            Point3::new(loa, hb, depth),
            Point3::new(0.0, hb, depth),
        ];
        let indices = vec![
            [0, 2, 1], [0, 3, 2],
            [4, 5, 6], [4, 6, 7],
            [0, 1, 5], [0, 5, 4],
            [2, 3, 7], [2, 7, 6],
            [0, 4, 7], [0, 7, 3],
            [1, 2, 6], [1, 6, 5],
        ];
        let mesh = TriMesh::new(vertices, indices).unwrap();
        Hull::from_mesh(mesh)
    }

    #[test]
    fn test_gz_at_zero_heel() {
        let hull = create_box_hull(10.0, 10.0, 10.0);
        let vessel = Vessel::new(hull);
        let calc = StabilityCalculator::new(&vessel, 1025.0);

        let cog = [5.0, 0.0, 2.0]; // Center of box, low VCG
        let displacement = 500.0 * 1025.0; // 500 m³ at 5m draft

        let curve = calc.calculate_gz_curve(displacement, cog, &[0.0]);

        // At zero heel for symmetric hull, GZ should be ~0
        assert!(curve.points[0].value.abs() < 0.01, "GZ at 0 heel = {}", curve.points[0].value);
    }
}
