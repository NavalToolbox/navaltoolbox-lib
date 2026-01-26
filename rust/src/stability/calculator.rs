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

use super::complete::{CompleteStabilityResult, WindHeelingData};
use super::{StabilityCurve, StabilityPoint};
use crate::hydrostatics::HydrostaticsCalculator;
use crate::mesh::{clip_at_waterline, transform_mesh, transform_point};
use crate::vessel::Vessel;
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
    ///
    /// This function uses parallel processing (Rayon) for improved performance.
    pub fn gz_curve(&self, displacement_mass: f64, cog: [f64; 3], heels: &[f64]) -> StabilityCurve {
        use crate::hull::Hull;
        use crate::vessel::Vessel;
        use rayon::prelude::*;

        // Configuration
        const SIMPLIFICATION_THRESHOLD: usize = 2000;
        const TARGET_TRIANGLES: usize = 1000;

        // Pre-calculate constant geometric properties
        let bounds = self.vessel.get_bounds();
        let center_x = (bounds.0 + bounds.1) / 2.0;
        let center_y = (bounds.2 + bounds.3) / 2.0;
        let z_min = bounds.4;
        let z_max = bounds.5;

        // Calculate total mass and volume
        let ship_mass = displacement_mass;
        let ship_cog = cog;
        let total_fluid_mass: f64 = self.vessel.tanks().iter().map(|t| t.fluid_mass()).sum();
        let total_mass = ship_mass + total_fluid_mass;
        let target_volume = total_mass / self.water_density;

        // Mesh Simplification Logic
        let total_triangles: usize = self.vessel.hulls().iter().map(|h| h.num_triangles()).sum();

        let simplified_vessel_storage = if total_triangles > SIMPLIFICATION_THRESHOLD {
            let simplified_hulls: Vec<Hull> = self
                .vessel
                .hulls()
                .iter()
                .map(|h| h.to_simplified(TARGET_TRIANGLES))
                .collect();
            // Note: new_multi creates a vessel without tanks/openings, which is fine for equilibrium search
            Vessel::new_multi(simplified_hulls).ok()
        } else {
            None
        };

        // Determine which calculator to use for the SEARCH phase.
        // We handle the 'borrow checking' by constructing the proxy calculator inside the scope
        // or using reference.
        // Since we need to use it inside parallel loop, we create a reference wrapper?
        // Actually, we can just decide inside the loop or create a struct.
        // Simplest: pass the relevant `&Vessel` to finding functions, but those are methods on Calculator.
        // So we create a `proxy_calc` here.

        let proxy_calc = if let Some(ref v) = simplified_vessel_storage {
            StabilityCalculator::new(v, self.water_density)
        } else {
            StabilityCalculator::new(self.vessel, self.water_density)
        };

        // Warm start: compute upright equilibrium using PROXY
        let upright_draft = proxy_calc.find_draft_for_volume(
            target_volume,
            0.0,
            0.0,
            center_x,
            center_y,
            z_min,
            z_max,
            None,
        );

        // Parallel processing
        let points: Vec<StabilityPoint> = heels
            .par_iter()
            .map(|&heel| {
                // 1. Calculate Effective COG (using tanks from ORIGINAL vessel)
                let mut total_moment_x = ship_mass * ship_cog[0];
                let mut total_moment_y = ship_mass * ship_cog[1];
                let mut total_moment_z = ship_mass * ship_cog[2];

                if total_fluid_mass > 0.0 {
                    for tank in self.vessel.tanks() {
                        let mass = tank.fluid_mass();
                        if mass > 0.0 {
                            let tank_cog = tank.center_of_gravity_at(heel, 0.0);
                            total_moment_x += mass * tank_cog[0];
                            total_moment_y += mass * tank_cog[1];
                            total_moment_z += mass * tank_cog[2];
                        }
                    }
                }

                let effective_cog = if total_mass > 1e-6 {
                    [
                        total_moment_x / total_mass,
                        total_moment_y / total_mass,
                        total_moment_z / total_mass,
                    ]
                } else {
                    ship_cog
                };

                // 2. Find Equilibrium (using PROXY calculator)
                // Note: We use proxy_calc which might use simplified mesh
                let (draft, trim, _approx_gz) = proxy_calc.find_equilibrium_at_heel(
                    target_volume,
                    effective_cog,
                    heel,
                    0.0,
                    Some(upright_draft),
                    center_x,
                    center_y,
                    z_min,
                    z_max,
                );

                // 3. Calculate Exact GZ (using ORIGINAL/SELF calculator)
                // If proxy was used, the GZ returned is approx. computation on simplified mesh.
                // We MUST recompute GZ on the full mesh for accuracy.
                // If NO proxy was used, proxy_calc IS self, so _approx_gz IS exact.
                // However, re-running one computation is negligible compared to the search.
                // Optimally:
                let gz = if simplified_vessel_storage.is_some() {
                    self.compute_gz_at_state(draft, trim, heel, effective_cog, center_x, center_y)
                } else {
                    _approx_gz // Already computed on full mesh
                };

                // 4. Check Downflooding (using ORIGINAL vessel)
                let pivot = [center_x, center_y, draft];
                let flooded_openings = crate::downflooding::check_openings_submerged(
                    self.vessel.downflooding_openings(),
                    heel,
                    trim,
                    pivot,
                    draft,
                );
                let is_flooding = !flooded_openings.is_empty();

                StabilityPoint {
                    heel,
                    draft,
                    trim,
                    value: gz,
                    is_flooding,
                    flooded_openings,
                }
            })
            .collect();

        StabilityCurve::new_gz(displacement_mass, cog, points)
    }

    /// Computes GZ at a specific state using the current (full) vessel geometry.
    ///
    /// Helper for refinement step after simplified search.
    fn compute_gz_at_state(
        &self,
        draft: f64,
        trim: f64,
        heel: f64,
        cog: [f64; 3],
        center_x: f64, // Bounds center for pivot
        center_y: f64,
    ) -> f64 {
        let pivot = Point3::new(center_x, center_y, draft);
        let mut total_volume = 0.0;
        let mut total_moment = [0.0f64; 3];

        for hull in self.vessel.hulls() {
            let transformed = transform_mesh(hull.mesh(), heel, trim, pivot);
            let (clipped, _) = clip_at_waterline(&transformed, draft);

            if let Some(mesh) = clipped {
                let mass_props = mesh.mass_properties(1.0);
                let vol = mass_props.mass();
                let com = mass_props.local_com;

                total_volume += vol;
                total_moment[0] += vol * com.x;
                total_moment[1] += vol * com.y;
                total_moment[2] += vol * com.z;
            }
        }

        if total_volume <= 0.0 {
            return 0.0;
        }

        // TCB
        let tcb = total_moment[1] / total_volume;

        // Transform Ship COG
        let g_ship = Point3::new(cog[0], cog[1], cog[2]);
        let g_transformed = transform_point(g_ship, heel, trim, pivot);

        // GZ = -(B_y - G_y)
        -(tcb - g_transformed.y)
    }

    /// Calculates KN curves (Righting Lever from Keel) for multiple displacements.
    ///
    /// This is equivalent to calculating GZ curves with VCG = 0.
    /// Returns one curve per displacement, useful for cross-curves of stability.
    ///
    /// # Arguments
    /// * `displacements` - List of target displacements in kg
    /// * `lcg` - Longitudinal Center of Gravity (m)
    /// * `tcg` - Transverse Center of Gravity (m)
    /// * `heels` - List of heel angles in degrees
    pub fn kn_curve(
        &self,
        displacements: &[f64],
        lcg: f64,
        tcg: f64,
        heels: &[f64],
    ) -> Vec<StabilityCurve> {
        // KN is GZ calculated with VCG = 0 (Keel as reference).
        let cog = [lcg, tcg, 0.0];
        displacements
            .iter()
            .map(|&disp| self.gz_curve(disp, cog, heels))
            .collect()
    }

    /// Find draft for target volume at given heel and trim.
    ///
    /// Uses warm start: if initial_draft is provided, uses it as starting point
    /// for faster convergence.
    #[allow(clippy::too_many_arguments)]
    fn find_draft_for_volume(
        &self,
        target_volume: f64,
        heel: f64,
        trim: f64,
        center_x: f64,
        center_y: f64,
        z_min: f64,
        z_max: f64,
        initial_draft: Option<f64>,
    ) -> f64 {
        let tolerance = target_volume * 1e-4;
        let max_iter = 50;

        // Warm start: use initial_draft if provided, otherwise use midpoint
        let (mut low, mut high) = if let Some(init) = initial_draft {
            // Start search around initial draft with a reasonable margin
            let margin = (z_max - z_min) * 0.2;
            ((init - margin).max(z_min), (init + margin).min(z_max))
        } else {
            (z_min, z_max)
        };

        // Initial guess
        let mut mid = if let Some(init) = initial_draft {
            init.clamp(low, high)
        } else {
            (low + high) / 2.0
        };

        for _ in 0..max_iter {
            let pivot = Point3::new(center_x, center_y, mid);

            let mut total_volume = 0.0;
            let mut total_aw = 0.0;

            for hull in self.vessel.hulls() {
                let transformed = transform_mesh(hull.mesh(), heel, trim, pivot);
                if let (Some(mesh), aw) = clip_at_waterline(&transformed, mid) {
                    total_volume += mesh.mass_properties(1.0).mass();
                    total_aw += aw;
                }
            }

            let diff = total_volume - target_volume;

            if diff.abs() < tolerance {
                return mid;
            }

            // Update bounds (Volume is monotonic)
            if diff > 0.0 {
                high = mid;
            } else {
                low = mid;
            }

            // Newton step: z_new = z - diff / Aw
            // Safe Newton: if step falls in bounds, use it. Else bisection.
            if total_aw > 1e-9 {
                let step = diff / total_aw;
                let z_new = mid - step;
                if z_new > low && z_new < high {
                    mid = z_new;
                    continue;
                }
            }

            // Fallback to bisection
            mid = (low + high) / 2.0;
        }

        (low + high) / 2.0
    }

    /// Find equilibrium state at a specific heel angle.
    ///
    /// Optimized: Combines draft search and property calculation in a single pass
    /// to avoid redundant mesh operations.
    #[allow(clippy::too_many_arguments)]
    fn find_equilibrium_at_heel(
        &self,
        target_volume: f64,
        cog: [f64; 3],
        heel: f64,
        initial_trim: f64,
        initial_draft: Option<f64>,
        center_x: f64,
        center_y: f64,
        z_min: f64,
        z_max: f64,
    ) -> (f64, f64, f64) {
        let lcb_tolerance = 0.5;
        let volume_tolerance = target_volume * 1e-4;
        let max_trim_iter = 15;
        let max_draft_iter = 50;

        let mut trim = initial_trim;
        let mut best_draft = initial_draft.unwrap_or((z_min + z_max) / 2.0);
        let mut best_trim = trim;
        let mut best_gz = 0.0;
        let mut best_error = f64::INFINITY;

        // Warm start bounds for draft search
        let mut draft_low = z_min;
        let mut draft_high = z_max;
        if let Some(init) = initial_draft {
            let margin = (z_max - z_min) * 0.2;
            draft_low = (init - margin).max(z_min);
            draft_high = (init + margin).min(z_max);
        }

        for _ in 0..max_trim_iter {
            // Find draft for target volume using bisection
            // Start draft search with current bounds
            let mut low = draft_low;
            let mut high = draft_high;

            // Use best_draft as initial guess if available and within bounds
            let mut mid = best_draft.clamp(low, high);

            let mut final_draft = mid;
            let mut final_volume = 0.0;
            let mut final_moment = [0.0f64; 3];

            for _ in 0..max_draft_iter {
                let pivot = Point3::new(center_x, center_y, mid);

                // Single-pass: compute volume, moment, AND waterplane area for Newton
                let mut total_volume = 0.0;
                let mut total_moment = [0.0f64; 3];
                let mut total_aw = 0.0;

                for hull in self.vessel.hulls() {
                    let transformed = transform_mesh(hull.mesh(), heel, trim, pivot);
                    if let (Some(mesh), aw) = clip_at_waterline(&transformed, mid) {
                        let mass_props = mesh.mass_properties(1.0);
                        let vol = mass_props.mass();
                        let com = mass_props.local_com;

                        total_volume += vol;
                        total_moment[0] += vol * com.x;
                        total_moment[1] += vol * com.y;
                        total_moment[2] += vol * com.z;
                        total_aw += aw;
                    }
                }

                let diff = total_volume - target_volume;

                // Store the latest values
                final_draft = mid;
                final_volume = total_volume;
                final_moment = total_moment;

                if diff.abs() < volume_tolerance {
                    break;
                }

                // Update Bounds
                if diff > 0.0 {
                    high = mid;
                } else {
                    low = mid;
                }

                // Safe Newton Step
                if total_aw > 1e-9 {
                    let step = diff / total_aw;
                    let z_new = mid - step;
                    if z_new > low && z_new < high {
                        mid = z_new;
                        continue;
                    }
                }

                // Fallback Bisection
                mid = (low + high) / 2.0;
            }

            // Use cached COB values from the converged draft (no recomputation!)
            if final_volume <= 0.0 {
                continue;
            }

            let lcb = final_moment[0] / final_volume;
            let tcb = final_moment[1] / final_volume;

            let pivot = Point3::new(center_x, center_y, final_draft);

            // Transform CoG
            let g_ship = Point3::new(cog[0], cog[1], cog[2]);
            let g_transformed = transform_point(g_ship, heel, trim, pivot);

            // GZ = -(B_y - G_y)
            let gz = -(tcb - g_transformed.y);

            // LCB error
            let lcb_error = (lcb - g_transformed.x).abs();

            if lcb_error < best_error {
                best_error = lcb_error;
                best_draft = final_draft;
                best_trim = trim;
                best_gz = gz;

                // Update warm start bounds for next trim iteration
                let margin = (z_max - z_min) * 0.1;
                draft_low = (final_draft - margin).max(z_min);
                draft_high = (final_draft + margin).min(z_max);
            }

            if lcb_error < lcb_tolerance {
                return (final_draft, trim, gz);
            }

            // Adjust trim
            let trim_gain = 0.05;
            trim += (lcb - g_transformed.x) * trim_gain;
            trim = trim.clamp(-10.0, 10.0);
        }

        (best_draft, best_trim, best_gz)
    }

    /// Calculate complete stability analysis for a loading condition.
    ///
    /// Combines hydrostatic calculations, GZ curve, and wind heeling data
    /// (if silhouettes are available) for a single loading condition.
    ///
    /// # Arguments
    /// * `displacement_mass` - Target displacement in kg (ship mass, tanks are added)
    /// * `cog` - Center of gravity (LCG, TCG, VCG) for the ship portion
    /// * `heels` - Heel angles for GZ curve calculation in degrees
    ///
    /// # Returns
    /// A `CompleteStabilityResult` containing:
    /// - Hydrostatic state at equilibrium (GM0, draft, trim, etc.)
    /// - GZ curve for the specified heel angles
    /// - Wind heeling data (if silhouettes exist)
    pub fn complete_stability(
        &self,
        displacement_mass: f64,
        cog: [f64; 3],
        heels: &[f64],
    ) -> CompleteStabilityResult {
        // Calculate hydrostatics at equilibrium
        let hydro_calc = HydrostaticsCalculator::new(self.vessel, self.water_density);
        let hydrostatics = hydro_calc
            .from_displacement(displacement_mass, None, Some(cog), None, None)
            .unwrap_or_default();

        // Calculate GZ curve
        let gz_curve = self.gz_curve(displacement_mass, cog, heels);

        // Calculate wind heeling data if silhouettes exist
        let wind_data = if self.vessel.has_silhouettes() {
            let waterline_z = hydrostatics.draft;
            let emerged_area = self.vessel.get_total_emerged_area(waterline_z);
            let emerged_centroid = self.vessel.get_combined_emerged_centroid(waterline_z);

            if emerged_area > 0.0 {
                Some(WindHeelingData::new(
                    emerged_area,
                    emerged_centroid,
                    waterline_z,
                ))
            } else {
                None
            }
        } else {
            None
        };

        CompleteStabilityResult::new(hydrostatics, gz_curve, wind_data, displacement_mass, cog)
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

        let curve = calc.gz_curve(displacement, cog, &[0.0]);

        // At zero heel for symmetric hull, GZ should be ~0
        assert!(
            curve.points[0].value.abs() < 0.01,
            "GZ at 0 heel = {}",
            curve.points[0].value
        );
    }
    #[test]
    fn test_fsc_gz_reduction() {
        use crate::tanks::Tank;

        let hull = create_box_hull(10.0, 10.0, 10.0);
        let mut vessel = Vessel::new(hull);

        // Add tank with free surface
        // 5x5x2 tank, 50% fill, water density inside
        let tank = Tank::from_box("FSC_Test", 0.0, 5.0, -2.5, 2.5, 0.0, 2.0, 1000.0);
        let mut tank = tank;
        tank.set_fill_percent(50.0);
        vessel.add_tank(tank.clone());

        let calc = StabilityCalculator::new(&vessel, 1025.0);
        let target_total_displacement = 500.0 * 1025.0; // 500m³ * 1.025
        let tank_mass = tank.fluid_mass();
        // Since calculator adds tank mass, we subtract it from input to keep total same
        let ship_mass = target_total_displacement - tank_mass;

        // Ship COG. We want the Total Upright COG to be [0,0,5] for comparison.
        // Total_Moment = Ship_M + Tank_M = Total_Mass * Total_COG
        // Ship_M = Total_M - Tank_M
        // Ship_COG = (Total_M * Target_COG - Tank_M * Tank_COG) / Ship_Mass
        // Tank is at z=0..2 (centered at z=1). COG_tank = [0, 0, 0.5] approx (for 50% full? 0..1m filled -> z=0.5)
        // Let's assume input cog is just the ship cog and we accept the resulting total cog
        // but for the verification logic (GG' reduction), we need to know the effective VCG.
        //
        // SIMPLIFICATION:
        // Let's just run the dry case with the same TOTAL properties (Mass, COG) as the wet case's UPRIGHT state.

        let ship_cog = [0.0, 0.0, 5.0];
        let heel: f64 = 10.0;

        // Calculate GZ with FSC (Wet)
        let curve_wet = calc.gz_curve(ship_mass, ship_cog, &[heel]);
        let gz_wet = curve_wet.points[0].value;

        // Calculate Dry reference
        // We need the exact total mass and exact upright COG of the wet vessel to match
        let total_mass = ship_mass + tank_mass;
        let tank_cog_upright = tank.center_of_gravity();
        let total_cog_z = (ship_mass * ship_cog[2] + tank_mass * tank_cog_upright[2]) / total_mass;
        // X and Y are 0.
        let total_cog = [0.0, 0.0, total_cog_z];

        vessel.remove_tank(0);
        let calc_dry = StabilityCalculator::new(&vessel, 1025.0);
        let curve_dry = calc_dry.gz_curve(total_mass, total_cog, &[heel]);
        let gz_dry = curve_dry.points[0].value;

        // Theoretical Reduction GG'
        // FSM * rho / Total_Mass
        let output_reduction = gz_dry - gz_wet;

        let fsm_inertia = 5.0 * 5.0f64.powi(3) / 12.0;
        let correction_gg = (fsm_inertia * 1000.0) / total_mass;
        let expected_reduction = correction_gg * heel.to_radians().sin();

        assert!(
            (output_reduction - expected_reduction).abs() < 0.02, // slightly looser tolerance for dynamic method
            "FSC Reduction mismatch. Expected: {:.4}, Actual: {:.4}, Dry: {:.4}, Wet: {:.4}",
            expected_reduction,
            output_reduction,
            gz_dry,
            gz_wet
        );
    }
}
