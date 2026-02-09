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

//! Hydrostatics calculator.
//!
//! Calculates hydrostatic properties for vessel geometries.

use super::HydrostaticState;
use crate::mesh::{clip_at_waterline, get_bounds, transform_mesh};
use crate::vessel::Vessel;
use nalgebra::{Point3, Rotation3, Vector3};
use parry3d_f64::shape::Shape;

/// Calculator for hydrostatic properties.
pub struct HydrostaticsCalculator<'a> {
    vessel: &'a Vessel,
    water_density: f64,
}

impl<'a> HydrostaticsCalculator<'a> {
    /// Creates a new hydrostatics calculator.
    ///
    /// # Arguments
    /// * `vessel` - The vessel to calculate hydrostatics for
    /// * `water_density` - Water density in kg/m³ (default: 1025 for seawater)
    pub fn new(vessel: &'a Vessel, water_density: f64) -> Self {
        Self {
            vessel,
            water_density,
        }
    }

    /// Calculates hydrostatics for a fixed draft, trim, and heel.
    ///
    /// # Arguments
    /// * `draft` - Draft at the reference point in meters
    /// * `trim` - Trim angle in degrees (positive = bow down)
    /// * `heel` - Heel angle in degrees (positive = starboard down)
    /// * `vcg` - Optional vertical center of gravity for GM calculation
    /// * `num_stations` - Number of stations for sectional area curve (default: 21)
    pub fn from_draft(
        &self,
        draft: f64,
        trim: f64,
        heel: f64,
        vcg: Option<f64>,
    ) -> Option<HydrostaticState> {
        self.from_draft_with_stations(draft, trim, heel, vcg, None)
    }

    /// Calculates hydrostatics with configurable number of stations for sectional area curve.
    ///
    /// # Arguments
    /// * `draft` - Draft at the reference point in meters
    /// * `trim` - Trim angle in degrees (positive = bow down)
    /// * `heel` - Heel angle in degrees (positive = starboard down)
    /// * `vcg` - Optional vertical center of gravity for GM calculation
    /// * `num_stations` - Number of stations for sectional area curve (default: 21)
    pub fn from_draft_with_stations(
        &self,
        draft: f64,
        trim: f64,
        heel: f64,
        vcg: Option<f64>,
        num_stations: Option<usize>,
    ) -> Option<HydrostaticState> {
        // Use AP/FP (defaults to bounds min/max if not set)
        let ap = self.vessel.ap();
        let fp = self.vessel.fp();
        let mp_x = (ap + fp) / 2.0;

        // Use bounds center for Y
        let bounds = self.vessel.get_bounds();
        let center_y = (bounds.2 + bounds.3) / 2.0;

        // Pivot is now at MP
        let pivot = Point3::new(mp_x, center_y, draft);

        // Calculate specific drafts
        // Trim is positive bow down -> draft increases forward
        let tan_trim = trim.to_radians().tan();
        let draft_mp = draft;
        let draft_fp = draft + (fp - mp_x) * tan_trim;
        let draft_ap = draft + (ap - mp_x) * tan_trim;

        let mut total_volume = 0.0;
        let mut total_moment = [0.0, 0.0, 0.0];

        let mut total_wetted_surface = 0.0;
        let mut total_midship_area = 0.0;

        // We will combine waterplane properties from all hulls
        let mut combined_wp_area = 0.0;
        let mut combined_wp_moment_x = 0.0;
        let mut combined_wp_moment_y = 0.0;
        let mut combined_i_trans = 0.0;
        let mut combined_i_long = 0.0;

        let mut min_x = f64::MAX;
        let mut max_x = f64::MIN;
        let mut min_y = f64::MAX;
        let mut max_y = f64::MIN;

        // Track submerged length (LOS)
        let mut min_x_submerged = f64::MAX;
        let mut max_x_submerged = f64::MIN;

        // Process each hull
        for hull in self.vessel.hulls() {
            // Transform hull
            let transformed = transform_mesh(hull.mesh(), heel, trim, pivot);
            let bounds = get_bounds(&transformed);

            // Clip at waterline
            if let Some(clipped) = clip_at_waterline(&transformed, draft).0 {
                let mass_props = clipped.mass_properties(1.0);
                let vol = mass_props.mass();
                let cob = mass_props.local_com;

                total_volume += vol;
                total_moment[0] += vol * cob.x;
                total_moment[1] += vol * cob.y;
                total_moment[2] += vol * cob.z;

                // Update LOS bounds from clipped mesh vertices
                for v in clipped.vertices() {
                    if v.x < min_x_submerged {
                        min_x_submerged = v.x;
                    }
                    if v.x > max_x_submerged {
                        max_x_submerged = v.x;
                    }
                }

                // Wetted Surface Area: Area(ClippedMesh) - Area(WaterplaneCap)
                let mesh_area = calculate_mesh_area(&clipped);

                // Waterplane Properties
                if let Some(wp) = crate::hydrostatics::waterplane::calculate_waterplane_properties(
                    &transformed,
                    draft,
                ) {
                    total_wetted_surface += (mesh_area - wp.area).max(0.0);

                    combined_wp_area += wp.area;
                    combined_wp_moment_x += wp.area * wp.centroid[0];
                    combined_wp_moment_y += wp.area * wp.centroid[1];

                    // Parallel axis theorem accumulation (relative to origin first)
                    let i_xx_origin = wp.i_transverse + wp.area * wp.centroid[1].powi(2);
                    let i_yy_origin = wp.i_longitudinal + wp.area * wp.centroid[0].powi(2);

                    combined_i_trans += i_xx_origin;
                    combined_i_long += i_yy_origin;

                    min_x = min_x.min(wp.min_x);
                    max_x = max_x.max(wp.max_x);
                    min_y = min_y.min(wp.min_y);
                    max_y = max_y.max(wp.max_y);
                } else {
                    total_wetted_surface += mesh_area;
                }

                // Midship Area: Slice at X = (bounds.0 + bounds.1) / 2.0
                let mid_x = (bounds.0 + bounds.1) / 2.0;
                let ma = calculate_section_area(&clipped, mid_x);
                total_midship_area += ma;
            }
        }

        if total_volume <= 1e-9 {
            return None;
        }

        let lcb = total_moment[0] / total_volume;
        let tcb = total_moment[1] / total_volume;
        let vcb = total_moment[2] / total_volume;
        let cob = [lcb, tcb, vcb];

        let displacement = total_volume * self.water_density;

        // Final Waterplane Properties (Combined)
        let (wp_area, lcf, bmt, bml, lwl, bwl) = if combined_wp_area > 1e-9 {
            let cx = combined_wp_moment_x / combined_wp_area;
            let cy = combined_wp_moment_y / combined_wp_area;

            // Convert inertias back to centroidal
            let i_trans = combined_i_trans - combined_wp_area * cy.powi(2);
            let i_long = combined_i_long - combined_wp_area * cx.powi(2);

            let bmt_val = i_trans / total_volume;
            let bml_val = i_long / total_volume;

            (
                combined_wp_area,
                cx,
                bmt_val,
                bml_val,
                max_x - min_x,
                max_y - min_y,
            )
        } else {
            (0.0, lcb, 0.0, 0.0, 0.0, 0.0)
        };

        // Coefficients
        let cb = if lwl * bwl * draft > 1e-6 {
            total_volume / (lwl * bwl * draft)
        } else {
            0.0
        };
        let cp = if total_midship_area * lwl > 1e-6 {
            total_volume / (total_midship_area * lwl)
        } else {
            0.0
        };
        let cm = if bwl * draft > 1e-6 {
            total_midship_area / (bwl * draft)
        } else {
            0.0
        };

        // Free Surface Correction
        let mut fsm_mass_moment_t = 0.0;
        let mut fsm_mass_moment_l = 0.0;
        for tank in self.vessel.tanks() {
            let rho = tank.fluid_density();
            fsm_mass_moment_t += tank.free_surface_moment_t() * rho;
            fsm_mass_moment_l += tank.free_surface_moment_l() * rho;
        }

        let fsc_t = if displacement > 0.0 {
            fsm_mass_moment_t / displacement
        } else {
            0.0
        };
        let fsc_l = if displacement > 0.0 {
            fsm_mass_moment_l / displacement
        } else {
            0.0
        };

        // GM Calculations
        let (gmt_dry, gml_dry, gmt_wet, gml_wet) = if let Some(vcg_val) = vcg {
            let kb = vcb;
            let kg = vcg_val;

            let gm_t_dry = kb + bmt - kg;
            let gm_l_dry = kb + bml - kg;

            (
                Some(gm_t_dry),
                Some(gm_l_dry),
                Some(gm_t_dry - fsc_t),
                Some(gm_l_dry - fsc_l),
            )
        } else {
            (None, None, None, None)
        };

        // Stiffness Matrix
        let g = 9.81;
        let rho_g = self.water_density * g;
        let mut k = [0.0; 36];

        // Heave (3,3)
        k[14] = rho_g * wp_area;

        // Pitch-Heave Coupling
        let c35 = -rho_g * wp_area * lcf;
        k[16] = c35; // Row 2, Col 4 (3,5)
        k[26] = c35; // Row 4, Col 2 (5,3)

        // Roll (4,4)
        if let Some(gmt) = gmt_wet {
            k[21] = displacement * g * gmt;
        }

        // Pitch (5,5)
        if let Some(gml) = gml_wet {
            k[28] = displacement * g * gml;
        }

        let cog_ret = vcg.map(|z| [lcb, tcb, z]);

        Some(HydrostaticState {
            draft,
            trim,
            heel,
            draft_ap,
            draft_fp,
            draft_mp,
            volume: total_volume,
            displacement,
            cob,
            cog: cog_ret,
            waterplane_area: wp_area,
            lcf,
            bmt,
            bml,
            gmt: gmt_wet,
            gml: gml_wet,
            gmt_dry,
            gml_dry,
            free_surface_correction_t: fsc_t,
            free_surface_correction_l: fsc_l,
            lwl,
            bwl,
            wetted_surface_area: total_wetted_surface,
            midship_area: total_midship_area,
            cm,
            cb,
            cp,
            stiffness_matrix: k,
            los: if max_x_submerged > min_x_submerged {
                max_x_submerged - min_x_submerged
            } else {
                0.0
            },
            sectional_areas: self.calculate_sectional_areas(draft, trim, heel, num_stations),
            freeboard: self.vessel.get_min_freeboard(heel, trim, draft),
        })
    }

    /// Calculates hydrostatics from drafts at Aft and Forward Perpendiculars.
    ///
    /// This is a convenience method that calculates the equivalent mean draft (at MP)
    /// and trim angle, then calls `from_draft`.
    ///
    /// # Arguments
    /// * `draft_ap` - Draft at Aft Perpendicular in meters
    /// * `draft_fp` - Draft at Forward Perpendicular in meters
    /// * `heel` - Heel angle in degrees
    /// * `vcg` - Optional vertical center of gravity for GM calculation
    pub fn from_drafts(
        &self,
        draft_ap: f64,
        draft_fp: f64,
        heel: f64,
        vcg: Option<f64>,
    ) -> Option<HydrostaticState> {
        let ap = self.vessel.ap();
        let fp = self.vessel.fp();
        let lpp = fp - ap;

        if lpp.abs() < 1e-4 {
            // Lpp is too small, assume zero trim
            return self.from_draft(draft_ap, 0.0, heel, vcg);
        }

        // Calculate trim: positive bow down (fp draft > ap draft)
        // tan(trim) = (T_fp - T_ap) / Lpp
        let trim_rad = ((draft_fp - draft_ap) / lpp).atan();
        let trim_deg = trim_rad.to_degrees();

        // Calculate draft at MP (midship)
        // T_mp = T_ap + (MP - AP) * tan(trim)
        // MP = (AP + FP) / 2
        // MP - AP = (FP - AP) / 2 = Lpp / 2
        // T_mp = T_ap + (Lpp/2) * (T_fp - T_ap) / Lpp
        // T_mp = T_ap + (T_fp - T_ap) / 2
        // T_mp = (T_ap + T_fp) / 2
        let draft_mp = (draft_ap + draft_fp) / 2.0;

        self.from_draft(draft_mp, trim_deg, heel, vcg)
    }

    /// Calculate hydrostatics for a given displacement with optional constraints.
    ///
    /// # Arguments
    /// * `displacement_mass` - Target displacement in kg
    /// * `cog` - Optional center of gravity [LCG, TCG, VCG]. Used for GM calculations and equilibrium.
    /// * `trim` - Optional fixed trim angle in degrees
    /// * `heel` - Optional fixed heel angle in degrees
    ///
    /// # Returns
    /// Complete HydrostaticState or error if constraints are invalid/unsatisfiable
    ///
    /// # Constraint Validation
    /// - Cannot specify both trim and LCG (conflicting longitudinal constraints)
    /// - Cannot specify both heel and TCG (conflicting transverse constraints)
    ///
    /// # Valid Constraint Combinations
    /// - Displacement only → finds draft, level trim/heel
    /// - Displacement + VCG only → finds draft, level, computes GMT/GML
    /// - Displacement + VCG + trim → finds draft with fixed trim, free heel
    /// - Displacement + VCG + heel → finds draft with fixed heel, free trim
    /// - Displacement + COG (full) → finds draft, level, full COG specified
    /// - Displacement + trim + heel → finds draft with fixed attitude
    ///
    /// # Arguments
    /// * `displacement_mass` - Target displacement in kg
    /// * `vcg` - Optional VCG only (m) for GM calculation
    /// * `cog` - Optional full COG [LCG, TCG, VCG] (overrides vcg if provided)
    /// * `trim` - Optional trim angle in degrees (default 0.0)
    /// * `heel` - Optional heel angle in degrees (default 0.0)
    pub fn from_displacement(
        &self,
        displacement_mass: f64,
        vcg: Option<f64>,
        cog: Option<[f64; 3]>,
        trim: Option<f64>,
        heel: Option<f64>,
    ) -> Result<HydrostaticState, String> {
        // Validate COG constraints (only if full COG is provided)
        if let Some(cog_val) = cog {
            if trim.is_some() && cog_val[0] != 0.0 {
                return Err(
                    "Cannot specify both trim and LCG: conflicting longitudinal constraints"
                        .to_string(),
                );
            }
            if heel.is_some() && cog_val[1] != 0.0 {
                return Err(
                    "Cannot specify both heel and TCG: conflicting transverse constraints"
                        .to_string(),
                );
            }
        }

        // Initialize state
        let solve_trim = trim.is_none() && cog.is_some();
        let solve_heel = heel.is_none() && cog.is_some();

        // Initial guess or user-provided values
        let mut fixed_trim = trim.unwrap_or(0.0);
        let mut fixed_heel = heel.unwrap_or(0.0);

        // Try Newton-Raphson first (fast)
        if let Some(state) = self.solve_equilibrium_newton(
            displacement_mass,
            vcg,
            cog,
            &mut fixed_trim,
            &mut fixed_heel,
            solve_trim,
            solve_heel,
        )? {
            return Ok(HydrostaticState {
                cog: cog.or(state.cog),
                ..state
            });
        }

        // Fallback to Robust Solver (coordinate descent with bisection)
        // Only if we need to solve for equilibrium
        if solve_trim || solve_heel {
            // Revert to initial guess if Newton diverged wildly
            fixed_trim = trim.unwrap_or(0.0);
            fixed_heel = heel.unwrap_or(0.0);

            if let Some(state) = self.solve_equilibrium_robust(
                displacement_mass,
                vcg,
                cog,
                &mut fixed_trim,
                &mut fixed_heel,
                solve_trim,
                solve_heel,
            ) {
                return Ok(HydrostaticState {
                    cog: cog.or(state.cog),
                    ..state
                });
            }
        }

        Err("Failed to find equilibrium state".to_string())
    }

    /// Primary equilibrium solver using Newton-Raphson (fast but can diverge)
    #[allow(clippy::too_many_arguments)]
    fn solve_equilibrium_newton(
        &self,
        displacement_mass: f64,
        vcg: Option<f64>,
        cog: Option<[f64; 3]>,
        fixed_trim: &mut f64,
        fixed_heel: &mut f64,
        solve_trim: bool,
        solve_heel: bool,
    ) -> Result<Option<HydrostaticState>, String> {
        let tolerance = displacement_mass / self.water_density * 1e-4;
        let target_volume = displacement_mass / self.water_density;

        let bounds = self.vessel.get_bounds();
        let z_min = bounds.4;
        let z_max = bounds.5;

        // Use adaptive iteration count: if we fallback, we don't need too many Newton steps
        let max_equi_iter = 20;
        let mut final_state = None;

        for iter in 0..max_equi_iter {
            // Determine VCG: COG takes precedence over vcg parameter
            let effective_vcg = if let Some(full_cog) = cog {
                Some(full_cog[2])
            } else {
                vcg
            };

            let mut low = z_min;
            let mut high = z_max;
            let mut found_draft = None;

            // Inner loop: Find draft for current heel/trim
            for _ in 0..50 {
                let mid = (low + high) / 2.0;

                if let Some(state) = self.from_draft(mid, *fixed_trim, *fixed_heel, effective_vcg) {
                    let diff = state.volume - target_volume;

                    if diff.abs() < tolerance {
                        found_draft = Some(state);
                        break;
                    }

                    if diff > 0.0 {
                        high = mid;
                    } else {
                        low = mid;
                    }
                } else {
                    low = mid;
                }
            }

            // If we failed to find a valid draft even once, try best estimate
            if found_draft.is_none() {
                let draft = (low + high) / 2.0;
                found_draft = self.from_draft(draft, *fixed_trim, *fixed_heel, effective_vcg);
            }

            match found_draft {
                Some(state) => {
                    // Check for equilibrium convergence
                    if !solve_trim && !solve_heel {
                        return Ok(Some(state));
                    }

                    // COG in ship frame
                    let cog_ship = Vector3::from(cog.unwrap());

                    // Reconstruct pivot used in from_draft
                    let ap = self.vessel.ap();
                    let fp = self.vessel.fp();
                    let mp_x = (ap + fp) / 2.0;
                    let bounds = self.vessel.get_bounds();
                    let center_y = (bounds.2 + bounds.3) / 2.0;
                    let pivot = Point3::new(mp_x, center_y, state.draft);

                    // COB is in global frame, transform to ship frame
                    // p_ship = pivot + Rot_inv * (p_global - pivot)
                    let heel_rad = fixed_heel.to_radians();
                    let trim_rad = fixed_trim.to_radians();
                    let rot_x = Rotation3::from_axis_angle(&Vector3::x_axis(), heel_rad);
                    let rot_y = Rotation3::from_axis_angle(&Vector3::y_axis(), trim_rad);
                    let rotation = rot_y * rot_x;
                    let inv_rotation = rotation.inverse();

                    let cob_global = Point3::new(state.lcb(), state.tcb(), state.vcb());
                    let cob_ship = pivot + inv_rotation * (cob_global - pivot);
                    let cob_ship_vec = Vector3::new(cob_ship.x, cob_ship.y, cob_ship.z);

                    // Difference in ship frame
                    let diff_ship = cog_ship - cob_ship_vec;
                    let mut converged = true;

                    if solve_heel {
                        let gmt = state.gmt.unwrap_or(1.0).max(0.1);
                        let d_heel = -(diff_ship.y / gmt).to_degrees();

                        if d_heel.abs() > 0.001 {
                            // Adaptive damping
                            let damping = if fixed_heel.abs() > 25.0 { 0.3 } else { 0.5 };
                            let step = (d_heel * damping).clamp(-5.0, 5.0);
                            *fixed_heel += step;
                            converged = false;
                        }
                    }

                    if solve_trim {
                        let gml = state.gml.unwrap_or(100.0).max(1.0);
                        let d_trim = (diff_ship.x / gml).to_degrees();

                        if d_trim.abs() > 0.001 {
                            // Adaptive damping
                            let damping = if fixed_heel.abs() > 15.0 { 0.3 } else { 0.5 };
                            let max_step = if fixed_heel.abs() > 15.0 { 1.0 } else { 2.0 };
                            let step = (d_trim * damping).clamp(-max_step, max_step);
                            *fixed_trim += step;
                            converged = false;
                        }
                    }

                    if converged {
                        return Ok(Some(state));
                    }

                    // Keep last state in case we hit max iter but result is "good enough"
                    if iter == max_equi_iter - 1 {
                        final_state = Some(state);
                    }
                }
                None => {
                    return Err(format!(
                        "Could not find draft for displacement {} kg",
                        displacement_mass
                    ));
                }
            }
        }

        // Check if final state is acceptable
        if let Some(state) = final_state {
            let cog_ship = Vector3::from(cog.unwrap());

            // Reconstruct pivot
            let ap = self.vessel.ap();
            let fp = self.vessel.fp();
            let mp_x = (ap + fp) / 2.0;
            let bounds = self.vessel.get_bounds();
            let center_y = (bounds.2 + bounds.3) / 2.0;
            let pivot = Point3::new(mp_x, center_y, state.draft);

            let heel_rad = state.heel.to_radians();
            let trim_rad = state.trim.to_radians();
            let rot = Rotation3::from_axis_angle(&Vector3::y_axis(), trim_rad)
                * Rotation3::from_axis_angle(&Vector3::x_axis(), heel_rad);

            let cob_global = Point3::new(state.lcb(), state.tcb(), state.vcb());
            let cob_ship = pivot + rot.inverse() * (cob_global - pivot);
            let cob_ship_vec = Vector3::new(cob_ship.x, cob_ship.y, cob_ship.z);

            let diff = cog_ship - cob_ship_vec;
            // If error is large (> 0.1m), signal failure to trigger fallback
            if diff.x.abs() > 0.1 || diff.y.abs() > 0.1 {
                return Ok(None);
            }
            return Ok(Some(state));
        }

        Ok(None)
    }

    /// Robust fallback solver using Coordinate Descent with Bisection
    #[allow(clippy::too_many_arguments)]
    fn solve_equilibrium_robust(
        &self,
        displacement_mass: f64,
        _vcg: Option<f64>,
        cog: Option<[f64; 3]>,
        fixed_trim: &mut f64,
        fixed_heel: &mut f64,
        solve_trim: bool,
        solve_heel: bool,
    ) -> Option<HydrostaticState> {
        let max_iter = 10;
        let cog_val = cog.unwrap();
        let target_lcg = cog_val[0];
        let target_tcg = cog_val[1];
        let effective_vcg = Some(cog_val[2]);

        for _ in 0..max_iter {
            let mut improved = false;

            // 1. Solve Heel (fixing trim)
            if solve_heel {
                let (new_heel, _err) = self.bisect_angle(
                    displacement_mass,
                    *fixed_trim,
                    *fixed_heel,
                    target_tcg,
                    true, // is_heel
                    effective_vcg,
                );
                if (new_heel - *fixed_heel).abs() > 0.01 {
                    *fixed_heel = new_heel;
                    improved = true;
                }
            }

            // 2. Solve Trim (fixing heel)
            if solve_trim {
                let (new_trim, _err) = self.bisect_angle(
                    displacement_mass,
                    *fixed_trim,
                    *fixed_heel,
                    target_lcg,
                    false, // is_heel (so it's trim)
                    effective_vcg,
                );
                if (new_trim - *fixed_trim).abs() > 0.01 {
                    *fixed_trim = new_trim;
                    improved = true;
                }
            }

            if !improved {
                break;
            }
        }

        // Return final state
        self.find_draft_for_displacement(displacement_mass, *fixed_trim, *fixed_heel, effective_vcg)
    }

    /// Helper for bisection search
    fn bisect_angle(
        &self,
        disp: f64,
        current_trim: f64,
        current_heel: f64,
        target_val: f64,
        is_heel: bool,
        vcg: Option<f64>,
    ) -> (f64, f64) {
        let (min, max) = if is_heel {
            (-89.0, 89.0)
        } else {
            (-85.0, 85.0)
        };

        // Try to bracket around current value first
        let center = if is_heel { current_heel } else { current_trim };
        let range = 20.0; // Increased range for robustness
        let mut b_min = (center - range).max(min);
        let mut b_max = (center + range).min(max);

        let mut final_diff = 0.0;

        // 20 iterations of bisection
        for _ in 0..20 {
            let mid = (b_min + b_max) / 2.0;
            let (t, h) = if is_heel {
                (current_trim, mid)
            } else {
                (mid, current_heel)
            };

            if let Some(state) = self.find_draft_for_displacement(disp, t, h, vcg) {
                // Transform COB to ship frame to compare with target TCG/LCG
                // Reconstruct pivot (needed for correct transform)
                let ap = self.vessel.ap();
                let fp = self.vessel.fp();
                let mp_x = (ap + fp) / 2.0;
                let bounds = self.vessel.get_bounds();
                let center_y = (bounds.2 + bounds.3) / 2.0;
                let pivot = Point3::new(mp_x, center_y, state.draft);

                let h_rad = h.to_radians();
                let t_rad = t.to_radians();
                let rot = Rotation3::from_axis_angle(&Vector3::y_axis(), t_rad)
                    * Rotation3::from_axis_angle(&Vector3::x_axis(), h_rad);

                let cob_global = Point3::new(state.lcb(), state.tcb(), state.vcb());
                let cob_ship = pivot + rot.inverse() * (cob_global - pivot);

                let val = if is_heel { cob_ship.y } else { cob_ship.x };
                let diff = val - target_val; // B_pos - G_pos
                final_diff = diff;

                if is_heel {
                    // If diff < 0 (B_y < G_y), B is stbd of G. Need to heel port down (negative).
                    // So we need smaller heel (towards min).
                    if diff < 0.0 {
                        b_max = mid;
                    } else {
                        b_min = mid;
                    }
                } else {
                    // Trim: if diff < 0 (B_x < G_x), B is aft of G.
                    // Need B to move fwd. Positive trim (bow down) moves B fwd.
                    // So we need larger trim (towards max).
                    if diff < 0.0 {
                        b_min = mid;
                    } else {
                        b_max = mid;
                    }
                }
            } else {
                // Draft finding failed, narrow range towards 0
                // If mid is positive, it might be too large, try smaller max
                if mid > 0.0 {
                    b_max = mid;
                } else {
                    b_min = mid;
                }
            }
        }

        let best = (b_min + b_max) / 2.0;
        (best, final_diff)
    }

    fn find_draft_for_displacement(
        &self,
        target_disp: f64,
        trim: f64,
        heel: f64,
        vcg: Option<f64>,
    ) -> Option<HydrostaticState> {
        let tolerance = target_disp / self.water_density * 1e-4;
        let target_vol = target_disp / self.water_density;
        let bounds = self.vessel.get_bounds();
        let z_min = bounds.4 - 2.0;
        let z_max = bounds.5 + 2.0;
        let mut low = z_min;
        let mut high = z_max;

        for _ in 0..50 {
            let mid = (low + high) / 2.0;
            if let Some(state) = self.from_draft(mid, trim, heel, vcg) {
                let diff = state.volume - target_vol;
                if diff.abs() < tolerance {
                    return Some(state);
                }
                if diff > 0.0 {
                    high = mid;
                } else {
                    low = mid;
                }
            } else {
                low = mid;
            }
        }
        self.from_draft((low + high) / 2.0, trim, heel, vcg)
    }

    /// Returns the water density.
    pub fn water_density(&self) -> f64 {
        self.water_density
    }

    /// Calculates sectional area curve at evenly spaced stations between AP and FP.
    ///
    /// # Arguments
    /// * `draft` - Draft at the reference point in meters
    /// * `trim` - Trim angle in degrees
    /// * `heel` - Heel angle in degrees
    /// * `num_stations` - Number of stations (default: 21)
    ///
    /// # Returns
    /// Vec of (x_position, area) tuples from AP to FP
    fn calculate_sectional_areas(
        &self,
        draft: f64,
        trim: f64,
        heel: f64,
        num_stations: Option<usize>,
    ) -> Vec<(f64, f64)> {
        let n = num_stations.unwrap_or(21).max(2);
        let ap = self.vessel.ap();
        let fp = self.vessel.fp();
        let lbp = fp - ap;

        if lbp <= 0.0 {
            return Vec::new();
        }

        // Calculate pivot and transform
        let mp_x = (ap + fp) / 2.0;
        let bounds = self.vessel.get_bounds();
        let center_y = (bounds.2 + bounds.3) / 2.0;
        let pivot = Point3::new(mp_x, center_y, draft);

        // Combine clipped meshes from all hulls
        let mut clipped_meshes = Vec::new();
        for hull in self.vessel.hulls() {
            let transformed = transform_mesh(hull.mesh(), heel, trim, pivot);
            if let Some(clipped) = clip_at_waterline(&transformed, draft).0 {
                clipped_meshes.push(clipped);
            }
        }

        // Calculate area at each station
        let mut areas = Vec::with_capacity(n);
        for i in 0..n {
            let x = ap + (i as f64) * lbp / ((n - 1) as f64);
            let mut total_area = 0.0;
            for mesh in &clipped_meshes {
                total_area += calculate_section_area(mesh, x);
            }
            areas.push((x, total_area));
        }

        areas
    }
}

/// Calculate total surface area of a mesh
fn calculate_mesh_area(mesh: &parry3d_f64::shape::TriMesh) -> f64 {
    let vertices = mesh.vertices();
    let indices = mesh.indices();
    let mut area = 0.0;

    for tri in indices {
        let v0 = vertices[tri[0] as usize];
        let v1 = vertices[tri[1] as usize];
        let v2 = vertices[tri[2] as usize];

        let ab = v1 - v0;
        let ac = v2 - v0;
        let cross = ab.cross(&ac);
        area += 0.5 * cross.norm();
    }
    area
}

/// Calculate the cross-sectional area of a mesh at a given X plane
fn calculate_section_area(mesh: &parry3d_f64::shape::TriMesh, x_plane: f64) -> f64 {
    let vertices = mesh.vertices();
    let indices = mesh.indices();

    // Find intersection segments with X plane
    let mut segments: Vec<(Point3<f64>, Point3<f64>)> = Vec::new();
    let tolerance = 1e-6;

    for tri in indices {
        let v0 = vertices[tri[0] as usize];
        let v1 = vertices[tri[1] as usize];
        let v2 = vertices[tri[2] as usize];

        // Calculate signed distances and signs
        let dists: [f64; 3] = [v0.x - x_plane, v1.x - x_plane, v2.x - x_plane];
        let signs: [i32; 3] = dists.map(|d| {
            if d.abs() < tolerance {
                0
            } else if d < 0.0 {
                -1
            } else {
                1
            }
        });

        // Helper to interpolate
        let interp = |i: usize, j: usize| -> Point3<f64> {
            let p_i = vertices[tri[i] as usize];
            let p_j = vertices[tri[j] as usize];
            let d_i = dists[i];
            let d_j = dists[j];
            let t = d_i / (d_i - d_j);
            Point3::new(
                x_plane,
                p_i.y + (p_j.y - p_i.y) * t,
                p_i.z + (p_j.z - p_i.z) * t,
            )
        };

        // Match cases based on signs (s0, s1, s2)
        match signs {
            // No intersection (all same side strict)
            [1, 1, 1] | [-1, -1, -1] => {}

            // Single vertex on plane, others same side (touching)
            [0, 1, 1] | [0, -1, -1] => {} // Point contact v0
            [1, 0, 1] | [-1, 0, -1] => {} // Point contact v1
            [1, 1, 0] | [-1, -1, 0] => {} // Point contact v2

            // Edge on plane, other same side (edge contact)
            [0, 0, 1] | [0, 0, -1] => {
                segments.push((vertices[tri[0] as usize], vertices[tri[1] as usize]));
            }
            [1, 0, 0] | [-1, 0, 0] => {
                segments.push((vertices[tri[1] as usize], vertices[tri[2] as usize]));
            }
            [0, 1, 0] | [0, -1, 0] => {
                segments.push((vertices[tri[2] as usize], vertices[tri[0] as usize]));
            }

            // Face on plane (degenerate)
            [0, 0, 0] => {
                let p0 = vertices[tri[0] as usize];
                let p1 = vertices[tri[1] as usize];
                let p2 = vertices[tri[2] as usize];
                segments.push((p0, p1));
                segments.push((p1, p2));
                segments.push((p2, p0));
            }

            // Standard Crossing (one distinct side) -> 2 intersections
            // v2 is alone
            [s0, s1, s2] if s0 == s1 && s2 != s0 => {
                let p_a = interp(1, 2);
                let p_b = interp(2, 0);
                segments.push((p_a, p_b));
            }
            // v0 is alone
            [s0, s1, s2] if s1 == s2 && s0 != s1 => {
                let p_a = interp(2, 0);
                let p_b = interp(0, 1);
                segments.push((p_a, p_b));
            }
            // v1 is alone
            [s0, s1, s2] if s2 == s0 && s1 != s2 => {
                let p_a = interp(0, 1);
                let p_b = interp(1, 2);
                segments.push((p_a, p_b));
            }

            // Single vertex on plane, others splitted
            [0, -1, 1] | [0, 1, -1] => {
                let p0 = vertices[tri[0] as usize];
                let p_cross = interp(1, 2);
                segments.push((p0, p_cross));
            }
            [-1, 0, 1] | [1, 0, -1] => {
                let p1 = vertices[tri[1] as usize];
                let p_cross = interp(2, 0);
                segments.push((p1, p_cross));
            }
            [-1, 1, 0] | [1, -1, 0] => {
                let p2 = vertices[tri[2] as usize];
                let p_cross = interp(0, 1);
                segments.push((p2, p_cross));
            }

            // impossible cases
            _ => {}
        }
    }

    // println!("Total segments found: {}", segments.len());

    if segments.is_empty() {
        // println!("No segments found at X={}", x_plane);
        return 0.0;
    }

    // Chain segments into contours
    // Naive O(N^2) chaining
    let mut contours: Vec<Vec<Point3<f64>>> = Vec::new();

    while let Some((start, mut current)) = segments.pop() {
        // Start a new contour with the last segment (pop for efficiency)
        let mut contour = vec![start, current];

        let mut loop_closed = false;

        while !loop_closed {
            // Find a segment starting at 'current'
            let mut found_idx = None;

            for (i, (s, e)) in segments.iter().enumerate() {
                if (s - current).norm_squared() < tolerance {
                    found_idx = Some((i, *e, false)); // Standard direction
                    break;
                } else if (e - current).norm_squared() < tolerance {
                    found_idx = Some((i, *s, true)); // Reversed direction
                    break;
                }
            }

            if let Some((idx, next_pt, _reversed)) = found_idx {
                segments.swap_remove(idx); // Remove found segment
                current = next_pt;
                contour.push(current);

                // Check if loop closed
                if (current - start).norm_squared() < tolerance {
                    loop_closed = true;
                }
            } else {
                // Broken loop or open chain
                break;
            }
        }

        if loop_closed {
            contours.push(contour);
        }
    }

    // Calculate area of contours (YZ plane)
    let mut total_area = 0.0;
    for contour in contours {
        let mut loop_area = 0.0;
        let n = contour.len();
        if n < 3 {
            continue;
        }

        for i in 0..n {
            let p1 = contour[i];
            let p2 = contour[(i + 1) % n];
            // Shoelace formula: (y1 + y2)(z1 - z2) / 2
            loop_area += (p1.y + p2.y) * (p1.z - p2.z);
        }
        total_area += 0.5 * loop_area.abs(); // Assume disjoint loops add up
    }

    total_area
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
    fn test_box_barge_volume() {
        // 10m x 10m x 10m box
        let hull = create_box_hull(10.0, 10.0, 10.0);
        let vessel = Vessel::new(hull);
        let calc = HydrostaticsCalculator::new(&vessel, 1025.0);

        // At draft 5m, volume should be 10 * 10 * 5 = 500 m³
        let state = calc.from_draft(5.0, 0.0, 0.0, None).unwrap();
        assert!(
            (state.volume - 500.0).abs() < 1.0,
            "Volume was {}",
            state.volume
        );
    }

    #[test]
    fn test_from_displacement_level() {
        let hull = create_box_hull(10.0, 10.0, 10.0);
        let vessel = Vessel::new(hull);
        let calc = HydrostaticsCalculator::new(&vessel, 1025.0);

        // Target displacement: 500 m³ * 1025 kg/m³ = 512500 kg
        let target_disp = 500.0 * 1025.0;

        // Calculate at displacement with no other constraints (level keel)
        let state = calc
            .from_displacement(target_disp, None, None, None, None)
            .expect("Calculation failed");

        assert!(
            (state.draft - 5.0).abs() < 0.01,
            "Draft should be ~5.0m, got {}",
            state.draft
        );
        assert!(
            (state.displacement - target_disp).abs() < 1.0,
            "Displacement mismatch"
        );
        assert_eq!(state.trim, 0.0);
        assert_eq!(state.heel, 0.0);
    }

    #[test]
    fn test_calculate_section_area() {
        // Create a 10x10x10 box
        let hull = create_box_hull(10.0, 10.0, 10.0);
        let mesh = hull.mesh();

        // Full slice at X=5.0
        let area = calculate_section_area(mesh, 5.0);
        assert!(
            (area - 100.0).abs() < 1.0,
            "Full area should be 100.0, got {}",
            area
        );

        // Clipped slice at X=5.0, Draft=5.0
        // Expect 10 (width) * 5 (draft) = 50 m²
        if let Some(clipped) = crate::mesh::clip_at_waterline(mesh, 5.0).0 {
            println!("Clipped Vertices: {:?}", clipped.vertices());
            println!("Clipped Indices: {:?}", clipped.indices());

            let area_clipped = calculate_section_area(&clipped, 5.0);
            assert!(
                (area_clipped - 50.0).abs() < 1.0,
                "Clipped area should be 50.0, got {}",
                area_clipped
            );
        } else {
            panic!("Clipping failed");
        }
    }

    #[test]
    fn test_from_displacement_with_vcg() {
        let hull = create_box_hull(10.0, 10.0, 10.0);
        let vessel = Vessel::new(hull);
        let calc = HydrostaticsCalculator::new(&vessel, 1025.0);
        let target_disp = 512500.0; // 5m draft condition

        // With VCG provided, should compute GMT/GML
        // Note: LCB/TCB assumed 0.0 for box hull, so just set VCG=7.0
        let state = calc
            .from_displacement(target_disp, Some(7.0), None, None, None)
            .expect("Calculation failed");

        assert!((state.draft - 5.0).abs() < 0.01);

        // Check that GMT is computed (VCG was provided)
        // For vcg-only mode, cog should now be populated with [LCB, TCB, VCG]
        // because from_draft populates it when vcg is passed, and we now preserve it.
        assert!(state.cog.is_some(), "COG should be Some for vcg-only mode");
        let cog = state.cog.unwrap();
        assert!((cog[2] - 7.0).abs() < 1e-6, "VCG should matches input");
        // For a box hull, LCB=5.0, TCB=0.0
        assert!((cog[0] - state.lcb()).abs() < 1e-6, "LCG should match LCB");
        assert!((cog[1] - state.tcb()).abs() < 1e-6, "TCG should match TCB");

        // Check Stability calculation
        // BM_t = 10²/60 = 1.667
        // VCB = 2.5
        // KM_t = 4.167
        // GMT_dry = 4.167 - 7.0 = -2.833
        assert!(state.gmt.is_some());
        assert!((state.gmt_dry.unwrap() - -2.833).abs() < 0.1);
    }

    #[test]
    fn test_constraints_validation() {
        let hull = create_box_hull(10.0, 10.0, 10.0);
        let vessel = Vessel::new(hull);
        let calc = HydrostaticsCalculator::new(&vessel, 1025.0);

        // Invalid: Trim provided but also LCG constrained (non-zero)
        let res = calc.from_displacement(100000.0, None, Some([5.0, 0.0, 0.0]), Some(0.0), None);
        assert!(res.is_err(), "Should fail for both LCG and Trim specified");

        // Invalid: Heel provided but also TCG constrained
        let res = calc.from_displacement(100000.0, None, Some([0.0, 5.0, 0.0]), None, Some(0.0));
        assert!(res.is_err(), "Should fail for both TCG and Heel specified");
    }
    #[test]
    fn test_from_drafts() {
        let hull = create_box_hull(100.0, 20.0, 10.0);
        let mut vessel = Vessel::new(hull);
        // Set AP/FP explicitly
        vessel.set_ap(0.0);
        vessel.set_fp(100.0);

        let calc = HydrostaticsCalculator::new(&vessel, 1025.0);

        // Case 1: Even Keel (Draft=5.0)
        let state1 = calc.from_drafts(5.0, 5.0, 0.0, None).unwrap();
        assert!((state1.draft - 5.0).abs() < 1e-6);
        assert!((state1.draft_ap - 5.0).abs() < 1e-6);
        assert!((state1.draft_fp - 5.0).abs() < 1e-6);
        assert!(state1.trim.abs() < 1e-6);

        // Case 2: Trimmed by stern (AP=6.0, FP=4.0)
        // MP draft should be 5.0
        // Trim = atan((4-6)/100) = atan(-0.02)
        let state2 = calc.from_drafts(6.0, 4.0, 0.0, None).unwrap();
        assert!((state2.draft_mp - 5.0).abs() < 1e-6);
        assert!((state2.draft_ap - 6.0).abs() < 1e-6);
        assert!((state2.draft_fp - 4.0).abs() < 1e-6);
        assert!(state2.trim < 0.0); // Stern down is negative trim
    }

    #[test]
    fn test_vcg_handling_scenarios() {
        let hull = create_box_hull(10.0, 10.0, 10.0);
        let vessel = Vessel::new(hull);
        let calc = HydrostaticsCalculator::new(&vessel, 1025.0);
        let target_disp = 512500.0; // 5m draft condition

        // Scenario 1: No VCG, No COG -> COG should be None
        let state1 = calc
            .from_displacement(target_disp, None, None, None, None)
            .expect("Calc failed S1");
        assert!(state1.cog.is_none(), "S1: COG should be None");

        // Scenario 2: VCG only -> COG should be [LCB, TCB, VCG]
        let state2 = calc
            .from_displacement(target_disp, Some(6.0), None, None, None)
            .expect("Calc failed S2");
        assert!(state2.cog.is_some(), "S2: COG should be Some");
        let cog2 = state2.cog.unwrap();
        assert!((cog2[2] - 6.0).abs() < 1e-6, "S2: VCG mismatch");
        assert!(
            (cog2[0] - state2.lcb()).abs() < 1e-6,
            "S2: LCG should match LCB"
        );
        assert!(
            (cog2[1] - state2.tcb()).abs() < 1e-6,
            "S2: TCG should match TCB"
        );

        // Scenario 3: Full COG provided -> COG should match input exactly (overriding VCG arg if any)
        // Note: from_displacement prefers 'cog' arg over 'vcg' arg for the computation
        let input_cog = [2.0, 1.0, 8.0];
        let state3 = calc
            .from_displacement(target_disp, Some(5.0), Some(input_cog), None, None)
            .expect("Calc failed S3");
        assert!(state3.cog.is_some(), "S3: COG should be Some");
        let cog3 = state3.cog.unwrap();
        assert_eq!(cog3, input_cog, "S3: COG should match input full COG");
        // Ensure it didn't use the '5.0' from vcg arg
        assert!((cog3[2] - 8.0).abs() < 1e-6);
    }

    #[test]
    fn test_equilibrium_heel_from_tcg_offset() {
        // Box hull 10x10x10, centered at y=0
        // Right-hand system: Y+ = port, positive heel = port UP (starboard down)
        // Weight on port (TcG > 0) → port goes DOWN → NEGATIVE heel
        let hull = create_box_hull(10.0, 10.0, 10.0);
        let vessel = Vessel::new(hull);
        let calc = HydrostaticsCalculator::new(&vessel, 1025.0);
        let target_disp = 512500.0; // 5m draft condition

        // COG offset to port (TcG = +1.0m)
        let cog = [5.0, 1.0, 5.0]; // LcG, TcG, VcG

        let state = calc
            .from_displacement(target_disp, None, Some(cog), None, None)
            .expect("Calculation failed");

        // Port weight → port down → negative heel (right-hand convention)
        assert!(
            state.heel < -0.1,
            "Heel should be negative for port TcG, got {}",
            state.heel
        );
    }

    #[test]
    fn test_equilibrium_trim_from_lcg_offset() {
        // Box hull 10x10x10, LCB at x=5.0
        // LcG < LCB (aft) → expect negative trim (stern down)
        // Note: Trim is rotation around Y. Positive trim = X+ (fwd) goes down (bow down).
        // If we put weight aft (small X), stern goes down, bow goes up. So negative trim.
        let hull = create_box_hull(10.0, 10.0, 10.0);
        let vessel = Vessel::new(hull);
        let calc = HydrostaticsCalculator::new(&vessel, 1025.0);
        let target_disp = 512500.0; // 5m draft

        // COG offset to aft (LcG = 3.0m, LCB = 5.0m)
        let cog = [3.0, 0.0, 5.0]; // LcG, TcG, VcG

        let state = calc
            .from_displacement(target_disp, None, Some(cog), None, None)
            .expect("Calculation failed");

        // Trim positive = bow down. Stern down = negative trim.
        assert!(
            state.trim < -0.1,
            "Trim should be negative for aft LcG, got {}",
            state.trim
        );
    }

    #[test]
    fn test_equilibrium_vertical_alignment() {
        // This test verifies the ACTUAL equilibrium condition:
        // COG and COB must be vertically aligned in the global (waterplane) frame.
        // This catches frame mixing bugs that simpler sign-only tests miss.

        let hull = create_box_hull(10.0, 10.0, 10.0);
        let vessel = Vessel::new(hull);
        let calc = HydrostaticsCalculator::new(&vessel, 1025.0);
        let target_disp = 512500.0; // 5m draft

        // Test combined offset: weight on port AND aft
        let cog = [3.0, 1.0, 5.0]; // LcG=3 (aft of LCB=5), TcG=1 (port), VcG=5

        let state = calc
            .from_displacement(target_disp, None, Some(cog), None, None)
            .expect("Calculation failed");

        // Now verify equilibrium: COG and COB should be aligned in SHIP frame
        // The solver transforms COB to ship frame using a pivot at (MP, 0, Draft)
        let heel_rad = state.heel.to_radians();
        let trim_rad = state.trim.to_radians();

        let rot_x = nalgebra::Rotation3::from_axis_angle(&nalgebra::Vector3::x_axis(), heel_rad);
        let rot_y = nalgebra::Rotation3::from_axis_angle(&nalgebra::Vector3::y_axis(), trim_rad);
        let rotation = rot_y * rot_x;
        let inv_rotation = rotation.inverse();

        // Reconstruct Pivot
        let ap = vessel.ap();
        let fp = vessel.fp();
        let mp_x = (ap + fp) / 2.0;
        let bounds = vessel.get_bounds();
        let center_y = (bounds.2 + bounds.3) / 2.0;
        let pivot = nalgebra::Point3::new(mp_x, center_y, state.draft);

        let cog_ship = nalgebra::Vector3::new(cog[0], cog[1], cog[2]);
        let cob_global = nalgebra::Point3::new(state.lcb(), state.tcb(), state.vcb());

        // p_ship = pivot + Rot_inv * (p_global - pivot)
        let cob_ship_point = pivot + inv_rotation * (cob_global - pivot);
        let cob_ship = nalgebra::Vector3::new(cob_ship_point.x, cob_ship_point.y, cob_ship_point.z);

        // Also compute for debug: COG in global
        // p_global = pivot + Rot * (p_ship - pivot)
        // Not strictly needed for the test assertion but useful debug
        // let cog_ship_point = pivot + (cog_ship - nalgebra::Vector3::new(pivot.x, pivot.y, pivot.z));
        // let _cog_global = pivot + rotation * (cog_ship_point - pivot);

        // At equilibrium, horizontal components in ship frame should match closely
        let diff_x = (cog_ship.x - cob_ship.x).abs();
        let diff_y = (cog_ship.y - cob_ship.y).abs();

        assert!(
            diff_x < 0.1,
            "Equilibrium failed: COG_x ({:.4}) != COB_x ({:.4}) in ship frame, diff={:.4}",
            cog_ship.x,
            cob_ship.x,
            diff_x
        );
        assert!(
            diff_y < 0.1,
            "Equilibrium failed: COG_y ({:.4}) != COB_y ({:.4}) in ship frame, diff={:.4}",
            cog_ship.y,
            cob_ship.y,
            diff_y
        );
    }
}
