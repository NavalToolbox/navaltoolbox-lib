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
use nalgebra::Point3;
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
    pub fn calculate_at_draft(
        &self,
        draft: f64,
        trim: f64,
        heel: f64,
        vcg: Option<f64>,
    ) -> Option<HydrostaticState> {
        let bounds = self.vessel.get_bounds();
        let center_x = (bounds.0 + bounds.1) / 2.0;
        let center_y = (bounds.2 + bounds.3) / 2.0;
        let pivot = Point3::new(center_x, center_y, draft);

        let mut total_volume = 0.0;
        let mut total_moment = [0.0, 0.0, 0.0];

        // Process each hull
        for hull in self.vessel.hulls() {
            // Transform hull
            let transformed = transform_mesh(hull.mesh(), heel, trim, pivot);
            let _trans_bounds = get_bounds(&transformed);

            // Clip at waterline
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
            return None;
        }

        let lcb = total_moment[0] / total_volume;
        let tcb = total_moment[1] / total_volume;
        let vcb = total_moment[2] / total_volume;
        
        let cob = [lcb, tcb, vcb];

        let displacement = total_volume * self.water_density;

        // Calculate waterplane properties for metacentric calculations
        // We need to calculate on the transformed mesh to account for heel/trim
        let (waterplane_area, lcf, bmt, bml) = if total_volume > 0.0 {
            // For each hull, calculate its waterplane contribution
            let mut total_wp_area = 0.0;
            let mut total_i_t = 0.0;
            let mut total_i_l = 0.0;
            let mut wp_moment_x = 0.0;
            
            for hull in self.vessel.hulls() {
                // Transform hull for heel/trim
                let transformed = transform_mesh(hull.mesh(), heel, trim, pivot);
                
                // Calculate waterplane properties at this draft
                if let Some(wp) = super::waterplane::calculate_waterplane_properties(&transformed, draft) {
                    total_wp_area += wp.area;
                    wp_moment_x += wp.area * wp.centroid[0];
                    total_i_t += wp.i_transverse + wp.area * wp.centroid[1].powi(2); // Parallel axis
                    total_i_l += wp.i_longitudinal + wp.area * wp.centroid[0].powi(2); // Parallel axis
                }
            }
            
            if total_wp_area > 0.0 {
                let lcf_calc = wp_moment_x / total_wp_area;
                
                // Transfer to centroidal axes
                let wp_tcf = tcb; // Assume waterplane TCF aligns with COB TCB
                let wp_lcf = lcf_calc;
                
                let i_t_centroidal = total_i_t - total_wp_area * wp_tcf.powi(2);
                let i_l_centroidal = total_i_l - total_wp_area * wp_lcf.powi(2);
                
                // Metacentric radii
                let bm_t = i_t_centroidal / total_volume;
                let bm_l = i_l_centroidal / total_volume;
                
                (total_wp_area, wp_lcf, bm_t.max(0.0), bm_l.max(0.0))
            } else {
                // Fallback if waterplane calculation fails
                (0.0, lcb, 0.0, 0.0)
            }
        } else {
            (0.0, lcb, 0.0, 0.0)
        };
        
        // Set COG with LCB, TCB from equilibrium and provided VCG
        // GMT/GML are only defined if VCG is provided
        let (cog, gmt, gml, gmt_dry, gml_dry) = if let Some(v) = vcg {
            let full_cog = [lcb, tcb, v];
            let km_t = vcb + bmt;
            let km_l = vcb + bml;
            
            // Calculate dry GM (without free surface correction)
            let gmt_dry_val = km_t - v;
            let gml_dry_val = km_l - v;
            
            // Calculate free surface correction from vessel tanks
            let (fsm_t, fsm_l) = self.vessel.get_total_free_surface_moment();
            
            // FSC = FSM / displacement (in meters)
            let fsc_t = if displacement > 1e-6 {
                fsm_t / displacement
            } else {
                0.0
            };
            
            let fsc_l = if displacement > 1e-6 {
                fsm_l / displacement
            } else {
                0.0
            };
            
            // Calculate wet GM (with free surface correction - more conservative)
            let gmt_wet_val = gmt_dry_val - fsc_t;
            let gml_wet_val = gml_dry_val - fsc_l;
            
            (
                Some(full_cog),
                Some(gmt_wet_val),  // gmt = wet (conservative)
                Some(gml_wet_val),  // gml = wet
                Some(gmt_dry_val),  // gmt_dry
                Some(gml_dry_val),  // gml_dry
            )
        } else {
            (None, None, None, None, None)
        };

        Some(HydrostaticState {
            draft,
            trim,
            heel,
            volume: total_volume,
            displacement,
            cob,
            cog,
            waterplane_area,
            lcf,
            bmt,
            bml,
            gmt,
            gml,
            gmt_dry,
            gml_dry,
            lwl: bounds.1 - bounds.0, // Simplified
            bwl: bounds.3 - bounds.2, // Simplified
        })
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
    /// - Displacement + LCG + TCG + VCG → finds draft, level, full COG specified
    /// - Displacement + trim + heel → finds draft with fixed attitude
    pub fn calculate_at_displacement(
        &self,
        displacement_mass: f64,
        cog: Option<[f64; 3]>,
        trim: Option<f64>,
        heel: Option<f64>,
    ) -> Result<HydrostaticState, String> {
        // Validate constraints
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

        let target_volume = displacement_mass / self.water_density;
        let bounds = self.vessel.get_bounds();
        let z_min = bounds.4;
        let z_max = bounds.5;

        let tolerance = target_volume * 1e-4;
        let max_iter = 50;

        let mut low = z_min;
        let mut high = z_max;

        let fixed_trim = trim.unwrap_or(0.0);
        let fixed_heel = heel.unwrap_or(0.0);
        
        // Extract VCG if provided via cog
        let vcg = cog.map(|c| c[2]);

        // Bisection search for draft
        for _ in 0..max_iter {
            let mid = (low + high) / 2.0;

            if let Some(state) = self.calculate_at_draft(mid, fixed_trim, fixed_heel, vcg) {
                let diff = state.volume - target_volume;

                if diff.abs() < tolerance {
                    // Found the draft! Now update COG if full COG was specified
                    let final_cog = if let Some(full_cog) = cog {
                        // User specified full COG, keep it
                        Some(full_cog)
                    } else {
                        // No COG specified or only VCG specified
                        state.cog
                    };
                    
                    return Ok(HydrostaticState {
                        cog: final_cog,
                        ..state
                    });
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

        // Convergence not achieved within max iterations
        // Return best estimate with a warning
        let final_draft = (low + high) / 2.0;
        self.calculate_at_draft(final_draft, fixed_trim, fixed_heel, vcg)
            .map(|state| {
                let final_cog = if let Some(full_cog) = cog {
                    Some(full_cog)
                } else {
                    state.cog
                };
                HydrostaticState {
                    cog: final_cog,
                    ..state
                }
            })
            .ok_or_else(|| {
                format!("Could not find draft for displacement {} kg", displacement_mass)
            })
    }

    /// Returns the water density.
    pub fn water_density(&self) -> f64 {
        self.water_density
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
    fn test_box_barge_volume() {
        // 10m x 10m x 10m box
        let hull = create_box_hull(10.0, 10.0, 10.0);
        let vessel = Vessel::new(hull);
        let calc = HydrostaticsCalculator::new(&vessel, 1025.0);

        // At draft 5m, volume should be 10 * 10 * 5 = 500 m³
        let state = calc.calculate_at_draft(5.0, 0.0, 0.0, None).unwrap();
        assert!(
            (state.volume - 500.0).abs() < 1.0,
            "Volume was {}",
            state.volume
        );
    }


    #[test]
    fn test_calculate_at_displacement_level() {
        let hull = create_box_hull(10.0, 10.0, 10.0);
        let vessel = Vessel::new(hull);
        let calc = HydrostaticsCalculator::new(&vessel, 1025.0);

        // Target displacement: 500 m³ * 1025 kg/m³ = 512500 kg
        let target_disp = 500.0 * 1025.0;
        
        // Calculate at displacement with no other constraints (level keel)
        let state = calc.calculate_at_displacement(
            target_disp,
            None, None, None
        ).expect("Calculation failed");

        assert!((state.draft - 5.0).abs() < 0.01, "Draft should be ~5.0m, got {}", state.draft);
        assert!((state.displacement - target_disp).abs() < 1.0, "Displacement mismatch");
        assert_eq!(state.trim, 0.0);
        assert_eq!(state.heel, 0.0);
    }

    #[test]
    fn test_calculate_at_displacement_with_vcg() {
        let hull = create_box_hull(10.0, 10.0, 10.0);
        let vessel = Vessel::new(hull);
        let calc = HydrostaticsCalculator::new(&vessel, 1025.0);
        let target_disp = 512500.0; // 5m draft condition

        // With VCG provided, should compute GMT/GML
        // Note: LCB/TCB assumed 0.0 for box hull, so just set VCG=7.0
        let state = calc.calculate_at_displacement(
            target_disp,
            Some([0.0, 0.0, 7.0]), None, None
        ).expect("Calculation failed");

        assert!((state.draft - 5.0).abs() < 0.01);
        
        // Check VCG is set
        let cog = state.cog.expect("COG should be set");
        assert_eq!(cog[2], 7.0);

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
        let res = calc.calculate_at_displacement(
             100000.0,
             Some([5.0, 0.0, 0.0]), Some(0.0), None
        );
        assert!(res.is_err(), "Should fail for both LCG and Trim specified");

        // Invalid: Heel provided but also TCG constrained
        let res = calc.calculate_at_displacement(
             100000.0,
             Some([0.0, 5.0, 0.0]), None, Some(0.0)
        );
        assert!(res.is_err(), "Should fail for both TCG and Heel specified");
    }
}
