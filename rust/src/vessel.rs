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

//! Vessel module.
//!
//! Provides the Vessel container for hull geometries, tanks, and vessel-level properties.

use crate::downflooding::DownfloodingOpening;
use crate::hull::Hull;
use crate::silhouette::Silhouette;
use crate::tanks::Tank;

/// Represents a vessel containing hull geometries, tanks, and vessel-level properties.
///
/// The Vessel class serves as a container for hull geometries and manages
/// vessel-level reference positions such as the forward and aft perpendiculars
/// (FP and AP). It supports both single-hull vessels (monohulls) and multi-hull
/// vessels (catamarans, trimarans).
#[derive(Clone)]
pub struct Vessel {
    /// List of hull geometries
    hulls: Vec<Hull>,
    /// List of tanks
    tanks: Vec<Tank>,
    /// Aft Perpendicular position (None = auto from bounds)
    ap: Option<f64>,
    /// Forward Perpendicular position (None = auto from bounds)
    fp: Option<f64>,
    /// Wind silhouette profiles (hull, superstructure, containers, etc.)
    silhouettes: Vec<Silhouette>,
    /// Downflooding openings for θf calculation
    downflooding_openings: Vec<DownfloodingOpening>,
}

impl Vessel {
    /// Creates a new Vessel with a single hull.
    pub fn new(hull: Hull) -> Self {
        Self {
            hulls: vec![hull],
            tanks: Vec::new(),
            ap: None,
            fp: None,
            silhouettes: Vec::new(),
            downflooding_openings: Vec::new(),
        }
    }

    /// Creates a new Vessel with multiple hulls (catamaran, trimaran).
    pub fn new_multi(hulls: Vec<Hull>) -> Result<Self, &'static str> {
        if hulls.is_empty() {
            return Err("At least one hull must be provided");
        }
        Ok(Self {
            hulls,
            tanks: Vec::new(),
            ap: None,
            fp: None,
            silhouettes: Vec::new(),
            downflooding_openings: Vec::new(),
        })
    }

    /// Creates a new Vessel with perpendicular positions.
    pub fn with_perpendiculars(hull: Hull, ap: f64, fp: f64) -> Self {
        Self {
            hulls: vec![hull],
            tanks: Vec::new(),
            ap: Some(ap),
            fp: Some(fp),
            silhouettes: Vec::new(),
            downflooding_openings: Vec::new(),
        }
    }

    /// Returns the list of hull geometries.
    pub fn hulls(&self) -> &[Hull] {
        &self.hulls
    }

    /// Returns a mutable reference to the hulls.
    pub fn hulls_mut(&mut self) -> &mut Vec<Hull> {
        &mut self.hulls
    }

    /// Returns true if this is a multi-hull vessel.
    pub fn is_multihull(&self) -> bool {
        self.hulls.len() > 1
    }

    /// Returns the list of tanks.
    pub fn tanks(&self) -> &[Tank] {
        &self.tanks
    }

    /// Returns a mutable reference to the tanks.
    pub fn tanks_mut(&mut self) -> &mut Vec<Tank> {
        &mut self.tanks
    }

    /// Returns the Aft Perpendicular position.
    ///
    /// If not explicitly set, returns the minimum X of the combined bounds.
    pub fn ap(&self) -> f64 {
        self.ap.unwrap_or_else(|| self.get_bounds().0)
    }

    /// Returns the Forward Perpendicular position.
    ///
    /// If not explicitly set, returns the maximum X of the combined bounds.
    pub fn fp(&self) -> f64 {
        self.fp.unwrap_or_else(|| self.get_bounds().1)
    }

    /// Sets the Aft Perpendicular position.
    pub fn set_ap(&mut self, ap: f64) {
        self.ap = Some(ap);
    }

    /// Sets the Forward Perpendicular position.
    pub fn set_fp(&mut self, fp: f64) {
        self.fp = Some(fp);
    }

    /// Returns the Length Between Perpendiculars (LBP).
    pub fn lbp(&self) -> f64 {
        self.fp() - self.ap()
    }

    /// Returns the bounding box of all hull geometries combined.
    ///
    /// Returns (xmin, xmax, ymin, ymax, zmin, zmax).
    pub fn get_bounds(&self) -> (f64, f64, f64, f64, f64, f64) {
        if self.hulls.len() == 1 {
            return self.hulls[0].get_bounds();
        }

        let all_bounds: Vec<_> = self.hulls.iter().map(|h| h.get_bounds()).collect();

        let xmin = all_bounds.iter().map(|b| b.0).fold(f64::INFINITY, f64::min);
        let xmax = all_bounds
            .iter()
            .map(|b| b.1)
            .fold(f64::NEG_INFINITY, f64::max);
        let ymin = all_bounds.iter().map(|b| b.2).fold(f64::INFINITY, f64::min);
        let ymax = all_bounds
            .iter()
            .map(|b| b.3)
            .fold(f64::NEG_INFINITY, f64::max);
        let zmin = all_bounds.iter().map(|b| b.4).fold(f64::INFINITY, f64::min);
        let zmax = all_bounds
            .iter()
            .map(|b| b.5)
            .fold(f64::NEG_INFINITY, f64::max);

        (xmin, xmax, ymin, ymax, zmin, zmax)
    }

    // =========================================================================
    // Tank Management
    // =========================================================================

    /// Adds a tank to the vessel.
    pub fn add_tank(&mut self, tank: Tank) {
        self.tanks.push(tank);
    }

    /// Removes a tank from the vessel by index.
    pub fn remove_tank(&mut self, index: usize) -> Option<Tank> {
        if index < self.tanks.len() {
            Some(self.tanks.remove(index))
        } else {
            None
        }
    }

    /// Finds a tank by its name.
    pub fn get_tank_by_name(&self, name: &str) -> Option<&Tank> {
        self.tanks.iter().find(|t| t.name() == name)
    }

    /// Finds a tank by its name (mutable).
    pub fn get_tank_by_name_mut(&mut self, name: &str) -> Option<&mut Tank> {
        self.tanks.iter_mut().find(|t| t.name() == name)
    }

    /// Calculates the total mass of all fluid in tanks.
    pub fn get_total_tanks_mass(&self) -> f64 {
        self.tanks.iter().map(|t| t.fluid_mass()).sum()
    }

    /// Calculates the combined center of gravity of all tank fluids.
    ///
    /// Returns the mass-weighted average of individual tank CoGs.
    pub fn get_tanks_center_of_gravity(&self) -> [f64; 3] {
        let total_mass = self.get_total_tanks_mass();
        if total_mass <= 0.0 {
            return [0.0, 0.0, 0.0];
        }

        let mut moment = [0.0, 0.0, 0.0];
        for tank in &self.tanks {
            if tank.fluid_mass() > 0.0 {
                let cog = tank.center_of_gravity();
                moment[0] += tank.fluid_mass() * cog[0];
                moment[1] += tank.fluid_mass() * cog[1];
                moment[2] += tank.fluid_mass() * cog[2];
            }
        }

        [
            moment[0] / total_mass,
            moment[1] / total_mass,
            moment[2] / total_mass,
        ]
    }

    /// Calculates the total free surface moment from all tanks.
    ///
    /// Returns (transverse_moment, longitudinal_moment) in m⁴.
    pub fn get_total_free_surface_moment(&self) -> (f64, f64) {
        let fsm_t: f64 = self.tanks.iter().map(|t| t.free_surface_moment_t()).sum();
        let fsm_l: f64 = self.tanks.iter().map(|t| t.free_surface_moment_l()).sum();
        (fsm_t, fsm_l)
    }

    /// Calculates the total free surface correction from all tanks.
    ///
    /// Returns (transverse_correction, longitudinal_correction) in m⁴.
    pub fn get_total_free_surface_correction(&self) -> (f64, f64) {
        let fsc_t: f64 = self
            .tanks
            .iter()
            .map(|t| t.free_surface_correction_t())
            .sum();
        let fsc_l: f64 = self
            .tanks
            .iter()
            .map(|t| t.free_surface_correction_l())
            .sum();
        (fsc_t, fsc_l)
    }

    // =========================================================================
    // Silhouette Management
    // =========================================================================

    /// Adds a wind silhouette profile to the vessel.
    pub fn add_silhouette(&mut self, silhouette: Silhouette) {
        self.silhouettes.push(silhouette);
    }

    /// Returns a reference to all wind silhouettes.
    pub fn silhouettes(&self) -> &[Silhouette] {
        &self.silhouettes
    }

    /// Returns a mutable reference to the silhouettes.
    pub fn silhouettes_mut(&mut self) -> &mut Vec<Silhouette> {
        &mut self.silhouettes
    }

    /// Returns the number of silhouettes.
    pub fn num_silhouettes(&self) -> usize {
        self.silhouettes.len()
    }

    /// Returns true if there are any silhouettes.
    pub fn has_silhouettes(&self) -> bool {
        !self.silhouettes.is_empty()
    }

    /// Finds a silhouette by its name.
    pub fn get_silhouette_by_name(&self, name: &str) -> Option<&Silhouette> {
        self.silhouettes.iter().find(|s| s.name() == name)
    }

    /// Removes a silhouette by index.
    pub fn remove_silhouette(&mut self, index: usize) -> Option<Silhouette> {
        if index < self.silhouettes.len() {
            Some(self.silhouettes.remove(index))
        } else {
            None
        }
    }

    /// Removes all silhouettes.
    pub fn clear_silhouettes(&mut self) {
        self.silhouettes.clear();
    }

    /// Calculates the total emerged area from all silhouettes.
    pub fn get_total_emerged_area(&self, waterline_z: f64) -> f64 {
        self.silhouettes
            .iter()
            .map(|s| s.get_emerged_area(waterline_z))
            .sum()
    }

    /// Calculates the combined centroid of all emerged areas.
    pub fn get_combined_emerged_centroid(&self, waterline_z: f64) -> [f64; 2] {
        let total_area = self.get_total_emerged_area(waterline_z);
        if total_area < 1e-9 {
            return [0.0, 0.0];
        }

        let mut cx = 0.0;
        let mut cz = 0.0;
        for s in &self.silhouettes {
            let area = s.get_emerged_area(waterline_z);
            if area > 1e-9 {
                let centroid = s.get_emerged_centroid(waterline_z);
                cx += centroid[0] * area;
                cz += centroid[1] * area;
            }
        }

        [cx / total_area, cz / total_area]
    }

    // =========================================================================
    // Downflooding Openings Management
    // =========================================================================

    /// Adds a downflooding opening to the vessel.
    pub fn add_downflooding_opening(&mut self, opening: DownfloodingOpening) {
        self.downflooding_openings.push(opening);
    }

    /// Returns a reference to all downflooding openings.
    pub fn downflooding_openings(&self) -> &[DownfloodingOpening] {
        &self.downflooding_openings
    }

    /// Returns a mutable reference to downflooding openings.
    pub fn downflooding_openings_mut(&mut self) -> &mut Vec<DownfloodingOpening> {
        &mut self.downflooding_openings
    }

    /// Returns the number of downflooding openings.
    pub fn num_downflooding_openings(&self) -> usize {
        self.downflooding_openings.len()
    }

    /// Returns true if any downflooding openings are defined.
    pub fn has_downflooding_openings(&self) -> bool {
        !self.downflooding_openings.is_empty()
    }

    /// Removes all downflooding openings from the vessel.
    pub fn clear_downflooding_openings(&mut self) {
        self.downflooding_openings.clear();
    }
}

impl std::fmt::Debug for Vessel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let bounds = self.get_bounds();
        f.debug_struct("Vessel")
            .field("hulls", &self.hulls.len())
            .field("tanks", &self.tanks.len())
            .field("ap", &self.ap())
            .field("fp", &self.fp())
            .field("lbp", &self.lbp())
            .field("bounds", &bounds)
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;
    use parry3d_f64::shape::TriMesh;

    fn create_test_hull() -> Hull {
        let vertices = vec![
            Point3::new(0.0, -5.0, 0.0),
            Point3::new(100.0, -5.0, 0.0),
            Point3::new(100.0, 5.0, 0.0),
            Point3::new(0.0, 5.0, 0.0),
            Point3::new(0.0, -5.0, 10.0),
            Point3::new(100.0, -5.0, 10.0),
            Point3::new(100.0, 5.0, 10.0),
            Point3::new(0.0, 5.0, 10.0),
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
    fn test_vessel_bounds() {
        let hull = create_test_hull();
        let vessel = Vessel::new(hull);

        let bounds = vessel.get_bounds();
        assert!((bounds.0 - 0.0).abs() < 1e-6);
        assert!((bounds.1 - 100.0).abs() < 1e-6);
    }

    #[test]
    fn test_vessel_perpendiculars() {
        let hull = create_test_hull();
        let vessel = Vessel::new(hull);

        assert!((vessel.ap() - 0.0).abs() < 1e-6);
        assert!((vessel.fp() - 100.0).abs() < 1e-6);
        assert!((vessel.lbp() - 100.0).abs() < 1e-6);
    }
}
