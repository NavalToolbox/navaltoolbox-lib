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

//! Python bindings for navaltoolbox.
//!
//! This module provides PyO3 bindings for the Rust library, exposing
//! Hull, Vessel, HydrostaticsCalculator, StabilityCalculator, and Tank.

use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;

use crate::hull::Hull as RustHull;
use crate::hydrostatics::{
    HydrostaticState as RustHydroState, HydrostaticsCalculator as RustHydroCalc,
};
use crate::scripting::{
    CriteriaContext as RustCriteriaContext, CriteriaResult as RustCriteriaResult,
    ScriptEngine as RustScriptEngine,
};
use crate::stability::{
    CompleteStabilityResult as RustCompleteStabilityResult, StabilityCalculator as RustStabCalc,
    StabilityCurve as RustStabCurve, WindHeelingData as RustWindHeelingData,
};
use crate::tanks::Tank as RustTank;
use crate::vessel::Vessel as RustVessel;

use std::path::Path;

// ============================================================================
// Hull Python Wrapper
// ============================================================================

/// A hull geometry loaded from an STL file.
#[pyclass(name = "Hull")]
pub struct PyHull {
    inner: RustHull,
}

#[pymethods]
impl PyHull {
    /// Load a hull from an STL file.
    #[new]
    fn new(file_path: &str) -> PyResult<Self> {
        let path = Path::new(file_path);
        let hull = RustHull::from_stl(path)
            .map_err(|e| PyIOError::new_err(format!("Failed to load STL: {}", e)))?;
        Ok(Self { inner: hull })
    }

    /// Create a box hull.
    ///
    /// Args:
    ///     length: Length of the box in meters
    ///     breadth: Breadth of the box in meters
    ///     depth: Depth of the box in meters
    #[staticmethod]
    fn from_box(length: f64, breadth: f64, depth: f64) -> Self {
        let hull = RustHull::from_box(length, breadth, depth);
        Self { inner: hull }
    }

    /// Returns the bounding box (xmin, xmax, ymin, ymax, zmin, zmax).
    fn get_bounds(&self) -> (f64, f64, f64, f64, f64, f64) {
        self.inner.get_bounds()
    }

    /// Returns the number of triangles.
    fn num_triangles(&self) -> usize {
        self.inner.num_triangles()
    }

    /// Returns the number of vertices.
    fn num_vertices(&self) -> usize {
        self.inner.num_vertices()
    }

    /// Applies a transformation to the hull.
    fn transform(
        &mut self,
        translation: (f64, f64, f64),
        rotation: (f64, f64, f64),
        pivot: (f64, f64, f64),
    ) {
        self.inner.transform(translation, rotation, pivot);
    }

    /// Scales the hull uniformly.
    fn scale(&mut self, factor: f64) {
        self.inner.scale(factor);
    }

    /// Scales the hull non-uniformly.
    fn scale_xyz(&mut self, sx: f64, sy: f64, sz: f64) {
        self.inner.scale_xyz(sx, sy, sz);
    }

    /// Simplifies the hull mesh to a target number of triangles.
    ///
    /// Args:
    ///     target_count: Target number of triangles for the simplified mesh.
    fn simplify(&mut self, target_count: usize) {
        self.inner.simplify(target_count);
    }

    /// Returns a simplified copy of the hull.
    ///
    /// Args:
    ///     target_count: Target number of triangles for the simplified mesh.
    fn to_simplified(&self, target_count: usize) -> Self {
        Self {
            inner: self.inner.to_simplified(target_count),
        }
    }

    /// Exports the hull to an STL file.
    fn export_stl(&self, file_path: &str) -> PyResult<()> {
        let path = Path::new(file_path);
        self.inner
            .export_stl(path)
            .map_err(|e| PyIOError::new_err(format!("Failed to export STL: {}", e)))
    }

    /// Returns vertices as list of tuples (x, y, z).
    fn get_vertices(&self) -> Vec<(f64, f64, f64)> {
        self.inner
            .mesh()
            .vertices()
            .iter()
            .map(|v| (v.x, v.y, v.z))
            .collect()
    }

    /// Returns faces as list of tuples (i, j, k).
    fn get_faces(&self) -> Vec<(u32, u32, u32)> {
        self.inner
            .mesh()
            .indices()
            .iter()
            .map(|idx| (idx[0], idx[1], idx[2]))
            .collect()
    }

    fn __repr__(&self) -> String {
        let bounds = self.inner.get_bounds();
        format!(
            "Hull(triangles={}, vertices={}, bounds=({:.2}, {:.2}, {:.2}, {:.2}, {:.2}, {:.2}))",
            self.inner.num_triangles(),
            self.inner.num_vertices(),
            bounds.0,
            bounds.1,
            bounds.2,
            bounds.3,
            bounds.4,
            bounds.5
        )
    }
}

// ============================================================================
// Vessel Python Wrapper
// ============================================================================

/// A vessel containing one or more hulls, tanks, and silhouettes.
#[pyclass(name = "Vessel")]
pub struct PyVessel {
    inner: RustVessel,
}

#[pymethods]
impl PyVessel {
    /// Create a vessel from a hull.
    #[new]
    fn new(hull: &PyHull) -> Self {
        Self {
            inner: RustVessel::new(hull.inner.clone()),
        }
    }

    /// Returns the bounding box of all hulls.
    fn get_bounds(&self) -> (f64, f64, f64, f64, f64, f64) {
        self.inner.get_bounds()
    }

    /// Returns the Aft Perpendicular position.
    #[getter]
    fn ap(&self) -> f64 {
        self.inner.ap()
    }

    /// Sets the Aft Perpendicular position.
    #[setter]
    fn set_ap(&mut self, ap: f64) {
        self.inner.set_ap(ap);
    }

    /// Returns the Forward Perpendicular position.
    #[getter]
    fn fp(&self) -> f64 {
        self.inner.fp()
    }

    /// Sets the Forward Perpendicular position.
    #[setter]
    fn set_fp(&mut self, fp: f64) {
        self.inner.set_fp(fp);
    }

    /// Returns the Length Between Perpendiculars.
    #[getter]
    fn lbp(&self) -> f64 {
        self.inner.lbp()
    }

    /// Returns the number of hulls.
    fn num_hulls(&self) -> usize {
        self.inner.hulls().len()
    }

    /// Returns the number of tanks.
    fn num_tanks(&self) -> usize {
        self.inner.tanks().len()
    }

    /// Add a tank to the vessel.
    fn add_tank(&mut self, tank: &PyTank) {
        self.inner.add_tank(tank.inner.clone());
    }

    /// Returns the total tanks mass in kg.
    fn get_total_tanks_mass(&self) -> f64 {
        self.inner.get_total_tanks_mass()
    }

    /// Returns the tanks center of gravity [x, y, z].
    fn get_tanks_center_of_gravity(&self) -> [f64; 3] {
        self.inner.get_tanks_center_of_gravity()
    }

    // =========================================================================
    // Silhouette methods
    // =========================================================================

    /// Add a silhouette profile to the vessel.
    fn add_silhouette(&mut self, silhouette: &PySilhouette) {
        self.inner.add_silhouette(silhouette.inner.clone());
    }

    /// Returns the number of silhouettes.
    fn num_silhouettes(&self) -> usize {
        self.inner.num_silhouettes()
    }

    /// Returns true if there are any silhouettes.
    fn has_silhouettes(&self) -> bool {
        self.inner.has_silhouettes()
    }

    /// Removes all silhouettes.
    fn clear_silhouettes(&mut self) {
        self.inner.clear_silhouettes();
    }

    /// Returns the total emerged area from all silhouettes (m²).
    fn get_total_emerged_area(&self, waterline_z: f64) -> f64 {
        self.inner.get_total_emerged_area(waterline_z)
    }

    /// Returns the combined emerged centroid [x, z].
    fn get_combined_emerged_centroid(&self, waterline_z: f64) -> [f64; 2] {
        self.inner.get_combined_emerged_centroid(waterline_z)
    }

    // =========================================================================
    // Downflooding Opening methods
    // =========================================================================

    /// Add a downflooding opening to the vessel.
    fn add_opening(&mut self, opening: &PyDownfloodingOpening) {
        self.inner.add_downflooding_opening(opening.inner.clone());
    }

    /// Returns the number of downflooding openings.
    fn num_openings(&self) -> usize {
        self.inner.num_downflooding_openings()
    }

    /// Removes all downflooding openings.
    fn clear_openings(&mut self) {
        self.inner.clear_downflooding_openings();
    }

    // =========================================================================
    // Component Getters for Visualization
    // =========================================================================

    /// Get all hulls.
    fn get_hulls(&self) -> Vec<PyHull> {
        self.inner
            .hulls()
            .iter()
            .map(|h| PyHull { inner: h.clone() })
            .collect()
    }

    /// Get all tanks.
    fn get_tanks(&self) -> Vec<PyTank> {
        self.inner
            .tanks()
            .iter()
            .map(|t| PyTank { inner: t.clone() })
            .collect()
    }

    /// Get all silhouettes.
    fn get_silhouettes(&self) -> Vec<PySilhouette> {
        self.inner
            .silhouettes()
            .iter()
            .map(|s| PySilhouette { inner: s.clone() })
            .collect()
    }

    /// Get all downflooding openings.
    fn get_openings(&self) -> Vec<PyDownfloodingOpening> {
        self.inner
            .downflooding_openings()
            .iter()
            .map(|o| PyDownfloodingOpening { inner: o.clone() })
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "Vessel(hulls={}, tanks={}, silhouettes={}, lbp={:.2}m)",
            self.inner.hulls().len(),
            self.inner.tanks().len(),
            self.inner.num_silhouettes(),
            self.inner.lbp()
        )
    }
}

// ============================================================================
// Silhouette Python Wrapper
// ============================================================================

use crate::silhouette::Silhouette as RustSilhouette;

/// A 2D silhouette profile in the X-Z plane for wind heeling calculations.
#[pyclass(name = "Silhouette")]
pub struct PySilhouette {
    inner: RustSilhouette,
}

#[pymethods]
impl PySilhouette {
    /// Load a silhouette from a file (DXF or VTK).
    #[new]
    fn new(file_path: &str) -> PyResult<Self> {
        let path = Path::new(file_path);
        let silhouette = RustSilhouette::from_file(path)
            .map_err(|e| PyIOError::new_err(format!("Failed to load silhouette: {}", e)))?;
        Ok(Self { inner: silhouette })
    }

    /// Create a silhouette from a list of points [(x, z), ...].
    #[staticmethod]
    fn from_points(points: Vec<(f64, f64)>, name: &str) -> Self {
        let pts: Vec<[f64; 3]> = points.iter().map(|(x, z)| [*x, 0.0, *z]).collect();
        Self {
            inner: RustSilhouette::new(pts, name.to_string()),
        }
    }

    /// Returns the silhouette name.
    #[getter]
    fn name(&self) -> &str {
        self.inner.name()
    }

    /// Returns the number of points.
    fn num_points(&self) -> usize {
        self.inner.num_points()
    }

    /// Returns true if the contour is closed.
    fn is_closed(&self) -> bool {
        self.inner.is_closed()
    }

    /// Returns the points as a list of tuples [(x, y, z), ...].
    fn get_points(&self) -> Vec<(f64, f64, f64)> {
        self.inner
            .points()
            .iter()
            .map(|p| (p[0], p[1], p[2]))
            .collect()
    }

    /// Returns the total lateral area in m².
    fn get_area(&self) -> f64 {
        self.inner.get_area()
    }

    /// Returns the centroid [x, z].
    fn get_centroid(&self) -> [f64; 2] {
        self.inner.get_centroid()
    }

    /// Returns the bounding box (x_min, x_max, z_min, z_max).
    fn get_bounds(&self) -> (f64, f64, f64, f64) {
        self.inner.get_bounds()
    }

    /// Returns the emerged area above waterline (m²).
    fn get_emerged_area(&self, waterline_z: f64) -> f64 {
        self.inner.get_emerged_area(waterline_z)
    }

    /// Returns the centroid of emerged area [x, z].
    fn get_emerged_centroid(&self, waterline_z: f64) -> [f64; 2] {
        self.inner.get_emerged_centroid(waterline_z)
    }

    fn __repr__(&self) -> String {
        format!(
            "Silhouette(name='{}', points={}, area={:.2}m²)",
            self.inner.name(),
            self.inner.num_points(),
            self.inner.get_area()
        )
    }
}

// ============================================================================
// DownfloodingOpening Python Wrapper
// ============================================================================

use crate::downflooding::{
    DownfloodingOpening as RustDownfloodingOpening, OpeningGeometry, OpeningType as RustOpeningType,
};

/// Type of opening that can cause downflooding.
#[pyclass(name = "OpeningType")]
#[derive(Clone)]
pub struct PyOpeningType {
    inner: RustOpeningType,
}

#[pymethods]
impl PyOpeningType {
    #[staticmethod]
    fn vent() -> Self {
        Self {
            inner: RustOpeningType::Vent,
        }
    }

    #[staticmethod]
    fn air_pipe() -> Self {
        Self {
            inner: RustOpeningType::AirPipe,
        }
    }

    #[staticmethod]
    fn hatch() -> Self {
        Self {
            inner: RustOpeningType::Hatch,
        }
    }

    #[staticmethod]
    fn door() -> Self {
        Self {
            inner: RustOpeningType::Door,
        }
    }

    #[staticmethod]
    fn window() -> Self {
        Self {
            inner: RustOpeningType::Window,
        }
    }

    #[staticmethod]
    fn other(name: &str) -> Self {
        Self {
            inner: RustOpeningType::Other(name.to_string()),
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

/// A downflooding opening point or contour.
#[pyclass(name = "DownfloodingOpening")]
pub struct PyDownfloodingOpening {
    pub(crate) inner: RustDownfloodingOpening,
}

#[pymethods]
impl PyDownfloodingOpening {
    /// Create a downflooding opening from a single point.
    #[staticmethod]
    fn from_point(name: &str, position: (f64, f64, f64), opening_type: &PyOpeningType) -> Self {
        Self {
            inner: RustDownfloodingOpening::from_point(
                name.to_string(),
                [position.0, position.1, position.2],
                opening_type.inner.clone(),
            ),
        }
    }

    /// Create a downflooding opening from a contour (polyline).
    #[staticmethod]
    fn from_contour(
        name: &str,
        points: Vec<(f64, f64, f64)>,
        opening_type: &PyOpeningType,
    ) -> Self {
        let pts: Vec<[f64; 3]> = points.iter().map(|(x, y, z)| [*x, *y, *z]).collect();
        Self {
            inner: RustDownfloodingOpening::from_contour(
                name.to_string(),
                pts,
                opening_type.inner.clone(),
            ),
        }
    }

    /// Load openings from a file (DXF or VTK).
    /// Returns a list of DownfloodingOpening objects.
    #[staticmethod]
    #[pyo3(signature = (file_path, default_type, name=None))]
    fn from_file(
        file_path: &str,
        default_type: &PyOpeningType,
        name: Option<String>,
    ) -> PyResult<Vec<Self>> {
        let path = Path::new(file_path);
        let mut openings = RustDownfloodingOpening::from_file(path, default_type.inner.clone())
            .map_err(|e| PyIOError::new_err(format!("Failed to load openings: {}", e)))?;

        // If name is provided, rename openings
        if let Some(base_name) = name {
            if openings.len() == 1 {
                openings[0].set_name(base_name);
            } else {
                for (i, opening) in openings.iter_mut().enumerate() {
                    opening.set_name(format!("{}_{}", base_name, i + 1));
                }
            }
        }

        Ok(openings.into_iter().map(|o| Self { inner: o }).collect())
    }

    /// Returns the opening name.
    #[getter]
    fn name(&self) -> &str {
        self.inner.name()
    }

    /// Check if opening is active.
    #[getter]
    fn is_active(&self) -> bool {
        self.inner.is_active()
    }

    /// Set opening active state.
    fn set_active(&mut self, active: bool) {
        self.inner.set_active(active);
    }

    /// Get number of points.
    fn num_points(&self) -> usize {
        self.inner.get_points().len()
    }

    /// Get all points as [(x, y, z), ...].
    fn get_points(&self) -> Vec<(f64, f64, f64)> {
        self.inner
            .get_points()
            .iter()
            .map(|p| (p[0], p[1], p[2]))
            .collect()
    }

    /// Check if submerged at given heel/trim/draft.
    fn is_submerged(&self, heel: f64, trim: f64, pivot: (f64, f64, f64), waterline_z: f64) -> bool {
        self.inner
            .is_submerged(heel, trim, [pivot.0, pivot.1, pivot.2], waterline_z)
    }

    fn __repr__(&self) -> String {
        let pts = self.inner.get_points();
        let geometry = match &self.inner.geometry() {
            OpeningGeometry::Point(_) => "Point",
            OpeningGeometry::Contour(_) => "Contour",
        };
        format!(
            "DownfloodingOpening(name='{}', type={:?}, geometry={}, points={})",
            self.inner.name(),
            self.inner.opening_type(),
            geometry,
            pts.len()
        )
    }
}

// ============================================================================
// HydrostaticState Python Wrapper
// ============================================================================

/// Result of hydrostatic calculations.
#[pyclass(name = "HydrostaticState")]
#[derive(Clone)]
pub struct PyHydrostaticState {
    #[pyo3(get)]
    pub draft: f64,
    #[pyo3(get)]
    pub trim: f64,
    #[pyo3(get)]
    pub heel: f64,
    #[pyo3(get)]
    pub draft_ap: f64,
    #[pyo3(get)]
    pub draft_fp: f64,
    #[pyo3(get)]
    pub draft_mp: f64,
    #[pyo3(get)]
    pub volume: f64,
    #[pyo3(get)]
    pub displacement: f64,

    // Internal storage as vectors
    cob_internal: [f64; 3],
    cog_internal: Option<[f64; 3]>,

    // Expose all hydrostatic properties
    #[pyo3(get)]
    pub waterplane_area: f64,
    #[pyo3(get)]
    pub lcf: f64,
    #[pyo3(get)]
    pub bmt: f64,
    #[pyo3(get)]
    pub bml: f64,

    // Optional GMT/GML (None if VCG not specified)
    gmt_internal: Option<f64>,     // wet (with FSC)
    gml_internal: Option<f64>,     // wet (with FSC)
    gmt_dry_internal: Option<f64>, // dry (without FSC)
    gml_dry_internal: Option<f64>, // dry (without FSC)

    #[pyo3(get)]
    pub lwl: f64,
    #[pyo3(get)]
    pub bwl: f64,
    #[pyo3(get)]
    pub los: f64,

    #[pyo3(get)]
    pub wetted_surface_area: f64,
    #[pyo3(get)]
    pub midship_area: f64,
    #[pyo3(get)]
    pub cm: f64,
    #[pyo3(get)]
    pub cb: f64,
    #[pyo3(get)]
    pub cp: f64,
    #[pyo3(get)]
    pub free_surface_correction_t: f64,
    #[pyo3(get)]
    pub free_surface_correction_l: f64,
    #[pyo3(get)]
    pub stiffness_matrix: Vec<f64>,
}

impl From<RustHydroState> for PyHydrostaticState {
    fn from(state: RustHydroState) -> Self {
        Self {
            draft: state.draft,
            trim: state.trim,
            heel: state.heel,
            draft_ap: state.draft_ap,
            draft_fp: state.draft_fp,
            draft_mp: state.draft_mp,
            volume: state.volume,
            displacement: state.displacement,
            cob_internal: state.cob,
            cog_internal: state.cog,
            waterplane_area: state.waterplane_area,
            lcf: state.lcf,
            bmt: state.bmt,
            bml: state.bml,
            gmt_internal: state.gmt,
            gml_internal: state.gml,
            gmt_dry_internal: state.gmt_dry,
            gml_dry_internal: state.gml_dry,
            lwl: state.lwl,
            bwl: state.bwl,
            los: state.los,
            wetted_surface_area: state.wetted_surface_area,
            midship_area: state.midship_area,
            cm: state.cm,
            cb: state.cb,
            cp: state.cp,
            free_surface_correction_t: state.free_surface_correction_t,
            free_surface_correction_l: state.free_surface_correction_l,
            stiffness_matrix: state.stiffness_matrix.to_vec(),
        }
    }
}

#[pymethods]
impl PyHydrostaticState {
    /// Returns center of buoyancy as tuple (lcb, tcb, vcb)
    #[getter]
    fn cob(&self) -> (f64, f64, f64) {
        (
            self.cob_internal[0],
            self.cob_internal[1],
            self.cob_internal[2],
        )
    }

    /// Returns center of gravity as tuple (lcg, tcg, vcg) if specified, None otherwise
    #[getter]
    fn cog(&self) -> Option<(f64, f64, f64)> {
        self.cog_internal.map(|c| (c[0], c[1], c[2]))
    }

    // Convenience getters for individual COB components
    #[getter]
    fn lcb(&self) -> f64 {
        self.cob_internal[0]
    }

    #[getter]
    fn tcb(&self) -> f64 {
        self.cob_internal[1]
    }

    #[getter]
    fn vcb(&self) -> f64 {
        self.cob_internal[2]
    }

    // Convenience getters for individual COG components
    #[getter]
    fn lcg(&self) -> Option<f64> {
        self.cog_internal.map(|c| c[0])
    }

    #[getter]
    fn tcg(&self) -> Option<f64> {
        self.cog_internal.map(|c| c[1])
    }

    #[getter]
    fn vcg(&self) -> Option<f64> {
        self.cog_internal.map(|c| c[2])
    }

    // GMT/GML getters (optional)
    /// GMT with free surface correction (wet - conservative)
    #[getter]
    fn gmt(&self) -> Option<f64> {
        self.gmt_internal
    }

    /// GML with free surface correction (wet - conservative)
    #[getter]
    fn gml(&self) -> Option<f64> {
        self.gml_internal
    }

    /// GMT without free surface correction (dry - reference)
    #[getter]
    fn gmt_dry(&self) -> Option<f64> {
        self.gmt_dry_internal
    }

    /// GML without free surface correction (dry - reference)
    #[getter]
    fn gml_dry(&self) -> Option<f64> {
        self.gml_dry_internal
    }

    fn __repr__(&self) -> String {
        let cog_str = if let Some(c) = self.cog_internal {
            format!("COG=({:.2}, {:.2}, {:.2})", c[0], c[1], c[2])
        } else {
            "COG=None".to_string()
        };

        format!(
            "HydrostaticState(draft_mp={:.3}m (AP={:.3}, FP={:.3}), volume={:.2}m³, displacement={:.0}kg, {})",
            self.draft_mp, self.draft_ap, self.draft_fp, self.volume, self.displacement, cog_str
        )
    }
}

// ============================================================================
// HydrostaticsCalculator Python Wrapper
// ============================================================================

/// Calculator for hydrostatic properties.
#[pyclass(name = "HydrostaticsCalculator")]
pub struct PyHydrostaticsCalculator {
    vessel: RustVessel,
    water_density: f64,
}

#[pymethods]
impl PyHydrostaticsCalculator {
    /// Create a hydrostatics calculator for a vessel.
    #[new]
    #[pyo3(signature = (vessel, water_density=1025.0))]
    fn new(vessel: &PyVessel, water_density: f64) -> Self {
        Self {
            vessel: vessel.inner.clone(),
            water_density,
        }
    }

    /// Calculate hydrostatics at a given draft, trim, and heel.
    ///
    /// Args:
    ///     draft: Draft at reference point in meters
    ///     trim: Trim angle in degrees (default 0.0)
    ///     heel: Heel angle in degrees (default 0.0)
    ///     vcg: Optional vertical center of gravity for GMT/GML calculation
    ///
    /// Returns:
    ///     HydrostaticState with all properties
    #[pyo3(signature = (draft, trim=0.0, heel=0.0, vcg=None), name = "from_draft")]
    #[allow(clippy::wrong_self_convention)]
    fn from_draft(
        &self,
        draft: f64,
        trim: f64,
        heel: f64,
        vcg: Option<f64>,
    ) -> PyResult<PyHydrostaticState> {
        let calc = RustHydroCalc::new(&self.vessel, self.water_density);
        calc.from_draft(draft, trim, heel, vcg)
            .map(|s| s.into())
            .ok_or_else(|| PyValueError::new_err("No submerged volume at this draft"))
    }

    /// Calculate hydrostatics from drafts at Aft and Forward Perpendiculars.
    ///
    /// Args:
    ///     draft_ap: Draft at Aft Perpendicular in meters.
    ///     draft_fp: Draft at Forward Perpendicular in meters.
    ///     heel: Heel angle in degrees (default 0.0).
    ///     vcg: Optional vertical center of gravity for GM calculation.
    ///
    /// Returns:
    ///     HydrostaticState with all properties
    ///
    /// Raises:
    ///     ValueError: If no submerged volume at this draft.
    #[pyo3(signature = (draft_ap, draft_fp, heel=0.0, vcg=None))]
    #[allow(clippy::wrong_self_convention)]
    fn from_drafts(
        &self,
        draft_ap: f64,
        draft_fp: f64,
        heel: f64,
        vcg: Option<f64>,
    ) -> PyResult<PyHydrostaticState> {
        let calc = RustHydroCalc::new(&self.vessel, self.water_density);
        calc.from_drafts(draft_ap, draft_fp, heel, vcg)
            .map(|s| s.into())
            .ok_or_else(|| PyValueError::new_err("No submerged volume at these drafts"))
    }

    /// Calculate hydrostatics for a given displacement with optional constraints.
    ///
    /// Args:
    ///     displacement_mass: Target displacement in kg
    ///     vcg: Optional vertical center of gravity (m) for GM calculations
    ///     cog: Optional (lcg, tcg, vcg) tuple in meters for full COG specification
    ///          (overrides vcg if both are provided)
    ///     trim: Optional trim angle in degrees (default 0.0)
    ///     heel: Optional heel angle in degrees (default 0.0)
    ///
    /// Returns:
    ///     Complete HydrostaticState
    ///
    /// Raises:
    ///     ValueError: If constraints are invalid or unsatisfiable
    ///
    /// Examples:
    ///     >>> # Basic: find draft for displacement
    ///     >>> state = calc.from_displacement(8635000.0)
    ///     
    ///     >>> # With VCG only: compute GMT/GML
    ///     >>> state = calc.from_displacement(8635000.0, vcg=7.555)
    ///     
    ///     >>> # With full COG: for trim optimization
    ///     >>> state = calc.from_displacement(8635000.0, cog=(71.67, 0.0, 7.555))
    ///     
    ///     >>> # With trim constraint
    ///     >>> state = calc.from_displacement(8635000.0, vcg=7.5, trim=2.0)
    #[pyo3(signature = (displacement_mass, vcg=None, cog=None, trim=None, heel=None), name = "from_displacement")]
    #[allow(clippy::wrong_self_convention)]
    fn from_displacement(
        &self,
        displacement_mass: f64,
        vcg: Option<f64>,
        cog: Option<(f64, f64, f64)>,
        trim: Option<f64>,
        heel: Option<f64>,
    ) -> PyResult<PyHydrostaticState> {
        let calc = RustHydroCalc::new(&self.vessel, self.water_density);

        // Convert cog tuple to array if provided
        let cog_array = cog.map(|(lcg, tcg, vcg_val)| [lcg, tcg, vcg_val]);

        calc.from_displacement(displacement_mass, vcg, cog_array, trim, heel)
            .map(|s| s.into())
            .map_err(PyValueError::new_err)
    }

    /// Returns the water density.\n    #[getter]
    fn water_density(&self) -> f64 {
        self.water_density
    }
}

// ============================================================================
// StabilityPoint and StabilityCurve Python Wrappers
// ============================================================================

/// A point on a stability curve.
#[pyclass(name = "StabilityPoint")]
#[derive(Clone)]
pub struct PyStabilityPoint {
    #[pyo3(get)]
    pub heel: f64,
    #[pyo3(get)]
    pub draft: f64,
    #[pyo3(get)]
    pub trim: f64,
    #[pyo3(get)]
    pub gz: f64,
    #[pyo3(get)]
    pub is_flooding: bool,
    #[pyo3(get)]
    pub flooded_openings: Vec<String>,
}

/// A complete GZ stability curve.
#[pyclass(name = "StabilityCurve")]
pub struct PyStabilityCurve {
    inner: RustStabCurve,
}

#[pymethods]
impl PyStabilityCurve {
    /// Returns the heel angles.
    fn heels(&self) -> Vec<f64> {
        self.inner.heels()
    }

    /// Returns the GZ values.
    fn values(&self) -> Vec<f64> {
        self.inner.values()
    }

    /// Returns the points as a list of tuples (heel, draft, trim, gz).
    fn points(&self) -> Vec<(f64, f64, f64, f64)> {
        self.inner
            .points
            .iter()
            .map(|p| (p.heel, p.draft, p.trim, p.value))
            .collect()
    }

    /// Returns the points as a list of StabilityPoint objects.
    fn get_stability_points(&self) -> Vec<PyStabilityPoint> {
        self.inner
            .points
            .iter()
            .map(|p| PyStabilityPoint {
                heel: p.heel,
                draft: p.draft,
                trim: p.trim,
                gz: p.value,
                is_flooding: p.is_flooding,
                flooded_openings: p.flooded_openings.clone(),
            })
            .collect()
    }

    /// Returns the displacement in kg.
    #[getter]
    fn displacement(&self) -> f64 {
        self.inner.displacement
    }

    fn __repr__(&self) -> String {
        format!(
            "StabilityCurve(displacement={:.0}kg, points={})",
            self.inner.displacement,
            self.inner.points.len()
        )
    }
}

// ============================================================================
// WindHeelingData Python Wrapper
// ============================================================================

/// Wind heeling data from silhouette calculations.
#[pyclass(name = "WindHeelingData")]
#[derive(Clone)]
pub struct PyWindHeelingData {
    #[pyo3(get)]
    pub emerged_area: f64,
    emerged_centroid_internal: [f64; 2],
    #[pyo3(get)]
    pub wind_lever_arm: f64,
    #[pyo3(get)]
    pub waterline_z: f64,
}

impl From<RustWindHeelingData> for PyWindHeelingData {
    fn from(data: RustWindHeelingData) -> Self {
        Self {
            emerged_area: data.emerged_area,
            emerged_centroid_internal: data.emerged_centroid,
            wind_lever_arm: data.wind_lever_arm,
            waterline_z: data.waterline_z,
        }
    }
}

#[pymethods]
impl PyWindHeelingData {
    /// Returns the centroid of emerged area [x, z].
    #[getter]
    fn emerged_centroid(&self) -> (f64, f64) {
        (
            self.emerged_centroid_internal[0],
            self.emerged_centroid_internal[1],
        )
    }

    fn __repr__(&self) -> String {
        format!(
            "WindHeelingData(emerged_area={:.2}m², lever_arm={:.2}m)",
            self.emerged_area, self.wind_lever_arm
        )
    }
}

// ============================================================================
// CompleteStabilityResult Python Wrapper
// ============================================================================

/// Complete stability calculation result.
///
/// Combines hydrostatic properties, GZ curve, and wind heeling data
/// for a single loading condition.
#[pyclass(name = "CompleteStabilityResult")]
pub struct PyCompleteStabilityResult {
    inner: RustCompleteStabilityResult,
}

impl From<RustCompleteStabilityResult> for PyCompleteStabilityResult {
    fn from(result: RustCompleteStabilityResult) -> Self {
        Self { inner: result }
    }
}

#[pymethods]
impl PyCompleteStabilityResult {
    /// Returns the hydrostatic state at equilibrium.
    #[getter]
    fn hydrostatics(&self) -> PyHydrostaticState {
        self.inner.hydrostatics.clone().into()
    }

    /// Returns the GZ stability curve.
    #[getter]
    fn gz_curve(&self) -> PyStabilityCurve {
        PyStabilityCurve {
            inner: self.inner.gz_curve.clone(),
        }
    }

    /// Returns the wind heeling data (if silhouettes are defined).
    #[getter]
    fn wind_data(&self) -> Option<PyWindHeelingData> {
        self.inner.wind_data.clone().map(|d| d.into())
    }

    /// Returns the displacement mass in kg.
    #[getter]
    fn displacement(&self) -> f64 {
        self.inner.displacement
    }

    /// Returns the center of gravity (lcg, tcg, vcg).
    #[getter]
    fn cog(&self) -> (f64, f64, f64) {
        (self.inner.cog[0], self.inner.cog[1], self.inner.cog[2])
    }

    /// Returns the initial transverse metacentric height (GM0) with FSC.
    #[getter]
    fn gm0(&self) -> Option<f64> {
        self.inner.gm0()
    }

    /// Returns the initial transverse metacentric height without FSC.
    #[getter]
    fn gm0_dry(&self) -> Option<f64> {
        self.inner.gm0_dry()
    }

    /// Returns the maximum GZ value.
    #[getter]
    fn max_gz(&self) -> Option<f64> {
        self.inner.max_gz()
    }

    /// Returns the heel angle at maximum GZ.
    #[getter]
    fn heel_at_max_gz(&self) -> Option<f64> {
        self.inner.heel_at_max_gz()
    }

    /// Returns true if wind heeling data is available.
    fn has_wind_data(&self) -> bool {
        self.inner.has_wind_data()
    }

    fn __repr__(&self) -> String {
        let gm_str = self
            .inner
            .gm0()
            .map(|gm| format!("{:.3}m", gm))
            .unwrap_or_else(|| "N/A".to_string());
        let max_gz_str = self
            .inner
            .max_gz()
            .map(|gz| format!("{:.3}m", gz))
            .unwrap_or_else(|| "N/A".to_string());
        format!(
            "CompleteStabilityResult(GM0={}, max_GZ={}, wind_data={})",
            gm_str,
            max_gz_str,
            self.inner.has_wind_data()
        )
    }
}

// ============================================================================
// StabilityCalculator Python Wrapper
// ============================================================================

/// Calculator for stability curves (GZ).
#[pyclass(name = "StabilityCalculator")]
pub struct PyStabilityCalculator {
    vessel: RustVessel,
    water_density: f64,
}

#[pymethods]
impl PyStabilityCalculator {
    /// Create a stability calculator for a vessel.
    #[new]
    #[pyo3(signature = (vessel, water_density=1025.0))]
    fn new(vessel: &PyVessel, water_density: f64) -> Self {
        Self {
            vessel: vessel.inner.clone(),
            water_density,
        }
    }

    /// Calculate the GZ curve for a given loading condition.
    fn gz_curve(
        &self,
        displacement_mass: f64,
        cog: (f64, f64, f64),
        heels: Vec<f64>,
    ) -> PyStabilityCurve {
        let calc = RustStabCalc::new(&self.vessel, self.water_density);
        let curve = calc.gz_curve(displacement_mass, [cog.0, cog.1, cog.2], &heels);
        PyStabilityCurve { inner: curve }
    }

    /// Calculate KN curves (Righting Lever from Keel) for multiple displacements.
    ///
    /// This is equivalent to calculating GZ curves with VCG = 0.
    ///
    /// Args:
    ///     displacements: List of target displacements in kg
    ///     lcg: Longitudinal Center of Gravity (m) (default 0.0)
    ///     tcg: Transverse Center of Gravity (m) (default 0.0)
    ///     heels: List of heel angles in degrees
    ///
    /// Returns:
    ///     List[StabilityCurve]: One curve per displacement
    #[pyo3(signature = (displacements, heels, lcg=0.0, tcg=0.0))]
    fn kn_curve(
        &self,
        displacements: Vec<f64>,
        heels: Vec<f64>,
        lcg: f64,
        tcg: f64,
    ) -> Vec<PyStabilityCurve> {
        let calc = RustStabCalc::new(&self.vessel, self.water_density);
        let curves = calc.kn_curve(&displacements, lcg, tcg, &heels);
        curves
            .into_iter()
            .map(|c| PyStabilityCurve { inner: c })
            .collect()
    }

    /// Calculate complete stability analysis for a loading condition.
    ///
    /// Combines hydrostatic calculations, GZ curve, and wind heeling data
    /// (if silhouettes are available) for a single loading condition.
    ///
    /// Args:
    ///     displacement_mass: Target displacement in kg
    ///     cog: Center of gravity (lcg, tcg, vcg) tuple
    ///     heels: List of heel angles for GZ curve in degrees
    ///
    /// Returns:
    ///     CompleteStabilityResult with hydrostatics, GZ curve, and wind data
    fn complete_stability(
        &self,
        displacement_mass: f64,
        cog: (f64, f64, f64),
        heels: Vec<f64>,
    ) -> PyCompleteStabilityResult {
        let calc = RustStabCalc::new(&self.vessel, self.water_density);
        let result = calc.complete_stability(displacement_mass, [cog.0, cog.1, cog.2], &heels);
        result.into()
    }
}

// ============================================================================
// Tank Python Wrapper
// ============================================================================

/// A tank with fluid management capabilities.
#[pyclass(name = "Tank")]
pub struct PyTank {
    inner: RustTank,
}

#[pymethods]
impl PyTank {
    /// Create a Tank from a file (STL or VTK).
    #[new]
    #[pyo3(signature = (file_path, fluid_density=1025.0, name=None))]
    fn new(file_path: &str, fluid_density: f64, name: Option<&str>) -> PyResult<Self> {
        let path = Path::new(file_path);
        let mut tank = RustTank::from_file(path, fluid_density)
            .map_err(|e| PyValueError::new_err(format!("Failed to load tank: {}", e)))?;

        if let Some(n) = name {
            tank.set_name(n);
        }
        Ok(Self { inner: tank })
    }

    /// Create a Tank as the intersection of a box with a hull geometry.
    #[staticmethod]
    #[pyo3(signature = (hull, x_min, x_max, y_min, y_max, z_min, z_max, fluid_density=1025.0, name="HullTank"))]
    #[allow(clippy::too_many_arguments)]
    fn from_box_hull_intersection(
        hull: &PyHull,
        x_min: f64,
        x_max: f64,
        y_min: f64,
        y_max: f64,
        z_min: f64,
        z_max: f64,
        fluid_density: f64,
        name: &str,
    ) -> PyResult<Self> {
        let tank = RustTank::from_box_hull_intersection(
            name,
            &hull.inner,
            x_min,
            x_max,
            y_min,
            y_max,
            z_min,
            z_max,
            fluid_density,
        )
        .map_err(|e| {
            PyValueError::new_err(format!("Failed to create tank from intersection: {}", e))
        })?;

        Ok(Self { inner: tank })
    }

    /// Create a box-shaped tank.
    #[staticmethod]
    #[allow(clippy::too_many_arguments)]
    fn from_box(
        name: &str,
        x_min: f64,
        x_max: f64,
        y_min: f64,
        y_max: f64,
        z_min: f64,
        z_max: f64,
        fluid_density: f64,
    ) -> Self {
        Self {
            inner: RustTank::from_box(
                name,
                x_min,
                x_max,
                y_min,
                y_max,
                z_min,
                z_max,
                fluid_density,
            ),
        }
    }

    /// Returns the tank name.
    #[getter]
    fn name(&self) -> &str {
        self.inner.name()
    }

    /// Returns the total volume in m³.
    #[getter]
    fn total_volume(&self) -> f64 {
        self.inner.total_volume()
    }

    /// Returns the fill level as a fraction (0-1).
    #[getter]
    fn fill_level(&self) -> f64 {
        self.inner.fill_level()
    }

    /// Sets the fill level as a fraction (0-1).
    #[setter]
    fn set_fill_level(&mut self, level: f64) {
        self.inner.set_fill_level(level);
    }

    /// Returns the fill level as a percentage (0-100).
    #[getter]
    fn fill_percent(&self) -> f64 {
        self.inner.fill_percent()
    }

    /// Sets the fill level as a percentage (0-100).
    #[setter]
    fn set_fill_percent(&mut self, percent: f64) {
        self.inner.set_fill_percent(percent);
    }

    /// Returns the filled volume in m³.
    #[getter]
    fn fill_volume(&self) -> f64 {
        self.inner.fill_volume()
    }

    /// Returns the fluid mass in kg.
    #[getter]
    fn fluid_mass(&self) -> f64 {
        self.inner.fluid_mass()
    }

    /// Returns the center of gravity [x, y, z].
    #[getter]
    fn center_of_gravity(&self) -> [f64; 3] {
        self.inner.center_of_gravity()
    }

    /// Returns the transverse free surface moment in m⁴.
    #[getter]
    fn free_surface_moment_t(&self) -> f64 {
        self.inner.free_surface_moment_t()
    }

    /// Returns the longitudinal free surface moment in m⁴.
    #[getter]
    fn free_surface_moment_l(&self) -> f64 {
        self.inner.free_surface_moment_l()
    }

    fn __repr__(&self) -> String {
        format!(
            "Tank(name='{}', volume={:.2}m³, fill={:.1}%)",
            self.inner.name(),
            self.inner.total_volume(),
            self.inner.fill_percent()
        )
    }

    /// Returns tank container vertices [(x,y,z)].
    fn get_vertices(&self) -> Vec<(f64, f64, f64)> {
        self.inner
            .mesh()
            .vertices()
            .iter()
            .map(|v| (v.x, v.y, v.z))
            .collect()
    }

    /// Returns tank container faces [(i,j,k)].
    fn get_faces(&self) -> Vec<(u32, u32, u32)> {
        self.inner
            .mesh()
            .indices()
            .iter()
            .map(|idx| (idx[0], idx[1], idx[2]))
            .collect()
    }

    /// Returns fluid mesh vertices [(x,y,z)] or empty list if empty.
    /// If heel/trim not specified, assumes 0.
    #[pyo3(signature = (heel=0.0, trim=0.0))]
    fn get_fluid_vertices(&self, heel: f64, trim: f64) -> Vec<(f64, f64, f64)> {
        self.inner
            .get_fluid_mesh_at(heel, trim)
            .map(|m| m.vertices().iter().map(|v| (v.x, v.y, v.z)).collect())
            .unwrap_or_default()
    }

    /// Returns fluid mesh faces [(i,j,k)] or empty list if empty.
    /// If heel/trim not specified, assumes 0.
    #[pyo3(signature = (heel=0.0, trim=0.0))]
    fn get_fluid_faces(&self, heel: f64, trim: f64) -> Vec<(u32, u32, u32)> {
        self.inner
            .get_fluid_mesh_at(heel, trim)
            .map(|m| {
                m.indices()
                    .iter()
                    .map(|idx| (idx[0], idx[1], idx[2]))
                    .collect()
            })
            .unwrap_or_default()
    }
}

// ============================================================================
// Scripting Python Wrappers
// ============================================================================

/// Result of a single criterion check.
#[pyclass(name = "CriterionResult")]
#[derive(Clone)]
pub struct PyCriterionResult {
    #[pyo3(get)]
    pub name: String,
    #[pyo3(get)]
    pub description: String,
    #[pyo3(get)]
    pub required_value: f64,
    #[pyo3(get)]
    pub actual_value: f64,
    #[pyo3(get)]
    pub unit: String,
    #[pyo3(get)]
    pub status: String,
    #[pyo3(get)]
    pub margin: f64,
    #[pyo3(get)]
    pub notes: Option<String>,
    #[pyo3(get)]
    pub plot_id: Option<String>,
}

#[pymethods]
impl PyCriterionResult {
    fn __repr__(&self) -> String {
        format!("<CriterionResult '{}': {}>", self.name, self.status)
    }
}

/// Result of a criteria verification script.
#[pyclass(name = "CriteriaResult")]
#[derive(Clone)]
pub struct PyCriteriaResult {
    #[pyo3(get)]
    pub regulation_name: String,
    #[pyo3(get)]
    pub regulation_reference: String,
    #[pyo3(get)]
    pub vessel_name: String,
    #[pyo3(get)]
    pub loading_condition: String,
    #[pyo3(get)]
    pub displacement: f64,
    #[pyo3(get)]
    pub overall_pass: bool,
    #[pyo3(get)]
    pub pass_count: usize,
    #[pyo3(get)]
    pub fail_count: usize,
    #[pyo3(get)]
    pub notes: String,
    #[pyo3(get)]
    pub criteria: Vec<PyCriterionResult>,
    /// Plots as JSON strings (list of serialized PlotData)
    #[pyo3(get)]
    pub plots: Vec<String>,
}

#[pymethods]
impl PyCriteriaResult {
    fn __repr__(&self) -> String {
        format!(
            "<CriteriaResult '{}': {} (Passed {}/{})>",
            self.regulation_name,
            if self.overall_pass { "PASS" } else { "FAIL" },
            self.pass_count,
            self.criteria.len()
        )
    }
}

/// Context for Rhai scripts.
#[pyclass(name = "CriteriaContext")]
#[derive(Clone)]
pub struct PyCriteriaContext {
    inner: RustCriteriaContext,
}

#[pymethods]
impl PyCriteriaContext {
    /// Create a context from a CompleteStabilityResult.
    #[staticmethod]
    fn from_result(
        result: &PyCompleteStabilityResult,
        vessel_name: String,
        loading_condition: String,
    ) -> Self {
        Self {
            inner: RustCriteriaContext::new(result.inner.clone(), vessel_name, loading_condition),
        }
    }

    /// Get first flooding angle (or None).
    fn get_first_flooding_angle(&self) -> Option<f64> {
        self.inner.get_first_flooding_angle().try_cast::<f64>()
    }

    /// Find equilibrium angle for a given heeling arm.
    fn find_equilibrium_angle(&self, heeling_arm: f64) -> Option<f64> {
        self.inner
            .find_equilibrium_angle(heeling_arm)
            .try_cast::<f64>()
    }

    /// Find second intercept angle for a given heeling arm.
    fn find_second_intercept(&self, heeling_arm: f64) -> Option<f64> {
        self.inner
            .find_second_intercept(heeling_arm)
            .try_cast::<f64>()
    }

    /// Set a parameter for the script.
    /// Supports string, float, and bool.
    fn set_param(&mut self, key: &str, value: &Bound<PyAny>) -> PyResult<()> {
        let val = if let Ok(s) = value.extract::<String>() {
            rhai::Dynamic::from(s)
        } else if let Ok(f) = value.extract::<f64>() {
            rhai::Dynamic::from(f)
        } else if let Ok(b) = value.extract::<bool>() {
            rhai::Dynamic::from(b)
        } else {
            return Err(PyValueError::new_err(
                "Unsupported parameter type. Use str, float, or bool.",
            ));
        };
        self.inner.set_param(key, val);
        Ok(())
    }
}

/// Script execution engine.
#[pyclass(name = "ScriptEngine")]
pub struct PyScriptEngine {
    inner: RustScriptEngine,
}

#[pymethods]
impl PyScriptEngine {
    #[new]
    fn new() -> Self {
        Self {
            inner: RustScriptEngine::new(),
        }
    }

    /// Run a script from a file path.
    fn run_script_file(
        &self,
        path: &str,
        context: &PyCriteriaContext,
    ) -> PyResult<PyCriteriaResult> {
        let result = self
            .inner
            .run_script_file(path, context.inner.clone())
            .map_err(|e| PyValueError::new_err(format!("Script error: {}", e)))?;

        Ok(map_result_to_py(result))
    }

    /// Run a script from a string.
    fn run_script(&self, script: &str, context: &PyCriteriaContext) -> PyResult<PyCriteriaResult> {
        let result = self
            .inner
            .run_script(script, context.inner.clone())
            .map_err(|e| PyValueError::new_err(format!("Script error: {}", e)))?;

        Ok(map_result_to_py(result))
    }
}

fn map_result_to_py(res: RustCriteriaResult) -> PyCriteriaResult {
    let criteria: Vec<PyCriterionResult> = res
        .criteria
        .iter()
        .map(|c| PyCriterionResult {
            name: c.name.clone(),
            description: c.description.clone(),
            required_value: c.required_value,
            actual_value: c.actual_value,
            unit: c.unit.clone(),
            status: c.status.to_string(),
            margin: c.margin,
            notes: c.notes.clone(),
            plot_id: c.plot_id.clone(),
        })
        .collect();

    // Serialize plots to JSON strings
    let plots: Vec<String> = res
        .plots
        .iter()
        .map(|p| serde_json::to_string(p).unwrap_or_default())
        .collect();

    PyCriteriaResult {
        regulation_name: res.regulation_name,
        regulation_reference: res.regulation_reference,
        vessel_name: res.vessel_name,
        loading_condition: res.loading_condition,
        displacement: res.displacement,
        overall_pass: res.overall_pass,
        pass_count: res.pass_count,
        fail_count: res.fail_count,
        notes: res.notes,
        criteria,
        plots,
    }
}

// ============================================================================
// Python Module Definition
// ============================================================================

/// Naval architecture library for hydrostatics, stability, and tank calculations.
#[pymodule]
fn navaltoolbox(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyHull>()?;
    m.add_class::<PyVessel>()?;
    m.add_class::<PySilhouette>()?;
    m.add_class::<PyOpeningType>()?;
    m.add_class::<PyDownfloodingOpening>()?;
    m.add_class::<PyHydrostaticState>()?;
    m.add_class::<PyHydrostaticsCalculator>()?;
    m.add_class::<PyStabilityPoint>()?;
    m.add_class::<PyStabilityCurve>()?;
    m.add_class::<PyWindHeelingData>()?;
    m.add_class::<PyCompleteStabilityResult>()?;
    m.add_class::<PyStabilityCalculator>()?;
    m.add_class::<PyTank>()?;

    // Scripting
    m.add_class::<PyCriterionResult>()?;
    m.add_class::<PyCriteriaResult>()?;
    m.add_class::<PyCriteriaContext>()?;
    m.add_class::<PyScriptEngine>()?;

    Ok(())
}
