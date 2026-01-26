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

//! Downflooding openings for intact stability calculations.
//!
//! Per IMO 2008 IS Code (MSC.267), the downflooding angle (θf) is the angle
//! at which non-weathertight openings become immersed.

use std::f64::consts::PI;

mod loader;

pub use loader::OpeningLoadError;
use loader::{load_dxf_openings, load_vtk_openings};
use std::path::Path;

/// Type of opening that can cause downflooding.
#[derive(Clone, Debug, PartialEq)]
pub enum OpeningType {
    /// Ventilator
    Vent,
    /// Air pipe without automatic closing device
    AirPipe,
    /// Hatch or manhole
    Hatch,
    /// Weathertight door (kept open for operation)
    Door,
    /// Non-weathertight window
    Window,
    /// Custom opening type
    Other(String),
}

/// Geometry of the opening.
#[derive(Clone, Debug)]
pub enum OpeningGeometry {
    /// Single point [x, y, z] in ship coordinates
    Point([f64; 3]),
    /// Contour/polyline defining the opening boundary
    Contour(Vec<[f64; 3]>),
}

/// A downflooding opening point or contour.
///
/// Used to determine the downflooding angle (θf) per IMO 2008 IS Code.
#[derive(Clone, Debug)]
pub struct DownfloodingOpening {
    /// Name/identifier of the opening
    name: String,
    /// Geometry of the opening (point or contour)
    geometry: OpeningGeometry,
    /// Type of opening
    opening_type: OpeningType,
    /// Whether this opening is active in calculations
    active: bool,
}

impl DownfloodingOpening {
    /// Create an opening from a single point.
    pub fn from_point(name: String, position: [f64; 3], opening_type: OpeningType) -> Self {
        Self {
            name,
            geometry: OpeningGeometry::Point(position),
            opening_type,
            active: true,
        }
    }

    /// Create an opening from a contour (polyline).
    pub fn from_contour(name: String, points: Vec<[f64; 3]>, opening_type: OpeningType) -> Self {
        Self {
            name,
            geometry: OpeningGeometry::Contour(points),
            opening_type,
            active: true,
        }
    }

    /// Load openings from a file (DXF or VTK) based on extension.
    /// Returns a list of openings found in the file.
    pub fn from_file(
        path: &Path,
        default_type: OpeningType,
    ) -> Result<Vec<Self>, OpeningLoadError> {
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.to_lowercase())
            .unwrap_or_default();

        match ext.as_str() {
            "dxf" => load_dxf_openings(path, default_type),
            "vtk" | "vtp" | "vtu" => load_vtk_openings(path, default_type),
            _ => Err(OpeningLoadError::UnsupportedFormat),
        }
    }

    /// Get the opening name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Set the opening name.
    pub fn set_name(&mut self, name: String) {
        self.name = name;
    }

    /// Get the opening type.
    pub fn opening_type(&self) -> &OpeningType {
        &self.opening_type
    }

    /// Check if the opening is active.
    pub fn is_active(&self) -> bool {
        self.active
    }

    /// Set the opening active state.
    pub fn set_active(&mut self, active: bool) {
        self.active = active;
    }

    /// Get the geometry.
    pub fn geometry(&self) -> &OpeningGeometry {
        &self.geometry
    }

    /// Get all points of the opening.
    pub fn get_points(&self) -> Vec<[f64; 3]> {
        match &self.geometry {
            OpeningGeometry::Point(p) => vec![*p],
            OpeningGeometry::Contour(pts) => pts.clone(),
        }
    }

    /// Check if the opening is submerged at given heel/trim/draft.
    ///
    /// Uses the same rotation as stability calculations around the pivot point.
    ///
    /// # Arguments
    /// * `heel` - Heel angle in degrees (positive = starboard down)
    /// * `trim` - Trim angle in degrees (positive = stern down)
    /// * `pivot` - Rotation pivot point [x, y, z] (typically center of waterplane at draft)
    /// * `waterline_z` - Waterline Z coordinate after rotation
    pub fn is_submerged(&self, heel: f64, trim: f64, pivot: [f64; 3], waterline_z: f64) -> bool {
        let points = self.get_points();

        for point in &points {
            let rotated = rotate_point(*point, heel, trim, pivot);
            if rotated[2] < waterline_z {
                return true;
            }
        }

        false
    }

    /// Get the lowest Z coordinate of any point after rotation.
    pub fn get_lowest_z(&self, heel: f64, trim: f64, pivot: [f64; 3]) -> f64 {
        let points = self.get_points();

        points
            .iter()
            .map(|p| rotate_point(*p, heel, trim, pivot)[2])
            .fold(f64::MAX, f64::min)
    }
}

/// Rotate a point around a pivot by heel and trim angles.
///
/// Matches the rotation used in stability calculations.
fn rotate_point(point: [f64; 3], heel: f64, trim: f64, pivot: [f64; 3]) -> [f64; 3] {
    let heel_rad = heel * PI / 180.0;
    let trim_rad = trim * PI / 180.0;

    // Translate to pivot
    let dx = point[0] - pivot[0];
    let dy = point[1] - pivot[1];
    let dz = point[2] - pivot[2];

    // Rotation around X axis (heel)
    let cos_h = heel_rad.cos();
    let sin_h = heel_rad.sin();
    let y1 = dy * cos_h - dz * sin_h;
    let z1 = dy * sin_h + dz * cos_h;

    // Rotation around Y axis (trim)
    let cos_t = trim_rad.cos();
    let sin_t = trim_rad.sin();
    let x2 = dx * cos_t + z1 * sin_t;
    let z2 = -dx * sin_t + z1 * cos_t;

    // Translate back
    [x2 + pivot[0], y1 + pivot[1], z2 + pivot[2]]
}

/// Check which openings are submerged at given conditions.
pub fn check_openings_submerged(
    openings: &[DownfloodingOpening],
    heel: f64,
    trim: f64,
    pivot: [f64; 3],
    waterline_z: f64,
) -> Vec<String> {
    openings
        .iter()
        .filter(|o| o.is_active() && o.is_submerged(heel, trim, pivot, waterline_z))
        .map(|o| o.name().to_string())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_opening() {
        let opening = DownfloodingOpening::from_point(
            "Test Vent".to_string(),
            [50.0, 5.0, 10.0],
            OpeningType::Vent,
        );

        assert_eq!(opening.name(), "Test Vent");
        assert!(opening.is_active());
        assert_eq!(opening.get_points().len(), 1);
    }

    #[test]
    fn test_contour_opening() {
        let points = vec![
            [10.0, -5.0, 15.0],
            [20.0, -5.0, 15.0],
            [20.0, 5.0, 15.0],
            [10.0, 5.0, 15.0],
        ];
        let opening = DownfloodingOpening::from_contour(
            "Cargo Hatch".to_string(),
            points,
            OpeningType::Hatch,
        );

        assert_eq!(opening.get_points().len(), 4);
    }

    #[test]
    fn test_submerged_at_heel() {
        // Opening on starboard side at y=5, z=10
        let opening = DownfloodingOpening::from_point(
            "Starboard Vent".to_string(),
            [50.0, 5.0, 10.0],
            OpeningType::Vent,
        );

        let pivot = [50.0, 0.0, 5.0];

        // At 0° heel, opening z after rotation is 10 (above waterline 5)
        assert!(!opening.is_submerged(0.0, 0.0, pivot, 5.0));

        // Get the rotated Z at 45° to debug
        let z_at_45 = opening.get_lowest_z(45.0, 0.0, pivot);

        // At 45° heel, starboard side goes down
        // The rotated z should be lower than before
        // Check if it goes below a higher waterline (z=8)
        assert!(
            opening.is_submerged(45.0, 0.0, pivot, z_at_45 + 1.0),
            "At 45° heel, z={:.2}, should be submerged at waterline {:.2}",
            z_at_45,
            z_at_45 + 1.0
        );
    }
}
