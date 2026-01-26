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

//! Silhouette struct for wind heeling moment calculations.

use super::loader::{load_dxf_silhouette, load_vtk_silhouette, SilhouetteLoadError};
use std::path::Path;

/// A 2D silhouette profile in the X-Z plane (ship lateral view).
///
/// Used for calculating wind heeling moments per IMO 2008 IS Code.
#[derive(Clone, Debug)]
pub struct Silhouette {
    /// Points defining the silhouette contour [x, y, z].
    /// Y is always 0 (X-Z plane). Stored as 3D for future extensibility.
    points: Vec<[f64; 3]>,
    /// Name of the silhouette (from filename or user-defined).
    name: String,
}

impl Silhouette {
    /// Create a new silhouette from a list of 3D points.
    pub fn new(points: Vec<[f64; 3]>, name: String) -> Self {
        Self { points, name }
    }

    /// Load a silhouette from a DXF file.
    ///
    /// Extracts the first LWPOLYLINE or POLYLINE entity.
    /// If Y coordinates are non-zero, they are set to 0 with a warning.
    pub fn from_dxf(path: &Path) -> Result<Self, SilhouetteLoadError> {
        let (points, name) = load_dxf_silhouette(path)?;
        Ok(Self { points, name })
    }

    /// Load a silhouette from a VTK file (.vtk or .vtp).
    ///
    /// Extracts the first polyline from the PolyData.
    /// If Y coordinates are non-zero, they are set to 0 with a warning.
    pub fn from_vtk(path: &Path) -> Result<Self, SilhouetteLoadError> {
        let (points, name) = load_vtk_silhouette(path)?;
        Ok(Self { points, name })
    }

    /// Load a silhouette from a file (DXF or VTK) based on extension.
    pub fn from_file(path: &Path) -> Result<Self, SilhouetteLoadError> {
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.to_lowercase())
            .unwrap_or_default();

        match ext.as_str() {
            "dxf" => Self::from_dxf(path),
            "vtk" | "vtp" | "vtu" => Self::from_vtk(path),
            _ => Err(SilhouetteLoadError::UnsupportedFormat),
        }
    }

    /// Get the silhouette name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Get the points defining the silhouette.
    pub fn points(&self) -> &[[f64; 3]] {
        &self.points
    }

    /// Get number of points.
    pub fn num_points(&self) -> usize {
        self.points.len()
    }

    /// Check if the contour is closed (first point == last point).
    pub fn is_closed(&self) -> bool {
        if self.points.len() < 2 {
            return false;
        }
        let first = self.points.first().unwrap();
        let last = self.points.last().unwrap();
        let eps = 1e-6;
        (first[0] - last[0]).abs() < eps
            && (first[1] - last[1]).abs() < eps
            && (first[2] - last[2]).abs() < eps
    }

    /// Calculate the total lateral area of the silhouette (m²).
    ///
    /// Uses the shoelace formula in the X-Z plane.
    pub fn get_area(&self) -> f64 {
        if self.points.len() < 3 {
            return 0.0;
        }

        let mut area = 0.0;
        let n = self.points.len();

        for i in 0..n {
            let j = (i + 1) % n;
            let xi = self.points[i][0];
            let zi = self.points[i][2];
            let xj = self.points[j][0];
            let zj = self.points[j][2];
            area += xi * zj - xj * zi;
        }

        (area / 2.0).abs()
    }

    /// Calculate the centroid of the silhouette in the X-Z plane.
    ///
    /// Returns [x_center, z_center].
    pub fn get_centroid(&self) -> [f64; 2] {
        if self.points.len() < 3 {
            return [0.0, 0.0];
        }

        let area = self.get_area();
        if area < 1e-9 {
            return [0.0, 0.0];
        }

        let mut cx = 0.0;
        let mut cz = 0.0;
        let n = self.points.len();

        for i in 0..n {
            let j = (i + 1) % n;
            let xi = self.points[i][0];
            let zi = self.points[i][2];
            let xj = self.points[j][0];
            let zj = self.points[j][2];
            let cross = xi * zj - xj * zi;
            cx += (xi + xj) * cross;
            cz += (zi + zj) * cross;
        }

        let factor = 1.0 / (6.0 * area);
        [cx.abs() * factor, cz.abs() * factor]
    }

    /// Get the bounding box of the silhouette.
    ///
    /// Returns (x_min, x_max, z_min, z_max).
    pub fn get_bounds(&self) -> (f64, f64, f64, f64) {
        if self.points.is_empty() {
            return (0.0, 0.0, 0.0, 0.0);
        }

        let mut x_min = f64::MAX;
        let mut x_max = f64::MIN;
        let mut z_min = f64::MAX;
        let mut z_max = f64::MIN;

        for p in &self.points {
            x_min = x_min.min(p[0]);
            x_max = x_max.max(p[0]);
            z_min = z_min.min(p[2]);
            z_max = z_max.max(p[2]);
        }

        (x_min, x_max, z_min, z_max)
    }

    /// Calculate the emerged area (above waterline) in m².
    ///
    /// Clips the silhouette at the given waterline Z and returns
    /// the area of the portion above the waterline.
    pub fn get_emerged_area(&self, waterline_z: f64) -> f64 {
        let clipped = self.clip_above(waterline_z);
        Self::polygon_area(&clipped)
    }

    /// Calculate the centroid of the emerged area.
    ///
    /// Returns [x_center, z_center] of the area above waterline.
    pub fn get_emerged_centroid(&self, waterline_z: f64) -> [f64; 2] {
        let clipped = self.clip_above(waterline_z);
        Self::polygon_centroid(&clipped)
    }

    /// Clip the silhouette to keep only points above the waterline.
    fn clip_above(&self, waterline_z: f64) -> Vec<[f64; 2]> {
        if self.points.len() < 2 {
            return Vec::new();
        }

        let mut result: Vec<[f64; 2]> = Vec::new();
        let n = self.points.len();

        for i in 0..n {
            let j = (i + 1) % n;
            let p1 = [self.points[i][0], self.points[i][2]];
            let p2 = [self.points[j][0], self.points[j][2]];

            let z1 = p1[1];
            let z2 = p2[1];

            // Both above
            if z1 >= waterline_z && z2 >= waterline_z {
                result.push(p1);
            }
            // p1 above, p2 below - add p1 and intersection
            else if z1 >= waterline_z && z2 < waterline_z {
                result.push(p1);
                let t = (waterline_z - z1) / (z2 - z1);
                let x_int = p1[0] + t * (p2[0] - p1[0]);
                result.push([x_int, waterline_z]);
            }
            // p1 below, p2 above - add intersection
            else if z1 < waterline_z && z2 >= waterline_z {
                let t = (waterline_z - z1) / (z2 - z1);
                let x_int = p1[0] + t * (p2[0] - p1[0]);
                result.push([x_int, waterline_z]);
            }
            // Both below - skip
        }

        result
    }

    /// Calculate area of a 2D polygon using shoelace formula.
    fn polygon_area(points: &[[f64; 2]]) -> f64 {
        if points.len() < 3 {
            return 0.0;
        }

        let mut area = 0.0;
        let n = points.len();

        for i in 0..n {
            let j = (i + 1) % n;
            area += points[i][0] * points[j][1] - points[j][0] * points[i][1];
        }

        (area / 2.0).abs()
    }

    /// Calculate centroid of a 2D polygon.
    fn polygon_centroid(points: &[[f64; 2]]) -> [f64; 2] {
        if points.len() < 3 {
            return [0.0, 0.0];
        }

        let area = Self::polygon_area(points);
        if area < 1e-9 {
            return [0.0, 0.0];
        }

        let mut cx = 0.0;
        let mut cy = 0.0;
        let n = points.len();

        for i in 0..n {
            let j = (i + 1) % n;
            let cross = points[i][0] * points[j][1] - points[j][0] * points[i][1];
            cx += (points[i][0] + points[j][0]) * cross;
            cy += (points[i][1] + points[j][1]) * cross;
        }

        let factor = 1.0 / (6.0 * area);
        [cx.abs() * factor, cy.abs() * factor]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_rectangle(x_min: f64, x_max: f64, z_min: f64, z_max: f64) -> Silhouette {
        let points = vec![
            [x_min, 0.0, z_min],
            [x_max, 0.0, z_min],
            [x_max, 0.0, z_max],
            [x_min, 0.0, z_max],
            [x_min, 0.0, z_min], // Close the polygon
        ];
        Silhouette::new(points, "rectangle".to_string())
    }

    #[test]
    fn test_rectangle_is_closed() {
        let rect = create_rectangle(0.0, 100.0, 0.0, 20.0);
        assert!(rect.is_closed());
    }

    #[test]
    fn test_rectangle_area() {
        let rect = create_rectangle(0.0, 100.0, 0.0, 20.0);
        let area = rect.get_area();
        assert!(
            (area - 2000.0).abs() < 1.0,
            "Area should be 2000 m², got {}",
            area
        );
    }

    #[test]
    fn test_emerged_area_full() {
        let rect = create_rectangle(0.0, 100.0, 0.0, 20.0);
        let emerged = rect.get_emerged_area(0.0);
        assert!(
            (emerged - 2000.0).abs() < 1.0,
            "Full emerged area should be 2000 m²"
        );
    }

    #[test]
    fn test_emerged_area_half() {
        let rect = create_rectangle(0.0, 100.0, 0.0, 20.0);
        let emerged = rect.get_emerged_area(10.0);
        assert!(
            (emerged - 1000.0).abs() < 10.0,
            "Half emerged area should be ~1000 m², got {}",
            emerged
        );
    }

    #[test]
    fn test_emerged_area_none() {
        let rect = create_rectangle(0.0, 100.0, 0.0, 20.0);
        let emerged = rect.get_emerged_area(25.0);
        assert!(emerged < 1.0, "No emerged area above z=25, got {}", emerged);
    }
}
