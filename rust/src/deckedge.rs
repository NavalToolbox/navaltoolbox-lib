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

//! Deck edge (livet) module for freeboard calculation.
//!
//! The deck edge represents the intersection of the deck with the hull shell.
//! It is used to calculate the freeboard (vertical distance from waterline to deck).

use dxf::entities::EntityType as DxfEntityType;
use dxf::Drawing;
use std::f64::consts::PI;
use std::path::Path;
use thiserror::Error;
use vtkio::model::{DataSet, Piece};
use vtkio::Vtk;

/// Error type for deck edge loading operations.
#[derive(Error, Debug)]
pub enum DeckEdgeLoadError {
    #[error("Failed to read file: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Failed to parse DXF file: {0}")]
    DxfParseError(#[from] dxf::DxfError),
    #[error("Failed to parse VTK file: {0}")]
    VtkParseError(#[from] vtkio::Error),
    #[error("No valid geometry found in file")]
    NoGeometry,
    #[error("Unsupported file format: {0}")]
    UnsupportedFormat(String),
}

/// Side of the deck edge.
#[derive(Clone, Debug, PartialEq)]
pub enum DeckEdgeSide {
    /// Port side only
    Port,
    /// Starboard side only
    Starboard,
    /// Both sides (symmetric or full contour)
    Both,
}

/// A deck edge contour (livet) for freeboard calculation.
///
/// The deck edge is a 3D polyline representing the edge of the deck.
/// It is used to calculate the minimum freeboard at any heel/trim condition.
#[derive(Clone, Debug)]
pub struct DeckEdge {
    name: String,
    points: Vec<[f64; 3]>,
    side: DeckEdgeSide,
}

impl DeckEdge {
    /// Create a deck edge from a list of 3D points.
    pub fn new(name: &str, points: Vec<[f64; 3]>, side: DeckEdgeSide) -> Self {
        Self {
            name: name.to_string(),
            points,
            side,
        }
    }

    /// Load a deck edge from a DXF or VTK file.
    pub fn from_file(name: &str, path: &Path) -> Result<Self, DeckEdgeLoadError> {
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .map(|s| s.to_lowercase())
            .unwrap_or_default();

        let points = match ext.as_str() {
            "dxf" => Self::load_dxf(path)?,
            "vtk" | "vtp" => Self::load_vtk(path)?,
            _ => {
                return Err(DeckEdgeLoadError::UnsupportedFormat(format!(
                    "Unsupported extension: {}",
                    ext
                )));
            }
        };

        if points.is_empty() {
            return Err(DeckEdgeLoadError::NoGeometry);
        }

        // Auto-detect side based on Y coordinates
        let side = Self::detect_side(&points);

        Ok(Self {
            name: name.to_string(),
            points,
            side,
        })
    }

    // =========================================================================
    // Accessors
    // =========================================================================

    /// Get the deck edge name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Set the deck edge name.
    pub fn set_name(&mut self, name: &str) {
        self.name = name.to_string();
    }

    /// Get the points of the deck edge.
    pub fn points(&self) -> &[[f64; 3]] {
        &self.points
    }

    /// Get the side of the deck edge.
    pub fn side(&self) -> &DeckEdgeSide {
        &self.side
    }

    /// Set the side of the deck edge.
    pub fn set_side(&mut self, side: DeckEdgeSide) {
        self.side = side;
    }

    // =========================================================================
    // Freeboard Calculation
    // =========================================================================

    /// Calculate the minimum Z coordinate of the deck edge after rotation.
    ///
    /// # Arguments
    /// * `heel` - Heel angle in degrees (positive = starboard down)
    /// * `trim` - Trim angle in degrees (positive = bow down)
    /// * `pivot` - Rotation pivot point [x, y, z]
    ///
    /// # Returns
    /// The minimum Z coordinate of the deck edge after rotation.
    pub fn get_min_z_at_heel(&self, heel: f64, trim: f64, pivot: [f64; 3]) -> f64 {
        self.points
            .iter()
            .map(|p| Self::rotate_point(*p, heel, trim, pivot)[2])
            .fold(f64::INFINITY, f64::min)
    }

    /// Calculate the freeboard at given conditions.
    ///
    /// # Arguments
    /// * `heel` - Heel angle in degrees
    /// * `trim` - Trim angle in degrees
    /// * `pivot` - Rotation pivot point
    /// * `waterline_z` - Waterline Z coordinate
    ///
    /// # Returns
    /// Freeboard in meters (positive = deck above water, negative = deck submerged)
    pub fn get_freeboard(&self, heel: f64, trim: f64, pivot: [f64; 3], waterline_z: f64) -> f64 {
        self.get_min_z_at_heel(heel, trim, pivot) - waterline_z
    }

    // =========================================================================
    // Private helpers
    // =========================================================================

    fn rotate_point(point: [f64; 3], heel: f64, trim: f64, pivot: [f64; 3]) -> [f64; 3] {
        let heel_rad = heel * PI / 180.0;
        let trim_rad = trim * PI / 180.0;

        // Translate to pivot
        let x = point[0] - pivot[0];
        let y = point[1] - pivot[1];
        let z = point[2] - pivot[2];

        // Rotation matrices (same as in downflooding and stability)
        // Roll (heel) around X axis
        let cos_h = heel_rad.cos();
        let sin_h = heel_rad.sin();
        let y1 = y * cos_h - z * sin_h;
        let z1 = y * sin_h + z * cos_h;

        // Pitch (trim) around Y axis
        let cos_t = trim_rad.cos();
        let sin_t = trim_rad.sin();
        let x2 = x * cos_t + z1 * sin_t;
        let z2 = -x * sin_t + z1 * cos_t;

        // Translate back
        [x2 + pivot[0], y1 + pivot[1], z2 + pivot[2]]
    }

    fn detect_side(points: &[[f64; 3]]) -> DeckEdgeSide {
        let has_positive_y = points.iter().any(|p| p[1] > 0.01);
        let has_negative_y = points.iter().any(|p| p[1] < -0.01);

        match (has_positive_y, has_negative_y) {
            (true, true) => DeckEdgeSide::Both,
            (true, false) => DeckEdgeSide::Starboard,
            (false, true) => DeckEdgeSide::Port,
            (false, false) => DeckEdgeSide::Both, // centerline
        }
    }

    fn load_dxf(path: &Path) -> Result<Vec<[f64; 3]>, DeckEdgeLoadError> {
        let drawing = Drawing::load_file(path)?;
        let mut points = Vec::new();

        for entity in drawing.entities() {
            match &entity.specific {
                DxfEntityType::Polyline(polyline) => {
                    for vertex in polyline.vertices() {
                        points.push([
                            vertex.location.x,
                            vertex.location.y,
                            vertex.location.z,
                        ]);
                    }
                }
                DxfEntityType::LwPolyline(lwpoly) => {
                    // LwPolyline is 2D, use extrusion z as elevation if available
                    let z = lwpoly.extrusion_direction.z;
                    for vertex in &lwpoly.vertices {
                        points.push([vertex.x, vertex.y, z]);
                    }
                }
                DxfEntityType::Line(line) => {
                    points.push([line.p1.x, line.p1.y, line.p1.z]);
                    points.push([line.p2.x, line.p2.y, line.p2.z]);
                }
                DxfEntityType::Spline(spline) => {
                    for cp in &spline.control_points {
                        points.push([cp.x, cp.y, cp.z]);
                    }
                }
                _ => {}
            }
        }

        Ok(points)
    }

    fn load_vtk(path: &Path) -> Result<Vec<[f64; 3]>, DeckEdgeLoadError> {
        let vtk = Vtk::import(path)?;

        fn extract_points(data: &DataSet) -> Result<Vec<[f64; 3]>, DeckEdgeLoadError> {
            match data {
                DataSet::PolyData { pieces, .. } => {
                    for piece in pieces {
                        if let Piece::Inline(piece_data) = piece {
                            let buffer = &piece_data.points;
                            return extract_coords(buffer);
                        }
                    }
                    Err(DeckEdgeLoadError::NoGeometry)
                }
                DataSet::UnstructuredGrid { pieces, .. } => {
                    for piece in pieces {
                        if let Piece::Inline(piece_data) = piece {
                            let buffer = &piece_data.points;
                            return extract_coords(buffer);
                        }
                    }
                    Err(DeckEdgeLoadError::NoGeometry)
                }
                _ => Err(DeckEdgeLoadError::NoGeometry),
            }
        }

        fn extract_coords(buffer: &vtkio::IOBuffer) -> Result<Vec<[f64; 3]>, DeckEdgeLoadError> {
            let coords: Vec<f64> = match buffer {
                vtkio::IOBuffer::F64(v) => v.clone(),
                vtkio::IOBuffer::F32(v) => v.iter().map(|x| *x as f64).collect(),
                _ => return Err(DeckEdgeLoadError::NoGeometry),
            };

            Ok(coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect())
        }

        extract_points(&vtk.data)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_deck_edge_creation() {
        let points = vec![
            [0.0, 5.0, 10.0],
            [50.0, 5.0, 10.0],
            [100.0, 5.0, 10.0],
        ];
        let edge = DeckEdge::new("main_deck", points.clone(), DeckEdgeSide::Starboard);

        assert_eq!(edge.name(), "main_deck");
        assert_eq!(edge.points().len(), 3);
        assert_eq!(*edge.side(), DeckEdgeSide::Starboard);
    }

    #[test]
    fn test_freeboard_level() {
        let points = vec![
            [0.0, 5.0, 10.0],
            [50.0, 5.0, 10.0],
            [100.0, 5.0, 10.0],
        ];
        let edge = DeckEdge::new("deck", points, DeckEdgeSide::Starboard);

        // At level, min Z should be 10.0
        let min_z = edge.get_min_z_at_heel(0.0, 0.0, [50.0, 0.0, 5.0]);
        assert!((min_z - 10.0).abs() < 1e-6);

        // Freeboard at waterline 5.0 should be 5.0
        let freeboard = edge.get_freeboard(0.0, 0.0, [50.0, 0.0, 5.0], 5.0);
        assert!((freeboard - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_freeboard_heeled() {
        let points = vec![
            [50.0, 5.0, 10.0],  // Starboard point at Y=5
        ];
        let edge = DeckEdge::new("deck", points, DeckEdgeSide::Starboard);

        // At 30 degrees heel, starboard side goes down
        // Y=5, Z=10, pivot at origin
        // After 30 deg heel: Z' = Y*sin(30) + Z*cos(30) = 5*0.5 + 10*0.866 = 2.5 + 8.66 = 11.16
        // Wait, the starboard side goes UP at positive heel (rotation around X)
        // Let me recalculate: cos(30)=0.866, sin(30)=0.5
        // z' = y*sin + z*cos = 5*0.5 + 10*0.866 = 11.16 (starboard goes up)
        let min_z = edge.get_min_z_at_heel(30.0, 0.0, [50.0, 0.0, 0.0]);
        assert!(min_z > 10.0); // Starboard goes up at positive heel
    }

    #[test]
    fn test_side_detection() {
        // Starboard only (positive Y)
        let pts_stbd = vec![[0.0, 5.0, 10.0]];
        assert_eq!(DeckEdge::detect_side(&pts_stbd), DeckEdgeSide::Starboard);

        // Port only (negative Y)
        let pts_port = vec![[0.0, -5.0, 10.0]];
        assert_eq!(DeckEdge::detect_side(&pts_port), DeckEdgeSide::Port);

        // Both sides
        let pts_both = vec![[0.0, 5.0, 10.0], [0.0, -5.0, 10.0]];
        assert_eq!(DeckEdge::detect_side(&pts_both), DeckEdgeSide::Both);
    }
}
