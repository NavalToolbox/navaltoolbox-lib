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

//! Appendage module for additional volume elements.
//!
//! Appendages represent additional volume elements attached to the hull,
//! such as keels, rudders, propellers, etc.

use parry3d_f64::shape::TriMesh;
use std::path::Path;
use thiserror::Error;

/// Error type for appendage loading operations.
#[derive(Error, Debug)]
pub enum AppendageLoadError {
    #[error("Failed to read file: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Failed to parse STL file: {0}")]
    StlParseError(String),
    #[error("Failed to parse VTK file: {0}")]
    VtkParseError(#[from] vtkio::Error),
    #[error("Unsupported file format: {0}")]
    UnsupportedFormat(String),
}

/// Geometry type for an appendage.
#[derive(Clone, Debug)]
pub enum AppendageGeometry {
    /// Point appendage with fixed volume at a position
    Point {
        center: [f64; 3],
        volume: f64,
    },
    /// Mesh appendage loaded from STL/VTK
    Mesh(TriMesh),
    /// Box appendage defined by bounds
    Box {
        bounds: (f64, f64, f64, f64, f64, f64), // xmin, xmax, ymin, ymax, zmin, zmax
    },
    /// Sphere appendage defined by center and volume
    Sphere {
        center: [f64; 3],
        volume: f64,
    },
    /// Cube appendage defined by center and volume
    Cube {
        center: [f64; 3],
        volume: f64,
    },
}

/// An appendage element attached to the vessel.
///
/// Appendages add volume and wetted surface to hydrostatic calculations.
#[derive(Clone, Debug)]
pub struct Appendage {
    name: String,
    geometry: AppendageGeometry,
    wetted_surface: Option<f64>,
}

impl Appendage {
    /// Create an appendage from a point (fixed volume at a position).
    pub fn from_point(name: &str, center: [f64; 3], volume: f64) -> Self {
        Self {
            name: name.to_string(),
            geometry: AppendageGeometry::Point { center, volume },
            wetted_surface: None,
        }
    }

    /// Create an appendage from an STL or VTK file.
    pub fn from_file(name: &str, path: &Path) -> Result<Self, AppendageLoadError> {
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .map(|s| s.to_lowercase())
            .unwrap_or_default();

        let mesh = match ext.as_str() {
            "stl" => {
                let file = std::fs::File::open(path)?;
                let mut reader = std::io::BufReader::new(file);
                let stl = stl_io::read_stl(&mut reader)
                    .map_err(|e| AppendageLoadError::StlParseError(e.to_string()))?;

                let vertices: Vec<nalgebra::Point3<f64>> = stl
                    .vertices
                    .iter()
                    .map(|v| nalgebra::Point3::new(v[0] as f64, v[1] as f64, v[2] as f64))
                    .collect();

                let indices: Vec<[u32; 3]> = stl
                    .faces
                    .iter()
                    .map(|f| [f.vertices[0] as u32, f.vertices[1] as u32, f.vertices[2] as u32])
                    .collect();

                TriMesh::new(vertices, indices)
                    .map_err(|e| AppendageLoadError::StlParseError(format!("{:?}", e)))?
            }
            "vtk" | "vtp" => {
                let vtk = vtkio::Vtk::import(path)?;
                Self::vtk_to_mesh(&vtk)?
            }
            _ => {
                return Err(AppendageLoadError::UnsupportedFormat(format!(
                    "Unsupported extension: {}",
                    ext
                )));
            }
        };

        Ok(Self {
            name: name.to_string(),
            geometry: AppendageGeometry::Mesh(mesh),
            wetted_surface: None,
        })
    }

    /// Create an appendage from a box (parallelepiped) defined by bounds.
    pub fn from_box(name: &str, bounds: (f64, f64, f64, f64, f64, f64)) -> Self {
        Self {
            name: name.to_string(),
            geometry: AppendageGeometry::Box { bounds },
            wetted_surface: None,
        }
    }

    /// Create an appendage from a cube defined by center and volume.
    pub fn from_cube(name: &str, center: [f64; 3], volume: f64) -> Self {
        Self {
            name: name.to_string(),
            geometry: AppendageGeometry::Cube { center, volume },
            wetted_surface: None,
        }
    }

    /// Create an appendage from a sphere defined by center and volume.
    pub fn from_sphere(name: &str, center: [f64; 3], volume: f64) -> Self {
        Self {
            name: name.to_string(),
            geometry: AppendageGeometry::Sphere { center, volume },
            wetted_surface: None,
        }
    }

    // =========================================================================
    // Accessors
    // =========================================================================

    /// Get the appendage name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Set the appendage name.
    pub fn set_name(&mut self, name: &str) {
        self.name = name.to_string();
    }

    /// Get the geometry type.
    pub fn geometry(&self) -> &AppendageGeometry {
        &self.geometry
    }

    /// Get the wetted surface area (if specified).
    pub fn wetted_surface(&self) -> Option<f64> {
        self.wetted_surface
    }

    /// Set the wetted surface area.
    pub fn set_wetted_surface(&mut self, surface: Option<f64>) {
        self.wetted_surface = surface;
    }

    /// Calculate the volume of the appendage.
    pub fn volume(&self) -> f64 {
        match &self.geometry {
            AppendageGeometry::Point { volume, .. } => *volume,
            AppendageGeometry::Mesh(mesh) => {
                use parry3d_f64::shape::Shape;
                mesh.mass_properties(1.0).mass()
            }
            AppendageGeometry::Box { bounds } => {
                let (xmin, xmax, ymin, ymax, zmin, zmax) = *bounds;
                (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
            }
            AppendageGeometry::Sphere { volume, .. } => *volume,
            AppendageGeometry::Cube { volume, .. } => *volume,
        }
    }

    /// Calculate the center of volume of the appendage.
    pub fn center(&self) -> [f64; 3] {
        match &self.geometry {
            AppendageGeometry::Point { center, .. } => *center,
            AppendageGeometry::Mesh(mesh) => {
                use parry3d_f64::shape::Shape;
                let com = mesh.mass_properties(1.0).local_com;
                [com.x, com.y, com.z]
            }
            AppendageGeometry::Box { bounds } => {
                let (xmin, xmax, ymin, ymax, zmin, zmax) = *bounds;
                [
                    (xmin + xmax) / 2.0,
                    (ymin + ymax) / 2.0,
                    (zmin + zmax) / 2.0,
                ]
            }
            AppendageGeometry::Sphere { center, .. } => *center,
            AppendageGeometry::Cube { center, .. } => *center,
        }
    }

    // =========================================================================
    // Point-specific mutators
    // =========================================================================

    /// Set the volume for a Point appendage.
    ///
    /// Returns an error if the geometry is not a Point.
    pub fn set_point_volume(&mut self, volume: f64) -> Result<(), &'static str> {
        match &mut self.geometry {
            AppendageGeometry::Point { volume: v, .. } => {
                *v = volume;
                Ok(())
            }
            _ => Err("Cannot set volume: appendage is not a Point"),
        }
    }

    /// Set the center for a Point appendage.
    ///
    /// Returns an error if the geometry is not a Point.
    pub fn set_point_center(&mut self, center: [f64; 3]) -> Result<(), &'static str> {
        match &mut self.geometry {
            AppendageGeometry::Point { center: c, .. } => {
                *c = center;
                Ok(())
            }
            _ => Err("Cannot set center: appendage is not a Point"),
        }
    }

    // =========================================================================
    // Private helpers
    // =========================================================================

    fn vtk_to_mesh(vtk: &vtkio::Vtk) -> Result<TriMesh, AppendageLoadError> {
        use vtkio::model::{DataSet, Piece};

        fn extract_mesh_inner(
            data: &DataSet,
        ) -> Result<TriMesh, AppendageLoadError> {
            match data {
                DataSet::PolyData { pieces, .. } => {
                    for piece in pieces {
                        if let Piece::Inline(piece_data) = piece {
                            let points = &piece_data.points;
                            let vertices: Vec<nalgebra::Point3<f64>> =
                                extract_points_from_iobuffer(points)?;

                            if let Some(polys) = &piece_data.polys {
                                let indices = extract_triangles_from_vertex_numbers(polys)?;
                                return TriMesh::new(vertices, indices).map_err(|e| {
                                    AppendageLoadError::StlParseError(format!("{:?}", e))
                                });
                            }
                        }
                    }
                    Err(AppendageLoadError::StlParseError(
                        "No valid polygon data found".to_string(),
                    ))
                }
                DataSet::UnstructuredGrid { pieces, .. } => {
                    for piece in pieces {
                        if let Piece::Inline(piece_data) = piece {
                            let points = &piece_data.points;
                            let vertices: Vec<nalgebra::Point3<f64>> =
                                extract_points_from_iobuffer(points)?;

                            let cells = &piece_data.cells;
                            let indices = extract_cells_triangles_inner(cells)?;
                            return TriMesh::new(vertices, indices).map_err(|e| {
                                AppendageLoadError::StlParseError(format!("{:?}", e))
                            });
                        }
                    }
                    Err(AppendageLoadError::StlParseError(
                        "No valid unstructured grid data found".to_string(),
                    ))
                }
                _ => Err(AppendageLoadError::StlParseError(
                    "Unsupported VTK data set type".to_string(),
                )),
            }
        }

        fn extract_points_from_iobuffer(
            buffer: &vtkio::IOBuffer,
        ) -> Result<Vec<nalgebra::Point3<f64>>, AppendageLoadError> {
            let coords: Vec<f64> = match buffer {
                vtkio::IOBuffer::F64(v) => v.clone(),
                vtkio::IOBuffer::F32(v) => v.iter().map(|x| *x as f64).collect(),
                _ => {
                    return Err(AppendageLoadError::StlParseError(
                        "Unsupported point data type".to_string(),
                    ))
                }
            };

            Ok(coords
                .chunks(3)
                .map(|c| nalgebra::Point3::new(c[0], c[1], c[2]))
                .collect())
        }

        fn extract_triangles_from_vertex_numbers(
            polys: &vtkio::model::VertexNumbers,
        ) -> Result<Vec<[u32; 3]>, AppendageLoadError> {
            let mut indices = Vec::new();
            match polys {
                vtkio::model::VertexNumbers::Legacy {
                    num_cells: _,
                    vertices,
                } => {
                    let mut i = 0;
                    while i < vertices.len() {
                        let n = vertices[i] as usize;
                        if n == 3 && i + 3 < vertices.len() {
                            indices.push([vertices[i + 1], vertices[i + 2], vertices[i + 3]]);
                        }
                        i += n + 1;
                    }
                }
                vtkio::model::VertexNumbers::XML {
                    connectivity,
                    offsets,
                } => {
                    let mut start = 0u64;
                    for &end in offsets.iter() {
                        let end_usize = end as usize;
                        if end_usize - start as usize == 3 {
                            indices.push([
                                connectivity[start as usize] as u32,
                                connectivity[start as usize + 1] as u32,
                                connectivity[start as usize + 2] as u32,
                            ]);
                        }
                        start = end as u64;
                    }
                }
            }
            Ok(indices)
        }

        fn extract_cells_triangles_inner(
            cells: &vtkio::model::Cells,
        ) -> Result<Vec<[u32; 3]>, AppendageLoadError> {
            let mut indices = Vec::new();

            let cell_verts = match &cells.cell_verts {
                vtkio::model::VertexNumbers::Legacy { vertices, .. } => vertices.clone(),
                vtkio::model::VertexNumbers::XML { connectivity, .. } => {
                    connectivity.iter().map(|x| *x as u32).collect()
                }
            };

            let mut i = 0;
            while i + 3 < cell_verts.len() {
                let n = cell_verts[i] as usize;
                if n == 3 {
                    indices.push([cell_verts[i + 1], cell_verts[i + 2], cell_verts[i + 3]]);
                }
                i += n + 1;
            }

            Ok(indices)
        }

        extract_mesh_inner(&vtk.data)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_appendage() {
        let app = Appendage::from_point("rudder", [50.0, 0.0, -2.0], 0.5);
        assert_eq!(app.name(), "rudder");
        assert!((app.volume() - 0.5).abs() < 1e-9);
        assert_eq!(app.center(), [50.0, 0.0, -2.0]);
    }

    #[test]
    fn test_box_appendage() {
        let app = Appendage::from_box("keel", (40.0, 60.0, -1.0, 1.0, -3.0, 0.0));
        // Volume = 20 * 2 * 3 = 120 m³
        assert!((app.volume() - 120.0).abs() < 1e-9);
        // Center = (50, 0, -1.5)
        let center = app.center();
        assert!((center[0] - 50.0).abs() < 1e-9);
        assert!((center[1] - 0.0).abs() < 1e-9);
        assert!((center[2] - -1.5).abs() < 1e-9);
    }

    #[test]
    fn test_sphere_appendage() {
        let volume = 10.0;
        let app = Appendage::from_sphere("bulb", [0.0, 0.0, 0.0], volume);
        assert!((app.volume() - volume).abs() < 1e-9);
    }

    #[test]
    fn test_cube_appendage() {
        let volume = 27.0; // 3³
        let app = Appendage::from_cube("block", [0.0, 0.0, 0.0], volume);
        assert!((app.volume() - volume).abs() < 1e-9);
    }

    #[test]
    fn test_set_point_volume() {
        let mut app = Appendage::from_point("test", [0.0, 0.0, 0.0], 1.0);
        assert!(app.set_point_volume(2.0).is_ok());
        assert!((app.volume() - 2.0).abs() < 1e-9);

        let mut box_app = Appendage::from_box("box", (0.0, 1.0, 0.0, 1.0, 0.0, 1.0));
        assert!(box_app.set_point_volume(2.0).is_err());
    }

    #[test]
    fn test_set_point_center() {
        let mut app = Appendage::from_point("test", [0.0, 0.0, 0.0], 1.0);
        assert!(app.set_point_center([1.0, 2.0, 3.0]).is_ok());
        assert_eq!(app.center(), [1.0, 2.0, 3.0]);
    }
}
