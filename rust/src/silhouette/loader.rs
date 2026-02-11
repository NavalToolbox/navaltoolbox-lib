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

//! File loaders for silhouette profiles (DXF and VTK).

use dxf::entities::EntityType;
use dxf::Drawing;
use std::fs;
use std::path::Path;
use thiserror::Error;
use vtkio::model::{DataSet, Piece};
use vtkio::Vtk;

#[derive(Error, Debug)]
pub enum SilhouetteLoadError {
    #[error("Failed to read file: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Failed to parse DXF file: {0}")]
    DxfParseError(#[from] dxf::DxfError),
    #[error("Failed to parse VTK file: {0}")]
    VtkParseError(#[from] vtkio::Error),
    #[error("No polyline found in file")]
    NoPolyline,
    #[error("Unsupported file format")]
    UnsupportedFormat,
    #[error("Failed to parse CSV: {0}")]
    CsvError(String),
}

// Legacy alias for DxfError

/// Load a silhouette from a DXF file.
///
/// Handles both LwPolyline and Polyline (AcDb2dPolyline) entities.
/// For 2D polylines, the OCS (Object Coordinate System) normal is used
/// to correctly map vertex coordinates to the X-Z world plane.
pub fn load_dxf_silhouette(path: &Path) -> Result<(Vec<[f64; 3]>, String), SilhouetteLoadError> {
    let drawing = Drawing::load_file(path)?;

    let name = path
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "silhouette".to_string());

    for entity in drawing.entities() {
        match &entity.specific {
            EntityType::LwPolyline(lwpoly) => {
                // LwPolyline vertices are 2D (x, y) in the OCS plane.
                // Default extrusion is (0,0,1) → x,y map to world X,Y.
                // For silhouettes we want X-Z plane, so we map (x,y) → (x, 0, y).
                // This is correct for the default case where the polyline
                // is drawn in a 2D view looking down the Z axis.
                let mut points: Vec<[f64; 3]> =
                    lwpoly.vertices.iter().map(|v| [v.x, 0.0, v.y]).collect();

                let is_closed = (lwpoly.flags & 1) != 0;
                if is_closed && !points.is_empty() && points.first() != points.last() {
                    points.push(points[0]);
                }

                if !points.is_empty() {
                    log::info!(
                        "Loaded LwPolyline silhouette '{}' with {} points",
                        name,
                        points.len()
                    );
                    return Ok((points, name));
                }
            }
            EntityType::Polyline(poly) => {
                // Polyline (AcDb2dPolyline) has an OCS normal vector.
                // We need to map vertex coordinates from OCS to world X-Z plane.
                let normal = &poly.normal;
                let is_2d_xz_plane = normal.y.abs() > 0.9;
                let is_2d_xy_plane = normal.z.abs() > 0.9;

                let mut points: Vec<[f64; 3]> = Vec::new();

                for vertex in poly.vertices() {
                    let (world_x, world_z) = if is_2d_xz_plane {
                        // OCS normal ≈ (0, ±1, 0): polyline is in world X-Z plane
                        // OCS (x, y) → world (x, z=y)
                        (vertex.location.x, vertex.location.y)
                    } else if is_2d_xy_plane {
                        // OCS normal ≈ (0, 0, ±1): polyline is in world X-Y plane
                        // OCS (x, y) → world (x, z=y) — same mapping for silhouettes
                        (vertex.location.x, vertex.location.y)
                    } else {
                        // Fallback: use x and z directly (3D polyline)
                        log::warn!(
                            "Polyline has unusual OCS normal ({:.2}, {:.2}, {:.2}). \
                             Using vertex.x and vertex.z directly.",
                            normal.x,
                            normal.y,
                            normal.z
                        );
                        (vertex.location.x, vertex.location.z)
                    };

                    points.push([world_x, 0.0, world_z]);
                }

                let is_closed = (poly.flags & 1) != 0;
                if is_closed && !points.is_empty() && points.first() != points.last() {
                    points.push(points[0]);
                }

                if !points.is_empty() {
                    log::info!(
                        "Loaded Polyline silhouette '{}' with {} points (OCS normal: {:.2}, {:.2}, {:.2})",
                        name,
                        points.len(),
                        normal.x,
                        normal.y,
                        normal.z
                    );
                    return Ok((points, name));
                }
            }
            _ => continue,
        }
    }

    Err(SilhouetteLoadError::NoPolyline)
}

/// Load a silhouette from a VTK file (.vtk or .vtp).
pub fn load_vtk_silhouette(path: &Path) -> Result<(Vec<[f64; 3]>, String), SilhouetteLoadError> {
    let vtk = Vtk::import(path)?;

    let name = path
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "silhouette".to_string());

    // Extract points based on data type
    match &vtk.data {
        DataSet::PolyData { pieces, .. } => {
            // Get the first inline piece with data
            for piece in pieces {
                if let Piece::Inline(polydata) = piece {
                    return extract_points_from_iobuffer(&polydata.points, &name);
                }
            }
            Err(SilhouetteLoadError::NoPolyline)
        }
        DataSet::UnstructuredGrid { pieces, .. } => {
            // Get the first inline piece with data
            for piece in pieces {
                if let Piece::Inline(grid) = piece {
                    return extract_points_from_iobuffer(&grid.points, &name);
                }
            }
            Err(SilhouetteLoadError::NoPolyline)
        }
        _ => Err(SilhouetteLoadError::UnsupportedFormat),
    }
}

/// Load a silhouette from a CSV or TXT file (X,Z or X,Y,Z points).
pub fn load_csv_silhouette(path: &Path) -> Result<(Vec<[f64; 3]>, String), SilhouetteLoadError> {
    let content = fs::read_to_string(path).map_err(SilhouetteLoadError::IoError)?;
    let mut points = Vec::new();

    let name = path
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "silhouette".to_string());

    for (i, line) in content.lines().enumerate() {
        let line = line.trim();
        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#') || line.starts_with("//") {
            continue;
        }

        // Split by comma, semicolon, space, or tab
        let parts: Vec<&str> = line
            .split(&[',', ';', ' ', '\t'][..])
            .filter(|s| !s.is_empty())
            .collect();

        // Skip headers (non-numeric first token)
        if let Some(first) = parts.first() {
            if first.parse::<f64>().is_err() {
                continue;
            }
        }

        if parts.len() >= 2 {
            let x = parts[0].parse::<f64>().map_err(|_| {
                SilhouetteLoadError::CsvError(format!("Line {}: invalid X coordinate", i + 1))
            })?;

            // Assume 2 columns = X, Z (Profile)
            // Assume 3 columns = X, Y, Z (3D point) -> Project to X-Z (use X and Z)
            let z = if parts.len() >= 3 {
                parts[2].parse::<f64>().map_err(|_| {
                    SilhouetteLoadError::CsvError(format!("Line {}: invalid Z coordinate", i + 1))
                })?
            } else {
                parts[1].parse::<f64>().map_err(|_| {
                    SilhouetteLoadError::CsvError(format!("Line {}: invalid Z coordinate", i + 1))
                })?
            };

            points.push([x, 0.0, z]);
        }
    }

    if points.is_empty() {
        return Err(SilhouetteLoadError::NoPolyline);
    }

    // Auto-close
    if points.first() != points.last() {
        points.push(points[0]);
    }

    log::info!(
        "Loaded CSV/TXT silhouette '{}' with {} points",
        name,
        points.len()
    );
    Ok((points, name))
}

fn extract_points_from_iobuffer(
    buffer: &vtkio::IOBuffer,
    name: &str,
) -> Result<(Vec<[f64; 3]>, String), SilhouetteLoadError> {
    let points_f64 = match buffer {
        vtkio::IOBuffer::F64(data) => data
            .chunks(3)
            .filter(|c| c.len() == 3)
            .map(|c| [c[0], c[1], c[2]])
            .collect(),
        vtkio::IOBuffer::F32(data) => data
            .chunks(3)
            .filter(|c| c.len() == 3)
            .map(|c| [c[0] as f64, c[1] as f64, c[2] as f64])
            .collect(),
        _ => Vec::new(),
    };

    if points_f64.is_empty() {
        return Err(SilhouetteLoadError::NoPolyline);
    }

    // Convert to X-Z plane (force Y=0)
    let mut has_nonzero_y = false;
    let mut result: Vec<[f64; 3]> = points_f64
        .iter()
        .map(|p| {
            if p[1].abs() > 1e-6 {
                has_nonzero_y = true;
            }
            [p[0], 0.0, p[2]]
        })
        .collect();

    if has_nonzero_y {
        log::warn!("VTK polyline has non-zero Y. Setting Y=0 for X-Z plane.");
    }

    // Close if not already closed
    if !result.is_empty() && result.first() != result.last() {
        result.push(result[0]);
    }

    Ok((result, name.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_csv(content: &str, ext: &str) -> NamedTempFile {
        let mut f = tempfile::Builder::new().suffix(ext).tempfile().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
    }

    #[test]
    fn test_csv_two_columns_comma() {
        let f = write_csv("0.0, 0.0\n10.0, 0.0\n10.0, 5.0\n0.0, 5.0\n", ".csv");
        let (points, _name) = load_csv_silhouette(f.path()).unwrap();
        // 4 data points + 1 auto-close = 5
        assert_eq!(points.len(), 5);
        assert_eq!(points[0], [0.0, 0.0, 0.0]);
        assert_eq!(points[1], [10.0, 0.0, 0.0]);
        assert_eq!(points[2], [10.0, 0.0, 5.0]);
        assert_eq!(points[3], [0.0, 0.0, 5.0]);
        assert_eq!(points[4], points[0]); // auto-closed
    }

    #[test]
    fn test_csv_three_columns_space() {
        // X Y Z format, Y should be ignored (set to 0)
        let f = write_csv(
            "0.0 1.0 0.0\n10.0 2.0 0.0\n10.0 3.0 5.0\n0.0 4.0 5.0\n",
            ".txt",
        );
        let (points, _name) = load_csv_silhouette(f.path()).unwrap();
        assert_eq!(points.len(), 5); // 4 + auto-close
                                     // Y should always be 0
        for p in &points {
            assert_eq!(p[1], 0.0);
        }
        // Z comes from third column
        assert_eq!(points[2], [10.0, 0.0, 5.0]);
    }

    #[test]
    fn test_csv_semicolons() {
        let f = write_csv("0.0; 0.0\n10.0; 0.0\n10.0; 5.0\n", ".csv");
        let (points, _name) = load_csv_silhouette(f.path()).unwrap();
        assert_eq!(points.len(), 4); // 3 + auto-close
        assert_eq!(points[0], [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_csv_tabs() {
        let f = write_csv("0.0\t0.0\n10.0\t0.0\n10.0\t5.0\n", ".txt");
        let (points, _name) = load_csv_silhouette(f.path()).unwrap();
        assert_eq!(points.len(), 4); // 3 + auto-close
    }

    #[test]
    fn test_csv_skips_comments_and_headers() {
        let content = "# This is a comment\n// Another comment\nX, Z\n\n0.0, 0.0\n10.0, 5.0\n";
        let f = write_csv(content, ".csv");
        let (points, _name) = load_csv_silhouette(f.path()).unwrap();
        assert_eq!(points.len(), 3); // 2 data + auto-close
    }

    #[test]
    fn test_csv_empty_file_fails() {
        let f = write_csv("# Only comments\n// nothing here\n", ".csv");
        let result = load_csv_silhouette(f.path());
        assert!(result.is_err());
    }

    #[test]
    fn test_csv_already_closed() {
        // First == Last, so no auto-close needed
        let f = write_csv("0.0, 0.0\n10.0, 0.0\n10.0, 5.0\n0.0, 0.0\n", ".csv");
        let (points, _name) = load_csv_silhouette(f.path()).unwrap();
        assert_eq!(points.len(), 4); // No extra point added
        assert_eq!(points[0], points[3]);
    }

    #[test]
    fn test_csv_name_from_filename() {
        let f = write_csv("0.0, 0.0\n10.0, 5.0\n", ".csv");
        let (_points, name) = load_csv_silhouette(f.path()).unwrap();
        // Name comes from file_stem, should not be empty
        assert!(!name.is_empty());
    }
}
