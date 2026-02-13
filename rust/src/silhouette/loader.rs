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
use nalgebra::{Matrix4, Point3, Vector3};
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

/// Compute the Object Coordinate System (OCS) to World Coordinate System (WCS) transformation matrix
/// using the Arbitrary Axis Algorithm.
fn get_ocs_to_wcs_matrix(normal: Vector3<f64>, elevation: f64) -> Matrix4<f64> {
    // Arbitrary Axis Algorithm
    let threshold = 1.0 / 64.0;
    let x_axis = if normal.x.abs() < threshold && normal.y.abs() < threshold {
        Vector3::y().cross(&normal)
    } else {
        Vector3::z().cross(&normal)
    };

    let x_axis = x_axis.normalize();
    let y_axis = normal.cross(&x_axis).normalize();
    let z_axis = normal.normalize();

    let mut matrix = Matrix4::identity();
    matrix.set_column(0, &x_axis.to_homogeneous());
    matrix.set_column(1, &y_axis.to_homogeneous());
    matrix.set_column(2, &z_axis.to_homogeneous());

    // The OCS origin in WCS is defined by the elevation along the Z axis (normal).
    let origin = z_axis * elevation;
    matrix[(0, 3)] = origin.x;
    matrix[(1, 3)] = origin.y;
    matrix[(2, 3)] = origin.z;

    matrix
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
                // Get OCS normal (extrusion direction)
                let normal = Vector3::new(
                    lwpoly.extrusion_direction.x,
                    lwpoly.extrusion_direction.y,
                    lwpoly.extrusion_direction.z,
                );

                // Note: dxf crate 0.5 LwPolyline struct seemingly doesn't expose elevation directly?
                // For now assuming 0.0 or TODO: check if it's in a different field or needs crate update.
                let elevation = 0.0;

                let transform = get_ocs_to_wcs_matrix(normal, elevation);
                let is_xy_plane = normal.z.abs() > 0.9; // Check if basically Top/Bottom view

                let mut points: Vec<[f64; 3]> = Vec::new();

                for v in &lwpoly.vertices {
                    // LwPolyline vertices are (x, y) in OCS. Z is elevation.
                    // Elevation is handled in the transform matrix (translation along Z axis).
                    let p_ocs = Point3::new(v.x, v.y, 0.0);
                    let p_wcs = transform.transform_point(&p_ocs);

                    // Mapping strategy:
                    // If the entity is in the XY plane (Top/Bottom view), we map WCS Y -> Output Z (Profile view convention).
                    // Otherwise (Front/Side/Arbitrary), we assume WCS Z -> Output Z.
                    // Output Y is always 0 (Centerline/Silhouette plane).

                    let x = p_wcs.x;
                    let z = if is_xy_plane { p_wcs.y } else { p_wcs.z };

                    points.push([x, 0.0, z]);
                }

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
                let normal = Vector3::new(poly.normal.x, poly.normal.y, poly.normal.z);
                // Note: Polyline struct also missing elevation in this crate version?
                let elevation = 0.0;

                let transform = get_ocs_to_wcs_matrix(normal, elevation);
                let is_xy_plane = normal.z.abs() > 0.9;

                let mut points: Vec<[f64; 3]> = Vec::new();

                for vertex in poly.vertices() {
                    // For AcDb2dPolyline, vertices are in OCS.
                    // For AcDb3dPolyline, vertices are in WCS (and normal is 0,0,1 usually).
                    // If normal is (0,1,0), it's definitely OCS-based 2D polyline.

                    // We treat vertex location as OCS param.
                    // Note: vertex.location.z should be 0 for 2D polyline, but if not, we include it.
                    // Elevation is handled by the matrix translation.
                    let p_ocs =
                        Point3::new(vertex.location.x, vertex.location.y, vertex.location.z);
                    let p_wcs = transform.transform_point(&p_ocs);

                    let x = p_wcs.x;
                    let z = if is_xy_plane { p_wcs.y } else { p_wcs.z };

                    points.push([x, 0.0, z]);
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
