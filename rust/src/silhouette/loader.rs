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
}

// Legacy alias for DxfError

/// Load a silhouette from a DXF file.
pub fn load_dxf_silhouette(path: &Path) -> Result<(Vec<[f64; 3]>, String), SilhouetteLoadError> {
    let drawing = Drawing::load_file(path)?;

    let name = path
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "silhouette".to_string());

    for entity in drawing.entities() {
        match &entity.specific {
            EntityType::LwPolyline(lwpoly) => {
                let mut points: Vec<[f64; 3]> =
                    lwpoly.vertices.iter().map(|v| [v.x, 0.0, v.y]).collect();

                let is_closed = (lwpoly.flags & 1) != 0;
                if is_closed && !points.is_empty() && points.first() != points.last() {
                    points.push(points[0]);
                }

                if !points.is_empty() {
                    return Ok((points, name));
                }
            }
            EntityType::Polyline(poly) => {
                let mut points: Vec<[f64; 3]> = Vec::new();
                let mut has_nonzero_y = false;

                for vertex in poly.vertices() {
                    let y = vertex.location.y;
                    if y.abs() > 1e-6 {
                        has_nonzero_y = true;
                    }
                    points.push([vertex.location.x, 0.0, vertex.location.z]);
                }

                if has_nonzero_y {
                    log::warn!("DXF polyline has non-zero Y. Setting Y=0 for X-Z plane.");
                }

                let is_closed = (poly.flags & 1) != 0;
                if is_closed && !points.is_empty() && points.first() != points.last() {
                    points.push(points[0]);
                }

                if !points.is_empty() {
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
    // Tests require actual DXF/VTK files
}
