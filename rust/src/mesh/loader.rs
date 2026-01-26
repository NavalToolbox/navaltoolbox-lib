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

//! Mesh loading from STL files.
//!
//! Note: VTK support can be added later via vtkio crate.

use nalgebra::Point3;
use parry3d_f64::shape::TriMesh;
use std::fs::File;
use std::io::{Cursor, Error, ErrorKind, Read};
use std::path::Path;

/// Loads a mesh from an STL file.
///
/// Supports both binary and ASCII STL formats.
/// Uses a robust strategy: loads entire file into memory before parsing
/// to avoid buffer underrun issues with corrupted binary STL headers.
pub fn load_stl(path: &Path) -> Result<TriMesh, Error> {
    // Read the entire file into memory first
    // This prevents issues with binary STL files that have incorrect triangle counts in the header
    let mut file = File::open(path)?;
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;

    if buffer.is_empty() {
        return Err(Error::new(ErrorKind::InvalidData, "STL file is empty"));
    }

    // Handle leading whitespace in ASCII STL files.
    // Some CAD software exports ASCII STL with leading spaces/newlines before "solid".
    // We trim leading whitespace to ensure correct parsing.
    let trimmed_start = buffer
        .iter()
        .position(|&c| !c.is_ascii_whitespace())
        .unwrap_or(0);
    let buffer = if trimmed_start > 0 {
        buffer[trimmed_start..].to_vec()
    } else {
        buffer
    };

    // Create a cursor over the buffer for parsing
    let mut cursor = Cursor::new(buffer);

    // Try to parse the STL
    let stl = stl_io::read_stl(&mut cursor).map_err(|e| {
        let error_msg = format!("{}", e);

        // Provide helpful error messages based on the error type
        if error_msg.contains("failed to fill whole buffer") {
            Error::new(
                ErrorKind::InvalidData,
                format!(
                    "STL parse error: The file appears to be a binary STL with an incorrect triangle count in the header, \
                    or the file is truncated/corrupted. Original error: {}",
                    error_msg
                ),
            )
        } else {
            Error::new(
                ErrorKind::InvalidData,
                format!("STL parse error: {}", error_msg),
            )
        }
    })?;

    if stl.vertices.is_empty() {
        return Err(Error::new(
            ErrorKind::InvalidData,
            "No geometry data in STL file",
        ));
    }

    let vertices: Vec<Point3<f64>> = stl
        .vertices
        .iter()
        .map(|v| Point3::new(v[0] as f64, v[1] as f64, v[2] as f64))
        .collect();

    let indices: Vec<[u32; 3]> = stl
        .faces
        .iter()
        .map(|f| {
            [
                f.vertices[0] as u32,
                f.vertices[1] as u32,
                f.vertices[2] as u32,
            ]
        })
        .collect();

    TriMesh::new(vertices, indices).map_err(|e| {
        Error::new(
            ErrorKind::InvalidData,
            format!("Mesh creation error: {:?}", e),
        )
    })
}

/// Placeholder for VTK loading (to be implemented later).
pub fn load_vtk(_path: &Path) -> Result<TriMesh, Error> {
    Err(Error::new(
        ErrorKind::Unsupported,
        "VTK loading not yet implemented. Use STL files.",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_load_stl_box() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/box_10x10.stl");

        if path.exists() {
            let mesh = load_stl(&path).expect("Failed to load STL");
            assert!(!mesh.vertices().is_empty());
            assert!(!mesh.indices().is_empty());
        }
    }
}
