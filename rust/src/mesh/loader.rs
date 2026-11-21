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

fn validate_mesh_topology(indices: &[[u32; 3]]) -> Result<(), String> {
    use std::collections::HashMap;
    let mut edge_count: HashMap<(u32, u32), u32> = HashMap::new();
    for tri in indices {
        for i in 0..3 {
            let a = tri[i];
            let b = tri[(i + 1) % 3];
            *edge_count.entry((a, b)).or_insert(0) += 1;
        }
    }
    let bad: u32 = edge_count
        .iter()
        .filter(|((a, b), &k)| k != 1 || edge_count.get(&(*b, *a)).copied().unwrap_or(0) != 1)
        .count() as u32;
    if bad > 0 {
        Err(format!(
            "Mesh topology error: {} directed edge(s) fail the pairing check. \
             The mesh is not a closed orientable manifold (open boundary or flipped triangles).",
            bad
        ))
    } else {
        Ok(())
    }
}

const ORIENTATION_EPSILON: f64 = 1e-10;

fn validate_mesh_orientation(vertices: &[Point3<f64>], indices: &[[u32; 3]]) -> Result<(), String> {
    let signed_vol: f64 = indices
        .iter()
        .map(|tri| {
            let [a, b, c] = *tri;
            let p = vertices[a as usize];
            let q = vertices[b as usize];
            let r = vertices[c as usize];
            p.x * (q.y * r.z - q.z * r.y) - p.y * (q.x * r.z - q.z * r.x)
                + p.z * (q.x * r.y - q.y * r.x)
        })
        .sum::<f64>()
        / 6.0;

    if signed_vol < -ORIENTATION_EPSILON {
        Err(format!(
            "Mesh orientation error: signed volume is {:.6e} m³. \
             The mesh is inside-out (all normals point inward). \
             Reverse triangle winding in your CAD export.",
            signed_vol
        ))
    } else {
        Ok(())
    }
}

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

    let mesh = TriMesh::new(vertices.clone(), indices.clone()).map_err(|e| {
        Error::new(
            ErrorKind::InvalidData,
            format!("Mesh creation error: {:?}", e),
        )
    })?;

    validate_mesh_topology(&indices).map_err(|msg| Error::new(ErrorKind::InvalidData, msg))?;
    validate_mesh_orientation(&vertices, &indices)
        .map_err(|msg| Error::new(ErrorKind::InvalidData, msg))?;

    Ok(mesh)
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

    #[test]
    fn test_validate_mesh_topology_inside_out() {
        // A minimal closed tetrahedron with 4 vertices and 4 faces
        let valid_tetra = vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]];
        assert!(validate_mesh_topology(&valid_tetra).is_ok());

        // A mesh with a flipped triangle (inconsistent winding / inside-out face)
        let flipped_tetra = vec![[2, 0, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]];
        let err = validate_mesh_topology(&flipped_tetra).unwrap_err();
        assert!(err.contains("Mesh topology error"));
        assert!(err.contains("directed edge(s) fail the pairing check"));

        // An open mesh (missing face)
        let open_tetra = vec![[0, 2, 1], [0, 1, 3], [1, 2, 3]];
        let err_open = validate_mesh_topology(&open_tetra).unwrap_err();
        assert!(err_open.contains("3 directed edge(s) fail the pairing check"));
    }

    #[test]
    fn test_validate_mesh_orientation_inside_out() {
        // Regular tetrahedron with vertices at:
        let vertices = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ];

        // Correct outward-pointing winding:
        // [0, 2, 1] (normal in -z direction for z=0 plane), [0, 1, 3], [1, 2, 3], [2, 0, 3]
        let valid_tetra = vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]];
        assert!(validate_mesh_topology(&valid_tetra).is_ok());
        assert!(validate_mesh_orientation(&vertices, &valid_tetra).is_ok());

        // Inside-out tetrahedron: every triangle reversed ([a, b, c] -> [a, c, b])
        // Edge-pairing is still valid (directed edge (a, c) pairs with (c, a))
        let inside_out_tetra = vec![[0, 1, 2], [0, 3, 1], [1, 3, 2], [2, 3, 0]];
        assert!(validate_mesh_topology(&inside_out_tetra).is_ok());

        let err = validate_mesh_orientation(&vertices, &inside_out_tetra).unwrap_err();
        assert!(err.contains("Mesh orientation error"));
        assert!(err.contains("signed volume is"));
        assert!(err.contains("The mesh is inside-out"));
        assert!(err.contains("Reverse triangle winding in your CAD export."));
    }
}
