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

//! Hull geometry module.
//!
//! Provides loading, transformation, and export for hull geometries.

use nalgebra::{Point3, Vector3};
use parry3d_f64::shape::TriMesh;
use std::path::Path;
use thiserror::Error;

use crate::mesh::{get_bounds, load_stl, load_vtk, transform_mesh};

/// Errors that can occur during hull operations.
#[derive(Error, Debug)]
pub enum HullError {
    #[error("File not found: {0}")]
    FileNotFound(String),

    #[error("Unsupported format: {0}. Use .stl or .vtk")]
    UnsupportedFormat(String),

    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("No geometry data in file")]
    EmptyGeometry,
}

/// Represents a hull geometry with loading, transformation, and export capabilities.
#[derive(Clone)]
pub struct Hull {
    /// The underlying triangle mesh
    mesh: TriMesh,
    /// Original file path (if loaded from file)
    file_path: Option<String>,
}

impl Hull {
    /// Creates a Hull from an existing TriMesh.
    pub fn from_mesh(mesh: TriMesh) -> Self {
        Self {
            mesh,
            file_path: None,
        }
    }

    /// Creates a box hull with given dimensions.
    pub fn from_box(length: f64, breadth: f64, depth: f64) -> Self {
        let hb = breadth / 2.0;
        let _hl = length / 2.0;
        // Note: Convention usually puts origin at AP or Midships.
        // Let's assume origin at 0,0,0 (AP if x positive?), centered transversely
        // Box vertices [0..L], [-B/2..B/2], [0..D]

        // Vertices
        let vertices = vec![
            Point3::new(0.0, -hb, 0.0),
            Point3::new(length, -hb, 0.0),
            Point3::new(length, hb, 0.0),
            Point3::new(0.0, hb, 0.0),
            Point3::new(0.0, -hb, depth),
            Point3::new(length, -hb, depth),
            Point3::new(length, hb, depth),
            Point3::new(0.0, hb, depth),
        ];

        // Indices (12 triangles)
        let indices = vec![
            [0, 2, 1],
            [0, 3, 2], // Bottom
            [4, 5, 6],
            [4, 6, 7], // Top
            [0, 1, 5],
            [0, 5, 4], // Side 1
            [2, 3, 7],
            [2, 7, 6], // Side 2
            [0, 4, 7],
            [0, 7, 3], // Back
            [1, 2, 6],
            [1, 6, 5], // Front
        ];

        let mesh = TriMesh::new(vertices, indices).expect("Failed to create box mesh");
        Self::from_mesh(mesh)
    }

    /// Loads a hull geometry from an STL file.
    pub fn from_stl<P: AsRef<Path>>(path: P) -> Result<Self, HullError> {
        let path = path.as_ref();
        if !path.exists() {
            return Err(HullError::FileNotFound(path.display().to_string()));
        }

        let mesh = load_stl(path)?;

        if mesh.vertices().is_empty() {
            return Err(HullError::EmptyGeometry);
        }

        Ok(Self {
            mesh,
            file_path: Some(path.display().to_string()),
        })
    }

    /// Loads a hull geometry from a VTK file.
    pub fn from_vtk<P: AsRef<Path>>(path: P) -> Result<Self, HullError> {
        let path = path.as_ref();
        if !path.exists() {
            return Err(HullError::FileNotFound(path.display().to_string()));
        }

        let mesh = load_vtk(path)?;

        if mesh.vertices().is_empty() {
            return Err(HullError::EmptyGeometry);
        }

        Ok(Self {
            mesh,
            file_path: Some(path.display().to_string()),
        })
    }

    /// Loads a hull from any supported format (determined by extension).
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, HullError> {
        let path = path.as_ref();
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.to_lowercase())
            .unwrap_or_default();

        match ext.as_str() {
            "stl" => Self::from_stl(path),
            "vtk" | "vtp" | "vtu" => Self::from_vtk(path),
            _ => Err(HullError::UnsupportedFormat(ext)),
        }
    }

    /// Returns the bounding box of the hull.
    ///
    /// Returns (xmin, xmax, ymin, ymax, zmin, zmax).
    pub fn get_bounds(&self) -> (f64, f64, f64, f64, f64, f64) {
        get_bounds(&self.mesh)
    }

    /// Returns the underlying TriMesh (for internal calculations).
    pub fn mesh(&self) -> &TriMesh {
        &self.mesh
    }

    /// Returns a mutable reference to the underlying TriMesh.
    pub fn mesh_mut(&mut self) -> &mut TriMesh {
        &mut self.mesh
    }

    /// Returns the number of triangles in the hull.
    pub fn num_triangles(&self) -> usize {
        self.mesh.indices().len()
    }

    /// Returns the number of vertices in the hull.
    pub fn num_vertices(&self) -> usize {
        self.mesh.vertices().len()
    }

    /// Applies a transformation to the hull geometry.
    ///
    /// # Arguments
    /// * `translation` - (dx, dy, dz) translation vector
    /// * `rotation` - (rx, ry, rz) rotation angles in degrees around X, Y, Z axes
    /// * `pivot` - (px, py, pz) point around which rotation occurs
    pub fn transform(
        &mut self,
        translation: (f64, f64, f64),
        rotation: (f64, f64, f64),
        pivot: (f64, f64, f64),
    ) {
        let pivot_point = Point3::new(pivot.0, pivot.1, pivot.2);

        // Apply rotation (heel=rx, trim=ry for now, rz not used in transform_mesh)
        // For full rotation support, we'd need to extend transform_mesh
        self.mesh = transform_mesh(&self.mesh, rotation.0, rotation.1, pivot_point);

        // Apply translation
        if translation.0 != 0.0 || translation.1 != 0.0 || translation.2 != 0.0 {
            let trans_vec = Vector3::new(translation.0, translation.1, translation.2);
            let new_vertices: Vec<Point3<f64>> =
                self.mesh.vertices().iter().map(|v| v + trans_vec).collect();

            let indices: Vec<[u32; 3]> = self
                .mesh
                .indices()
                .iter()
                .map(|idx| [idx[0], idx[1], idx[2]])
                .collect();

            self.mesh =
                TriMesh::new(new_vertices, indices).expect("Failed to create transformed mesh");
        }
    }

    /// Scales the hull geometry uniformly.
    pub fn scale(&mut self, factor: f64) {
        self.scale_xyz(factor, factor, factor);
    }

    /// Scales the hull geometry non-uniformly.
    pub fn scale_xyz(&mut self, sx: f64, sy: f64, sz: f64) {
        let new_vertices: Vec<Point3<f64>> = self
            .mesh
            .vertices()
            .iter()
            .map(|v| Point3::new(v.x * sx, v.y * sy, v.z * sz))
            .collect();

        let indices: Vec<[u32; 3]> = self
            .mesh
            .indices()
            .iter()
            .map(|idx| [idx[0], idx[1], idx[2]])
            .collect();

        self.mesh = TriMesh::new(new_vertices, indices).expect("Failed to create scaled mesh");
    }

    /// Scales the hull to fit within target bounds.
    pub fn scale_to_bounds(&mut self, target: (f64, f64, f64, f64, f64, f64)) {
        let current = self.get_bounds();

        let x_range = current.1 - current.0;
        let y_range = current.3 - current.2;
        let z_range = current.5 - current.4;

        let target_x_range = target.1 - target.0;
        let target_y_range = target.3 - target.2;
        let target_z_range = target.5 - target.4;

        let sx = if x_range > 1e-9 {
            target_x_range / x_range
        } else {
            1.0
        };
        let sy = if y_range > 1e-9 {
            target_y_range / y_range
        } else {
            1.0
        };
        let sz = if z_range > 1e-9 {
            target_z_range / z_range
        } else {
            1.0
        };

        self.scale_xyz(sx, sy, sz);
    }

    /// Simplifies the hull mesh by reducing the number of triangles.
    ///
    /// # Arguments
    /// * `target_count` - Target number of triangles
    pub fn simplify(&mut self, target_count: usize) {
        if target_count >= self.num_triangles() {
            return;
        }

        let vertices = self.mesh.vertices();
        let indices = self.mesh.indices();

        // 1. Flatten vertices for meshopt (f32)
        // meshopt uses f32. We need to convert.
        let vertices_f32: Vec<f32> = vertices
            .iter()
            .flat_map(|v| vec![v.x as f32, v.y as f32, v.z as f32])
            .collect();

        // 2. Flatten indices
        let indices_u32: Vec<u32> = indices
            .iter()
            .flat_map(|tri| vec![tri[0], tri[1], tri[2]])
            .collect();

        // 3. Simplify
        // Target index count = target_triangles * 3
        let target_index_count = target_count * 3;
        let target_error = 1e-3; // Initial error tolerance

        // meshopt::simplify requires VertexDataAdapter with &[u8]
        // Stride is in bytes for f32 (3 * 4 = 12)
        let vertices_bytes = unsafe {
            std::slice::from_raw_parts(
                vertices_f32.as_ptr() as *const u8,
                vertices_f32.len() * std::mem::size_of::<f32>(),
            )
        };

        let vertex_data = meshopt::VertexDataAdapter::new(vertices_bytes, 12, 0)
            .expect("Failed to create vertex data adapter");

        let simplified_indices = meshopt::simplify(
            &indices_u32,
            &vertex_data,
            target_index_count,
            target_error as f32,                  // meshopt uses f32
            meshopt::SimplifyOptions::LockBorder, // Lock borders to avoid gaps between hulls? Or just default?
            None,                                 // result_error
        );

        // 4. Reconstruct mesh
        // Re-indexing is needed because simplify returns indices into original vertex buffer.
        // We want to compact the vertex buffer to remove unused vertices.

        let mut new_vertices = Vec::new();
        let mut new_indices = Vec::new();
        let mut vertex_map = std::collections::HashMap::new();

        for i in (0..simplified_indices.len()).step_by(3) {
            let idx0 = simplified_indices[i];
            let idx1 = simplified_indices[i + 1];
            let idx2 = simplified_indices[i + 2];

            let v0 = vertices[idx0 as usize];
            let v1 = vertices[idx1 as usize];
            let v2 = vertices[idx2 as usize];

            let map_vertex = |v: Point3<f64>,
                              list: &mut Vec<Point3<f64>>,
                              map: &mut std::collections::HashMap<u32, u32>,
                              original_idx: u32|
             -> u32 {
                if let Some(&new_idx) = map.get(&original_idx) {
                    new_idx
                } else {
                    let new_idx = list.len() as u32;
                    list.push(v);
                    map.insert(original_idx, new_idx);
                    new_idx
                }
            };

            let n0 = map_vertex(v0, &mut new_vertices, &mut vertex_map, idx0);
            let n1 = map_vertex(v1, &mut new_vertices, &mut vertex_map, idx1);
            let n2 = map_vertex(v2, &mut new_vertices, &mut vertex_map, idx2);

            new_indices.push([n0, n1, n2]);
        }

        // Update mesh
        if let Ok(new_mesh) = TriMesh::new(new_vertices, new_indices) {
            self.mesh = new_mesh;
        } else {
            eprintln!("Failed to construct simplified mesh");
        }
    }

    /// Returns a simplified copy of the hull.
    pub fn to_simplified(&self, target_count: usize) -> Self {
        let mut clone = self.clone();
        clone.simplify(target_count);
        clone
    }

    /// Exports the hull to an STL file.
    pub fn export_stl<P: AsRef<Path>>(&self, path: P) -> Result<(), HullError> {
        use std::fs::File;
        use std::io::BufWriter;

        let path = path.as_ref();
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Convert to stl_io format
        let vertices: Vec<stl_io::Vertex> = self
            .mesh
            .vertices()
            .iter()
            .map(|v| stl_io::Vertex::new([v.x as f32, v.y as f32, v.z as f32]))
            .collect();

        let triangles: Vec<stl_io::Triangle> = self
            .mesh
            .indices()
            .iter()
            .map(|idx| {
                stl_io::Triangle {
                    normal: stl_io::Normal::new([0.0, 0.0, 1.0]), // Placeholder normal
                    vertices: [
                        vertices[idx[0] as usize],
                        vertices[idx[1] as usize],
                        vertices[idx[2] as usize],
                    ],
                }
            })
            .collect();

        stl_io::write_stl(&mut writer, triangles.iter()).map_err(|e| {
            HullError::IoError(std::io::Error::other(format!("STL write error: {}", e)))
        })?;

        Ok(())
    }

    /// Returns the original file path if the hull was loaded from a file.
    pub fn file_path(&self) -> Option<&str> {
        self.file_path.as_deref()
    }
}

impl std::fmt::Debug for Hull {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let bounds = self.get_bounds();
        f.debug_struct("Hull")
            .field("triangles", &self.num_triangles())
            .field("vertices", &self.num_vertices())
            .field("bounds", &bounds)
            .field("file_path", &self.file_path)
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hull_bounds() {
        // Create a simple box hull
        let vertices = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.0, 5.0, 0.0),
            Point3::new(0.0, 5.0, 0.0),
            Point3::new(0.0, 0.0, 3.0),
            Point3::new(10.0, 0.0, 3.0),
            Point3::new(10.0, 5.0, 3.0),
            Point3::new(0.0, 5.0, 3.0),
        ];

        let indices = vec![
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [2, 3, 7],
            [2, 7, 6],
            [0, 4, 7],
            [0, 7, 3],
            [1, 2, 6],
            [1, 6, 5],
        ];

        let mesh = TriMesh::new(vertices, indices).unwrap();
        let hull = Hull::from_mesh(mesh);

        let bounds = hull.get_bounds();
        assert!((bounds.0 - 0.0).abs() < 1e-6);
        assert!((bounds.1 - 10.0).abs() < 1e-6);
        assert!((bounds.2 - 0.0).abs() < 1e-6);
        assert!((bounds.3 - 5.0).abs() < 1e-6);
        assert!((bounds.4 - 0.0).abs() < 1e-6);
        assert!((bounds.5 - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_hull_scale() {
        let vertices = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
        ];
        let indices = vec![[0, 1, 2]];
        let mesh = TriMesh::new(vertices, indices).unwrap();
        let mut hull = Hull::from_mesh(mesh);

        hull.scale(2.0);

        let bounds = hull.get_bounds();
        assert!((bounds.1 - 2.0).abs() < 1e-6);
        assert!((bounds.3 - 2.0).abs() < 1e-6);
    }
}
