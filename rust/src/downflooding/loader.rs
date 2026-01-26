// Copyright (C) 2026 Antoine ANCEAU
//
// This file is part of navaltoolbox.

use dxf::entities::EntityType as DxfEntityType;
use dxf::Drawing;
use std::path::Path;
use thiserror::Error;
use vtkio::model::{DataSet, Piece};
use vtkio::Vtk;

use super::{DownfloodingOpening, OpeningType};

#[derive(Error, Debug)]
pub enum OpeningLoadError {
    #[error("Failed to read file: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Failed to parse DXF file: {0}")]
    DxfParseError(#[from] dxf::DxfError),
    #[error("Failed to parse VTK file: {0}")]
    VtkParseError(#[from] vtkio::Error),
    #[error("No valid geometry found in file (Points or Polylines)")]
    NoGeometry,
    #[error("Unsupported file format")]
    UnsupportedFormat,
}

/// Load openings from a DXF file.
/// Points become Point openings.
/// Polylines become Contour openings.
pub fn load_dxf_openings(
    path: &Path,
    default_type: OpeningType,
) -> Result<Vec<DownfloodingOpening>, OpeningLoadError> {
    let drawing = Drawing::load_file(path)?;
    let mut openings = Vec::new();

    let base_name = path
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "opening".to_string());

    let mut counter = 0;

    for entity in drawing.entities() {
        match &entity.specific {
            DxfEntityType::ModelPoint(model_point) => {
                counter += 1;
                let name = format!("{}_{}", base_name, counter);
                openings.push(DownfloodingOpening::from_point(
                    name,
                    [
                        model_point.location.x,
                        model_point.location.y,
                        model_point.location.z,
                    ],
                    default_type.clone(),
                ));
            }
            DxfEntityType::LwPolyline(lwpoly) => {
                let points: Vec<[f64; 3]> = lwpoly
                    .vertices
                    .iter()
                    .map(|v| [v.x, v.y, entity.common.elevation]) // LWPolyline is 2D + elevation
                    .collect();

                if !points.is_empty() {
                    counter += 1;
                    let name = format!("{}_{}", base_name, counter);
                    openings.push(DownfloodingOpening::from_contour(
                        name,
                        points,
                        default_type.clone(),
                    ));
                }
            }
            DxfEntityType::Polyline(poly) => {
                let points: Vec<[f64; 3]> = poly
                    .vertices()
                    .map(|v| [v.location.x, v.location.y, v.location.z])
                    .collect();

                if !points.is_empty() {
                    counter += 1;
                    let name = format!("{}_{}", base_name, counter);
                    openings.push(DownfloodingOpening::from_contour(
                        name,
                        points,
                        default_type.clone(),
                    ));
                }
            }
            _ => continue,
        }
    }

    if openings.is_empty() {
        return Err(OpeningLoadError::NoGeometry);
    }

    Ok(openings)
}

/// Load openings from a VTK file (.vtk, .vtp).
/// Assumes generic Points/Polys.
pub fn load_vtk_openings(
    path: &Path,
    default_type: OpeningType,
) -> Result<Vec<DownfloodingOpening>, OpeningLoadError> {
    let vtk = Vtk::import(path)?;
    let mut openings = Vec::new();

    let base_name = path
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "opening".to_string());

    // Generic extractor
    fn process_dataset(
        data: &DataSet,
        base_name: &str,
        default_type: &OpeningType,
    ) -> Result<Vec<DownfloodingOpening>, OpeningLoadError> {
        let mut results = Vec::new();
        match data {
            DataSet::PolyData { pieces, .. } => {
                for piece in pieces {
                    if let Piece::Inline(polydata) = piece {
                        // Check for Polylines/Lines
                        // polydata.lines and polydata.polys are Option<Box<Topology>>.
                        // Topology has num_cells().
                        let has_lines = polydata.lines.as_ref().map_or(0, |t| t.num_cells()) > 0
                            || polydata.polys.as_ref().map_or(0, |t| t.num_cells()) > 0;

                        let points = extract_points(&polydata.points)?;

                        if has_lines {
                            // Treat entire point set as one contour for now
                            results.push(DownfloodingOpening::from_contour(
                                base_name.to_string(),
                                points,
                                default_type.clone(),
                            ));
                        } else {
                            // Treat as individual points
                            for (i, p) in points.iter().enumerate() {
                                results.push(DownfloodingOpening::from_point(
                                    format!("{}_{}", base_name, i + 1),
                                    *p,
                                    default_type.clone(),
                                ));
                            }
                        }
                    }
                }
            }
            DataSet::UnstructuredGrid { pieces, .. } => {
                for piece in pieces {
                    if let Piece::Inline(grid) = piece {
                        let points = extract_points(&grid.points)?;
                        // Assume points
                        for (i, p) in points.iter().enumerate() {
                            results.push(DownfloodingOpening::from_point(
                                format!("{}_{}", base_name, i + 1),
                                *p,
                                default_type.clone(),
                            ));
                        }
                    }
                }
            }
            _ => return Err(OpeningLoadError::UnsupportedFormat),
        }
        Ok(results)
    }

    match &vtk.data {
        DataSet::PolyData { .. } | DataSet::UnstructuredGrid { .. } => {
            openings.extend(process_dataset(&vtk.data, &base_name, &default_type)?);
        }
        _ => return Err(OpeningLoadError::UnsupportedFormat),
    }

    if openings.is_empty() {
        return Err(OpeningLoadError::NoGeometry);
    }

    Ok(openings)
}

fn extract_points(buffer: &vtkio::IOBuffer) -> Result<Vec<[f64; 3]>, OpeningLoadError> {
    match buffer {
        vtkio::IOBuffer::F64(data) => Ok(data
            .chunks(3)
            .filter(|c| c.len() == 3)
            .map(|c| [c[0], c[1], c[2]])
            .collect()),
        vtkio::IOBuffer::F32(data) => Ok(data
            .chunks(3)
            .filter(|c| c.len() == 3)
            .map(|c| [c[0] as f64, c[1] as f64, c[2] as f64])
            .collect()),
        _ => Err(OpeningLoadError::NoGeometry),
    }
}

#[cfg(test)]
mod tests {
    use super::super::OpeningGeometry;
    use super::*;
    use dxf::entities::{Entity, EntityType, ModelPoint, Polyline, Vertex};
    use dxf::Point;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn test_load_dxf_points() -> Result<(), Box<dyn std::error::Error>> {
        let dir = tempdir()?;
        let file_path = dir.path().join("points.dxf");

        let mut drawing = Drawing::new();

        // Add two points representing 2 openings
        let p1 = ModelPoint::new(Point::new(10.0, 5.0, 2.0));
        let p2 = ModelPoint::new(Point::new(20.0, 6.0, 3.0));

        drawing.add_entity(Entity::new(EntityType::ModelPoint(p1)));
        drawing.add_entity(Entity::new(EntityType::ModelPoint(p2)));

        drawing.save_file(&file_path)?;

        let openings = load_dxf_openings(&file_path, OpeningType::Vent)?;

        assert_eq!(openings.len(), 2);

        // Validate geometry
        let mut found_p1 = false;
        let mut found_p2 = false;

        for op in openings {
            if let OpeningGeometry::Point(p) = op.geometry {
                if (p[0] - 10.0).abs() < 1e-6
                    && (p[1] - 5.0).abs() < 1e-6
                    && (p[2] - 2.0).abs() < 1e-6
                {
                    found_p1 = true;
                } else if (p[0] - 20.0).abs() < 1e-6
                    && (p[1] - 6.0).abs() < 1e-6
                    && (p[2] - 3.0).abs() < 1e-6
                {
                    found_p2 = true;
                }
            }
        }
        assert!(found_p1, "Point 1 not found");
        assert!(found_p2, "Point 2 not found");

        Ok(())
    }

    #[test]
    fn test_load_dxf_polyline() -> Result<(), Box<dyn std::error::Error>> {
        let dir = tempdir()?;
        let file_path = dir.path().join("polyline.dxf");

        let mut drawing = Drawing::new();

        let mut poly = Polyline::default();
        poly.add_vertex(&mut drawing, Vertex::new(Point::new(0.0, 0.0, 0.0)));
        poly.add_vertex(&mut drawing, Vertex::new(Point::new(10.0, 0.0, 0.0)));
        poly.add_vertex(&mut drawing, Vertex::new(Point::new(10.0, 10.0, 2.0))); // Elevation 2.0 via Z

        drawing.add_entity(Entity::new(EntityType::Polyline(poly)));
        drawing.save_file(&file_path)?;

        let openings = load_dxf_openings(&file_path, OpeningType::Hatch)?;

        assert_eq!(openings.len(), 1);
        match &openings[0].geometry {
            OpeningGeometry::Contour(points) => {
                assert_eq!(points.len(), 3);
                assert!((points[2][0] - 10.0).abs() < 1e-6);
                assert!((points[2][1] - 10.0).abs() < 1e-6);
                assert!((points[2][2] - 2.0).abs() < 1e-6);
            }
            _ => panic!("Expected Contour geometry"),
        }

        Ok(())
    }

    #[test]
    fn test_load_vtk_points() -> Result<(), Box<dyn std::error::Error>> {
        let dir = tempdir()?;
        let file_path = dir.path().join("points.vtk");

        // Simple ASCII Legacy VTK
        let vtk_content = "\
# vtk DataFile Version 3.0
Test Points
ASCII
DATASET POLYDATA
POINTS 2 double
0.0 0.0 1.0
10.0 5.0 2.0
VERTICES 2 4
1 0
1 1
";
        {
            let mut file = File::create(&file_path)?;
            file.write_all(vtk_content.as_bytes())?;
        }

        let openings = load_vtk_openings(&file_path, OpeningType::Other("ChainLocker".into()))?;

        assert_eq!(openings.len(), 2);
        // Verify one of them
        let mut found = false;
        for op in openings {
            if let OpeningGeometry::Point(p) = op.geometry {
                if (p[0] - 10.0).abs() < 1e-6
                    && (p[1] - 5.0).abs() < 1e-6
                    && (p[2] - 2.0).abs() < 1e-6
                {
                    found = true;
                }
            }
        }
        assert!(found, "Point (10,5,2) not found in VTK");

        Ok(())
    }

    #[test]
    fn test_load_vtk_contour() -> Result<(), Box<dyn std::error::Error>> {
        let dir = tempdir()?;
        let file_path = dir.path().join("contour.vtk");

        // Polyline in VTK (LINES)
        let vtk_content = "\
# vtk DataFile Version 3.0
Test Contour
ASCII
DATASET POLYDATA
POINTS 3 double
0.0 0.0 0.0
5.0 0.0 0.0
5.0 5.0 1.0
LINES 1 4
3 0 1 2
";
        {
            let mut file = File::create(&file_path)?;
            file.write_all(vtk_content.as_bytes())?;
        }

        let openings = load_vtk_openings(&file_path, OpeningType::Other("FwTurn".into()))?;

        // Should retrieve 1 contour
        assert_eq!(openings.len(), 1);
        match &openings[0].geometry {
            OpeningGeometry::Contour(points) => {
                assert_eq!(points.len(), 3);
                assert!((points[2][0] - 5.0).abs() < 1e-6);
                assert!((points[2][1] - 5.0).abs() < 1e-6);
                assert!((points[2][2] - 1.0).abs() < 1e-6);
            }
            _ => panic!("Expected Contour geometry"),
        }

        Ok(())
    }
}
