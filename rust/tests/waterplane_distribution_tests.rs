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

use navaltoolbox::hull::Hull;
use navaltoolbox::hydrostatics::HydrostaticsCalculator;
use navaltoolbox::mesh::calculate_waterplane_properties;
use navaltoolbox::vessel::Vessel;
use parry3d_f64::math::Point;
use parry3d_f64::shape::TriMesh;

fn create_box_hull(loa: f64, boa: f64, depth: f64) -> Hull {
    let hb = boa / 2.0;
    let vertices = vec![
        Point::new(0.0, -hb, 0.0),
        Point::new(loa, -hb, 0.0),
        Point::new(loa, hb, 0.0),
        Point::new(0.0, hb, 0.0),
        Point::new(0.0, -hb, depth),
        Point::new(loa, -hb, depth),
        Point::new(loa, hb, depth),
        Point::new(0.0, hb, depth),
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
    Hull::from_mesh(mesh)
}

fn create_u_shape_mesh(loa: f64, boa: f64, depth: f64, thickness: f64) -> TriMesh {
    let _arm_width = thickness;

    let hb = boa / 2.0;
    let th = thickness;

    let z0 = 0.0;
    let z1 = depth;

    let profile = vec![
        Point::new(0.0, -hb, 0.0),      // 0
        Point::new(loa, -hb, 0.0),      // 1
        Point::new(loa, -hb + th, 0.0), // 2
        Point::new(th, -hb + th, 0.0),  // 3
        Point::new(th, hb - th, 0.0),   // 4
        Point::new(loa, hb - th, 0.0),  // 5
        Point::new(loa, hb, 0.0),       // 6
        Point::new(0.0, hb, 0.0),       // 7
    ];

    let mut vertices = Vec::new();
    // Bottom vertices (0-7)
    for p in &profile {
        vertices.push(Point::new(p.x, p.y, z0));
    }
    // Top vertices (8-15)
    for p in &profile {
        vertices.push(Point::new(p.x, p.y, z1));
    }

    let mut indices = vec![];

    // Bottom (reverse winding for -Z normal)
    indices.extend_from_slice(&[[0, 2, 1], [0, 3, 2]]);
    indices.extend_from_slice(&[[0, 4, 3], [0, 7, 4]]);
    indices.extend_from_slice(&[[4, 6, 5], [4, 7, 6]]);

    // Top (forward winding for +Z normal) offset by 8
    indices.extend_from_slice(&[[8, 9, 10], [8, 10, 11]]);
    indices.extend_from_slice(&[[8, 11, 12], [8, 12, 15]]);
    indices.extend_from_slice(&[[12, 13, 14], [12, 14, 15]]);

    // Sides
    let n = 8;
    for i in 0..n {
        let next = (i + 1) % n;
        indices.push([i as u32, next as u32, (next + 8) as u32]);
        indices.push([i as u32, (next + 8) as u32, (i + 8) as u32]);
    }

    TriMesh::new(vertices, indices).unwrap()
}

#[test]
fn test_u_shape_waterplane() {
    // U-shape: LOA 10, BOA 10, Depth 10, Thickness 2
    let mesh = create_u_shape_mesh(10.0, 10.0, 10.0, 2.0);
    let draft = 5.0;
    let wp = calculate_waterplane_properties(&mesh, draft).unwrap();

    // Area = 10*2 + 2*6 + 10*2 = 20 + 12 + 20 = 52.
    let expected_area = 52.0;
    assert!(
        (wp.area - expected_area).abs() < 1.0,
        "U-shape Area: got {:.2}, expected {:.2}",
        wp.area,
        expected_area
    );

    // LCF check.
    // Arm 1 Center: (5, -4), Area 20. Moment X = 100
    // Base Center: (1, 0), Area 12. Moment X = 12
    // Arm 2 Center: (5, 4), Area 20. Moment X = 100
    // Total Moment X = 212.
    // LCF = 212 / 52 = 4.077
    let expected_lcf = 212.0 / 52.0;
    assert!(
        (wp.centroid[0] - expected_lcf).abs() < 0.1,
        "U-shape LCF: got {:.3}, expected {:.3}",
        wp.centroid[0],
        expected_lcf
    );
}

#[test]
fn test_box_waterplane_area() {
    // 10m × 10m × 10m box
    let hull = create_box_hull(10.0, 10.0, 10.0);
    let mesh = hull.mesh();

    let draft = 5.0;
    let wp = calculate_waterplane_properties(mesh, draft).unwrap();

    // Expected area = L × B = 100 m²
    let expected_area = 100.0;
    assert!(
        (wp.area - expected_area).abs() < 1.0,
        "Area: got {:.2}, expected {:.2}",
        wp.area,
        expected_area
    );
}

#[test]
fn test_box_waterplane_centroid() {
    let hull = create_box_hull(10.0, 10.0, 10.0);
    let mesh = hull.mesh();

    let draft = 5.0;
    let wp = calculate_waterplane_properties(mesh, draft).unwrap();

    // Expected LCF = 5.0 (center at x=5)
    // Expected TCF = 0.0 (symmetric)
    assert!(
        (wp.centroid[0] - 5.0).abs() < 0.1,
        "LCF: got {:.2}, expected 5.0",
        wp.centroid[0]
    );
    assert!(
        wp.centroid[1].abs() < 0.1,
        "TCF: got {:.2}, expected 0.0",
        wp.centroid[1]
    );
}

#[test]
fn test_box_metacentric_heights() {
    let hull = create_box_hull(10.0, 10.0, 10.0);
    let vessel = Vessel::new(hull);
    let calc = HydrostaticsCalculator::new(&vessel, 1025.0);

    let draft = 5.0;
    let vcg = 7.0;

    let state = calc
        .from_draft(draft, 0.0, 0.0, Some(vcg), None, None, None, None)
        .unwrap();

    // Theoretical BM_t = B² / (12×T) = 100 / 60 = 1.667 m
    let expected_bmt = 10.0 * 10.0 / (12.0 * 5.0);
    assert!(
        (state.bmt - expected_bmt).abs() < 0.1,
        "BMt: got {:.3}, expected {:.3}",
        state.bmt,
        expected_bmt
    );

    // BM_l = L² / (12×T) = 100 / 60 = 1.667 m (same for square)
    assert!(
        (state.bml - expected_bmt).abs() < 0.1,
        "BMl: got {:.3}, expected {:.3}",
        state.bml,
        expected_bmt
    );

    // VCB = T/2 = 2.5 m
    let expected_vcb = 2.5;
    assert!(
        (state.vcb() - expected_vcb).abs() < 0.1,
        "VCB: got {:.3}, expected {:.3}",
        state.vcb(),
        expected_vcb
    );

    // KM_t = VCB + BM_t = 2.5 + 1.667 = 4.167 m
    // GMT = KM_t - VCG = 4.167 - 7.0 = -2.833 m (unstable)
    let expected_kmt = expected_vcb + expected_bmt;
    let expected_gmt = expected_kmt - vcg;

    assert!(
        (state.gmt.unwrap() - expected_gmt).abs() < 0.1,
        "GMT: got {:.3}, expected {:.3}",
        state.gmt.unwrap(),
        expected_gmt
    );
}

#[test]
fn test_catamaran_two_disjoint_boxes() {
    // Two 10m × 5m × 10m boxes separated by 5m
    // Total waterplane area = 2 × (10 × 5) = 100 m²
    let hull1 = create_box_hull(10.0, 5.0, 10.0);
    let mut hull2 = create_box_hull(10.0, 5.0, 10.0);

    // Shift hull2 by 10m in y direction (5m gap between)
    hull2.transform((0.0, 10.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0));

    let vessel = Vessel::new_multi(vec![hull1, hull2]).unwrap();

    let calc = HydrostaticsCalculator::new(&vessel, 1025.0);
    let state = calc
        .from_draft(5.0, 0.0, 0.0, None, None, None, None, None)
        .unwrap();

    // Expected total area = 100 m²
    assert!(
        (state.waterplane_area - 100.0).abs() < 2.0,
        "Catamaran waterplane area: got {:.2}, expected 100.0",
        state.waterplane_area
    );

    // Expected TCF = 5.0 (midpoint between two hulls)
    // Hull1 centroid: (5.0, 0.0)
    // Hull2 centroid: (5.0, 10.0)
    // Combined: (5.0, 5.0)
    assert!(
        (state.cob[1] - 5.0).abs() < 0.5,
        "Catamaran TCB: got {:.2}, expected 5.0",
        state.cob[1]
    );
}

#[test]
fn test_two_joined_boxes() {
    // Two 5m × 10m × 10m boxes joined side-by-side
    // Forms a 10m × 10m waterplane
    let hull1 = create_box_hull(5.0, 10.0, 10.0);
    let mut hull2 = create_box_hull(5.0, 10.0, 10.0);

    // Shift hull2 by 5m in x direction (joined at x=5)
    hull2.transform((5.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0));

    let vessel = Vessel::new_multi(vec![hull1, hull2]).unwrap();

    let calc = HydrostaticsCalculator::new(&vessel, 1025.0);
    let state = calc
        .from_draft(5.0, 0.0, 0.0, None, None, None, None, None)
        .unwrap();

    // Expected area = 10 × 10 = 100 m²
    assert!(
        (state.waterplane_area - 100.0).abs() < 2.0,
        "Joined boxes waterplane area: got {:.2}, expected 100.0",
        state.waterplane_area
    );

    // Expected LCF = 5.0 (midpoint)
    assert!(
        (state.cob[0] - 5.0).abs() < 0.3,
        "Joined boxes LCB: got {:.2}, expected 5.0",
        state.cob[0]
    );
}

#[test]
fn test_tank_fsm_modes() {
    use navaltoolbox::tanks::{FSMMode, Tank};

    // Create a 10x10x10 box tank
    let hull = create_box_hull(10.0, 10.0, 10.0);
    let mesh = hull.mesh();
    let mut tank = Tank::new("TestTank", mesh.clone(), 1025.0);

    // 1. Test Actual Mode (Default)
    tank.set_fill_level(0.5);
    // I_t for 10x10 box = 10*10^3/12 = 833.333
    let fsm_actual = tank.free_surface_moment_t();
    assert!(
        (fsm_actual - 833.333).abs() < 1.0,
        "Actual FSM incorrect: {}",
        fsm_actual
    );

    // Fill 100% -> 0
    tank.set_fill_level(1.0);
    assert_eq!(tank.free_surface_moment_t(), 0.0);

    // 2. Test Fixed Mode
    tank.set_fsm_mode(FSMMode::Fixed {
        t: 500.0,
        l: 1000.0,
    });
    tank.set_fill_level(0.5);
    assert_eq!(tank.free_surface_moment_t(), 500.0);
    assert_eq!(tank.free_surface_moment_l(), 1000.0);

    // 3. Test Maximum Mode
    tank.set_fsm_mode(FSMMode::Maximum);
    // For a box, Max FSM is constant 833.33.
    let fsm_max = tank.free_surface_moment_t();
    assert!(
        (fsm_max - 833.333).abs() < 10.0,
        "Max FSM incorrect: {}",
        fsm_max
    );

    // At very low fill, Actual is small, Max is large
    tank.set_fill_level(0.01);
    let fsm_low = tank.free_surface_moment_t();
    assert!(
        (fsm_low - 833.333).abs() < 10.0,
        "Max FSM at low fill incorrect: {}",
        fsm_low
    );
}
