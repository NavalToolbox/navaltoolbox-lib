use navaltoolbox::hull::Hull;
use navaltoolbox::hydrostatics::HydrostaticsCalculator;
use navaltoolbox::vessel::Vessel;

#[test]
fn test_box_hydrostatics() {
    // 10x10x10 box
    let hull = Hull::from_box(10.0, 10.0, 10.0);
    let vessel = Vessel::new(hull);
    let calc = HydrostaticsCalculator::new(&vessel, 1025.0);

    let draft = 5.0;

    // We calculate at draft 5.0, level trim/heel
    // Since VCG is optional, we pass None here for basic geo checks,
    // but the python test passed None explicitly or implicitly.
    // Python test:
    // state = self.calc.from_draft(draft, 0.0, 0.0)
    // Implicit VCG=None.

    let state = calc
        .from_draft(draft, 0.0, 0.0, None)
        .expect("Failed to calculate box hydrostatics");

    println!("Box State at T=5.0: {:#?}", state);

    // 1. Wetted Surface Area
    // Bottom: 10 * 10 = 100
    // Sides: 4 sides * (10 * 5) = 200
    // Total = 300
    assert!(
        (state.wetted_surface_area - 300.0).abs() < 1.0,
        "Wetted Surface mismatch: got {}, expected 300.0",
        state.wetted_surface_area
    );

    // 2. Midship Area
    // 10 width * 5 draft = 50
    assert!(
        (state.midship_area - 50.0).abs() < 0.5,
        "Midship Area mismatch: got {}, expected 50.0",
        state.midship_area
    );

    // 3. Coefficients
    // Volume = 10 * 10 * 5 = 500
    // Lwl = 10
    // Bwl = 10
    // T = 5
    // Cb = Vol / (L * B * T) = 500 / 500 = 1.0
    // Cm = Am / (B * T) = 50 / 50 = 1.0
    // Cp = Vol / (Am * L) = 500 / (50 * 10) = 1.0

    assert!(
        (state.cb - 1.0).abs() < 0.01,
        "Cb mismatch: got {}, expected 1.0",
        state.cb
    );
    assert!(
        (state.cm - 1.0).abs() < 0.01,
        "Cm mismatch: got {}, expected 1.0",
        state.cm
    );
    assert!(
        (state.cp - 1.0).abs() < 0.01,
        "Cp mismatch: got {}, expected 1.0",
        state.cp
    );

    // 4. Stiffness Matrix
    // C33 (Heave) = rho * g * Awp
    // Awp = 10 * 10 = 100
    // C33 = 1025 * 9.81 * 100 = 1005525.0
    let expected_c33 = 1025.0 * 9.81 * 100.0;
    let c33 = state.stiffness_matrix[14]; // Row 2, Col 2 (index 14 flattened?)
                                          // Matrix is 6x6 flattened.
                                          // Row 0, 1, 2 (Heave). Col 0, 1, 2.
                                          // Index = row * 6 + col = 2 * 6 + 2 = 14. Correct.

    assert!(
        (c33 - expected_c33).abs() < 100.0,
        "C33 stiffness mismatch: got {}, expected {}",
        c33,
        expected_c33
    );

    // 5. LOS
    // Should be 10.0
    assert!(
        (state.los - 10.0).abs() < 0.1,
        "LOS mismatch: got {}, expected 10.0",
        state.los
    );
}

#[test]
fn test_box_displacement_calculation() {
    let hull = Hull::from_box(10.0, 10.0, 10.0);
    let vessel = Vessel::new(hull);
    let calc = HydrostaticsCalculator::new(&vessel, 1025.0);

    // Target displacement for draft 5.0
    // Vol = 500
    let target_disp = 500.0 * 1025.0;

    // Calculate at displacement
    let state = calc
        .from_displacement(target_disp, None, None, None, None)
        .expect("Failed to calculate at displacement");

    // Check draft
    assert!(
        (state.draft - 5.0).abs() < 0.01,
        "Draft mismatch for displacement: got {}, expected 5.0",
        state.draft
    );

    // Check properties propagate
    assert!(
        (state.wetted_surface_area - 300.0).abs() < 1.0,
        "WSA mismatch in displacement calc: got {}",
        state.wetted_surface_area
    );
    assert!(
        (state.cb - 1.0).abs() < 0.01,
        "Cb mismatch in displacement calc: got {}",
        state.cb
    );
}
