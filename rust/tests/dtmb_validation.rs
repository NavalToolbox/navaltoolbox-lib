use navaltoolbox::hull::Hull;
use navaltoolbox::hydrostatics::HydrostaticsCalculator;
use navaltoolbox::vessel::Vessel;
use std::path::PathBuf;

#[test]
fn test_dtmb5415_simman_validation() {
    let mut stl_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    stl_path.push("tests/data/dtmb5415.stl");

    let hull = Hull::from_stl(stl_path).expect("Failed to load DTMB 5415 STL");
    let vessel = Vessel::new(hull);
    let calc = HydrostaticsCalculator::new(&vessel, 1025.0);

    // Target conditions
    let draft = 6.15;
    let vcg = 7.555;

    let state = calc
        .from_draft(draft, 0.0, 0.0, Some(vcg))
        .expect("Failed to calculate hydrostatics");

    println!("Hydrostatic State at T={}: {:#?}", draft, state);

    // 1. Beam at Waterline (Bwl) - Ref: 19.06 m
    assert!(
        (state.bwl - 19.06).abs() < 0.01,
        "Bwl mismatch: got {}, expected ~19.06",
        state.bwl
    );

    // 2. Length at Waterline (Lwl) - Ref: 142.18 m
    assert!(
        (state.lwl - 142.18).abs() < 0.1,
        "Lwl mismatch: got {}, expected ~142.18",
        state.lwl
    );

    // 2b. Length Overall Submerged (LOS)
    // Should be >= Lwl. The bulbous bow might extend slightly.
    println!("LOS: {}", state.los);
    assert!(state.los >= state.lwl, "LOS should be >= Lwl");
    assert!((state.los - 142.0).abs() < 2.0, "LOS reasonable check");

    // 3. Displacement Volume - Ref: 8424 m3 (Bare) / 8635 (Appended? ~8424 * 1.025?)
    // Actually 8635 MT / 1.025 = 8424 m3.
    assert!(
        (state.volume - 8424.0).abs() < 80.0,
        "Volume mismatch: got {}, expected ~8424",
        state.volume
    );

    // 4. Block Coefficient (CB) - Ref: 0.506
    assert!(
        (state.cb - 0.506).abs() < 0.02,
        "CB mismatch: got {}, expected ~0.506",
        state.cb
    );

    // 5. Midship Area Coefficient (Cm) - Ref: 0.816
    // Note: Python test checks Cm at x=72.0 specifically. Rust state.cm uses auto-detected mid-point.
    // Bounds of DTMB roughly 0 to 142. Mid ~71.
    assert!(
        (state.cm - 0.816).abs() < 0.05,
        "Cm mismatch: got {}, expected ~0.816",
        state.cm
    );

    // 6. Wetted Surface Area (S) - Ref: 2972.6 m2
    assert!(
        (state.wetted_surface_area - 2972.6).abs() < 20.0,
        "S mismatch: got {}, expected ~2973",
        state.wetted_surface_area
    );

    // 7. GMt - Ref: 1.95 m
    // GMt = KB + BMt - KG
    // KB = VCB
    // BMt = I_t / Vol
    // KG = 7.555
    if let Some(gmt) = state.gmt {
        assert!(
            (gmt - 1.95).abs() < 0.02,
            "GMt mismatch: got {}, expected ~1.95",
            gmt
        );
    } else {
        panic!("GMt was None");
    }
}
