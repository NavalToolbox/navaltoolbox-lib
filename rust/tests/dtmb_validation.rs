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
    }
}

#[test]
fn test_dtmb_equilibrium_trim_from_lcg_offset() {
    let mut stl_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    stl_path.push("tests/data/dtmb5415.stl");

    let hull = Hull::from_stl(stl_path).expect("Failed to load DTMB 5415 STL");
    let vessel = Vessel::new(hull);
    let calc = HydrostaticsCalculator::new(&vessel, 1025.0);

    // Initial state at even keel
    let draft = 6.15;
    let vcg = 7.555;
    let initial_state = calc
        .from_draft(draft, 0.0, 0.0, Some(vcg))
        .expect("Failed to calculate initial hydrostatics");

    // LCB at T=6.15
    let lcb = initial_state.lcb();
    let target_disp = initial_state.volume * 1025.0;

    // LCG slightly aft of LCB (e.g. -0.5m)
    // Positive trim = Bow down (check convention)
    // If LCG is aft of LCB, stern goes down -> Negative trim?
    // Let's check calculator conventions:
    // "test_equilibrium_trim_from_lcg_offset":
    // // LcG < LCB (aft) → expect negative trim (stern down)
    // // Trim positive = bow down. Stern down = negative trim.

    let lcg_offset = -0.5;
    let lcg = lcb + lcg_offset;
    let cog = [lcg, 0.0, vcg];

    let state = calc
        .from_displacement(target_disp, None, Some(cog), None, None)
        .expect("Calculation failed");

    println!("LCB: {:.3}, LCG: {:.3}, Trim: {:.3}", lcb, lcg, state.trim);

    assert!(
        state.trim < -0.01,
        "Trim should be negative (stern down) for aft LCG, got {}",
        state.trim
    );
}

#[test]
fn test_dtmb_equilibrium_combined() {
    let mut stl_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    stl_path.push("tests/data/dtmb5415.stl");

    let hull = Hull::from_stl(stl_path).expect("Failed to load DTMB 5415 STL");
    let vessel = Vessel::new(hull);
    let calc = HydrostaticsCalculator::new(&vessel, 1025.0);

    // Initial state at even keel
    let draft = 6.15;
    let vcg = 7.555;
    let initial_state = calc
        .from_draft(draft, 0.0, 0.0, Some(vcg))
        .expect("Failed to calculate initial hydrostatics");

    // LCB at T=6.15
    let lcb = initial_state.lcb();
    let tcb = initial_state.tcb(); // Should be ~0.0
    let target_disp = initial_state.volume * 1025.0;

    // Apply combined offsets:
    // 1. LCG slighty aft (-0.5m) -> Expect negative trim (stern down)
    // 2. TCG slightly to Port (+0.1m) -> Expect negative heel (port down)?
    //    Recall: Y+ is Port. Heel+ is Starboard Down (Port Up).
    //    So Port weight -> Port down -> Negative heel.

    let lcg_offset = -0.1;
    let tcg_offset = 0.01;

    let lcg = lcb + lcg_offset;
    let tcg = tcb + tcg_offset;
    let cog = [lcg, tcg, vcg];

    let state = calc
        .from_displacement(target_disp, None, Some(cog), None, None)
        .expect("Calculation failed");

    println!(
        "LCB: {:.3}, TCB: {:.3}, Trim: {:.3}, Heel: {:.3}",
        lcb, tcb, state.trim, state.heel
    );

    assert!(
        state.trim < -0.01,
        "Trim should be negative (stern down) for aft LCG, got {}",
        state.trim
    );
    assert!(
        state.heel < -0.01,
        "Heel should be negative (port down) for port TCG, got {}",
        state.heel
    );
}
