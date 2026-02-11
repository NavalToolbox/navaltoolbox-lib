use navaltoolbox::hull::Hull;
use navaltoolbox::stability::StabilityCalculator;
use navaltoolbox::tanks::Tank;
use navaltoolbox::vessel::Vessel;

#[test]
fn test_complete_stability_with_tanks() {
    // 1. Create Vessel (10x10x10 box)
    let hull = Hull::from_box(10.0, 10.0, 10.0);
    let mut vessel = Vessel::new(hull);

    // 2. Add Tank (centered at x=2.5, y=0, z=2.5)
    // 5x5x5 tank.
    // Mass = 5*5*5 * 1025 = 125 * 1025 = 128125 kg
    let mut tank = Tank::from_box(
        "TestTank", 0.0, 5.0, -2.5, 2.5, 0.0, 5.0, 1025.0, // Seawater density inside
    );
    tank.set_fill_percent(100.0);
    vessel.add_tank(tank);

    // 3. Setup Stability Calculator
    let calc = StabilityCalculator::new(&vessel, 1025.0);

    // 4. Define Ship Mass (Lightship + Deadweight excluding tanks)
    // Let's say ship mass is such that total draft without tanks would be 2.0m
    // Disp_ship = 10 * 10 * 2.0 * 1025 = 200 * 1025 = 205000 kg
    let ship_mass = 205000.0;
    let ship_cog = [5.0, 0.0, 2.0]; // Centered

    // 5. Calculate Complete Stability
    let heels = vec![0.0, 10.0, 20.0];
    let tank_opts = Some(navaltoolbox::hydrostatics::TankOptions::all());
    let result = calc.complete_stability(ship_mass, ship_cog, &heels, tank_opts);

    // 6. Verify Hydrostatics
    println!(
        "Hydrostatic State from Complete Stability: {:#?}",
        result.hydrostatics
    );

    let total_mass = ship_mass + 128125.0;
    assert!(
        (result.hydrostatics.displacement - total_mass).abs() < 20.0,
        "Hydrostatic displacement ({}) should match total mass ({})",
        result.hydrostatics.displacement,
        total_mass
    );

    // Check Draft
    assert!(
        (result.hydrostatics.draft - 3.25).abs() < 0.02,
        "Hydrostatic draft ({}) should match expected draft (3.25)",
        result.hydrostatics.draft
    );

    // Check Trim
    assert!(
        result.hydrostatics.trim < -10.0,
        "Vessel should trim heavily by stern, got {}",
        result.hydrostatics.trim
    );

    // Check GZ Curve properties
    assert_eq!(result.gz_curve.points.len(), 3);
}

#[test]
fn test_complete_stability_equilibrium_trim() {
    let hull = Hull::from_box(20.0, 10.0, 10.0);
    let mut vessel = Vessel::new(hull);

    let mut tank = Tank::from_box("FwdTank", 15.0, 20.0, -2.5, 2.5, 0.0, 5.0, 1025.0);
    tank.set_fill_percent(100.0);
    vessel.add_tank(tank);

    let calc = StabilityCalculator::new(&vessel, 1025.0);
    let ship_mass = 205000.0;
    let ship_cog = [10.0, 0.0, 2.0];

    let tank_opts = Some(navaltoolbox::hydrostatics::TankOptions::all());
    let result = calc.complete_stability(ship_mass, ship_cog, &vec![0.0], tank_opts);

    assert!(
        result.hydrostatics.trim > 0.1,
        "Vessel should trim by bow (positive), got {}",
        result.hydrostatics.trim
    );
}

#[test]
fn test_displacement_components() {
    let hull = Hull::from_box(10.0, 10.0, 10.0);
    let mut vessel = Vessel::new(hull);

    // Tank ~125m3 -> ~128t
    let mut tank = Tank::from_box("TestTank", 0.0, 5.0, -2.5, 2.5, 0.0, 5.0, 1025.0);
    tank.set_fill_percent(100.0);
    let tank_mass = tank.fluid_mass();
    vessel.add_tank(tank);

    let calc = StabilityCalculator::new(&vessel, 1025.0);
    let ship_mass = 205000.0; // ~200t
    let ship_cog = [5.0, 0.0, 2.0];

    // CASE 1: Include Mass = TRUE
    let opts_all = Some(navaltoolbox::hydrostatics::TankOptions::all());
    let res_all = calc.complete_stability(ship_mass, ship_cog, &[0.0], opts_all);
    let h_all = res_all.hydrostatics;

    println!(
        "All: Disp={}, Vessel={}, Tank={}",
        h_all.displacement, h_all.vessel_displacement, h_all.tank_displacement
    );

    assert!(
        (h_all.tank_displacement - tank_mass).abs() < 1e-3,
        "Tank displacement incorrect (all)"
    );
    // Vessel displacement should be roughly ship_mass input (though exact match depends on solver convergence)
    // Actually, vessel_displacement is derived: Total - Tank.
    // Total is found by solver to match target (Ship+Tank). So Total ~= Ship+Tank.
    // Thus Vessel = Total - Tank ~= Ship.
    assert!(
        (h_all.vessel_displacement - ship_mass).abs() < 50.0,
        "Vessel displacement mismatch (all): {} vs {}",
        h_all.vessel_displacement,
        ship_mass
    );
    assert!(
        (h_all.displacement - (ship_mass + tank_mass)).abs() < 50.0,
        "Total displacement mismatch (all)"
    );

    // CASE 2: Include Mass = FALSE
    let opts_nomass = Some(navaltoolbox::hydrostatics::TankOptions::fsm_only());
    let res_nm = calc.complete_stability(ship_mass, ship_cog, &[0.0], opts_nomass);
    let h_nm = res_nm.hydrostatics;

    println!(
        "NoMass: Disp={}, Vessel={}, Tank={}",
        h_nm.displacement, h_nm.vessel_displacement, h_nm.tank_displacement
    );

    assert_eq!(
        h_nm.tank_displacement, 0.0,
        "Tank displacement should be 0 when include_mass=false"
    );
    // Vessel displacement = Total - 0 = Total.
    assert!((h_nm.displacement - h_nm.vessel_displacement).abs() < 1e-9);
    // Total should match just ship_mass
    assert!(
        (h_nm.displacement - ship_mass).abs() < 50.0,
        "Total should match ship_mass when tanks excluded"
    );

    // CASE 3: Tank Options = None
    // Should behave like Case 2 (No Mass) regarding displacement, or just no tanks at all?
    // Current logic: None -> tanks ignored completely.
    let res_none = calc.complete_stability(ship_mass, ship_cog, &[0.0], None);
    let h_none = res_none.hydrostatics;

    assert_eq!(h_none.tank_displacement, 0.0);
    assert!((h_none.displacement - ship_mass).abs() < 50.0);
}
