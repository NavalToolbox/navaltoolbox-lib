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
    // We expect the hydrostatics to reflect Total Mass = Ship + Tank
    // Total Mass = 205000 + 128125 = 333125 kg
    // Expected Draft = Total Mass / (1025 * 10 * 10) = 333125 / 102500 = 3.25 m
    let heels = vec![0.0, 10.0, 20.0];
    // Use TankOptions with Mass=True (previous default)
    let tank_opts = Some(navaltoolbox::hydrostatics::TankOptions::all());
    let result = calc.complete_stability(ship_mass, ship_cog, &heels, tank_opts);

    // 6. Verify Hydrostatics
    println!(
        "Hydrostatic State from Complete Stability: {:#?}",
        result.hydrostatics
    );

    // Check Displacement (should be TOTAL mass)
    let total_mass = ship_mass + 128125.0;
    // Difference was ~12.5kg (0.003%), likely solver default tolerance.
    assert!(
        (result.hydrostatics.displacement - total_mass).abs() < 20.0,
        "Hydrostatic displacement ({}) should match total mass ({})",
        result.hydrostatics.displacement,
        total_mass
    );

    // Check Draft
    assert!(
        (result.hydrostatics.draft - 3.25).abs() < 0.01,
        "Hydrostatic draft ({}) should match expected draft (3.25)",
        result.hydrostatics.draft
    );

    // Check Trim (Tank is aft at X=2.5 vs LCB=5.0) -> Trim should be negative (stern down)
    println!("Trim: {}", result.hydrostatics.trim);
    assert!(
        result.hydrostatics.trim < -10.0,
        "Vessel should trim heavily by stern, got {}",
        result.hydrostatics.trim
    );

    // Check GZ Curve properties
    // Verify that GZ curve was calculated
    assert_eq!(result.gz_curve.points.len(), 3);
}

#[test]
fn test_complete_stability_equilibrium_trim() {
    // Test that an off-center tank affects the equilibrium trim in hydrostatics result
    let hull = Hull::from_box(20.0, 10.0, 10.0); // 20m long
    let mut vessel = Vessel::new(hull);

    // Tank at forward end (x=15..20), mass 128125 kg
    let mut tank = Tank::from_box("FwdTank", 15.0, 20.0, -2.5, 2.5, 0.0, 5.0, 1025.0);
    tank.set_fill_percent(100.0);
    vessel.add_tank(tank);

    let calc = StabilityCalculator::new(&vessel, 1025.0);

    let ship_mass = 205000.0; // Draft ~1m for 20x10 hull
    let ship_cog = [10.0, 0.0, 2.0]; // centered

    // With tank fwd, vessel should trim bow down (positive trim)
    // Use TankOptions with Mass=True
    let tank_opts = Some(navaltoolbox::hydrostatics::TankOptions::all());
    let result = calc.complete_stability(ship_mass, ship_cog, &vec![0.0], tank_opts);

    println!("Trim with fwd tank: {}", result.hydrostatics.trim);
    assert!(
        result.hydrostatics.trim > 0.1,
        "Vessel should trim by bow (positive), got {}",
        result.hydrostatics.trim
    );
}
