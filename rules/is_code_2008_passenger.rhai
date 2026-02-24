// ============================================================================
// IS Code 2008 Passenger Ships Additional Criteria
// ============================================================================
//
// Reference: IS Code 2008 (MSC.267(85)), Part A, Chapter 3
//
// This script checks the additional criteria for passenger ships:
// - Angle of heel due to passenger crowding ≤ 10°
// - Angle of heel due to turning ≤ 10°
//
// Parameters that must be provided via `ctx.set_param()`:
// - "passenger_heeling_lever" (in meters)
// - "turning_heeling_lever" (in meters)
//
// Usage from Python:
//   ctx = CriteriaContext.from_result(result, vessel_name="MV Example")
//   ctx.set_param("passenger_heeling_lever", 0.15)
//   ctx.set_param("turning_heeling_lever", 0.08)
//   engine = ScriptEngine()
//   result = engine.run_script_file("rules/is_code_2008_passenger.rhai", ctx)
//
// ============================================================================

fn check(ctx) {
    const MAX_HEEL_PASSENGER = 10.0;
    const MAX_HEEL_TURNING = 10.0;

    let results = [];
    let l_pass = ctx.get_param("passenger_heeling_lever");
    let l_turn = ctx.get_param("turning_heeling_lever");

    // Helper: find intersection
    // This is a naive intersection to find angle where GZ = lever
    fn find_intersection_angle(ctx, lever) {
        let angles = ctx.get_heels();
        let gz = ctx.get_gz_values();
        let len = angles.len();
        
        for i in 0..len-1 {
            let a1 = angles[i];
            let a2 = angles[i+1];
            let gz1 = gz[i];
            let gz2 = gz[i+1];
            
            // Check if lever crosses between these two points
            if (gz1 <= lever && gz2 >= lever) || (gz1 >= lever && gz2 <= lever) {
                // Linear interpolation to find angle
                if gz2 == gz1 { return a1; }
                let fraction = (lever - gz1) / (gz2 - gz1);
                return a1 + fraction * (a2 - a1);
            }
        }
        return ();
    }

    if l_pass != () {
        let angle_pass = find_intersection_angle(ctx, l_pass);
        if angle_pass != () {
            results.push(criterion(
                "Passenger Crowding Heel",
                "Angle of heel due to passenger crowding",
                MAX_HEEL_PASSENGER,
                angle_pass,
                "°"
            ));
        } else {
            results.push(#{
                name: "Passenger Crowding Heel",
                status: "FAIL",
                notes: "No equilibrium found for passenger crowding lever"
            });
        }
    } else {
        results.push(#{
            name: "Passenger Crowding Heel",
            status: "FAIL",
            notes: "Missing 'passenger_heeling_lever' parameter"
        });
    }

    if l_turn != () {
        let angle_turn = find_intersection_angle(ctx, l_turn);
        if angle_turn != () {
            results.push(criterion(
                "Turning Heel",
                "Angle of heel due to turning",
                MAX_HEEL_TURNING,
                angle_turn,
                "°"
            ));
        } else {
            results.push(#{
                name: "Turning Heel",
                status: "FAIL",
                notes: "No equilibrium found for turning lever"
            });
        }
    } else {
        results.push(#{
            name: "Turning Heel",
            status: "FAIL",
            notes: "Missing 'turning_heeling_lever' parameter"
        });
    }

    let fail_count = 0;
    for r in results {
        if r.status == "FAIL" { fail_count += 1; }
    }

    // Generate response map
    #{
        regulation_name: "IS Code 2008 Passenger Ships Stability",
        regulation_reference: "IS Code 2008, Part A, Chapter 3",
        vessel_name: ctx.get_vessel_name(),
        loading_condition: ctx.get_loading_condition(),
        displacement: ctx.get_displacement(),
        cog: ctx.get_cog(),
        criteria: results,
        overall_pass: fail_count == 0,
        notes: "",
        plots: []
    }
}
