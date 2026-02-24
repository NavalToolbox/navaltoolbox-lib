// ============================================================================
// HSC Code 2000 — Monohull Intact Stability (Annex 8)
// ============================================================================
//
// Reference: High-Speed Craft Code 2000, Annex 8
//
// Parameters that can be provided via `ctx.set_param()`:
// - "passenger_heeling_lever" (in meters)
// - "wind_heeling_lever" (in meters)
//
// Usage from Python:
//   ctx = CriteriaContext.from_result(result, vessel_name="Fast Ferry 1")
//   ctx.set_param("passenger_heeling_lever", 0.08)
//   ctx.set_param("wind_heeling_lever", 0.05)
//   engine = ScriptEngine()
//   result = engine.run_script_file("rules/hsc_2000_monohull_intact.rhai", ctx)
//
// ============================================================================

fn check(ctx) {
    let results = [];
    let plot_id = "gz_plot";

    // =========================================================================
    // Base Parameters
    // =========================================================================
    let gm0_min = 0.15;
    let gz_min_at_30 = 0.20;

    // =========================================================================
    // Limiting angle (min of 40° or flooding angle)
    // =========================================================================
    let limiting_angle = ctx.get_limiting_angle(40.0);

    // =========================================================================
    // 1. Initial GM (GM0)
    // =========================================================================
    let gm0 = ctx.get_gm0();
    if gm0 != () {
        results.push(criterion(
            "GM₀ initial",
            "Initial metacentric height ≥ 0.15 m",
            gm0_min,
            gm0,
            "m"
        ));
    } else {
        results.push(#{
            name: "GM₀ initial",
            status: "FAIL",
            notes: "GM₀ not available"
        });
    }

    // =========================================================================
    // 2. Maximum GZ Angle and GZ Curve Areas
    // =========================================================================
    let max_info = ctx.find_max_gz();
    let max_angle = max_info.angle;
    let max_gz = max_info.value;

    if max_angle != () {
        // Recommend ≥ 25°, required ≥ 15°
        let angle_status = if max_angle >= 15.0 { "PASS" } else { "FAIL" };
        let notes = if max_angle >= 15.0 && max_angle < 25.0 { 
            "Warning: Angle < 25° (preferably ≥ 25°)" 
        } else { "" };
        
        results.push(#{
            name: "Angle of GZmax",
            description: "Angle of maximum righting lever ≥ 15°",
            required: 15.0,
            actual: max_angle,
            unit: "°",
            status: angle_status,
            margin: max_angle - 15.0,
            notes: notes
        });

        // Area criteria depends on GZ max angle
        if max_angle <= 15.0 {
            // Case 1: max angle <= 15°
            let area_15 = ctx.area_under_curve(0.0, 15.0);
            results.push(criterion(
                "Area 0-15°",
                "Area under GZ curve up to 15° ≥ 0.070 m·rad (GZmax ≤ 15°)",
                0.070,
                area_15,
                "m·rad"
            ));
        } else if max_angle >= 30.0 {
            // Case 2: max angle >= 30°
            let area_30 = ctx.area_under_curve(0.0, 30.0);
            results.push(criterion(
                "Area 0-30°",
                "Area under GZ curve up to 30° ≥ 0.055 m·rad",
                0.055,
                area_30,
                "m·rad"
            ));

            let area_lim = ctx.area_under_curve(0.0, limiting_angle);
            results.push(criterion(
                "Area 0-40°(or θf)",
                "Area under GZ curve up to 40° or θf ≥ 0.090 m·rad",
                0.090,
                area_lim,
                "m·rad"
            ));

            if limiting_angle > 30.0 {
                let area_30_lim = ctx.area_under_curve(30.0, limiting_angle);
                results.push(criterion(
                    "Area 30-40°(or θf)",
                    "Area under GZ curve between 30° and 40° or θf ≥ 0.030 m·rad",
                    0.030,
                    area_30_lim,
                    "m·rad"
                ));
            } else {
                results.push(#{
                    name: "Area 30-40°(or θf)",
                    status: "FAIL",
                    notes: "Not applicable - limiting angle ≤ 30°"
                });
            }
        } else {
            // Case 3: 15° < max angle < 30°
            // Interpolate area requirement: 0.055 + 0.001 * (30 - max_angle)
            let req_area = 0.055 + 0.001 * (30.0 - max_angle);
            let area_max = ctx.area_under_curve(0.0, max_angle);
            results.push(criterion(
                "Area 0-θmax",
                "Interpolated area requirement for GZmax angle between 15° and 30°",
                req_area,
                area_max,
                "m·rad"
            ));
        }
    } else {
        results.push(#{
            name: "Angle of GZmax",
            status: "FAIL",
            notes: "Could not find maximum GZ"
        });
    }

    // =========================================================================
    // 3. GZ ≥ 0.20m at θ ≥ 30°
    // =========================================================================
    let gz_30 = ctx.gz_at_angle(30.0);
    results.push(criterion(
        "GZ at 30°+",
        "Righting lever ≥ 0.20 m at an angle of heel ≥ 30°",
        gz_min_at_30,
        gz_30,
        "m"
    ));

    // =========================================================================
    // 4. Heeling Levers (Wind & Passenger Crowding)
    // =========================================================================
    let l_pass = ctx.get_param("passenger_heeling_lever");
    let l_wind = ctx.get_param("wind_heeling_lever");
    
    // Find intersection for heeling angle
    fn find_intersection_angle(ctx, lever) {
        let angles = ctx.get_heels();
        let gz = ctx.get_gz_values();
        let len = angles.len();
        
        for i in 0..len-1 {
            let a1 = angles[i];
            let a2 = angles[i+1];
            let gz1 = gz[i];
            let gz2 = gz[i+1];
            
            if (gz1 <= lever && gz2 >= lever) || (gz1 >= lever && gz2 <= lever) {
                if gz2 == gz1 { return a1; }
                let fraction = (lever - gz1) / (gz2 - gz1);
                return a1 + fraction * (a2 - a1);
            }
        }
        return ();
    }

    // Max combined heeling angle limit = 10° by default 
    // (can be 16° if specific deck features exist, but we assume 10° for standard passenger)
    let max_heel_due_to_lever = 10.0;
    let dei = ctx.get_deck_edge_immersion_angle();
    if dei != () && dei > 0.0 {
        let dei_limit = 0.8 * dei;
        if dei_limit < max_heel_due_to_lever {
            max_heel_due_to_lever = dei_limit;
        }
    }

    if l_pass != () {
        let angle_pass = find_intersection_angle(ctx, l_pass);
        if angle_pass != () {
            results.push(criterion(
                "Passenger Crowding Heel",
                "Equilibrium angle due to passenger crowding ≤ " + max_heel_due_to_lever.toFixed(1) + "°",
                max_heel_due_to_lever,
                angle_pass,
                "°"
            ));
            // Actual margin should be inverted here since angle_pass must be <= max
            let r_len = results.len();
            results[r_len - 1].status = if angle_pass <= max_heel_due_to_lever { "PASS" } else { "FAIL" };
            results[r_len - 1].margin = max_heel_due_to_lever - angle_pass;
        } else {
            results.push(#{
                name: "Passenger Crowding Heel",
                status: "FAIL",
                notes: "No equilibrium found for passenger crowding lever"
            });
        }
    }

    if l_wind != () {
        let angle_wind = find_intersection_angle(ctx, l_wind);
        if angle_wind != () {
            results.push(criterion(
                "Wind Heeling Angle",
                "Equilibrium angle due to wind ≤ " + max_heel_due_to_lever.toFixed(1) + "°",
                max_heel_due_to_lever,
                angle_wind,
                "°"
            ));
            // Actual margin should be inverted
            let r_len = results.len();
            results[r_len - 1].status = if angle_wind <= max_heel_due_to_lever { "PASS" } else { "FAIL" };
            results[r_len - 1].margin = max_heel_due_to_lever - angle_wind;
        } else {
            results.push(#{
                name: "Wind Heeling Angle",
                status: "FAIL",
                notes: "No equilibrium found for wind heeling lever"
            });
        }
    }

    // =========================================================================
    // Overall result
    // =========================================================================
    let fail_count = 0;
    for r in results {
        if r.status == "FAIL" { fail_count += 1; }
        r.plot_id = plot_id;
    }

    // =========================================================================
    // Plot
    // =========================================================================
    let heels = ctx.get_heels();
    let gz_values = ctx.get_gz_values();

    let elements = [
        #{
            type: "Curve",
            name: "GZ",
            x: heels,
            y: gz_values,
            color: "blue"
        }
    ];

    if max_angle != () && max_gz != () {
        elements.push(#{
            type: "Point",
            name: "GZmax",
            x: max_angle,
            y: max_gz,
            marker: "o",
            color: "green"
        });
    }

    elements.push(#{
        type: "HorizontalLine",
        name: "GZ min (0.20m at 30°)",
        y: gz_min_at_30,
        x_min: 30.0,
        color: "red",
        style: "dashed"
    });

    if l_pass != () {
        elements.push(#{
            type: "Curve",
            name: "Passenger Lever",
            x: [heels[0], heels[heels.len()-1]],
            y: [l_pass, l_pass],
            color: "magenta",
            style: "dotted"
        });
    }
    
    if l_wind != () {
        elements.push(#{
            type: "Curve",
            name: "Wind Lever",
            x: [heels[0], heels[heels.len()-1]],
            y: [l_wind, l_wind],
            color: "cyan",
            style: "dotted"
        });
    }

    if limiting_angle < 40.0 {
        elements.push(#{
            type: "VerticalLine",
            name: "Limiting Angle (" + limiting_angle.toFixed(0) + "°)",
            x: limiting_angle,
            color: "orange",
            style: "dashed"
        });
    }

    let plot = #{
        id: plot_id,
        title: "Intact Stability — HSC Code 2000",
        x_label: "Heel (°)",
        y_label: "Righting Lever (m)",
        elements: elements
    };

    #{
        regulation_name: "HSC Code 2000 — Monohull Intact Stability",
        regulation_reference: "High-Speed Craft Code 2000, Annex 8",
        vessel_name: ctx.get_vessel_name(),
        loading_condition: ctx.get_loading_condition(),
        displacement: ctx.get_displacement(),
        cog: ctx.get_cog(),
        criteria: results,
        overall_pass: fail_count == 0,
        notes: "",
        plots: [plot]
    }
}
