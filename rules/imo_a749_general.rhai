// ============================================================================
// IMO A.749(18) General Intact Stability Criteria
// ============================================================================
//
// Reference: IMO Resolution A.749(18), Section 3.1
//
// This script checks the following criteria:
// - Area 0-30°: ≥ 0.055 m·rad
// - Area 0-40° (or θf): ≥ 0.090 m·rad
// - Area 30-40° (or θf): ≥ 0.030 m·rad
// - GZ at 30°+: ≥ 0.200 m
// - Angle of GZ max: ≥ 25° (preferably > 30°)
// - GM₀: ≥ 0.15 m
//
// Usage from Python:
//   ctx = CriteriaContext.from_result(result, vessel_name="MV Example")
//   engine = ScriptEngine()
//   result = engine.run_script_file("rules/imo_a749_general.rhai", ctx)
//
// ============================================================================

// Threshold values
// Threshold values moved inside check()




// Main check function - called by ScriptEngine
fn check(ctx) {
    // Threshold values
    const AREA_0_30_MIN = 0.055;      // m·rad
    const AREA_0_40_MIN = 0.090;      // m·rad
    const AREA_30_40_MIN = 0.030;     // m·rad
    const GZ_AT_30_MIN = 0.200;       // m
    const GZ_MAX_ANGLE_MIN = 25.0;    // degrees
    const GZ_MAX_ANGLE_PREF = 30.0;   // degrees (preferred)
    const GM0_MIN = 0.15;             // m
    let results = [];
    
    // Get limiting angle (40° or flooding angle, whichever is smaller)
    let limiting_angle = ctx.get_limiting_angle(40.0);
    
    // ========================================================================
    // Criterion 1: Area under GZ curve from 0° to 30°
    // ========================================================================
    let area_0_30 = ctx.area_under_curve(0.0, 30.0);
    results.push(criterion(
        "Area 0-30°",
        "Area under the righting lever curve up to 30°",
        AREA_0_30_MIN,
        area_0_30,
        "m·rad"
    ));
    
    // ========================================================================
    // Criterion 2: Area under GZ curve from 0° to 40° (or θf)
    // ========================================================================
    let to_angle = if limiting_angle < 40.0 { limiting_angle } else { 40.0 };
    let area_0_40 = ctx.area_under_curve(0.0, to_angle);
    let crit_0_40 = criterion(
        "Area 0-" + to_angle + "°",
        "Area under curve up to " + to_angle + "° or θf",
        AREA_0_40_MIN,
        area_0_40,
        "m·rad"
    );
    if limiting_angle < 40.0 {
        crit_0_40.notes = "Calculated to θf = " + limiting_angle + "°";
    }
    results.push(crit_0_40);
    
    // ========================================================================
    // Criterion 3: Area under GZ curve from 30° to 40° (or θf)
    // ========================================================================
    if to_angle > 30.0 {
        let area_30_40 = ctx.area_under_curve(30.0, to_angle);
        let crit_30_40 = criterion(
            "Area 30-" + to_angle + "°",
            "Area from 30° to " + to_angle + "°",
            AREA_30_40_MIN,
            area_30_40,
            "m·rad"
        );
        if limiting_angle < 40.0 {
            crit_30_40.notes = "Calculated to θf = " + limiting_angle + "°";
        }
        results.push(crit_30_40);
    } else {
        results.push(#{
            name: "Area 30-40°",
            description: "Not applicable - θf ≤ 30°",
            required: AREA_30_40_MIN,
            actual: 0.0,
            unit: "m·rad",
            status: "N/A",
            margin: 0.0,
            notes: "θf (" + limiting_angle + "°) ≤ 30°"
        });
    }
    
    // ========================================================================
    // Criterion 4: GZ at 30° or above ≥ 0.20 m
    // ========================================================================
    let gz_at_30 = ctx.gz_at_angle(30.0);
    let max_gz_info = ctx.find_max_gz();
    let gz_30_plus = if max_gz_info.angle >= 30.0 && max_gz_info.value > gz_at_30 {
        max_gz_info.value
    } else {
        gz_at_30
    };
    results.push(criterion(
        "GZ at 30°+",
        "Maximum GZ at angle ≥ 30°",
        GZ_AT_30_MIN,
        gz_30_plus,
        "m"
    ));
    
    // ========================================================================
    // Criterion 5: Angle of maximum GZ ≥ 25°
    // ========================================================================
    let max_angle = max_gz_info.angle;
    let max_value = max_gz_info.value;
    let crit_angle = criterion(
        "GZ max angle",
        "Angle of maximum GZ (= " + max_value + " m)",
        GZ_MAX_ANGLE_MIN,
        max_angle,
        "°"
    );
    if max_angle >= GZ_MAX_ANGLE_PREF {
        crit_angle.notes = "Exceeds preferred 30°";
    } else if max_angle >= GZ_MAX_ANGLE_MIN {
        crit_angle.notes = "Meets 25° but below preferred 30°";
    }
    results.push(crit_angle);
    
    // ========================================================================
    // Criterion 6: GM₀ ≥ 0.15 m
    // ========================================================================
    let gm0 = ctx.get_gm0();
    if gm0 != () {
        results.push(criterion(
            "Initial GM₀",
            "Initial metacentric height",
            GM0_MIN,
            gm0,
            "m"
        ));
    } else {
        results.push(#{
            name: "Initial GM₀",
            description: "GM₀ not available (no VCG specified)",
            required: GM0_MIN,
            actual: 0.0,
            unit: "m",
            status: "N/A",
            margin: 0.0
        });
    }
    
    // ========================================================================
    // Build final result
    // ========================================================================
    // Build GZ curve plot data
    let plot_id = "gz_curve";
    let pass_count = 0;
    let fail_count = 0;
    for r in results {
        if r.status == "PASS" { pass_count += 1; }
        if r.status == "FAIL" { fail_count += 1; }
        // Link to plot
        r.plot_id = plot_id;
    }
    
    
    let plot = #{
        id: plot_id,
        title: "GZ Curve - IMO A.749",
        x_label: "Heel (°)",
        y_label: "GZ (m)",
        elements: [
            #{
                type: "Curve",
                name: "GZ",
                x: ctx.get_heels(),
                y: ctx.get_gz_values(),
                color: "blue",
                style: "solid"
            },
            #{
                type: "Point",
                name: "GZ max",
                x: max_angle,
                y: max_value,
                marker: "circle",
                color: "red"
            }
        ]
    };
    
    // Add flooding line if applicable
    let flooding = ctx.get_first_flooding_angle();
    if flooding != () {
        plot.elements.push(#{
            type: "VerticalLine",
            name: "θf",
            x: flooding,
            color: "red",
            style: "dashed"
        });
    }
    
    #{
        regulation_name: "IMO A.749(18) General Intact Stability Criteria",
        regulation_reference: "IMO Resolution A.749(18), Section 3.1",
        vessel_name: ctx.get_vessel_name(),
        loading_condition: ctx.get_loading_condition(),
        displacement: ctx.get_displacement(),
        cog: ctx.get_cog(),
        criteria: results,
        overall_pass: fail_count == 0,
        notes: if limiting_angle < 40.0 { "Limiting angle θf = " + limiting_angle + "°" } else { "" },
        plots: [plot]
    }
}
