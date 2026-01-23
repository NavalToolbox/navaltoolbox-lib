// ============================================================================
// Custom Criteria Template
// ============================================================================
//
// Copy and modify this template to create your own stability criteria scripts.
//
// Usage:
//   ctx = CriteriaContext.from_result(result, vessel_name="MV Example")
//   ctx.set_param("my_param", 123.0)  // Set custom parameters
//   engine = ScriptEngine()
//   result = engine.run_script_file("rules/my_criteria.rhai", ctx)
//
// ============================================================================

// === YOUR THRESHOLD VALUES ===
const MY_AREA_MIN = 0.055;       // m·rad
const MY_GZ_MIN = 0.20;          // m
const MY_GM_MIN = 0.15;          // m



// === MAIN CHECK FUNCTION ===
fn check(ctx) {
    let results = [];
    
    // -----------------------------------------------------------------
    // Example: Get user-defined parameters
    // -----------------------------------------------------------------
    let my_param = ctx.get_param("my_param");
    if my_param == () {
        my_param = 1.0;  // Default value if not provided
    }
    
    // -----------------------------------------------------------------
    // Example: Access hydrostatic data
    // -----------------------------------------------------------------
    let draft = ctx.get_draft();
    let bwl = ctx.get_bwl();
    let cb = ctx.get_cb();
    let gm0 = ctx.get_gm0();
    
    // -----------------------------------------------------------------
    // Example: Check area under GZ curve
    // -----------------------------------------------------------------
    let area = ctx.area_under_curve(0.0, 30.0);
    results.push(criterion(
        "Area 0-30°",
        "Your custom area check",
        MY_AREA_MIN,
        area,
        "m·rad"
    ));
    
    // -----------------------------------------------------------------
    // Example: Check GZ at specific angle
    // -----------------------------------------------------------------
    let gz_at_30 = ctx.gz_at_angle(30.0);
    results.push(criterion(
        "GZ at 30°",
        "Custom GZ threshold check",
        MY_GZ_MIN,
        gz_at_30,
        "m"
    ));
    
    // -----------------------------------------------------------------
    // Example: Check GM0
    // -----------------------------------------------------------------
    if gm0 != () {
        results.push(criterion(
            "GM₀",
            "Metacentric height",
            MY_GM_MIN,
            gm0,
            "m"
        ));
    }
    
    // -----------------------------------------------------------------
    // Example: Find max GZ
    // -----------------------------------------------------------------
    let max_info = ctx.find_max_gz();
    let max_angle = max_info.angle;
    let max_gz = max_info.value;
    
    // -----------------------------------------------------------------
    // Compute overall result
    // -----------------------------------------------------------------
    let fail_count = 0;
    for r in results {
        if r.status == "FAIL" { fail_count += 1; }
    }
    
    #{
        regulation_name: "My Custom Criteria",
        regulation_reference: "CUSTOM-001",
        vessel_name: ctx.get_vessel_name(),
        loading_condition: ctx.get_loading_condition(),
        displacement: ctx.get_displacement(),
        cog: ctx.get_cog(),
        criteria: results,
        overall_pass: fail_count == 0,
        notes: "",
    // -----------------------------------------------------------------
    // Example: Optional Plot
    // -----------------------------------------------------------------
    let plot_id = "main_plot";
    // Associate criteria with plot
    for r in results {
        r.plot_id = plot_id;
    }

    let plot = #{
        id: plot_id,
        title: "GZ Curve",
        x_label: "Heel (°)",
        y_label: "GZ (m)",
        elements: [
            #{
                type: "Curve",
                name: "GZ",
                x: ctx.get_heels(),
                y: ctx.get_gz_values(),
                color: "blue"
            }
        ]
    };

    #{
        regulation_name: "My Custom Criteria",
        regulation_reference: "CUSTOM-001",
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
