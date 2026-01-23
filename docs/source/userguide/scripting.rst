Scriptable Verification
=======================

NavalToolbox includes a powerful scripting engine based on `Rhai <https://rhai.rs>`_ that allows you to define and verify custom stability criteria without recompiling the library.

Overview
--------

The verification process involves three steps:

1. **Calculate Stability**: Perform standard stability analysis (GZ curve, etc.).
2. **Create Context**: Wrap the result into a ``CriteriaContext``.
3. **Run Script**: Execute a Rhai script using ``ScriptEngine``.

Basic Example
-------------

.. code-block:: python

    from navaltoolbox import CriteriaContext, ScriptEngine, StabilityCalculator
    
    # 1. Calculate stability
    calc = StabilityCalculator(vessel)
    result = calc.complete_stability(...)
    
    # 2. Create context
    ctx = CriteriaContext.from_result(result, "MV Example", "Departure")
    
    # 3. Run script
    engine = ScriptEngine()
    criteria = engine.run_script_file("rules/imo_a749_general.rhai", ctx)
    
    # Check results
    if criteria.overall_pass:
        print("Vessel Compliant")
    else:
        print(f"Failed {criteria.fail_count} criteria")
        
    for c in criteria.criteria:
        print(f"{c.name}: {c.status} (Actual: {c.actual_value} {c.unit})")

Writing Custom Scripts
----------------------

Scripts are written in Rhai. You have access to the ``ctx`` object which provides helper methods.

.. code-block:: rust

    // Custom criterion script
    
    // 1. Get data
    let area_30 = ctx.area_under_curve(0.0, 30.0);
    let gz_30 = ctx.gz_at_angle(30.0);
    
    let results = [];
    
    // 2. Define criteria
    results.push(criterion(
        "Area 0-30",
        "Minimum area under GZ curve up to 30 deg",
        0.055,      // Required
        area_30,    // Actual
        "m.rad"     // Unit
    ));
    
    results.push(criterion(
        "GZ at 30",
        "Minimum GZ at 30 deg",
        0.20,
        gz_30,
        "m"
    ));
    
    // 3. Define Plot (Optional)
    let plot = #{
        id: "main_plot",
        title: "GZ Curve",
        x_label: "Heel (deg)",
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
    
    // Link criteria to plot
    for r in results { r.plot_id = "main_plot"; }
    
    // 4. Return result map
    let fail_count = 0;
    for r in results { if r.status == "FAIL" { fail_count += 1; } }
    
    #{
        regulation_name: "My Custom Rules",
        regulation_reference: "CUSTOM-01",
        vessel_name: ctx.get_vessel_name(),
        loading_condition: ctx.get_loading_condition(),
        displacement: ctx.get_displacement(),
        overall_pass: fail_count == 0,
        criteria: results,
        plots: [plot],
        notes: ""
    }

Plotting Results
----------------

You can visualize the verification results using the built-in plotting utility (requires Matplotlib).

.. code-block:: python

    from navaltoolbox import plotting
    
    # Display interactively
    plotting.plot_criteria_result(criteria)
    
    # Save to file
    plotting.plot_criteria_result(criteria, show=False, save_to="report.png")

This generates a graph showing the GZ curve along with visual indications of the criteria (threshold lines, areas, etc.).
