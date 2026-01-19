Hydrostatics Tutorial
=====================

This tutorial covers detailed hydrostatic calculations with NavalToolbox.

Loading the Hull
----------------

First, load a hull geometry:

.. code-block:: python

    from navaltoolbox import Hull, Vessel, HydrostaticsCalculator

    # Load hull from STL file
    hull = Hull("path/to/dtmb5415.stl")

    # Check dimensions
    bounds = hull.get_bounds()
    print(f"LOA: {bounds[1] - bounds[0]:.2f}m")
    print(f"BOA: {bounds[3] - bounds[2]:.2f}m")
    print(f"Depth: {bounds[5] - bounds[4]:.2f}m")

    # Create vessel
    vessel = Vessel(hull)

Basic Hydrostatics
------------------

Calculate hydrostatics at a specific draft:

.. code-block:: python

    # Create calculator (seawater density)
    calc = HydrostaticsCalculator(vessel, water_density=1025.0)

    # Calculate at draft 6.0m, even keel
    state = calc.calculate_at_draft(
        draft=6.0,
        trim=0.0,
        heel=0.0,
        vcg=7.5  # for GM calculation
    )

    print(f"Volume: {state.volume:.1f} m³")
    print(f"Displacement: {state.displacement:.0f} kg")
    print(f"LCB: {state.lcb:.2f}m from AP")
    print(f"TCB: {state.tcb:.3f}m (should be ~0)")
    print(f"VCB: {state.vcb:.2f}m (KB)")

Waterplane & Stability
----------------------

The calculator computes accurate waterplane properties and applies Free Surface Correction (FSC) automatically if tanks are present:

.. code-block:: python

    # Waterplane Properties
    print(f"Waterplane Area: {state.waterplane_area:.1f} m²")
    print(f"LCF: {state.lcf:.2f} m from AP")
    print(f"BMt: {state.bmt:.2f} m")
    print(f"BMl: {state.bml:.2f} m")

    # Metacentric Heights (requires VCG)
    # KM = KB + BM
    # GM_dry = KM - VCG
    # GM_wet = GM_dry - FSC
    print(f"KMt: {state.vcb + state.bmt:.2f} m")
    print(f"GMT (corrected/wet): {state.gmt:.3f} m")
    print(f"GMT (solid/dry):     {state.gmt_dry:.3f} m")

Finding Equilibrium Draft
-------------------------

Find the draft for a known displacement:

.. code-block:: python

    # DTMB 5415 reference displacement: 8635 tonnes
    target_displacement = 8635000  # kg

    # Find draft for displacement (level keel)
    state = calc.calculate_at_displacement(target_displacement)
    print(f"Equilibrium draft: {state.draft:.3f}m")
    print(f"Actual displacement: {state.displacement:.0f} kg")

    # Or with constraints (e.g. fixed VCG)
    state_stab = calc.calculate_at_displacement(
        target_displacement, 
        vcg=7.555
    )
    print(f"Draft with VCG: {state_stab.draft:.3f}m")
    print(f"GMT: {state_stab.gmt:.3f}m")

    # Verify
    state = calc.calculate_at_draft(draft, 0.0, 0.0, 7.555)
    print(f"Actual displacement: {state.displacement:.0f} kg")

Calculating at Heel
-------------------

Hydrostatics change with heel angle:

.. code-block:: python

    # Compare upright vs heeled
    for heel in [0, 10, 20, 30]:
        state = calc.calculate_at_draft(6.15, 0.0, heel, 7.555)
        print(f"Heel {heel:2d}°: Vol={state.volume:.1f}m³, VCB={state.vcb:.2f}m")

Fresh vs Salt Water
-------------------

Density affects displacement:

.. code-block:: python

    # Fresh water (rivers, lakes)
    calc_fresh = HydrostaticsCalculator(vessel, water_density=1000.0)

    # Salt water (sea)
    calc_salt = HydrostaticsCalculator(vessel, water_density=1025.0)

    draft = 6.0

    state_fresh = calc_fresh.calculate_at_draft(draft, 0.0, 0.0, 7.5)
    state_salt = calc_salt.calculate_at_draft(draft, 0.0, 0.0, 7.5)

    print(f"Fresh water: {state_fresh.displacement:.0f} kg")
    print(f"Salt water:  {state_salt.displacement:.0f} kg")
    print(f"Difference:  {state_salt.displacement - state_fresh.displacement:.0f} kg")

Reference Values
----------------

For DTMB 5415 at T=6.15m (SIMMAN 2008):

+------------------+--------------+
| Property         | Value        |
+==================+==============+
| Volume           | ~8424 m³     |
+------------------+--------------+
| CB               | 0.506        |
+------------------+--------------+
| CM               | 0.816        |
+------------------+--------------+
| LCB (from mid)   | -0.686% Lpp  |
+------------------+--------------+
