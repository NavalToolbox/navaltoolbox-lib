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

    # Alternatively, create a simple box hull for testing
    # hull = Hull.from_box(length=50.0, breadth=10.0, depth=5.0)

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
    state = calc.from_draft(
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

    # Drafts at perpendiculars
    print(f"Draft MP: {state.draft_mp:.3f}m")
    print(f"Draft AP: {state.draft_ap:.3f}m")
    print(f"Draft FP: {state.draft_fp:.3f}m")

Calculating from Drafts at Perpendiculars
-----------------------------------------

You can also calculate hydrostatics by specifying drafts at the Aft Perpendicular (AP) and Forward Perpendicular (FP):

.. code-block:: python

    # Calculate from AP/FP drafts
    # This automatically computes the trim and midship draft
    state = calc.from_drafts(
        draft_ap=6.5,
        draft_fp=5.5,
        heel=0.0,
        vcg=7.5
    )

    print(f"Calculated MP Draft: {state.draft_mp:.3f}m")
    print(f"Calculated Trim: {state.trim:.2f} degrees")
    print(f"Displacement: {state.displacement:.0f} kg")

Note on Reference System
------------------------

In NavalToolbox, the input draft for calculations corresponds to the draft at the **Midship Section (MP)**.
The longitudinal position of MP is calculated as the average of the Aft Perpendicular (AP) and Forward Perpendicular (FP).
If AP and FP are not explicitly set on the vessel, they default to the minimum and maximum X bounds of the geometry.

When an ongoing trim is applied (positive bow down), the drafts at perpendiculars are calculated as:

- **Draft MP**: Equal to the input draft at the pivot point
- **Draft AP**: Calculated from MP draft and trim angle
- **Draft FP**: Calculated from MP draft and trim angle

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
    state = calc.from_displacement(target_displacement)
    print(f"Equilibrium draft: {state.draft:.3f}m")
    print(f"Actual displacement: {state.displacement:.0f} kg")

    # With VCG for GM calculation
    state = calc.from_displacement(target_displacement, vcg=7.555)
    print(f"Draft: {state.draft:.3f}m")
    print(f"GMT: {state.gmt:.3f}m")

    # Or with full COG (LCG, TCG, VCG) for advanced use
    state = calc.from_displacement(
        target_displacement, 
        cog=(71.67, 0.0, 7.555)
    )

Calculating at Heel
-------------------

Hydrostatics change with heel angle:

.. code-block:: python

    # Compare upright vs heeled
    for heel in [0, 10, 20, 30]:
        state = calc.from_draft(6.15, 0.0, heel, 7.555)
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

    state_fresh = calc_fresh.from_draft(draft, 0.0, 0.0, 7.5)
    state_salt = calc_salt.from_draft(draft, 0.0, 0.0, 7.5)

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
