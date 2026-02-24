Stability Tutorial
==================

This tutorial covers GZ curve calculations for intact stability analysis.

Introduction
------------

The GZ (righting arm) curve shows how a vessel's stability changes with heel angle.
Positive GZ values indicate the vessel will return to upright.

Basic GZ Calculation
--------------------

.. code-block:: python

    from navaltoolbox import Hull, Vessel, StabilityCalculator

    # Load hull
    hull = Hull("dtmb5415.stl")
    vessel = Vessel(hull)

    # Create stability calculator
    calc = StabilityCalculator(vessel, water_density=1025.0)

    # Define loading condition
    # Note: This is the vessel mass EXCLUDING dynamic tank fluids.
    # The calculator will add the mass of any active tanks to this value.
    displacement = 8635000  # kg (8635 tonnes)
    cog = (71.67, 0.0, 7.555)  # LCG, TCG, VCG in meters (Base/Lightship COG)

    # Heel angles to calculate
    heels = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]

    # Calculate GZ curve
    curve = calc.gz_curve(displacement, cog, heels)

    # Print results
    print("Heel (°)   GZ (m)")
    print("-" * 20)
    for heel, gz in zip(curve.heels(), curve.values()):
        print(f"{heel:6.1f}    {gz:+.3f}")

Fixed Trim Option
-----------------

By default, the calculator finds the free-trim equilibrium for each heel angle. However, you can force a constant trim angle using the `fixed_trim` parameter:

.. code-block:: python

    # Calculate GZ curve maintaining exactly 0.0 degrees of trim
    fixed_trim_curve = calc.gz_curve(displacement, cog, heels, fixed_trim=0.0)

    # For cross curves (KN)
    kn_curves_fixed = calc.kn_curve(displacements, heels, fixed_trim=0.0)

Analyzing the Curve
-------------------

Extract key stability metrics:

.. code-block:: python

    # Get GZ values
    heels = curve.heels()
    gz_values = curve.values()

    # Find maximum GZ
    max_gz = max(gz_values)
    max_gz_idx = gz_values.index(max_gz)
    max_gz_heel = heels[max_gz_idx]

    print(f"Maximum GZ: {max_gz:.3f}m at {max_gz_heel:.1f}°")

    # Check positive range
    positive_range = [h for h, gz in zip(heels, gz_values) if gz > 0]
    print(f"Positive GZ range: 0° to {max(positive_range):.0f}°")

    # Read freeboard and draft at specific heels
    points = curve.get_stability_points()
    for pt in points:
        if pt.freeboard is not None and pt.freeboard < 0.0:
            print(f"Deck edge submerged at {pt.heel}°")

    # Initial GM (from slope at 0°)
    # GZ ≈ GM × sin(φ) for small angles
    if len(gz_values) > 1 and heels[1] > 0:
        gm_approx = gz_values[1] / math.sin(math.radians(heels[1]))
        print(f"Approximate GM: {gm_approx:.2f}m")

Effect of VCG
-------------

VCG (vertical center of gravity) significantly affects stability:

.. code-block:: python

    # Compare different VCG values
    vcg_values = [6.5, 7.0, 7.5, 8.0, 8.5]

    print("VCG (m)    Max GZ (m)    Max GZ Angle")
    print("-" * 45)

    for vcg in vcg_values:
        cog = (71.67, 0.0, vcg)
        curve = calc.gz_curve(displacement, cog, [0, 10, 20, 30, 40, 50])

        gz_vals = curve.values()
        max_gz = max(gz_vals)
        max_idx = gz_vals.index(max_gz)

        print(f"{vcg:6.1f}      {max_gz:.3f}        {curve.heels()[max_idx]:.0f}°")

DTMB 5415 Reference Values
--------------------------

NavalToolbox has been validated against reference data from Ariffin (2017).

+--------+-------------+---------------+
| Heel   | Reference   | NavalToolbox  |
+========+=============+===============+
| 10°    | 0.339 m     | 0.333 m       |
+--------+-------------+---------------+
| 20°    | 0.674 m     | 0.669 m       |
+--------+-------------+---------------+
| 30°    | 0.993 m     | 0.982 m       |
+--------+-------------+---------------+
| 40°    | 1.077 m     | 1.052 m       |
+--------+-------------+---------------+

**Maximum error: 2.5cm**

Plotting Results
----------------

Visualize the GZ curve using matplotlib:

.. code-block:: python

    import matplotlib.pyplot as plt

    heels = curve.heels()
    gz_values = curve.values()

    plt.figure(figsize=(10, 6))
    plt.plot(heels, gz_values, 'b-o', linewidth=2, markersize=6)
    plt.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    plt.xlabel('Heel Angle (degrees)')
    plt.ylabel('GZ (m)')
    plt.title('Stability Curve (GZ)')
    plt.grid(True, alpha=0.3)
    plt.xlim(0, max(heels))
    plt.show()

Including Tanks
---------------

Free surface effects reduce effective GM:

.. code-block:: python

    from navaltoolbox import Tank

    # Add a fuel tank
    fuel_tank = Tank.from_box(
        "Fuel Tank",
        x_min=50, x_max=70,
        y_min=-5, y_max=5,
        z_min=0, z_max=4,
        fluid_density=850.0
    )
    
    # Alternatively, load from file (STL/VTK):
    # fuel_tank = Tank("fuel_tank.stl", fluid_density=850.0, name="FO_1P")

    # Or create from hull intersection (e.g., double bottom):
    # fuel_tank = Tank.from_box_hull_intersection(
    #     hull, 
    #     x_min=20, x_max=50, y_min=-10, y_max=10, z_min=0, z_max=1.5,
    #     fluid_density=850.0, name="DB_1C"
    # )
    fuel_tank.fill_percent = 50  # Partial fill = free surface

    vessel.add_tank(fuel_tank)

    # Free surface correction
    fsm = fuel_tank.free_surface_moment_t
    fs_correction = (fsm * fuel_tank.fluid_mass / fuel_tank.fill_volume) / displacement

    print(f"Free surface moment: {fsm:.1f} m⁴")
    print(f"GM reduction (virtual VCG rise): {fs_correction:.3f}m")

    print(f"GM reduction (virtual VCG rise): {fs_correction:.3f}m")

    # Note: Both HydrostaticsCalculator and StabilityCalculator apply this
    # correction automatically if tanks are present in the vessel.
    # The GZ curve computed above ALREADY includes this reduction.
    # No manual adjustment of VCG is required.

Downflooding Points
-------------------

Downflooding points are critical openings (vents, hatches) through which water can enter the hull. They are used to calculate the downflooding angle (:math:`\theta_f`).

.. code-block:: python

    from navaltoolbox import DownfloodingOpening, OpeningType

    # Create opening from a point
    er_vent = DownfloodingOpening.from_point(
        "ER Vent Stbd",
        position=(60.0, 8.0, 12.0),
        opening_type=OpeningType.vent()
    )
    vessel.add_opening(er_vent)

    # Or load from file
    # openings = DownfloodingOpening.from_file("openings.dxf", OpeningType.hatch())
    # for op in openings:
    #     vessel.add_opening(op)

    # Check for flooding
    # Access detailed point data for specific heel angles
    points = curve.get_stability_points()
    for pt in points:
        if pt.is_flooding:
            print(f"Flooding detected at {pt.heel}° through: {pt.flooded_openings}")
            break

Wind Heeling
------------

Wind heeling moments are calculated based on the lateral projected area of the vessel's superstructure (silhouette).

.. code-block:: python

    from navaltoolbox import Silhouette

    # Load silhouette from DXF
    # silhouette = Silhouette("superstructure.dxf")
    
    # Or create manually
    silhouette = Silhouette.from_points([
        (0, 0), (100, 0), (100, 10), (80, 20), (20, 20), (0, 10)
    ], "Superstructure")
    
    vessel.add_silhouette(silhouette)

    # Calculate complete stability
    result = calc.complete_stability(displacement, cog, heels)
    
    if result.has_wind_data():
        wd = result.wind_data
        print(f"Wind Area: {wd.emerged_area:.1f} m²")
        print(f"Wind Moment Lever: {wd.wind_lever_arm:.2f} m")

Cross Curves (KN)
-----------------

Cross curves of stability (KN curves) describe the geometric stability of the hull form independent of the center of gravity height. They are typically plotted as **KN vs Displacement** for fixed heel angles.

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt

    # Define range of displacements
    # Example: 80% to 120% of design displacement
    base_disp = 8635000.0
    displacements = np.linspace(base_disp * 0.8, base_disp * 1.2, 5).tolist()

    # Fixed heel angles for cross curves (typically every 10°)
    heels = [10, 20, 30, 40, 50, 60]

    # Calculate KN curves for all displacements at once
    # Note: Returns a list of StabilityCurve objects, one per displacement
    kn_curves = calc.kn_curve(displacements, heels)

    # Plotting Cross Curves (KN vs Displacement)
    plt.figure(figsize=(10, 6))

    # We want to plot one line per heel angle
    # X-axis: Displacement
    # Y-axis: KN value

    # Extract data: matrix [displacement_idx][heel_idx]
    kn_values = [[curve.values()[h_idx] for curve in kn_curves] for h_idx in range(len(heels))]

    for h_idx, heel in enumerate(heels):
        plt.plot(displacements, kn_values[h_idx], 'o-', label=f'{heel}°')

    plt.xlabel('Displacement (kg)')
    plt.ylabel('KN (m)')
    plt.title('Cross Curves of Stability (KN)')
    plt.legend(title="Heel Angle")
    plt.grid(True, alpha=0.3)
    plt.show()

.. tip::
    The `kn_curve` method is much faster than looping over displacements manually because it can optimize geometry re-use internally.

Verification via Scripting
--------------------------

Once you have calculated stability, you can verify it against regulations using the built-in scripting engine.

.. code-block:: python

    from navaltoolbox import CriteriaContext, ScriptEngine, plotting

    # 1. Create context from result
    # (Assuming 'result' is the output of complete_stability())
    ctx = CriteriaContext.from_result(result, "My Vessel", "Departure")

    # 2. Run verification script
    engine = ScriptEngine()
    # Assuming rule file is available
    # criteria = engine.run_script_file("rules/is_code_2008_general.rhai", ctx)

    # 3. Plot results
    # plotting.plot_criteria_result(criteria)

For more details, see the :doc:`Scriptable Verification <../userguide/scripting>` guide.
