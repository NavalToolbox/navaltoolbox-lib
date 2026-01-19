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
    curve = calc.calculate_gz_curve(displacement, cog, heels)

    # Print results
    print("Heel (°)   GZ (m)")
    print("-" * 20)
    for heel, gz in zip(curve.heels(), curve.values()):
        print(f"{heel:6.1f}    {gz:+.3f}")

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
        curve = calc.calculate_gz_curve(displacement, cog, [0, 10, 20, 30, 40, 50])

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
