Quickstart
==========

This guide will help you get started with NavalToolbox in Python.

Loading a Hull
--------------

Load a hull geometry from an STL file:

.. code-block:: python

    from navaltoolbox import Hull

    # Load hull from STL file
    hull = Hull("ship.stl")

    # Check hull dimensions
    bounds = hull.get_bounds()
    loa = bounds[1] - bounds[0]  # Length overall
    boa = bounds[3] - bounds[2]  # Beam overall

    print(f"LOA: {loa:.2f}m, BOA: {boa:.2f}m")
    print(f"Triangles: {hull.num_triangles()}")

Creating a Vessel
-----------------

Wrap the hull in a Vessel container:

.. code-block:: python

    from navaltoolbox import Hull, Vessel

    hull = Hull("ship.stl")
    vessel = Vessel(hull)

    # Get perpendiculars (auto-calculated from bounds)
    print(f"AP: {vessel.ap:.2f}m")
    print(f"FP: {vessel.fp:.2f}m")
    print(f"LBP: {vessel.lbp:.2f}m")

Calculating Hydrostatics
------------------------

Calculate hydrostatic properties at a given draft:

.. code-block:: python

    from navaltoolbox import Hull, Vessel, HydrostaticsCalculator

    hull = Hull("ship.stl")
    vessel = Vessel(hull)
    calc = HydrostaticsCalculator(vessel, water_density=1025.0)

    # Calculate at specific draft
    state = calc.from_draft(
        draft=6.0,  # meters
        trim=0.0,   # degrees
        heel=0.0,   # degrees
        vcg=7.5     # VCG for GM calculation
    )

    print(f"Volume: {state.volume:.1f} m³")
    print(f"Displacement: {state.displacement:.0f} kg")
    print(f"Lwl: {state.lwl:.2f}m")
    print(f"Bwl: {state.bwl:.2f}m")
    print(f"Wetted Surface: {state.wetted_surface_area:.1f}m²")
    print(f"Cb: {state.cb:.3f}")
    print(f"LCB: {state.lcb:.2f}m")
    print(f"VCB: {state.vcb:.2f}m")
    print(f"GMT (corrected): {state.gmt:.3f}m")

Calculating GZ Curve
--------------------

Calculate the stability (GZ) curve:

.. code-block:: python

    from navaltoolbox import Hull, Vessel, StabilityCalculator

    hull = Hull("ship.stl")
    vessel = Vessel(hull)
    calc = StabilityCalculator(vessel, water_density=1025.0)

    # Define loading condition
    displacement = 8635000  # kg
    cog = (71.67, 0.0, 7.555)  # LCG, TCG, VCG

    # Heel angles to calculate
    heels = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]

    # Calculate GZ curve
    curve = calc.calculate_gz_curve(displacement, cog, heels)

    # Print results
    for heel, gz in zip(curve.heels(), curve.values()):
        print(f"Heel {heel:5.1f}°: GZ = {gz:+.3f}m")

Complete Stability Analysis
---------------------------

Calculate hydrostatics, GZ curve, and wind data in one call:

.. code-block:: python

    from navaltoolbox import Hull, Vessel, StabilityCalculator, Silhouette

    hull = Hull("ship.stl")
    vessel = Vessel(hull)

    # Optionally add silhouette for wind heeling data
    silhouette = Silhouette.from_points([
        (0, 0), (150, 0), (150, 15), (0, 15), (0, 0)
    ], "hull_profile")
    vessel.add_silhouette(silhouette)

    calc = StabilityCalculator(vessel, water_density=1025.0)

    # Complete stability analysis
    result = calc.calculate_complete_stability(
        displacement_mass=8635000,          # kg
        cog=(71.67, 0.0, 7.555),            # LCG, TCG, VCG
        heels=[0, 10, 20, 30, 40, 50, 60]   # degrees
    )

    # Access hydrostatics
    print(f"Draft: {result.hydrostatics.draft:.3f} m")
    print(f"GM0: {result.gm0:.3f} m")

    # Access GZ curve
    print(f"Max GZ: {result.max_gz:.3f} m at {result.heel_at_max_gz}°")

    # Access wind data (if silhouettes defined)
    if result.has_wind_data():
        wind = result.wind_data
        print(f"Emerged area: {wind.emerged_area:.1f} m²")
        print(f"Wind lever arm: {wind.wind_lever_arm:.2f} m")

Working with Tanks
------------------

Create and manage tanks:

.. code-block:: python

    from navaltoolbox import Tank

    # Create a box-shaped tank
    tank = Tank.from_box(
        name="Fuel Tank 1",
        x_min=40.0, x_max=60.0,  # longitudinal
        y_min=-8.0, y_max=8.0,   # transverse
        z_min=0.0, z_max=5.0,    # vertical
        fluid_density=850.0      # kg/m³ (fuel oil)
    )

    # Set fill level
    tank.fill_percent = 75.0

    # Get tank properties
    print(f"Total volume: {tank.total_volume:.1f} m³")
    print(f"Fill volume: {tank.fill_volume:.1f} m³")
    print(f"Fluid mass: {tank.fluid_mass:.0f} kg")
    print(f"CoG: {tank.center_of_gravity}")
    print(f"FSM (transverse): {tank.free_surface_moment_t:.1f} m⁴")

Working with Silhouettes
------------------------

Create silhouettes for wind heeling calculations:

.. code-block:: python

    from navaltoolbox import Hull, Vessel, Silhouette

    hull = Hull("ship.stl")
    vessel = Vessel(hull)

    # Load from DXF file
    hull_profile = Silhouette("hull_profile.dxf")

    # Or load from VTK file
    superstructure = Silhouette.from_vtk("superstructure.vtk")

    # Or create from points (x, z coordinates)
    container = Silhouette.from_points([
        (40, 15), (80, 15), (80, 25), (40, 25), (40, 15)
    ], "container")

    # Add to vessel
    vessel.add_silhouette(hull_profile)
    vessel.add_silhouette(superstructure)
    vessel.add_silhouette(container)

    print(f"Number of silhouettes: {vessel.num_silhouettes()}")

    # Calculate emerged area at waterline
    waterline_z = 6.15  # meters
    area = vessel.get_total_emerged_area(waterline_z)
    centroid = vessel.get_combined_emerged_centroid(waterline_z)

    print(f"Emerged area: {area:.1f} m²")
    print(f"Centroid: x={centroid[0]:.1f}m, z={centroid[1]:.1f}m")

Next Steps
----------

- See :doc:`../tutorials/hydrostatics` for detailed hydrostatics calculations
- See :doc:`../tutorials/stability` for complete stability analysis
- See :doc:`../api/python/index` for full API reference

