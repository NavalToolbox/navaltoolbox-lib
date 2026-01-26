Visualization Tutorial
======================

This tutorial demonstrates how to use the interactive 3D visualization capabilities of NavalToolbox.

The visualization module uses `Plotly <https://plotly.com/python/>`_ to create interactive HTML files that can be viewed in any web browser.

Setup
-----

First, we load a hull and create a vessel instance, as usual:

.. code-block:: python

    from navaltoolbox import Hull, Vessel, Tank, Silhouette, DownfloodingOpening, OpeningType
    from navaltoolbox.visualization import plot_vessel_3d, plot_hydrostatic_condition

    # 1. Load Hull
    hull = Hull("path/to/dtmb5415.stl")
    vessel = Vessel(hull)

Adding Components
-----------------

To make the visualization interesting, let's add some tanks, a silhouette (superstructure), and an opening.

.. code-block:: python

    # Add a Fuel Tank
    tank = Tank.from_box(
        name="Fuel Tank 1P",
        x_min=50.0, x_max=60.0,
        y_min=2.0, y_max=8.0,
        z_min=0.0, z_max=5.0,
        fluid_density=850.0
    )
    tank.fill_percent = 50.0 # Partially filled to show fluid surface
    vessel.add_tank(tank)

    # Add a Silhouette (Wind profile)
    # Ideally load from DXF, but here's a manual one for demo consistency
    # (Assuming you have a DXF file or points)
    # sil = Silhouette("superstructure.dxf")
    # vessel.add_silhouette(sil)

    # Add a Downflooding Opening
    opening = DownfloodingOpening.from_point(
        name="Vent 1",
        position=(75.0, 9.8, 12.0),
        opening_type=OpeningType.vent()
    )
    vessel.add_opening(opening)

General Arrangement Plot
------------------------

Use :func:`~navaltoolbox.visualization.plot_vessel_3d` to visualize the vessel and its internal components in the upright condition.

.. code-block:: python

    # Create the figure
    fig = plot_vessel_3d(
        vessel,
        show_hulls=True,
        show_tanks=True,
        show_silhouettes=True,
        show_openings=True,
        opacity_hull=0.3, # Semi-transparent to see tanks inside
        title="General Arrangement"
    )

    # Save to HTML
    fig.write_html("vessel_ga.html")
    # fig.show() # To open immediately in browser

The resulting plot includes:

*   **View Controls**: Buttons to switch between **ISO**, **STBD** (Profile), **TOP** (Plan), etc.
*   **Projection Toggle**: Switch between **Perspective** (realistic) and **Orthographic** (engineering) views.
*   **Hull Opacity Slider**: A slider to adjust hull transparency dynamically.

Hydrostatic Condition Plot
--------------------------

Use :func:`~navaltoolbox.visualization.plot_hydrostatic_condition` to visualize the vessel floating at a specific draft, trim, and heel.
This visualization:

1.   Fixes the waterplane at Z=0.
2.   Rotates and translates the vessel to match the floating condition.
3.   Shows the **Center of Gravity (COG)**.
4.  Shows the fluid level inside tanks, correctly inclined.

.. code-block:: python

    # Define floating condition
    draft_mid = 6.15
    heel = 15.0 # degrees stability test
    trim = 1.0  # degrees by stern

    # Assume we know the COG position
    cog_pos = (70.0, 0.0, 7.5) 

    # Create the hydrostatic plot
    fig_hydro = plot_hydrostatic_condition(
        vessel,
        draft=draft_mid,
        heel=heel,
        trim=trim,
        title=f"Hydrostatic Condition (T={draft_mid}m, Heel={heel}deg)",
        opacity_hull=0.5,
        cog=cog_pos,
        show_axes=True # Set to False for a cleaner look
    )

    # Save
    fig_hydro.write_html("hydrostatics.html")

Features:

*   The **Waterplane** is visualized as a blue surface at Z=0.
*   In **Orthographic Side Views** (e.g. STBD), a blue line represents the waterline level.
*   **Tanks** show the fluid surface parallel to the external waterplane (taking heel/trim into account).

.. note::
   The vessel is transformed (moved/rotated) relative to a fixed world frame where the water surface is always at Z=0. This makes it intuitive to check immersion and downflooding points relative to the water.
