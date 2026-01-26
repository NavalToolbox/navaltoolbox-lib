Visualization
=============

Interactive 3D visualization tools using Plotly.

.. py:module:: navaltoolbox.visualization

.. py:function:: plot_vessel_3d(vessel, show_hulls=True, show_tanks=True, show_silhouettes=True, show_openings=True, opacity_hull=0.5, opacity_tank=0.3, title="Vessel Visualization", enable_opacity_slider=True, show_axes=True)

   Create an interactive 3D plot of the vessel using Plotly.

   :param vessel: The vessel to visualize.
   :type vessel: Vessel
   :param show_hulls: Whether to show hull meshes.
   :type show_hulls: bool
   :param show_tanks: Whether to show tank meshes.
   :type show_tanks: bool
   :param show_silhouettes: Whether to show silhouette profiles.
   :type show_silhouettes: bool
   :param show_openings: Whether to show downflooding openings.
   :type show_openings: bool
   :param opacity_hull: Opacity for hull meshes (0.0 to 1.0).
   :type opacity_hull: float
   :param opacity_tank: Opacity for tank meshes (0.0 to 1.0).
   :type opacity_tank: float
   :param title: Title of the plot.
   :type title: str
   :param enable_opacity_slider: Whether to add a slider to control hull opacity.
   :type enable_opacity_slider: bool
   :param show_axes: Whether to show axes, grid, and background.
   :type show_axes: bool
   :returns: A Plotly Figure object.
   :rtype: plotly.graph_objects.Figure

.. py:function:: plot_hydrostatic_condition(vessel, draft, trim=0.0, heel=0.0, show_hulls=True, show_tanks=True, show_silhouettes=True, show_openings=True, opacity_hull=0.5, opacity_tank=0.3, title="Hydrostatic Condition", enable_opacity_slider=True, cog=None, show_axes=True)

   Visualize the vessel at a specific hydrostatic condition (floating in water).
   The waterplane is fixed at Z=0. The vessel is transformed (translated and rotated) to match the draft, trim, and heel.

   :param vessel: The vessel to visualize.
   :type vessel: Vessel
   :param draft: Draft at midship (or reference point) in meters.
   :type draft: float
   :param trim: Trim angle in degrees (positive = stern down usually).
   :type trim: float
   :param heel: Heel angle in degrees (positive = starboard down usually).
   :type heel: float
   :param show_hulls: Whether to show hull meshes.
   :type show_hulls: bool
   :param show_tanks: Whether to show tank meshes.
   :type show_tanks: bool
   :param show_silhouettes: Whether to show silhouette profiles.
   :type show_silhouettes: bool
   :param show_openings: Whether to show downflooding openings.
   :type show_openings: bool
   :param opacity_hull: Opacity for hull meshes (0.0 to 1.0).
   :type opacity_hull: float
   :param opacity_tank: Opacity for tank meshes (0.0 to 1.0).
   :type opacity_tank: float
   :param title: Title of the plot.
   :type title: str
   :param enable_opacity_slider: Whether to add a slider to control hull opacity.
   :type enable_opacity_slider: bool
   :param cog: Optional Center of Gravity (x, y, z) to display.
   :type cog: tuple[float, float, float] or None
   :param show_axes: Whether to show axes, grid, and background.
   :type show_axes: bool
   :returns: A Plotly Figure object.
   :rtype: plotly.graph_objects.Figure
