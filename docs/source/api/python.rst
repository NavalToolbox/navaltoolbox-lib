Python API Reference
====================

This page documents all Python classes and methods in NavalToolbox.

Hull
----

.. py:class:: Hull(file_path)

   A hull geometry loaded from an STL file.

   :param file_path: Path to the STL file
   :type file_path: str

   .. py:method:: get_bounds()

      Returns the bounding box of the hull.

      :returns: (xmin, xmax, ymin, ymax, zmin, zmax)
      :rtype: tuple[float, float, float, float, float, float]

   .. py:method:: num_triangles()

      Returns the number of triangles in the mesh.

      :rtype: int

   .. py:method:: num_vertices()

      Returns the number of vertices in the mesh.

      :rtype: int

   .. py:method:: transform(translation, rotation, pivot)

      Applies a transformation to the hull.

      :param translation: (dx, dy, dz) translation vector
      :param rotation: (rx, ry, rz) rotation angles in degrees
      :param pivot: (px, py, pz) pivot point for rotation

   .. py:method:: scale(factor)

      Scales the hull uniformly.

      :param factor: Scale factor

   .. py:method:: export_stl(file_path)

      Exports the hull to an STL file.

      :param file_path: Output file path

Vessel
------

.. py:class:: Vessel(hull)

   A vessel containing one or more hulls, tanks, and silhouettes.

   :param hull: The hull geometry
   :type hull: Hull

   .. py:attribute:: ap
      :type: float

      Aft Perpendicular position (x-coordinate).

   .. py:attribute:: fp
      :type: float

      Forward Perpendicular position (x-coordinate).

   .. py:attribute:: lbp
      :type: float

      Length Between Perpendiculars.

   .. py:method:: get_bounds()

      Returns the combined bounding box of all hulls.

   .. py:method:: num_hulls()

      Returns the number of hulls.

   .. py:method:: num_tanks()

      Returns the number of tanks.

   .. py:method:: add_tank(tank)

      Adds a tank to the vessel.

      :param tank: Tank to add
      :type tank: Tank

   .. py:method:: get_total_tanks_mass()

      Returns the total mass of all tank fluids in kg.

   .. py:method:: get_tanks_center_of_gravity()

      Returns the combined CoG of all tank fluids [x, y, z].

   .. py:method:: add_silhouette(silhouette)

      Adds a silhouette profile to the vessel.

      :param silhouette: Silhouette to add
      :type silhouette: Silhouette

   .. py:method:: num_silhouettes()

      Returns the number of silhouettes.

      :rtype: int

   .. py:method:: has_silhouettes()

      Returns True if any silhouettes are defined.

      :rtype: bool

   .. py:method:: clear_silhouettes()

      Removes all silhouettes from the vessel.

   .. py:method:: get_total_emerged_area(waterline_z)

      Returns the total emerged area from all silhouettes (m²).

      :param waterline_z: Waterline Z coordinate
      :type waterline_z: float
      :rtype: float

   .. py:method:: get_combined_emerged_centroid(waterline_z)

      Returns the combined centroid of all emerged areas [x, z].

      :param waterline_z: Waterline Z coordinate
      :type waterline_z: float
      :rtype: list[float]

Silhouette
----------

.. py:class:: Silhouette(file_path)

   A 2D silhouette profile in the X-Z plane for wind heeling calculations.
   
   Used for calculating wind heeling moments per IMO 2008 IS Code (MSC.267).

   :param file_path: Path to a DXF file
   :type file_path: str

   .. py:staticmethod:: from_vtk(file_path)

      Load a silhouette from a VTK file (.vtk or .vtp polyline).

      :param file_path: Path to the VTK file
      :type file_path: str
      :returns: Silhouette object
      :rtype: Silhouette

   .. py:staticmethod:: from_points(points, name)

      Create a silhouette from a list of (x, z) points.

      :param points: List of (x, z) tuples defining the contour
      :type points: list[tuple[float, float]]
      :param name: Silhouette name
      :type name: str
      :returns: Silhouette object
      :rtype: Silhouette

   .. py:attribute:: name
      :type: str

      Silhouette name (from filename or user-defined).

   .. py:method:: num_points()

      Returns the number of points in the contour.

      :rtype: int

   .. py:method:: is_closed()

      Returns True if the contour is closed (first == last point).

      :rtype: bool

   .. py:method:: get_points()

      Returns the points as a list of (x, y, z) tuples.

      :rtype: list[tuple[float, float, float]]

   .. py:method:: get_area()

      Returns the total lateral area (m²).

      :rtype: float

   .. py:method:: get_centroid()

      Returns the centroid [x, z].

      :rtype: list[float]

   .. py:method:: get_bounds()

      Returns the bounding box (x_min, x_max, z_min, z_max).

      :rtype: tuple[float, float, float, float]

   .. py:method:: get_emerged_area(waterline_z)

      Returns the emerged area above waterline (m²).

      :param waterline_z: Waterline Z coordinate
      :type waterline_z: float
      :rtype: float

   .. py:method:: get_emerged_centroid(waterline_z)

      Returns the centroid of the emerged area [x, z].

      :param waterline_z: Waterline Z coordinate
      :type waterline_z: float
      :rtype: list[float]

OpeningType
-----------

.. py:class:: OpeningType

   Type of opening that can cause downflooding.

   .. py:staticmethod:: vent()

      Ventilator opening.

   .. py:staticmethod:: air_pipe()

      Air pipe without automatic closing device.

   .. py:staticmethod:: hatch()

      Hatch or manhole.

   .. py:staticmethod:: door()

      Weathertight door (kept open for operation).

   .. py:staticmethod:: window()

      Non-weathertight window.

   .. py:staticmethod:: other(name)

      Custom opening type.

      :param name: Opening type name
      :type name: str

DownfloodingOpening
-------------------

.. py:class:: DownfloodingOpening

   A downflooding opening point or contour for θf calculation.
   
   Per IMO 2008 IS Code, used to determine when non-weathertight 
   openings become submerged during heel.

   .. py:staticmethod:: from_point(name, position, opening_type)

      Create opening from a single point.

      :param name: Opening name
      :type name: str
      :param position: (x, y, z) position in ship coordinates
      :type position: tuple[float, float, float]
      :param opening_type: Type of opening
      :type opening_type: OpeningType
      :returns: DownfloodingOpening object
      :rtype: DownfloodingOpening

   .. py:staticmethod:: from_contour(name, points, opening_type)

      Create opening from a contour (polyline boundary).

      :param name: Opening name
      :type name: str
      :param points: List of (x, y, z) tuples defining the boundary
      :type points: list[tuple[float, float, float]]
      :param opening_type: Type of opening
      :type opening_type: OpeningType
      :returns: DownfloodingOpening object
      :rtype: DownfloodingOpening

   .. py:attribute:: name
      :type: str

      Opening name.

   .. py:attribute:: is_active
      :type: bool

      Whether opening is active in calculations.

   .. py:method:: set_active(active)

      Set opening active state.

      :param active: True to include in calculations
      :type active: bool

   .. py:method:: num_points()

      Returns the number of points defining the opening.

      :rtype: int

   .. py:method:: get_points()

      Returns points as list of (x, y, z) tuples.

      :rtype: list[tuple[float, float, float]]

   .. py:method:: is_submerged(heel, trim, pivot, waterline_z)

      Check if opening is submerged at given conditions.

      :param heel: Heel angle in degrees
      :type heel: float
      :param trim: Trim angle in degrees
      :type trim: float
      :param pivot: Rotation pivot (x, y, z)
      :type pivot: tuple[float, float, float]
      :param waterline_z: Waterline Z coordinate
      :type waterline_z: float
      :rtype: bool


HydrostaticsCalculator
----------------------

.. py:class:: HydrostaticsCalculator(vessel, water_density=1025.0)

   Calculator for hydrostatic properties.

   :param vessel: The vessel to calculate hydrostatics for
   :param water_density: Water density in kg/m³ (default: 1025.0 for seawater)

   .. py:method:: calculate_at_draft(draft, trim=0.0, heel=0.0, vcg=None)

      Calculates hydrostatics at a given draft, trim, and heel.

      :param draft: Draft at reference point (m)
      :param trim: Trim angle (degrees)
      :param heel: Heel angle (degrees)
      :param vcg: Vertical center of gravity for GM calculation (optional)
      :type vcg: float or None
      :returns: HydrostaticState with all properties

   .. py:method:: calculate_at_displacement(displacement_mass, cog=None, trim=None, heel=None)

      Calculates hydrostatics for a given displacement, finding the required draft.

      :param displacement_mass: Target displacement in kg
      :param cog: Center of gravity (LCG, TCG, VCG) (optional)
      :type cog: tuple[float, float, float] or None
      :param trim: Fixed trim angle (optional)
      :param heel: Fixed heel angle (optional)
      :returns: HydrostaticState with all properties

HydrostaticState
----------------

.. py:class:: HydrostaticState

   Result of hydrostatic calculations.

   .. py:attribute:: draft
      :type: float

      Draft at reference point (m).

   .. py:attribute:: volume
      :type: float

      Submerged volume (m³).

   .. py:attribute:: displacement
      :type: float

      Displacement mass (kg).

   .. py:attribute:: cob
      :type: tuple[float, float, float]

      Center of Buoyancy (LCB, TCB, VCB) vector.

   .. py:attribute:: cog
      :type: tuple[float, float, float] or None

      Center of Gravity (LCG, TCG, VCG) vector, if provided.

   .. py:attribute:: waterplane_area
      :type: float

      Waterplane area (m²).

   .. py:attribute:: lcf
      :type: float

      Longitudinal Center of Flotation (m).

   .. py:attribute:: bmt
      :type: float

      Transverse metacentric radius (m).

   .. py:attribute:: bml
      :type: float

      Longitudinal metacentric radius (m).

   .. py:attribute:: gmt
      :type: float or None

      Transverse metacentric height (Wet/Corrected) (m).

   .. py:attribute:: gmt_dry
      :type: float or None

      Transverse metacentric height (Dry/Uncorrected) (m).

   .. py:attribute:: gml
      :type: float or None

      Longitudinal metacentric height (Wet/Corrected) (m).

StabilityCalculator
-------------------

.. py:class:: StabilityCalculator(vessel, water_density=1025.0)

   Calculator for stability analysis (GZ curves).

   :param vessel: The vessel to calculate stability for
   :param water_density: Water density in kg/m³

   .. py:method:: calculate_gz_curve(displacement_mass, cog, heels)

      Calculates the GZ curve for a given loading condition.

      :param displacement_mass: Base vessel displacement (excluding dynamic tanks) in kg
      :param cog: Base center of gravity (LCG, TCG, VCG) excluding dynamic tanks
      :param heels: List of heel angles in degrees
      :returns: StabilityCurve object with corrected GZ values

StabilityCurve
--------------

.. py:class:: StabilityCurve

   A complete GZ stability curve.

   .. py:method:: heels()

      Returns the heel angles.

      :rtype: list[float]

   .. py:method:: values()

      Returns the GZ values.

      :rtype: list[float]

   .. py:method:: points()

      Returns a list of StabilityPoint with (heel, draft, trim, gz, is_flooding, flooded_openings).
      
      Each point has:
      - ``heel``: Heel angle in degrees
      - ``draft``: Draft at equilibrium (m)
      - ``trim``: Trim angle in degrees
      - ``gz``: GZ value (m)
      - ``is_flooding``: True if any opening is submerged
      - ``flooded_openings``: List of flooded opening names

   .. py:attribute:: displacement
      :type: float

      Displacement in kg.

Tank
----

.. py:class:: Tank

   A tank with fluid management capabilities.

   .. py:staticmethod:: from_box(name, x_min, x_max, y_min, y_max, z_min, z_max, fluid_density)

      Creates a box-shaped tank.

      :param name: Tank name
      :param x_min: Minimum x coordinate
      :param x_max: Maximum x coordinate
      :param y_min: Minimum y coordinate
      :param y_max: Maximum y coordinate
      :param z_min: Minimum z coordinate
      :param z_max: Maximum z coordinate
      :param fluid_density: Fluid density in kg/m³

   .. py:attribute:: name
      :type: str

      Tank name.

   .. py:attribute:: total_volume
      :type: float

      Total tank volume (m³).

   .. py:attribute:: fill_level
      :type: float

      Fill level as fraction (0.0 to 1.0).

   .. py:attribute:: fill_percent
      :type: float

      Fill level as percentage (0 to 100).

   .. py:attribute:: fill_volume
      :type: float

      Filled volume (m³).

   .. py:attribute:: fluid_mass
      :type: float

      Fluid mass (kg).

   .. py:attribute:: center_of_gravity
      :type: list[float]

      Center of gravity [x, y, z].

   .. py:attribute:: free_surface_moment_t
      :type: float

      Transverse free surface moment (m⁴).

   .. py:attribute:: free_surface_moment_l
      :type: float

      Longitudinal free surface moment (m⁴).
