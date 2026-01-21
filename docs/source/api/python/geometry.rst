Geometry
========

Classes for defining hull geometry and vessel structure.

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
