Geometry
========

Classes for defining hull geometry and vessel structure.

Hull
----

.. py:class:: Hull(file_path)

   A hull geometry loaded from an STL file.

   :param file_path: Path to the STL file
   :type file_path: str

   .. py:staticmethod:: from_box(length, breadth, depth)

      Create a box hull.

      :param length: Length of the box in meters
      :type length: float
      :param breadth: Breadth of the box in meters
      :type breadth: float
      :param depth: Depth of the box in meters
      :type depth: float
      :returns: Hull object
      :rtype: Hull

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

   .. py:method:: scale_xyz(sx, sy, sz)

      Scales the hull non-uniformly along each axis.

      :param sx: Scale factor along X axis
      :param sy: Scale factor along Y axis
      :param sz: Scale factor along Z axis

   .. py:method:: simplify(target_count)

      Simplifies the hull mesh to a target number of triangles (in-place).

      :param target_count: Target number of triangles
      :type target_count: int

   .. py:method:: to_simplified(target_count)

      Returns a simplified copy of the hull.

      :param target_count: Target number of triangles
      :type target_count: int
      :returns: Simplified Hull object
      :rtype: Hull


   .. py:method:: export_stl(file_path)

      Exports the hull to an STL file.

      :param file_path: Output file path

   .. py:method:: get_vertices()

      Returns the vertices of the hull mesh.

      :returns: List of (x, y, z) tuples
      :rtype: list[tuple[float, float, float]]

   .. py:method:: get_faces()

      Returns the faces (triangles) of the hull mesh.

      :returns: List of (i, j, k) indices tuples
      :rtype: list[tuple[int, int, int]]

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

   .. py:method:: add_opening(opening)

      Adds a downflooding opening to the vessel.

      :param opening: Opening to add
      :type opening: DownfloodingOpening

   .. py:method:: num_openings()

      Returns the number of downflooding openings.

      :rtype: int

   .. py:method:: clear_openings()

      Removes all downflooding openings from the vessel.

   .. py:method:: get_hulls()

      Get all hulls in the vessel.

      :returns: List of Hull objects
      :rtype: list[Hull]

   .. py:method:: get_tanks()

      Get all tanks in the vessel.

      :returns: List of Tank objects
      :rtype: list[Tank]

   .. py:method:: get_silhouettes()

      Get all silhouettes in the vessel.

      :returns: List of Silhouette objects
      :rtype: list[Silhouette]

   .. py:method:: get_openings()

      Get all downflooding openings in the vessel.

      :returns: List of DownfloodingOpening objects
      :rtype: list[DownfloodingOpening]

Silhouette
----------

.. py:class:: Silhouette(file_path)

   A 2D silhouette profile in the X-Z plane for wind heeling calculations.
   
   Used for calculating wind heeling moments per IMO 2008 IS Code (MSC.267).

   :param file_path: Path to the geometry file (DXF, VTK, VTP)
   :type file_path: str

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

   .. py:staticmethod:: from_file(file_path, default_type, name=None)

      Load openings from a file (DXF or VTK).

      If ``name`` is provided:
      - Single opening: sets logic name to ``name``
      - Multiple openings: sets names to ``{name}_{i+1}``

      :param file_path: Path to the geometry file (DXF, VTK, VTP)
      :type file_path: str
      :param default_type: Default OpeningType for loaded openings
      :type default_type: OpeningType
      :param name: Optional base name for loaded openings
      :type name: str or None

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

Appendage
---------

.. py:class:: Appendage

   An appendage (additional volume element) attached to the vessel.
   
   Appendages represent volume contributions from items like keels, rudders,
   bulbous bows, etc. that are not part of the main hull geometry.

   .. py:staticmethod:: from_point(name, center, volume)

      Create an appendage from a point (fixed volume at position).

      :param name: Appendage name
      :type name: str
      :param center: (x, y, z) center position
      :type center: tuple[float, float, float]
      :param volume: Volume in m³
      :type volume: float
      :returns: Appendage object
      :rtype: Appendage

   .. py:staticmethod:: from_file(name, file_path)

      Create an appendage from an STL or VTK file.

      :param name: Appendage name
      :type name: str
      :param file_path: Path to the geometry file
      :type file_path: str
      :returns: Appendage object
      :rtype: Appendage

   .. py:staticmethod:: from_box(name, xmin, xmax, ymin, ymax, zmin, zmax)

      Create an appendage from a box (parallelepiped).

      :param name: Appendage name
      :type name: str
      :param xmin: Minimum X coordinate
      :param xmax: Maximum X coordinate
      :param ymin: Minimum Y coordinate
      :param ymax: Maximum Y coordinate
      :param zmin: Minimum Z coordinate
      :param zmax: Maximum Z coordinate
      :returns: Appendage object
      :rtype: Appendage

   .. py:staticmethod:: from_cube(name, center, volume)

      Create an appendage from a cube (center and volume).

      :param name: Appendage name
      :type name: str
      :param center: (x, y, z) center position
      :type center: tuple[float, float, float]
      :param volume: Volume in m³
      :type volume: float
      :returns: Appendage object
      :rtype: Appendage

   .. py:staticmethod:: from_sphere(name, center, volume)

      Create an appendage from a sphere (center and volume).

      :param name: Appendage name
      :type name: str
      :param center: (x, y, z) center position
      :type center: tuple[float, float, float]
      :param volume: Volume in m³
      :type volume: float
      :returns: Appendage object
      :rtype: Appendage

   .. py:attribute:: name
      :type: str

      The appendage name.

   .. py:attribute:: volume
      :type: float

      Volume in m³.

   .. py:attribute:: center
      :type: tuple[float, float, float]

      Center of volume (x, y, z) in meters.

   .. py:attribute:: wetted_surface
      :type: float or None

      Wetted surface area in m², or None if not set.

   .. py:attribute:: bounds
      :type: tuple[float, float, float, float, float, float] or None
      
      Returns bounds (xmin, xmax, ymin, ymax, zmin, zmax) or None if not applicable (Point).

   .. py:method:: geometry_type()

      Returns the geometry type (Point, Mesh, Box, etc.).

      :rtype: str

   .. py:method:: get_mesh_data()

      Returns mesh data (vertices, faces) if geometry is a mesh.

      :returns: Tuple of (vertices, faces), or None if not a mesh
      :rtype: tuple[list[tuple[float, float, float]], list[tuple[int, int, int]]] or None

DeckEdge
--------

.. py:class:: DeckEdgeSide

   Side of the deck edge (Port, Starboard, or Both).

   .. py:staticmethod:: port()

      Port side.

   .. py:staticmethod:: starboard()

      Starboard side.

   .. py:staticmethod:: both()

      Both sides (mirrored).

.. py:class:: DeckEdge

   A deck edge contour (livet) for freeboard calculation.

   .. py:staticmethod:: from_points(name, points, side)

      Create a deck edge from a list of 3D points.

      :param name: Deck edge name
      :type name: str
      :param points: List of (x, y, z) points
      :type points: list[tuple[float, float, float]]
      :param side: Side of the deck edge
      :type side: DeckEdgeSide
      :returns: DeckEdge object
      :rtype: DeckEdge

   .. py:staticmethod:: from_file(name, file_path)

      Load a deck edge from a DXF or VTK file.

      :param name: Deck edge name
      :type name: str
      :param file_path: Path to the geometry file
      :type file_path: str
      :returns: DeckEdge object
      :rtype: DeckEdge

   .. py:attribute:: name
      :type: str

      The deck edge name.

   .. py:method:: num_points()

      Returns the number of points.

      :rtype: int

   .. py:method:: get_points()

      Returns points as list of (x, y, z) tuples.

      :rtype: list[tuple[float, float, float]]

   .. py:method:: get_side()

      Returns the side as a string.

      :rtype: str

   .. py:method:: get_freeboard(heel, trim, pivot, waterline_z)

      Calculate freeboard at given conditions.

      :param heel: Heel angle in degrees
      :type heel: float
      :param trim: Trim angle in degrees
      :type trim: float
      :param pivot: Rotation pivot (x, y, z)
      :type pivot: tuple[float, float, float]
      :param waterline_z: Waterline Z coordinate
      :type waterline_z: float
      :returns: Minimum freeboard distance in meters
      :rtype: float
