Tanks
=====

Classes for tank management and free surface effects.

Tank
----

   .. py:class:: Tank(file_path, fluid_density=1025.0, name=None)

      A tank with fluid management capabilities.

      :param file_path: Path to the geometry file (STL or VTK)
      :param fluid_density: Fluid density in kg/m³
      :param name: Optional name for the tank

      .. py:method:: __init__(file_path, fluid_density=1025.0, name=None)

         Creates a Tank from a file.

   .. py:staticmethod:: from_box_hull_intersection(hull, x_min, x_max, y_min, y_max, z_min, z_max, fluid_density=1025.0, name="HullTank")

      Creates a Tank as the intersection of a box with a hull geometry.
      Useful for double bottom tanks or side tanks that follow the hull shape.

      :param hull: Hull object to intersect with
      :param x_min: Minimum x coordinate (aft)
      :param x_max: Maximum x coordinate (fwd)
      :param y_min: Minimum y coordinate (stbd)
      :param y_max: Maximum y coordinate (port)
      :param z_min: Minimum z coordinate (bottom)
      :param z_max: Maximum z coordinate (top)
      :param fluid_density: Fluid density in kg/m³
      :param name: Tank name

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

   **Properties**

   .. py:attribute:: name
      :type: str

      Tank name.

   .. py:attribute:: total_volume
      :type: float

      Total tank volume (m³).

   **Fill State**

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

   **Center of Gravity**

   .. py:attribute:: center_of_gravity
      :type: list[float]

      Center of gravity [x, y, z].

   **Free Surface Moments**

   .. py:attribute:: free_surface_moment_t
      :type: float

      Transverse free surface moment (m⁴).

   .. py:attribute:: free_surface_moment_l
      :type: float

      Longitudinal free surface moment (m⁴).

   **Geometry**

   .. py:method:: get_vertices()

      Returns the vertices of the tank shell mesh.

      :returns: List of (x, y, z) tuples
      :rtype: list[tuple[float, float, float]]

   .. py:method:: get_faces()

      Returns the faces (triangles) of the tank shell mesh.

      :returns: List of (i, j, k) indices tuples
      :rtype: list[tuple[int, int, int]]

   .. py:method:: get_fluid_vertices(heel=0.0, trim=0.0)

      Returns the vertices of the fluid mesh at specified attitude.

      :param heel: Heel angle in degrees
      :type heel: float
      :param trim: Trim angle in degrees
      :type trim: float
      :returns: List of (x, y, z) tuples
      :rtype: list[tuple[float, float, float]]

   .. py:method:: get_fluid_faces(heel=0.0, trim=0.0)

      Returns the faces of the fluid mesh at specified attitude.

      :param heel: Heel angle in degrees
      :type heel: float
      :param trim: Trim angle in degrees
      :type trim: float
      :returns: List of (i, j, k) indices tuples
      :rtype: list[tuple[int, int, int]]
