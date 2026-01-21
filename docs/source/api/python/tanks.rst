Tanks
=====

Classes for tank management and free surface effects.

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
