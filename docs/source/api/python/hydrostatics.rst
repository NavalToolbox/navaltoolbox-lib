Hydrostatics
============

Classes for hydrostatic calculations.

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

   .. py:method:: calculate_at_displacement(displacement_mass, vcg=None, cog=None, trim=None, heel=None)

      Calculates hydrostatics for a given displacement, finding the required draft.

      :param displacement_mass: Target displacement in kg
      :param vcg: Vertical center of gravity (m) for GM calculation (optional)
      :param cog: Full center of gravity (LCG, TCG, VCG) tuple (optional, overrides vcg)
      :param trim: Fixed trim angle in degrees (optional, default 0.0)
      :param heel: Fixed heel angle in degrees (optional, default 0.0)
      :returns: HydrostaticState with all properties


HydrostaticState
----------------

.. py:class:: HydrostaticState

   Result of hydrostatic calculations.

   **Basic Properties**

   .. py:attribute:: draft
      :type: float

      Draft at reference point (m).

   .. py:attribute:: volume
      :type: float

      Submerged volume (m³).

   .. py:attribute:: displacement
      :type: float

      Displacement mass (kg).

   **Centers**

   .. py:attribute:: cob
      :type: tuple[float, float, float]

      Center of Buoyancy (LCB, TCB, VCB) vector.

   .. py:attribute:: cog
      :type: tuple[float, float, float] or None

      Center of Gravity (LCG, TCG, VCG) vector, if provided.

   **Waterplane Properties**

   .. py:attribute:: waterplane_area
      :type: float

      Waterplane area (m²).

   .. py:attribute:: lcf
      :type: float

      Longitudinal Center of Flotation (m).

   **Metacentric Radii**

   .. py:attribute:: bmt
      :type: float

      Transverse metacentric radius (m).

   .. py:attribute:: bml
      :type: float

      Longitudinal metacentric radius (m).

   **Metacentric Heights (Corrected for Free Surface)**

   .. py:attribute:: gmt
      :type: float or None

      Transverse metacentric height (Wet/Corrected) (m).

   .. py:attribute:: gml
      :type: float or None

      Longitudinal metacentric height (Wet/Corrected) (m).

   **Metacentric Heights (Uncorrected)**

   .. py:attribute:: gmt_dry
      :type: float or None

      Transverse metacentric height (Dry/Uncorrected) (m).

   .. py:attribute:: gml_dry
      :type: float or None

      Longitudinal metacentric height (Dry/Uncorrected) (m).

   **Dimensions**

   .. py:attribute:: lwl
      :type: float

      Length at Waterline (m).

   .. py:attribute:: bwl
      :type: float

      Beam at Waterline (m).

   .. py:attribute:: los
      :type: float

      Length Overall Submerged (m).

   .. py:attribute:: wetted_surface_area
      :type: float

      Wetted Surface Area (m²).

   .. py:attribute:: midship_area
      :type: float

      Midship Section Area (m²).

   **Form Coefficients**

   .. py:attribute:: cb
      :type: float

      Block Coefficient.

      .. math::

          C_b = \frac{\nabla}{L_{wl} \cdot B_{wl} \cdot T}

   .. py:attribute:: cm
      :type: float

      Midship Section Coefficient.

      .. math::

          C_m = \frac{A_m}{B_{wl} \cdot T}

   .. py:attribute:: cp
      :type: float

      Prismatic Coefficient.

      .. math::

          C_p = \frac{\nabla}{A_m \cdot L_{wl}}

   **Free Surface Corrections**

   .. py:attribute:: free_surface_correction_t
      :type: float

      Transverse Free Surface Correction (m).

   .. py:attribute:: free_surface_correction_l
      :type: float

      Longitudinal Free Surface Correction (m).

   **Stiffness Matrix**

   .. py:attribute:: stiffness_matrix
      :type: list[float]

      6x6 Hydrostatic Stiffness Matrix (flattened row-major).
