Hydrostatics
============


TankOptions
-----------

.. py:class:: TankOptions(include_mass=False, include_fsm=True)

   Options for tank handling in hydrostatic calculations.

   :param include_mass: Include tank fluid mass in displacement (default: False)
   :param include_fsm: Apply Free Surface Moment correction to GM (default: True)

   .. py:staticmethod:: none()

      Create options with no tank effects (mass=False, fsm=False).

   .. py:staticmethod:: all()

      Create options with all tank effects (mass=True, fsm=True).

   .. py:staticmethod:: mass_only()

      Create options with only mass included (mass=True, fsm=False).

   .. py:staticmethod:: fsm_only()

      Create options with only FSM correction (mass=False, fsm=True).

   .. py:attribute:: include_mass
      :type: bool

      Whether tank mass is included.

   .. py:attribute:: include_fsm
      :type: bool

      Whether FSM correction is applied.


HydrostaticsCalculator
----------------------

.. py:class:: HydrostaticsCalculator(vessel, water_density=1025.0)

   Calculator for hydrostatic properties.

   :param vessel: The vessel to calculate hydrostatics for
   :param water_density: Water density in kg/m³ (default: 1025.0 for seawater)

   .. py:method:: from_draft(draft, trim=0.0, heel=0.0, vcg=None, num_stations=None, tank_options=None)

      Calculates hydrostatics at a given draft, trim, and heel.

      :param draft: Draft at reference point (m)
      :param trim: Trim angle (degrees)
      :param heel: Heel angle (degrees)
      :param vcg: Vertical center of gravity for GM calculation (optional)
      :type vcg: float or None
      :param num_stations: Number of stations (optional)
      :param tank_options: TankOptions object (optional)
      :returns: HydrostaticState with all properties

   .. py:method:: from_drafts(draft_ap, draft_fp, heel=0.0, vcg=None, tank_options=None)

      Calculates hydrostatics from drafts at Aft and Forward Perpendiculars.

      :param draft_ap: Draft at Aft Perpendicular (m)
      :param draft_fp: Draft at Forward Perpendicular (m)
      :param heel: Heel angle (degrees)
      :param vcg: Vertical center of gravity for GM calculation (optional)
      :type vcg: float or None
      :param tank_options: TankOptions object (optional)
      :returns: HydrostaticState with all properties

   .. py:method:: from_displacement(displacement_mass, vcg=None, cog=None, trim=None, heel=None, tank_options=None)

      Calculates hydrostatics for a given displacement, finding the required draft.

      :param displacement_mass: Target displacement in kg
      :param vcg: Vertical center of gravity (m) for GM calculation (optional)
      :param cog: Full center of gravity (LCG, TCG, VCG) tuple (optional, overrides vcg)
      :param trim: Fixed trim angle in degrees (optional, default 0.0)
      :param heel: Fixed heel angle in degrees (optional, default 0.0)
      :param tank_options: TankOptions object (optional)
      :returns: HydrostaticState with all properties

      .. note::
         When ``tank_options`` is provided with ``include_mass=True``, ``displacement_mass``
         is interpreted as the mass of the vessel **excluding** the fluid in the tanks.
         The total displacement used for calculation will be ``displacement_mass + tank_fluid_mass``.


HydrostaticState
----------------

.. py:class:: HydrostaticState

   Result of hydrostatic calculations.

   **Basic Properties**

   .. py:attribute:: draft
      :type: float

      Draft at Midship Section (MP) in meters.

   .. py:attribute:: draft_ap
      :type: float

      Draft at Aft Perpendicular (AP) in meters.

   .. py:attribute:: draft_fp
      :type: float

      Draft at Forward Perpendicular (FP) in meters.

   .. py:attribute:: draft_mp
      :type: float

      Draft at Midship Perpendicular (MP) in meters (alias for draft).


   .. py:attribute:: volume
      :type: float

      Submerged volume (m³).

   .. py:attribute:: displacement
      :type: float

      Total displacement mass in kg (Volume * Water Density). Represents the total weight of the floating system (Vessel + Tanks).

   .. py:attribute:: vessel_displacement
      :type: float

      Vessel displacement mass in kg (Total - Tank Mass). Represents the input displacement (Lightship + Deadweight excluding tanks).

   .. py:attribute:: tank_displacement
      :type: float

      Tank fluid mass in kg (sum of all tank fluid masses). 0.0 if ``include_mass`` is False.

   **Centers**

   .. py:attribute:: cob
      :type: tuple[float, float, float]

      Center of Buoyancy (LCB, TCB, VCB) vector.

   .. py:attribute:: cog
      :type: tuple[float, float, float] or None

      Center of Gravity (LCG, TCG, VCG) vector, if provided. Matches ``displacement``.

   .. py:attribute:: vessel_cog
      :type: tuple[float, float, float] or None

      Vessel Center of Gravity (LCG, TCG, VCG) vector (Ship-only). Matches ``vessel_displacement``.

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

   **Advanced Properties**

   .. py:attribute:: sectional_areas
      :type: list[tuple[float, float]]

      Sectional Area Curve as list of (x, area) pairs.

   .. py:attribute:: freeboard
      :type: float or None

      Minimum freeboard in meters, or None if no deck edges defined.
