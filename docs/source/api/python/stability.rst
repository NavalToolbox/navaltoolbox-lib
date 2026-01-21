Stability
=========

Classes for stability analysis and GZ curve calculations.

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

   .. py:method:: calculate_complete_stability(displacement_mass, cog, heels)

      Calculates complete stability analysis combining hydrostatics, GZ curve, and wind data.

      :param displacement_mass: Target displacement in kg
      :param cog: Center of gravity (LCG, TCG, VCG) tuple
      :param heels: List of heel angles for GZ curve in degrees
      :returns: CompleteStabilityResult with hydrostatics, GZ curve, and wind data

CompleteStabilityResult
-----------------------

.. py:class:: CompleteStabilityResult

   Complete stability calculation result combining hydrostatics, GZ curve, and wind data.

   .. py:attribute:: hydrostatics
      :type: HydrostaticState

      Hydrostatic state at equilibrium (draft, trim, GM0, etc.).

   .. py:attribute:: gz_curve
      :type: StabilityCurve

      GZ stability curve for the loading condition.

   .. py:attribute:: wind_data
      :type: WindHeelingData or None

      Wind heeling data (if silhouettes are defined).

   .. py:attribute:: displacement
      :type: float

      Displacement mass in kg.

   .. py:attribute:: cog
      :type: tuple[float, float, float]

      Center of gravity (LCG, TCG, VCG).

   .. py:attribute:: gm0
      :type: float or None

      Initial transverse metacentric height (GM0) with free surface correction.

   .. py:attribute:: gm0_dry
      :type: float or None

      Initial transverse metacentric height without free surface correction.

   .. py:attribute:: max_gz
      :type: float or None

      Maximum GZ value on the curve.

   .. py:attribute:: heel_at_max_gz
      :type: float or None

      Heel angle at maximum GZ.

   .. py:method:: has_wind_data()

      Returns True if wind heeling data is available.

      :rtype: bool

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

WindHeelingData
---------------

.. py:class:: WindHeelingData

   Wind heeling data from silhouette calculations.

   Used for wind heeling moment calculations per IMO 2008 IS Code (MSC.267).

   .. py:attribute:: emerged_area
      :type: float

      Total emerged lateral area above waterline (m²).

   .. py:attribute:: emerged_centroid
      :type: tuple[float, float]

      Centroid of emerged area (x, z) in meters.

   .. py:attribute:: wind_lever_arm
      :type: float

      Lever arm from waterline to centroid z coordinate (m).

   .. py:attribute:: waterline_z
      :type: float

      Waterline Z at which calculations were performed.
