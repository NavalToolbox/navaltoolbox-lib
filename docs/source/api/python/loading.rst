Loading Conditions
==================

Classes for managing mass inventories and tank fill overrides to define operational profiles.

MassCategory
------------

   .. py:class:: MassCategory

      Category of a mass item.

      .. py:staticmethod:: lightship()

         Creates a Lightship category.

      .. py:staticmethod:: deadweight()

         Creates a Deadweight category.

      .. py:staticmethod:: other()

         Creates an Other category (default).


MassItem
--------

   .. py:class:: MassItem(name, mass, cog, category=None)

      A single mass item with name, mass, position, and optional category.

      :param name: Identifier for the mass item.
      :param mass: Mass in kg.
      :param cog: Center of gravity (lcg, tcg, vcg) in meters.
      :param category: Optional MassCategory (default: Other).

      .. py:method:: __init__(name, mass, cog, category=None)

         Creates a MassItem.

   **Properties**

   .. py:attribute:: name
      :type: str
      :no-index:

   .. py:attribute:: mass
      :type: float

   .. py:attribute:: cog
      :type: tuple[float, float, float]

   .. py:attribute:: category
      :type: MassCategory

LoadingCondition
----------------

   .. py:class:: LoadingCondition(name)

      A complete loading condition with mass items and tank fill overrides.

      :param name: Name of the loading condition.

      .. py:method:: __init__(name)

         Creates a LoadingCondition.

   **Properties**

   .. py:attribute:: name
      :type: str
      :no-index:

      Name of the loading condition.

   **Mass Management**

   .. py:method:: add_mass(item)

      Add a MassItem.

      :param item: A MassItem instance.

   .. py:method:: add_mass_simple(name, mass, cog, category=None)

      Add a mass item by parameters (convenience).

      :param name: Identifier for the mass item.
      :param mass: Mass in kg.
      :param cog: Center of gravity (lcg, tcg, vcg) in meters.
      :param category: Optional MassCategory (default: Other).

   .. py:method:: remove_mass(name)

      Remove a mass item by name. 

      :param name: Name of the mass.
      :returns: True if found and removed.
      :rtype: bool

   .. py:method:: get_masses()

      Returns all mass items.

      :rtype: list[MassItem]

   .. py:method:: num_masses()

      Returns the number of mass items.

      :rtype: int

   **Tank Fill Overrides**

   .. py:method:: set_tank_fill(tank_name, fill_level)

      Set a tank fill override by fill level (0.0 to 1.0).

   .. py:method:: set_tank_fill_percent(tank_name, fill_percent)

      Set a tank fill override by percentage (0 to 100).

   .. py:method:: remove_tank_fill(tank_name)

      Remove a tank fill override.

      :returns: True if found and removed.
      :rtype: bool

   .. py:method:: num_tank_overrides()

      Returns the number of tank fill overrides.

      :rtype: int

   .. py:method:: get_tank_fills()

      Returns tank fill overrides as a dictionary.

      :returns: {name: fill_level}
      :rtype: dict[str, float]

   **Application & Calculation**

   .. py:method:: apply(vessel)

      Apply tank fill overrides to the vessel's tanks. Only tanks listed in overrides are modified. Other tanks keep their current fill level.

      :param vessel: Vessel instance.

   .. py:method:: total_displacement(vessel)

      Returns the total displacement (masses + tank fluid masses) in kg. Must be called after apply().

      :param vessel: Vessel instance.
      :rtype: float

   .. py:method:: total_cog(vessel)

      Returns the combined center of gravity in meters. Must be called after apply().

      :param vessel: Vessel instance.
      :returns: (lcg, tcg, vcg)
      :rtype: tuple[float, float, float]

   .. py:method:: resolve(vessel)

      Returns (total_displacement, (lcg, tcg, vcg)) in a single call. Must be called after apply().

      :param vessel: Vessel instance.
      :rtype: tuple[float, tuple[float, float, float]]

   .. py:method:: item_displacement()

      Returns the displacement of mass items only (excluding tank fluids) in kg.

      :rtype: float

   .. py:method:: item_cog()

      Returns the center of gravity of mass items only (lcg, tcg, vcg) in meters.

      :rtype: tuple[float, float, float]

   .. py:method:: resolve_items()

      Returns (item_displacement, (lcg, tcg, vcg)) in a single call. Use this for stability calculations to avoid double-counting tank masses.
      
      .. note::
         It is highly recommended to use the convenience methods directly on the calculators: :py:meth:`~navaltoolbox.StabilityCalculator.gz_curve_from_loading` and :py:meth:`~navaltoolbox.HydrostaticsCalculator.from_loading`. These methods automatically handle the save/restore process for tank configurations.

      :rtype: tuple[float, tuple[float, float, float]]

   **Serialization and Copying**

   .. py:method:: to_json()

      Serialize to JSON string.

      :rtype: str

   .. py:staticmethod:: from_json(json_str)

      Deserialize from JSON string.

      :rtype: LoadingCondition

   .. py:method:: save_json(path)

      Save to JSON file.

   .. py:staticmethod:: load_json(path)

      Load from JSON file.

   .. py:method:: copy(name=None)

      Create a copy, optionally with a new name.

      :param name: New name
      :rtype: LoadingCondition
