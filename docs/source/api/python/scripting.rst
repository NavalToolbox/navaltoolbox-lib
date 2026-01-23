Scripting
=========

Classes for Rhai scripting integration and criteria verification.

Script Engine
-------------

.. py:class:: ScriptEngine

   Rhai script execution engine.
   
   .. py:method:: run_script_file(path, context)
   
      Run a Rhai script from file.
      
      :param path: Path to .rhai script file
      :type path: str
      :param context: Data context for the script
      :type context: CriteriaContext
      :returns: Verification result
      :rtype: CriteriaResult
      :raises ValueError: If script execution fails
      
   .. py:method:: run_script(script, context)
   
      Run a Rhai script from string.
      
      :param script: Rhai script content
      :type script: str
      :param context: Data context for the script
      :type context: CriteriaContext
      :returns: Verification result
      :rtype: CriteriaResult
      
Criteria Context
----------------

.. py:class:: CriteriaContext

   Context for Rhai scripts, wrapping stability results.
   
   .. py:staticmethod:: from_result(result, vessel_name, loading_condition)
   
      Create a context from a CompleteStabilityResult.
      
      :param result: The stability calculation result
      :type result: CompleteStabilityResult
      :param vessel_name: Name of the vessel
      :type vessel_name: str
      :param loading_condition: Description of loading condition
      :type loading_condition: str
      :rtype: CriteriaContext
      
   .. py:method:: set_param(key, value)
   
      Set a parameter accessible to the script.
      
      :param key: Parameter name
      :type key: str
      :param value: Parameter value (str, float, or bool)
      :type value: Union[str, float, bool]

   .. py:method:: get_first_flooding_angle()
   
      Get the first angle where downflooding occurs.
      
      :returns: Angle in degrees or None
      :rtype: Optional[float]

   .. py:method:: find_equilibrium_angle(heeling_arm)
   
      Find the first stable equilibrium angle (where GZ = heeling_arm).
      
      :param heeling_arm: Heeling arm lever in meters
      :type heeling_arm: float
      :returns: Angle in degrees or None
      :rtype: Optional[float]

   .. py:method:: find_second_intercept(heeling_arm)
   
      Find the second intercept angle (unstable equilibrium).
      
      :param heeling_arm: Heeling arm lever in meters
      :type heeling_arm: float
      :returns: Angle in degrees or None
      :rtype: Optional[float]

Criteria Result
---------------

.. py:class:: CriteriaResult

   Result of a criteria verification script.
   
   .. py:attribute:: overall_pass
      :type: bool
      
      Whether all criteria passed.
      
   .. py:attribute:: pass_count
      :type: int
      
      Number of passing criteria.
      
   .. py:attribute:: fail_count
      :type: int
      
      Number of failing criteria.
      
   .. py:attribute:: criteria
      :type: list[CriterionResult]
      
      List of individual criterion results.
      
   .. py:attribute:: plots
      :type: list[str]
      
      List of plots data as JSON strings.

Criterion Result
----------------

.. py:class:: CriterionResult

   Result of a single criterion check.
   
   .. py:attribute:: name
      :type: str
      
   .. py:attribute:: status
      :type: str
      
      'PASS', 'FAIL', 'WARN', or 'N/A'.
      
   .. py:attribute:: actual_value
      :type: float
      
   .. py:attribute:: required_value
      :type: float
      
   .. py:attribute:: unit
      :type: str
      
   .. py:attribute:: margin
      :type: float
      
   .. py:attribute:: notes
      :type: str
      
   .. py:attribute:: plot_id
      :type: str
      
      ID of the associated plot, if any.
