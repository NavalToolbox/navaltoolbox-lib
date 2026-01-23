Plotting
========

Utilities for visualizing criteria results.

.. py:module:: navaltoolbox.plotting

.. py:function:: plot_criteria_result(result, show=True, plot_id=None, save_to=None)

   Plot the graphs defined in a CriteriaResult using Matplotlib.

   :param result: The CriteriaResult object from a script execution.
   :type result: CriteriaResult
   :param show: If True, calls plt.show() after creating figures.
   :type show: bool
   :param plot_id: If provided, only plots the graph with this specific ID.
   :type plot_id: str, optional
   :param save_to: Path to save the plot. If multiple plots are generated, the plot ID will be appended to the filename.
   :type save_to: str, optional
   :returns: List of created matplotlib Figure objects.
   :rtype: list[matplotlib.figure.Figure]
   :raises ImportError: If matplotlib is not installed.
