# Copyright (C) 2026 Antoine ANCEAU
#
# This file is part of navaltoolbox.

import json
import logging
from typing import List, Optional, Dict, Any, TYPE_CHECKING

try:
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

if TYPE_CHECKING:
    from . import CriteriaResult


def plot_criteria_result(
    result: "CriteriaResult",
    show: bool = True,
    plot_id: Optional[str] = None,
    save_to: Optional[str] = None
) -> List["Figure"]:
    """
    Plot the graphs defined in a CriteriaResult using Matplotlib.

    Args:
        result: The CriteriaResult object from a script execution.
        show: If True, calls plt.show() after creating figures.
        plot_id: If provided, only plots the graph with this specific ID.
        save_to: Path to save the plot. If multiple plots are generated,
                 the plot ID will be appended to the filename.

    Returns:
        List of created matplotlib Figure objects.
    """
    if not HAS_MATPLOTLIB:
        raise ImportError(
            "Matplotlib is required for plotting. "
            "Install it with 'pip install matplotlib'."
        )

    figures = []
    generated_ids = []

    # Identify relevant criteria for each plot
    criteria_by_plot: Dict[str, List] = {}
    for crit in result.criteria:
        pid = crit.plot_id or "main"
        if pid not in criteria_by_plot:
            criteria_by_plot[pid] = []
        criteria_by_plot[pid].append(crit)

    # Parse and plot each graph
    for plot_json in result.plots:
        try:
            plot_data = json.loads(plot_json)
        except json.JSONDecodeError:
            logging.error("Failed to decode plot JSON")
            continue

        pid = plot_data.get("id", "main")

        # Filter if requested
        if plot_id and pid != plot_id:
            continue

        fig = _render_single_plot(plot_data, criteria_by_plot.get(pid, []))
        figures.append(fig)
        generated_ids.append(pid)

    if save_to and figures:
        import os
        base, ext = os.path.splitext(save_to)
        if not ext:
            ext = ".png"

        if len(figures) == 1:
            figures[0].savefig(save_to)
            logging.info(f"Saved plot to {save_to}")
        else:
            for fig, pid in zip(figures, generated_ids):
                filename = f"{base}_{pid}{ext}"
                fig.savefig(filename)
                logging.info(f"Saved plot {pid} to {filename}")

    if show and figures:
        plt.show()

    return figures


def _render_single_plot(plot_data: Dict[str, Any], criteria: List) -> "Figure":
    """Render a single plot data dictionary to a Matplotlib figure."""

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))
    if fig.canvas.manager:
        fig.canvas.manager.set_window_title(plot_data.get("title", "Plot"))

    ax.set_title(plot_data.get("title", ""))
    ax.set_xlabel(plot_data.get("x_label", ""))
    ax.set_ylabel(plot_data.get("y_label", ""))
    ax.grid(True, linestyle='--', alpha=0.7)

    # Render elements
    elements = plot_data.get("elements", [])
    for elem in elements:
        etype = elem.get("type")
        name = elem.get("name", "")
        color = elem.get("color", "blue")
        style = elem.get("style", "-")

        # Map style strings to matplotlib
        mpl_style = "-"
        if style == "dashed":
            mpl_style = "--"
        elif style == "dotted":
            mpl_style = ":"

        if etype == "Curve":
            x = elem.get("x", [])
            y = elem.get("y", [])
            ax.plot(x, y, label=name, color=color, linestyle=mpl_style)

        elif etype == "HorizontalLine":
            y = elem.get("y", 0.0)
            xmin = elem.get("x_min")
            xmax = elem.get("x_max")
            if xmin is not None and xmax is not None:
                ax.hlines(
                    y, xmin, xmax, label=name, color=color, linestyle=mpl_style
                )
            else:
                ax.axhline(y, label=name, color=color, linestyle=mpl_style)

        elif etype == "VerticalLine":
            x = elem.get("x", 0.0)
            ymin = elem.get("y_min")
            ymax = elem.get("y_max")
            if ymin is not None and ymax is not None:
                ax.vlines(
                    x, ymin, ymax, label=name, color=color, linestyle=mpl_style
                )
            else:
                ax.axvline(x, label=name, color=color, linestyle=mpl_style)

        elif etype == "Point":
            x = elem.get("x", 0.0)
            y = elem.get("y", 0.0)
            marker = elem.get("marker", "o")
            # Convert marker name if needed
            if marker == "circle":
                marker = "o"
            ax.plot(
                x, y, marker=marker, color=color, label=name, linestyle='None'
            )

        elif etype == "FilledArea":
            x = elem.get("x", [])
            y1 = elem.get("y_lower", [])
            y2 = elem.get("y_upper", [])
            alpha = elem.get("alpha", 0.3)
            ax.fill_between(x, y1, y2, color=color, alpha=alpha, label=name)

    # Add legend
    ax.legend()

    # Criteria summary could be added here
    # pass
    return fig
