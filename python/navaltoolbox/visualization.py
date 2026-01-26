# Copyright (C) 2026 Antoine ANCEAU
#
# This file is part of navaltoolbox.
#
# navaltoolbox is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""
Visualization module for navaltoolbox using Plotly.
"""

from typing import Optional, List, Tuple, Any
import numpy as np
import plotly.graph_objects as go  # type: ignore

from . import Vessel, Hull


def plot_vessel_3d(
    vessel: Vessel,
    show_hulls: bool = True,
    show_tanks: bool = True,
    show_silhouettes: bool = True,
    show_openings: bool = True,
    opacity_hull: float = 0.5,
    opacity_tank: float = 0.3,
    title: str = "Vessel Visualization",
    enable_opacity_slider: bool = True,
    show_axes: bool = True,
) -> go.Figure:
    """
    Create an interactive 3D plot of the vessel.

    Args:
        vessel: The vessel to visualize.
        show_hulls: Whether to show hull meshes.
        show_tanks: Whether to show tank meshes.
        show_silhouettes: Whether to show silhouette profiles.
        show_openings: Whether to show downflooding openings.
        opacity_hull: Opacity for hull meshes (0.0 to 1.0).
        opacity_tank: Opacity for tank meshes (0.0 to 1.0).
        title: Title of the plot.
        enable_opacity_slider: Whether to add a slider to control hull opacity.
        show_axes: Whether to show axes, grid, and background.

    Returns:
        A Plotly Figure object.
    """
    fig = go.Figure()

    hull_indices = []

    # 1. Hulls
    if show_hulls:
        for i, hull in enumerate(vessel.get_hulls()):
            _add_hull_trace(fig, hull, f"Hull {i+1}", opacity_hull)
            # Track the index of the added trace
            hull_indices.append(len(fig.data) - 1)

    # 2. Tanks
    if show_tanks:
        for tank in vessel.get_tanks():
            # Tank Container (Wireframe or transparent mesh)
            if hasattr(tank, "get_vertices") and hasattr(tank, "get_faces"):
                verts = tank.get_vertices()
                faces = tank.get_faces()

                if verts and faces:
                    x, y, z = zip(*verts)
                    i_idx, j_idx, k_idx = zip(*faces)

                    fig.add_trace(
                        go.Mesh3d(
                            x=x,
                            y=y,
                            z=z,
                            i=i_idx,
                            j=j_idx,
                            k=k_idx,
                            color="lightgray",
                            opacity=opacity_tank,
                            name=f"{tank.name} (Tank)",
                            showscale=False,
                            alphahull=0,  # Use exact mesh
                        )
                    )

            # Fluid Volume (Higher opacity)
            if hasattr(tank, "get_fluid_vertices") and hasattr(
                tank, "get_fluid_faces"
            ):
                f_verts = tank.get_fluid_vertices()
                f_faces = tank.get_fluid_faces()

                if f_verts and f_faces:
                    fx, fy, fz = zip(*f_verts)
                    fi, fj, fk = zip(*f_faces)

                    fig.add_trace(
                        go.Mesh3d(
                            x=fx,
                            y=fy,
                            z=fz,
                            i=fi,
                            j=fj,
                            k=fk,
                            color="blue",
                            opacity=min(opacity_tank * 2.0, 0.9),
                            name=f"{tank.name} (Fluid)",
                            showscale=False,
                        )
                    )

            # Keep COG marker for reference
            cog = tank.center_of_gravity
            fig.add_trace(
                go.Scatter3d(
                    x=[cog[0]],
                    y=[cog[1]],
                    z=[cog[2]],
                    mode="markers",
                    marker=dict(size=3, color="black"),
                    name=f"{tank.name} COG",
                    showlegend=False,
                )
            )

    # 3. Silhouettes (Lines)
    if show_silhouettes:
        for sil in vessel.get_silhouettes():
            points = sil.get_points()
            if not points:
                continue

            x, y, z = zip(*points)
            # Close the loop if needed
            if sil.is_closed() and (x[0] != x[-1] or z[0] != z[-1]):
                x = x + (x[0],)
                y = y + (y[0],)
                z = z + (z[0],)

            fig.add_trace(
                go.Scatter3d(
                    x=x,
                    y=y,
                    z=z,
                    mode="lines",
                    line=dict(color="green", width=2),
                    name=f"Silhouette: {sil.name}",
                )
            )

    # 4. Openings (Markers)
    if show_openings:
        for opening in vessel.get_openings():
            points = opening.get_points()
            if not points:
                continue

            x, y, z = zip(*points)

            # If contour, draw lines
            if len(points) > 1:
                # Close loop
                x = x + (x[0],)
                y = y + (y[0],)
                z = z + (z[0],)
                mode = "lines+markers"
            else:
                mode = "markers"

            fig.add_trace(
                go.Scatter3d(
                    x=x,
                    y=y,
                    z=z,
                    mode=mode,
                    marker=dict(size=4, color="red"),
                    line=dict(color="red", width=2),
                    name=f"Opening: {opening.name}",
                )
            )

    # 5. Reference Points (AP/FP)
    # Using Y=0, Z=0 for reference
    ap_x = vessel.ap
    fp_x = vessel.fp

    fig.add_trace(
        go.Scatter3d(
            x=[ap_x, fp_x],
            y=[0, 0],
            z=[0, 0],
            mode="markers+text",
            marker=dict(size=5, color="black", symbol="diamond"),
            text=["AP", "FP"],
            textposition="top center",
            name="Perpendiculars",
        )
    )

    # Set axis labels and equal aspect ratio
    bounds = vessel.get_bounds()
    max_range = max(
        bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]
    )
    center_x = (bounds[0] + bounds[1]) / 2
    center_y = (bounds[2] + bounds[3]) / 2
    center_z = (bounds[4] + bounds[5]) / 2

    # Hack to force aspect ratio in Plotly
    # We add invisible points to set the range
    fig.add_trace(
        go.Scatter3d(
            x=[center_x - max_range / 2, center_x + max_range / 2],
            y=[center_y - max_range / 2, center_y + max_range / 2],
            z=[center_z - max_range / 2, center_z + max_range / 2],
            mode="markers",
            marker=dict(size=0, opacity=0),
            showlegend=False,
        )
    )

    # Sliders
    sliders = []
    if enable_opacity_slider and hull_indices:
        steps = []
        for op in np.arange(0, 1.1, 0.1):
            step = dict(
                method="restyle",
                args=[{"opacity": op}, hull_indices],
                label=f"{op:.1f}",
            )
            steps.append(step)

        sliders.append(
            dict(
                active=int(opacity_hull * 10),
                currentvalue={"prefix": "Hull Opacity: "},
                pad={"t": 50},
                steps=steps,
            )
        )

    scene_args: dict[str, Any] = dict(aspectmode="cube")

    if show_axes:
        scene_args.update(
            dict(
                xaxis_title="X (Longitudinal)",
                yaxis_title="Y (Transverse)",
                zaxis_title="Z (Vertical)",
            )
        )
    else:
        scene_args.update(
            dict(
                xaxis=dict(
                    visible=False, showgrid=False, showbackground=False
                ),
                yaxis=dict(
                    visible=False, showgrid=False, showbackground=False
                ),
                zaxis=dict(
                    visible=False, showgrid=False, showbackground=False
                ),
            )
        )

    fig.update_layout(
        title=title,
        scene=scene_args,
        sliders=sliders,
        updatemenus=[
            # 1. Projection (Left)
            dict(
                type="buttons",
                direction="right",
                buttons=[
                    dict(
                        label="Persp 👁️",
                        method="relayout",
                        args=["scene.camera.projection.type", "perspective"],
                    ),
                    dict(
                        label="Ortho 📐",
                        method="relayout",
                        args=["scene.camera.projection.type", "orthographic"],
                    ),
                ],
                pad={"r": 5, "t": 5},
                showactive=True,
                x=0.0,
                xanchor="left",
                y=1.1,
                yanchor="top",
                bgcolor="rgba(255, 255, 255, 0.7)",
            ),
            # 2. Views (Right, with gap)
            dict(
                type="buttons",
                direction="right",
                buttons=_get_camera_buttons(),
                pad={"r": 5, "t": 5},
                showactive=True,
                x=0.25,  # Gap from projection buttons
                xanchor="left",
                y=1.1,
                yanchor="top",
                bgcolor="rgba(255, 255, 255, 0.7)",
            ),
        ],
    )

    return fig


def plot_hydrostatic_condition(
    vessel: Vessel,
    draft: float,
    trim: float = 0.0,
    heel: float = 0.0,
    show_hulls: bool = True,
    show_tanks: bool = True,
    show_silhouettes: bool = True,
    show_openings: bool = True,
    opacity_hull: float = 0.5,
    opacity_tank: float = 0.3,
    title: str = "Hydrostatic Condition",
    enable_opacity_slider: bool = True,
    cog: Optional[Tuple[float, float, float]] = None,
    show_axes: bool = True,
) -> go.Figure:
    """
    Visualize the vessel at a specific hydrostatic condition (floating in
    water).
    The waterplane is fixed at Z=0.
    Vessel is transformed to match drafts.

    Args:
        vessel: The vessel to visualize.
        draft: Draft at midship (or reference point) in meters.
        trim: Trim angle in degrees (positive = stern down usually).
        heel: Heel angle in degrees (positive = starboard down usually).
        show_hulls: Whether to show hull meshes.
        show_tanks: Whether to show tank meshes.
        show_silhouettes: Whether to show silhouette profiles.
        show_openings: Whether to show downflooding openings.
        opacity_hull: Opacity for hull meshes (0.0 to 1.0).
        opacity_tank: Opacity for tank meshes (0.0 to 1.0).
        title: Title of the plot.
        enable_opacity_slider: Whether to add a slider to control hull opacity.
        cog: Optional Center of Gravity (x, y, z) to display.
        show_axes: Whether to show axes, grid, and background.
    """
    fig = go.Figure()
    hull_indices = []

    # Calculate Transformation Matrix
    # We want to move Ship -> World such that water is at Z=0.
    # Rotation Order: Pitch (Trim) then Roll (Heel).
    # R = Ry(trim) * Rx(heel)

    h_rad = np.deg2rad(heel)
    t_rad = np.deg2rad(trim)

    ch, sh = np.cos(h_rad), np.sin(h_rad)
    ct, st = np.cos(t_rad), np.sin(t_rad)

    # Rx (Heel around X)
    Rx = np.array([[1, 0, 0], [0, ch, -sh], [0, sh, ch]])

    # Ry (Trim around Y)
    Ry = np.array([[ct, 0, st], [0, 1, 0], [-st, 0, ct]])

    # Combined Rotation R = Rx @ Ry
    R_total = Rx @ Ry

    # Transform function
    def transform_points(points):
        if len(points) == 0:
            return []
        pts = np.array(points)
        # Apply rotation (pts is Nx3, R_total is 3x3)
        # (R @ v)^T = v^T @ R^T
        rotated = pts @ R_total.T

        # Apply translation
        # Z translation: -draft.
        # Ideally, we should pivot around the flotation center,
        # but here we just translate vertically.
        translated = rotated + [0, 0, -draft]
        return translated

    # 0. Waterplane
    bounds = vessel.get_bounds()
    # Big plane
    L = bounds[1] - bounds[0]
    B = bounds[3] - bounds[2]
    plane_size = max(L, B) * 1.5

    # Grid of points for water
    wx = np.linspace(-plane_size / 2, plane_size / 2, 2) + (
        bounds[0] + bounds[1]
    ) / 2
    wy = np.linspace(-plane_size / 2, plane_size / 2, 2)
    WX, WY = np.meshgrid(wx, wy)
    WZ = np.zeros_like(WX)

    fig.add_trace(
        go.Surface(
            x=WX,
            y=WY,
            z=WZ,
            opacity=0.3,
            colorscale=[[0, "blue"], [1, "blue"]],
            showscale=False,
            name="Waterplane",
        )
    )

    # Add waterplane outline (so it's visible in orthographic side views)
    cx = (bounds[0] + bounds[1]) / 2
    cy = 0
    half_size = plane_size / 2

    # Square aligned with axes
    outline_x = [
        cx - half_size,
        cx + half_size,
        cx + half_size,
        cx - half_size,
        cx - half_size,
    ]
    outline_y = [
        cy - half_size,
        cy - half_size,
        cy + half_size,
        cy + half_size,
        cy - half_size,
    ]
    outline_z = [0, 0, 0, 0, 0]

    fig.add_trace(
        go.Scatter3d(
            x=outline_x,
            y=outline_y,
            z=outline_z,
            mode="lines",
            line=dict(color="blue", width=4),
            name="Waterline Level",
        )
    )

    # 1. Hulls
    if show_hulls:
        for i, hull in enumerate(vessel.get_hulls()):
            if hasattr(hull, "get_vertices") and hasattr(hull, "get_faces"):
                verts = hull.get_vertices()
                faces = hull.get_faces()
                if verts and faces:
                    t_verts = transform_points(verts)
                    x, y, z = t_verts[:, 0], t_verts[:, 1], t_verts[:, 2]
                    i_idx, j_idx, k_idx = zip(*faces)

                    fig.add_trace(
                        go.Mesh3d(
                            x=x,
                            y=y,
                            z=z,
                            i=i_idx,
                            j=j_idx,
                            k=k_idx,
                            color="gray",
                            opacity=opacity_hull,
                            name=f"Hull {i+1}",
                            showscale=False,
                            flatshading=True,
                        )
                    )
                    hull_indices.append(len(fig.data) - 1)

    # 2. Tanks
    if show_tanks:
        for tank in vessel.get_tanks():
            # Container
            if hasattr(tank, "get_vertices"):
                verts = tank.get_vertices()
                faces = tank.get_faces()
                if verts:
                    t_verts = transform_points(verts)
                    x, y, z = t_verts[:, 0], t_verts[:, 1], t_verts[:, 2]
                    i_idx, j_idx, k_idx = zip(*faces)
                    fig.add_trace(
                        go.Mesh3d(
                            x=x,
                            y=y,
                            z=z,
                            i=i_idx,
                            j=j_idx,
                            k=k_idx,
                            color="lightgray",
                            opacity=opacity_tank,
                            name=f"{tank.name} (Tank)",
                            showscale=False,
                        )
                    )

            # Fluid (Use heel/trim-aware getter)
            if hasattr(tank, "get_fluid_vertices"):
                # Pass vessel heel/trim to get the fluid mesh that is
                # 'horizontal' in world but 'tilted' in ship frame.
                verts = tank.get_fluid_vertices(heel=heel, trim=trim)
                faces = tank.get_fluid_faces(heel=heel, trim=trim)
                if verts:
                    # Apply SHIP transformation to place it in world
                    t_verts = transform_points(verts)
                    x, y, z = t_verts[:, 0], t_verts[:, 1], t_verts[:, 2]
                    i_idx, j_idx, k_idx = zip(*faces)
                    fig.add_trace(
                        go.Mesh3d(
                            x=x,
                            y=y,
                            z=z,
                            i=i_idx,
                            j=j_idx,
                            k=k_idx,
                            color="blue",
                            opacity=0.8,
                            name=f"{tank.name} (Fluid)",
                            showscale=False,
                        )
                    )

    # 3. Silhouettes (Lines)
    if show_silhouettes:
        for sil in vessel.get_silhouettes():
            points = sil.get_points()
            if not points:
                continue

            # Close the loop if needed (in ship frame)
            pts_nx3 = np.array(points)
            if sil.is_closed():
                # Check distance between first and last point
                if np.linalg.norm(pts_nx3[0] - pts_nx3[-1]) > 1e-6:
                    pts_nx3 = np.vstack([pts_nx3, pts_nx3[0]])

            # Transform to world
            t_points = transform_points(pts_nx3)
            x, y, z = t_points[:, 0], t_points[:, 1], t_points[:, 2]

            fig.add_trace(
                go.Scatter3d(
                    x=x,
                    y=y,
                    z=z,
                    mode="lines",
                    line=dict(color="green", width=2),
                    name=f"Silhouette: {sil.name}",
                )
            )

    # 4. Openings (Markers)
    if show_openings:
        for opening in vessel.get_openings():
            points = opening.get_points()
            if not points:
                continue

            pts_nx3 = np.array(points)

            # If contour, draw lines
            if len(points) > 1:
                # Close loop
                if np.linalg.norm(pts_nx3[0] - pts_nx3[-1]) > 1e-6:
                    pts_nx3 = np.vstack([pts_nx3, pts_nx3[0]])

                t_points = transform_points(pts_nx3)
                x, y, z = t_points[:, 0], t_points[:, 1], t_points[:, 2]

                mode = "lines+markers"
            else:
                t_points = transform_points(pts_nx3)
                x, y, z = t_points[:, 0], t_points[:, 1], t_points[:, 2]
                mode = "markers"

            fig.add_trace(
                go.Scatter3d(
                    x=x,
                    y=y,
                    z=z,
                    mode=mode,
                    marker=dict(size=4, color="red"),
                    line=dict(color="red", width=2),
                    name=f"Opening: {opening.name}",
                )
            )

    # 5. Reference Points
    ap_pt = transform_points([(vessel.ap, 0, 0)])[0]
    fp_pt = transform_points([(vessel.fp, 0, 0)])[0]

    fig.add_trace(
        go.Scatter3d(
            x=[ap_pt[0], fp_pt[0]],
            y=[ap_pt[1], fp_pt[1]],
            z=[ap_pt[2], fp_pt[2]],
            mode="markers+text",
            marker=dict(size=5, color="black", symbol="diamond"),
            text=["AP", "FP"],
            name="Perpendiculars",
        )
    )

    # 6. Center of Gravity
    if cog:
        cog_pt = transform_points([cog])[0]
        fig.add_trace(
            go.Scatter3d(
                x=[cog_pt[0]],
                y=[cog_pt[1]],
                z=[cog_pt[2]],
                mode="markers+text",
                marker=dict(size=8, color="orange", symbol="diamond"),
                text=["COG"],
                textposition="top center",
                name="Center of Gravity",
            )
        )

    # Opacity Slider
    sliders = []
    if enable_opacity_slider and hull_indices:
        steps = []
        for op in np.arange(0, 1.1, 0.1):
            step = dict(
                method="restyle",
                args=[{"opacity": op}, hull_indices],
                label=f"{op:.1f}",
            )
            steps.append(step)

        sliders.append(
            dict(
                active=int(opacity_hull * 10),
                currentvalue={"prefix": "Hull Opacity: "},
                pad={"t": 50},
                steps=steps,
            )
        )

    scene_args: dict[str, Any] = dict(aspectmode="data")

    if show_axes:
        scene_args.update(
            dict(
                xaxis_title="X",
                yaxis_title="Y",
                zaxis_title="Z",
            )
        )
    else:
        scene_args.update(
            dict(
                xaxis=dict(
                    visible=False, showgrid=False, showbackground=False
                ),
                yaxis=dict(
                    visible=False, showgrid=False, showbackground=False
                ),
                zaxis=dict(
                    visible=False, showgrid=False, showbackground=False
                ),
            )
        )

    fig.update_layout(
        title=title,
        scene=scene_args,
        sliders=sliders,
        updatemenus=[
            # 1. Projection (Left)
            dict(
                type="buttons",
                direction="right",
                buttons=[
                    dict(
                        label="Persp 👁️",
                        method="relayout",
                        args=["scene.camera.projection.type", "perspective"],
                    ),
                    dict(
                        label="Ortho 📐",
                        method="relayout",
                        args=["scene.camera.projection.type", "orthographic"],
                    ),
                ],
                pad={"r": 5, "t": 5},
                showactive=True,
                x=0.0,
                xanchor="left",
                y=1.1,
                yanchor="top",
                bgcolor="rgba(255, 255, 255, 0.7)",
            ),
            # 2. Views (Right, with gap)
            dict(
                type="buttons",
                direction="right",
                buttons=_get_camera_buttons(),
                pad={"r": 5, "t": 5},
                showactive=True,
                x=0.25,  # Gap from projection buttons
                xanchor="left",
                y=1.1,
                yanchor="top",
                bgcolor="rgba(255, 255, 255, 0.7)",
            ),
        ],
    )
    return fig


def _get_camera_buttons() -> List[dict]:
    """Returns a list of camera view buttons for Plotly."""
    return [
        dict(
            label="ISO",
            method="relayout",
            args=[
                {
                    "scene.camera.eye": dict(x=1.25, y=1.25, z=1.25),
                    "scene.camera.up": dict(x=0, y=0, z=1),
                    "scene.camera.projection.type": "perspective",
                }
            ],
        ),
        dict(
            label="STBD",  # Starboard
            method="relayout",
            args=[
                {
                    "scene.camera.eye": dict(x=0, y=-2.5, z=0),
                    "scene.camera.up": dict(x=0, y=0, z=1),
                    "scene.camera.projection.type": "orthographic",
                }
            ],
        ),
        dict(
            label="PORT",  # Port
            method="relayout",
            args=[
                {
                    "scene.camera.eye": dict(x=0, y=2.5, z=0),
                    "scene.camera.up": dict(x=0, y=0, z=1),
                    "scene.camera.projection.type": "orthographic",
                }
            ],
        ),
        dict(
            label="TOP",
            method="relayout",
            args=[
                {
                    "scene.camera.eye": dict(x=0, y=0, z=2.5),
                    "scene.camera.up": dict(x=1, y=0, z=0),
                    "scene.camera.projection.type": "orthographic",
                }
            ],
        ),
        dict(
            label="BOW",
            method="relayout",
            args=[
                {
                    "scene.camera.eye": dict(x=2.5, y=0, z=0),
                    "scene.camera.up": dict(x=0, y=0, z=1),
                    "scene.camera.projection.type": "orthographic",
                }
            ],
        ),
        dict(
            label="AFT",  # Stern
            method="relayout",
            args=[
                {
                    "scene.camera.eye": dict(x=-2.5, y=0, z=0),
                    "scene.camera.up": dict(x=0, y=0, z=1),
                    "scene.camera.projection.type": "orthographic",
                }
            ],
        ),
    ]


def _add_hull_trace(fig: go.Figure, hull: Hull, name: str, opacity: float):
    """Helper to add hull mesh to figure."""
    if hasattr(hull, "get_vertices") and hasattr(hull, "get_faces"):
        verts = hull.get_vertices()
        faces = hull.get_faces()

        if not verts or not faces:
            return

        x, y, z = zip(*verts)
        i, j, k = zip(*faces)

        fig.add_trace(
            go.Mesh3d(
                x=x,
                y=y,
                z=z,
                i=i,
                j=j,
                k=k,
                color="gray",
                opacity=opacity,
                name=name,
                showscale=False,
                flatshading=True,
            )
        )
