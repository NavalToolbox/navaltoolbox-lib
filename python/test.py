import os
import sys
import navaltoolbox as nv
from navaltoolbox.visualization import plot_hydrostatic_condition

# Water density (kg/m³)
rho = 1025.0

# Main dimensions (m)
length = 20.0
breadth = 5.0
depth = 2.0
draft = 1.0

# Vessel
hull = nv.Hull.from_box(length, breadth, depth)
vessel = nv.Vessel(hull)

# Displacement (kg) and center of gravity (m)
displacement = length * breadth * draft * rho
cg = (0.7*length, 0.0, draft)

# Hydrostatic condition
hydro_calc = nv.HydrostaticsCalculator(vessel, water_density=rho)

state = hydro_calc.from_displacement(displacement, cog=cg)

fig_hydro = plot_hydrostatic_condition(
    vessel,
    draft=state.draft,  # I get the right result with state.draft_ap
    trim=state.trim,
    heel=state.heel,
    enable_opacity_slider=False,
    show_axes=False,
)

output_hydro = "dtmb_hydro.html"
fig_hydro.write_html(output_hydro)
if sys.platform == "darwin":  # macOS
    os.system(f"open {output_hydro}")
