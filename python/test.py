import navaltoolbox as nv

# Water density (kg/m³)
rho = 1025.0

# Main dimensions (m)
length = 20.0
breadth = 5.0
depth = 2.0
draft = 1.6  # draft = 1.5 -> submerged_centroid = (0, 0)

# Displacement (kg) and center of gravity (m)
displacement = length * breadth * draft * rho
cg = (0.5*length, 0.0, draft)

# Vessel
hull = nv.Hull.from_box(length, breadth, depth)
vessel = nv.Vessel(hull)

# Silhouette
silhouette = nv.Silhouette.from_points([
    (0.0, 0.0),
    (length, 0.0),
    (length, depth),
    (0.0, depth),
    (0.0, 0.0),
], "Windage")

vessel.add_silhouette(silhouette)

# Complete stability
stab_calc = nv.StabilityCalculator(vessel, water_density=rho)

result = stab_calc.complete_stability(
    displacement_mass=displacement,
    cog=cg,
    heels=[0.0]
)

# Submerged centroid
expected_submerged_centroid = (0.5 * length, 0.5 * draft)

print(f"Expected submerged centroid   = {expected_submerged_centroid}")
print(f"Calculated submerged centroid = {result.wind_data.submerged_centroid}")
