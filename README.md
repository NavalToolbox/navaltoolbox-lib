# NavalToolbox

High-performance naval architecture library for hydrostatics, stability, and tank calculations.

## Features

- **Hull geometry**: Load STL files, transform, scale, export
- **Multi-hull support**: Catamarans, trimarans
- **Hydrostatics**: Volume, center of buoyancy, waterplane properties
- **Stability**: GZ curve calculation with trim optimization
- **Tanks**: Fill level management, free surface effects
- **Silhouettes**: Wind heeling calculations (DXF/VTK support)

## Installation

```bash
pip install navaltoolbox
```

## Quick Start

```python
from navaltoolbox import Hull, Vessel, HydrostaticsCalculator, StabilityCalculator

# Load a hull
hull = Hull("ship.stl")
print(f"Bounds: {hull.get_bounds()}")

# Create a vessel
vessel = Vessel(hull)

# Calculate hydrostatics
calc = HydrostaticsCalculator(vessel, water_density=1025.0)
state = calc.calculate_at_draft(5.0)
print(f"Volume: {state.volume} m³")

# Calculate GZ curve
stab = StabilityCalculator(vessel, water_density=1025.0)
heels = [0, 10, 20, 30, 40, 50, 60]
curve = stab.calculate_gz_curve(
    displacement_mass=1000000,
    cog=(50.0, 0.0, 5.0),
    heels=heels
)
for heel, gz in zip(curve.heels(), curve.values()):
    print(f"Heel: {heel}°, GZ: {gz:.3f}m")
```

## License

AGPL-3.0-or-later

## Author

[Antoine ANCEAU](https://github.com/antoineanceau) · [Website](https://antoine.anceau.fr)
