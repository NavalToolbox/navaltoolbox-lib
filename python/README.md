# NavalToolbox

High-performance naval architecture library for Python, powered by Rust.

[![PyPI version](https://badge.fury.io/py/navaltoolbox.svg)](https://pypi.org/project/navaltoolbox/)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

## Overview

NavalToolbox provides fast and accurate naval architecture calculations through a Python API, with a high-performance Rust backend that handles the heavy computational work.

## Features

- ⚓ **Hull Geometry**: Load and manipulate ship hulls from STL/VTK files
- 🚢 **Multi-hull Support**: Catamarans, trimarans, and arbitrary configurations
- 📊 **Hydrostatics**: Volume, COB, Waterplane ($A_{wp}$, LCF, $BM_t$, $BM_l$), Free Surface Correction ($GM_{dry}/GM_{wet}$)
- ⚖️ **Stability Analysis**: GZ curve calculation with automatic trim optimization
- 🌊 **Downflooding Detection**: Automatic detection of submerged openings
- 🛢️ **Tank Management**: Fill levels, free surface effects, sounding tables
- ⚖️ **Loading Conditions**: Compose mass inventories and tank overrides to build operational profiles
- 💨 **Wind Heeling**: Silhouette-based wind calculations (DXF/VTK support)
- 📝 **Scriptable Verification**: Rhai scripting engine for custom stability criteria
- 🧊 **3D Visualization**: Interactive vessel and hydrostatic visualization with Plotly
- 📈 **Plotting**: Built-in 2D plotting utilities with Matplotlib integration
- ⚡ **High Performance**: Rust backend with Python convenience

## Installation

```bash
pip install navaltoolbox
```

**Requirements:**
- Python 3.9 or higher
- No additional dependencies required (all native code included in wheels)

## Quick Start

### Loading a Hull

```python
from navaltoolbox import Hull

# Load hull from STL file
hull = Hull("ship.stl")

# Get hull dimensions
bounds = hull.get_bounds()
loa = bounds[1] - bounds[0]  # Length overall
boa = bounds[3] - bounds[2]  # Beam overall

print(f"LOA: {loa:.2f}m, BOA: {boa:.2f}m")
print(f"Triangles: {hull.num_triangles()}")
```

### Hydrostatics Calculation

```python
from navaltoolbox import Hull, Vessel, HydrostaticsCalculator

# Create vessel
hull = Hull("ship.stl")
vessel = Vessel(hull)

# Calculate hydrostatics at a given draft
calc = HydrostaticsCalculator(vessel, water_density=1025.0)

# Option 1: At draft with VCG (computes stability)
state = calc.from_draft(draft=5.0, vcg=7.0)

print(f"Volume: {state.volume:.1f} m³")
print(f"Displacement: {state.displacement:.0f} kg")
print(f"COB: ({state.cob[0]:.2f}, {state.cob[1]:.2f}, {state.cob[2]:.2f})")
print(f"Waterplane Area: {state.waterplane_area:.1f} m²")
print(f"LCF: {state.lcf:.2f} m")
print(f"GMT (wet): {state.gmt:.3f} m")
print(f"GMT (dry): {state.gmt_dry:.3f} m")

# Option 2: Find draft for displacement
state_disp = calc.from_displacement(512500.0)
print(f"Draft for {state_disp.displacement:.0f}kg: {state_disp.draft:.3f} m")
```

### Stability Analysis (GZ Curve)

```python
from navaltoolbox import Hull, Vessel, StabilityCalculator

# Create vessel and calculator
hull = Hull("ship.stl")
vessel = Vessel(hull)
calc = StabilityCalculator(vessel, water_density=1025.0)

# Calculate GZ curve
displacement_mass = 8635000.0  # kg
cog = (71.67, 0.0, 7.555)      # LCG, TCG, VCG in meters
heels = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]

curve = calc.gz_curve(displacement_mass, cog, heels)

# Display results
print("Heel (°)  |  GZ (m)")
print("-" * 25)
for heel, gz in zip(curve.heels(), curve.values()):
    print(f"{heel:7.1f}  |  {gz:7.3f}")

# Get max GZ
gz_values = curve.values()
max_gz = max(gz_values)
max_idx = gz_values.index(max_gz)
max_heel = heels[max_idx]
print(f"\nMax GZ: {max_gz:.3f}m at {max_heel}°")
```

### Loading Conditions

```python
from navaltoolbox import Hull, Vessel, Tank, LoadingCondition, MassCategory

# Create vessel with a tank
vessel = Vessel(Hull("ship.stl"))
vessel.add_tank(Tank.from_box("FO_1", 20.0, 30.0, -5.0, 5.0, 0.0, 2.0, 1000.0))

# Define loading condition
lc = LoadingCondition("Arrival")
lc.add_mass_simple("Lightship", 5000000.0, (40.0, 0.0, 5.0), MassCategory.lightship())
lc.set_tank_fill_percent("FO_1", 50.0)

# Use convenience methods directly to seamlessly merge the solid masses and tank configurations:
from navaltoolbox import StabilityCalculator, HydrostaticsCalculator

stab_calc = StabilityCalculator(vessel)
curve = stab_calc.gz_curve_from_loading(lc, heels=[0, 10, 20, 30])
print(f"Calculated GZ max: {max(curve.values()):.3f} m")

hydro_calc = HydrostaticsCalculator(vessel)
state = hydro_calc.from_loading(lc)
print(f"Equilibrium Draft: {state.draft:.3f} m")
```

## Documentation

For more detailed documentation, examples, and API reference, visit:
- **Documentation**: [GitHub Pages](https://navaltoolbox.github.io/navaltoolbox-lib/)
- **GitHub Repository**: [NavalToolbox/navaltoolbox-lib](https://github.com/NavalToolbox/navaltoolbox-lib)
- **Issue Tracker**: [GitHub Issues](https://github.com/NavalToolbox/navaltoolbox-lib/issues)

## Performance

NavalToolbox is built with performance in mind:

- Written in **Rust** for maximum speed and memory safety
- Efficient mesh operations using `parry3d`
- Parallel processing where applicable
- Zero-copy data transfer between Python and Rust

Example benchmark (DTMB 5415 hull, 3436 triangles):
- Load STL: ~10ms
- Hydrostatics calculation: ~50ms
- GZ curve (13 points): ~650ms

## Use Cases

NavalToolbox is suitable for:

- 🎓 **Naval architecture education**: Teaching hydrostatics and stability
- 🔬 **Research**: Rapid prototyping of new methods and algorithms
- 🏭 **Engineering**: Production stability calculations and analysis
- 🤖 **Optimization**: Integration with optimization frameworks
- 📊 **Batch processing**: Analyzing multiple design variants

## Requirements and Compatibility

- **Python**: 3.9, 3.10, 3.11, 3.12, 3.13
- **Operating Systems**: 
  - macOS (Intel and Apple Silicon)
  - Linux (x86_64, aarch64)
  - Windows (x86_64)
- **File Formats**:
  - Hull geometry: STL (binary and ASCII), VTK
  - Silhouettes: DXF, VTK

## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests on GitHub.

## License

This project is licensed under the **GNU Affero General Public License v3.0 or later (AGPL-3.0-or-later)**.

This means:
- ✅ Free to use, modify, and distribute
- ✅ Can be used in commercial projects
- ⚠️ Source code must be made available under AGPL
- ⚠️ Network use is considered distribution (AGPL provision)

See the [LICENSE](https://github.com/NavalToolbox/navaltoolbox-lib/blob/main/LICENSE) file for details.

## Disclaimer

NavalToolbox has been developed with care to ensure that all models and methods are correct, and that calculations reflect the most accurate results achievable with the implemented algorithms. However, **results must not be considered as a guarantee of performance**. The author cannot be held responsible for any inaccuracies in the calculations or for any consequences arising from the use of this software. Users are advised to independently verify critical calculations and to use this software as a tool to support, not replace, professional engineering judgment.

## Citation

If you use NavalToolbox in your research, please cite:

```bibtex
@software{navaltoolbox2026,
  author = {Anceau, Antoine},
  title = {NavalToolbox: High-performance naval architecture library},
  year = {2026},
  url = {https://github.com/NavalToolbox/navaltoolbox-lib}
}
```

## Author

**Antoine ANCEAU**
- GitHub: [@antoineanceau](https://github.com/antoineanceau)
- Website: [antoine.anceau.fr](https://antoine.anceau.fr)

## Support

- 📖 **Documentation**: [GitHub Pages](https://navaltoolbox.github.io/navaltoolbox-lib)
- 🐛 **Bug Reports**: [Open an issue](https://github.com/NavalToolbox/navaltoolbox-lib/issues)
- 💬 **Discussions**: [GitHub Discussions](https://github.com/NavalToolbox/navaltoolbox-lib/discussions)
- 📧 **Email**: For private inquiries

---

**Note**: This package uses a Rust backend for high performance. Pre-built wheels are provided for common platforms. If a wheel is not available for your platform, the package will attempt to build from source (requires Rust toolchain).
