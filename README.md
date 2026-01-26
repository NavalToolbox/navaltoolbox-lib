# NavalToolbox

High-performance naval architecture library written in **Rust** with **Python bindings**.

## Architecture

NavalToolbox is built as a **Rust library** (`navaltoolbox`) with optional Python bindings via PyO3/Maturin. This architecture provides:

- ⚡ **High performance**: Rust's zero-cost abstractions and memory safety
- 🐍 **Python convenience**: Easy-to-use Python API for rapid prototyping
- 🔒 **Type safety**: Compile-time guarantees in Rust
- 🚀 **Production ready**: Deploy as Rust library or Python package

## Features

- **Hull geometry**: Load STL/VTK files, transform, scale, export
- **Multi-hull support**: Catamarans, trimarans, arbitrary configurations
- **Hydrostatics**: Volume, COB vector, Waterplane ($A_{wp}$, LCF, $BM_{t/l}$), Wetted Surface, Midship Area, LOS, Coefficients ($C_b, C_m, C_p$), Stiffness Matrix, Free Surface ($GM_{dry/wet}$)
- **Stability**: GZ curve calculation with trim optimization and downflooding detection
- **Complete stability analysis**: Combines hydrostatics, GZ curve, and wind heeling data
- **Tanks**: Fill level management, free surface effects
- **Silhouettes**: Wind heeling calculations (DXF/VTK support)
- **Verification**: Rhai scripting engine for custom stability criteria (IMO, localized rules)
- **Visualization**: Interactive 3D visualization using Plotly (Vessel, Tanks, Opening, Hydrostatics)
- **Plotting**: Built-in 2D plotting utilities with Matplotlib integration

## Installation

### Python Package

```bash
pip install navaltoolbox
```

### Rust Library

Add to your `Cargo.toml`:

```toml
[dependencies]
navaltoolbox = "0.2"
```

## Quick Start

### Python

```python
from navaltoolbox import Hull, Vessel, HydrostaticsCalculator, StabilityCalculator

# Load a hull
hull = Hull("ship.stl")
print(f"Bounds: {hull.get_bounds()}")

# Create a vessel
vessel = Vessel(hull)

# Calculate hydrostatics
calc = HydrostaticsCalculator(vessel, water_density=1025.0)

# Option 1: At draft with VCG
state = calc.from_draft(5.0, vcg=6.0)
print(f"Volume: {state.volume:.1f} m³")
print(f"Waterplane Area: {state.waterplane_area:.1f} m²")
print(f"GMT (wet): {state.gmt:.3f} m")

# Option 2: Find draft for displacement
state_disp = calc.from_displacement(512500.0)
print(f"Draft: {state_disp.draft:.3f} m")

# Option 3: From perpendicular drafts (AP/FP)
state_drafts = calc.from_drafts(draft_ap=6.0, draft_fp=4.0)
print(f"Trim: {state_drafts.trim:.2f}°")

# Calculate GZ curve
stab = StabilityCalculator(vessel, water_density=1025.0)
heels = [0, 10, 20, 30, 40, 50, 60]
curve = stab.gz_curve(
    displacement_mass=1000000,
    cog=(50.0, 0.0, 5.0),
    heels=heels
)
for heel, gz in zip(curve.heels(), curve.values()):
    print(f"Heel: {heel}°, GZ: {gz:.3f}m")

# Complete stability analysis (hydrostatics + GZ + wind data)
result = stab.complete_stability(
    displacement_mass=1000000,
    cog=(50.0, 0.0, 5.0),
    heels=heels
)
print(f"GM0: {result.gm0:.3f}m")
print(f"Max GZ: {result.max_gz:.3f}m at {result.heel_at_max_gz}°")

# Scripting & Verification
from navaltoolbox import CriteriaContext, ScriptEngine, plotting

# Create context
ctx = CriteriaContext.from_result(result, "MV Example", "Departure")

# Run verification script (IMO A.749)
engine = ScriptEngine()
# Assuming the rule file is available in 'rules/'
# criteria = engine.run_script_file("rules/imo_a749_general.rhai", ctx)

# Plot results
# plotting.plot_criteria_result(criteria, show=True)
```

### Rust

```rust
use navaltoolbox::{Hull, Vessel, HydrostaticsCalculator, StabilityCalculator};

// Load a hull
let hull = Hull::from_stl("ship.stl")?;
println!("Bounds: {:?}", hull.get_bounds());

// Create a vessel
let vessel = Vessel::new(hull);

// Calculate hydrostatics
let calc = HydrostaticsCalculator::new(&vessel, 1025.0);
let state = calc.from_draft(5.0, 0.0, 0.0, None)?;
println!("Volume: {} m³", state.volume);

// Calculate from AP/FP drafts
let state_drafts = calc.from_drafts(6.0, 4.0, 0.0, None).unwrap();
println!("Trim: {:.2}°", state_drafts.trim);

// Calculate GZ curve
let stab = StabilityCalculator::new(&vessel, 1025.0);
let heels = vec![0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0];
let curve = stab.gz_curve(1000000.0, [50.0, 0.0, 5.0], &heels);
for point in &curve.points {
    println!("Heel: {}°, GZ: {:.3}m", point.heel, point.value);
}
```

## Development

### Building from Source

```bash
# Build Rust library
cd rust
cargo build --release

# Build Python package
cd python
maturin develop --release
```

## License

AGPL-3.0-or-later

## Disclaimer

NavalToolbox has been developed with care to ensure that all models and methods are correct, and that calculations reflect the most accurate results achievable with the implemented algorithms. However, **results must not be considered as a guarantee of performance**. The author cannot be held responsible for any inaccuracies in the calculations or for any consequences arising from the use of this software. Users are advised to independently verify critical calculations and to use this software as a tool to support, not replace, professional engineering judgment.

## Author

[Antoine ANCEAU](https://github.com/antoineanceau) · [Website](https://antoine.anceau.fr)
