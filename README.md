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
- **Tanks**: Fill level management, free surface effects
- **Silhouettes**: Wind heeling calculations (DXF/VTK support)

## Installation

### Python Package

```bash
pip install navaltoolbox
```

### Rust Library

Add to your `Cargo.toml`:

```toml
[dependencies]
navaltoolbox = "0.1"
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
# Calculate hydrostatics
calc = HydrostaticsCalculator(vessel, water_density=1025.0)

# Option 1: At draft with VCG
state = calc.calculate_at_draft(5.0, vcg=6.0)
print(f"Volume: {state.volume:.1f} m³")
print(f"Waterplane Area: {state.waterplane_area:.1f} m²")
print(f"GMT (wet): {state.gmt:.3f} m")

# Option 2: Find draft for displacement
state_disp = calc.calculate_at_displacement(512500.0)
print(f"Draft: {state_disp.draft:.3f} m")

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
let state = calc.calculate_at_draft(5.0, 0.0, 0.0, 0.0)?;
println!("Volume: {} m³", state.volume);

// Calculate GZ curve
let stab = StabilityCalculator::new(&vessel, 1025.0);
let heels = vec![0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0];
let curve = stab.calculate_gz_curve(1000000.0, [50.0, 0.0, 5.0], &heels);
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

## Author

[Antoine ANCEAU](https://github.com/antoineanceau) · [Website](https://antoine.anceau.fr)
