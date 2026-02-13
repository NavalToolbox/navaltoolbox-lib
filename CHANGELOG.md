# Changelog

All notable changes to NavalToolbox will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.2] - 2026-02-13

### Added
- **hydrostatics**: Added `cog` (Total COG) and `vessel_cog` (Vessel COG) fields to `HydrostaticState`
- **stability**: Added `cog` and `vessel_cog` to `StabilityPoint` for detailed center of gravity tracking (Total vs Vessel mass)

### Fixed
- **loading**: Correctly handle OCS transformation in DXF silhouettes for inverted extrusion directions (fixes coordinate inversion issue)

## [0.6.1] - 2026-02-11

### Added
- **testing**: Comprehensive integrity test suite for silhouette profiles in `python/tests`
- **stability**: Warning notification when wind data is skipped during stability analysis
- **loading**: Support for loading silhouettes from CSV/TXT point files

### Fixed
- **loading**: Correctly load DXF `AcDb2dPolyline` entities with OCS normal (0,1,0) (fixes "no wind data" error)
- **stability**: Resolve issue where `StabilityCalculator` used a stale Vessel copy, ignoring subsequently added silhouettes
- **validation**: Add warnings for invalid silhouette geometries (zero area, open loops)

### Changed
- **python**: Refactor `StabilityCalculator` binding to hold a reference to the `Vessel` object, ensuring dynamic updates
- **docs**: Updated API documentation and type stubs for silhouette components

## [0.6.0] - 2026-02-11

### Added
- **stability**: Support for `tank_options` in `gz_curve` and `complete_stability` methods
- **hydrostatics**: Added `vessel_displacement` and `tank_displacement` fields to `HydrostaticState`

### Changed
- **hydrostatics**: **BREAKING** Removed `hull_displacement` field from `HydrostaticState` (replaced by `vessel_displacement`)
- **hydrostatics**: `displacement` field now explicitly represents the Total Displacement (Vessel + Tanks)
- **stability**: `complete_stability` now correctly delegates displacement calculation to avoid double-counting tank mass

## [0.5.1] - 2026-02-09

### Fixed
- **python**: Fix type inference error in `from_displacement` with `tank_options`

## [0.5.0] - 2026-02-09

### Added
- **hydrostatics**: Introduce `TankOptions` for fine-grained control over tank mass and FSM inclusion
- **hydrostatics**: Add `hull_displacement` (buoyancy) and `tank_mass` fields to `HydrostaticState`
- **tanks**: Implement `FSMMode` (Actual, Maximum, Fixed) for flexible free surface moment calculations
- **tanks**: Exact sorting of waterplane vertices for robust FSM calculation on complex hull shapes

### Changed
- **hydrostatics**: `from_draft`, `from_drafts`, and `from_displacement` now accept optional `tank_options` argument

## [0.4.2] - 2026-02-07

### Fixed
- **hydrostatics**: Implement robust fallback solver with Coordinate Descent + Bisection for extreme equilibrium cases
- **hydrostatics**: Restore VCG parameter handling in equilibrium solvers

## [0.4.1] - 2026-02-06

### Fixed
- **hydrostatics**: Correct heel sign convention (TcG>0 → negative heel)
- **hydrostatics**: Compute equilibrium heel/trim from off-center COG
- **hydrostatics**: Preserve COG in `from_displacement` when VCG provided

### Documentation
- Add coordinate system conventions to userguide
- Fix Y axis convention in Tank docstrings

## [0.4.0] - 2026-02-04

### Added
- Complete hydrostatic properties: LOS (Length Overall Submerged), Wetted Surface Area
- Sectional area curve calculation
- Freeboard calculation with deck edges
- Appendages support for additional volume/center corrections

### Fixed
- Python and Rust lint errors

### Documentation
- Update visualization tutorial for Appendages and Deck Edges
- Add sectional_areas and freeboard to HydrostaticState documentation

## [0.3.0] - 2026-01-19

### Added
- GZ curve calculation with trim optimization
- Downflooding openings detection during stability calculations
- Silhouette wind heeling calculations (DXF/VTK support)
- Rhai scripting engine for custom stability criteria

### Changed
- Refactored hydrostatics fields to use public properties

## [0.2.0] - 2026-01-08

### Added
- Multi-hull support (catamarans, trimarans)
- Tank management with free surface corrections
- Complete stability analysis combining hydrostatics, GZ curve, and wind data

## [0.1.0] - 2026-01-01

### Added
- Initial release
- Hull geometry loading (STL/VTK files)
- Basic hydrostatics: Volume, COB, Waterplane properties
- Python bindings via PyO3/Maturin

---
*Generated with [git-cliff](https://git-cliff.org/) and manually curated*
