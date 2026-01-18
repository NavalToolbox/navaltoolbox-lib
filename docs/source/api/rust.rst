Rust API Reference
==================

The Rust API is documented on `docs.rs <https://docs.rs/navaltoolbox>`_.

Quick Reference
---------------

Core Types
^^^^^^^^^^

- ``Hull`` - Hull geometry with STL loading and transformations
- ``Vessel`` - Container for hulls and tanks
- ``HydrostaticsCalculator`` - Hydrostatic property calculations
- ``StabilityCalculator`` - GZ curve calculations
- ``Tank`` - Tank with fluid management

Usage Example
^^^^^^^^^^^^^

.. code-block:: rust

    use navaltoolbox::{Hull, Vessel, StabilityCalculator};
    use std::path::Path;

    fn main() -> Result<(), Box<dyn std::error::Error>> {
        // Load hull
        let hull = Hull::from_stl(Path::new("ship.stl"))?;
        let vessel = Vessel::new(hull);

        // Calculate GZ curve
        let calc = StabilityCalculator::new(&vessel, 1025.0);
        let heels = vec![0.0, 10.0, 20.0, 30.0, 40.0];
        let cog = [71.67, 0.0, 7.555];

        let curve = calc.calculate_gz_curve(8635000.0, cog, &heels);

        for point in &curve.points {
            println!("Heel {:5.1}°: GZ = {:.3}m", point.heel, point.value);
        }

        Ok(())
    }

Modules
-------

hull
^^^^
Hull geometry loading and manipulation.

vessel
^^^^^^
Vessel container with multi-hull and tank support.

hydrostatics
^^^^^^^^^^^^
Hydrostatic calculations (volume, center of buoyancy, etc.).

stability
^^^^^^^^^
Stability calculations (GZ curves with trim optimization).

tanks
^^^^^
Tank management with free surface effects.

mesh
^^^^
Mesh utilities (STL loading, clipping, transformations).

Full Documentation
------------------

See the complete Rust API documentation at:

- **docs.rs**: https://docs.rs/navaltoolbox
- **Source**: https://github.com/NavalToolbox/navaltoolbox-lib
