NavalToolbox Documentation
===========================

**NavalToolbox** is a high-performance naval architecture library for hydrostatics, 
stability, and tank calculations. Built in Rust with Python bindings.

.. grid:: 1 1 2 2
    :gutter: 2

    .. grid-item-card:: 🚀 Getting Started
        :link: quickstart
        :link-type: doc

        New to NavalToolbox? Start here to learn the basics.

    .. grid-item-card:: 📖 API Reference
        :link: api/python
        :link-type: doc

        Complete reference for all Python classes and methods.

    .. grid-item-card:: 🎓 Tutorials
        :link: tutorials/hydrostatics
        :link-type: doc

        In-depth tutorials for hydrostatics and stability calculations.

    .. grid-item-card:: ⚙️ Installation
        :link: installation
        :link-type: doc

        How to install NavalToolbox for Python or Rust.

Features
--------

- **Fast**: Rust core with Python bindings via PyO3
- **Accurate**: Validated against DTMB 5415 reference hull (< 3.5cm GZ error)
- **Complete**: Hull, Vessel, Hydrostatics, Stability, Tanks modules
- **Dual API**: Use from Python or as a native Rust library

Quick Example
-------------

.. code-block:: python

    from navaltoolbox import Hull, Vessel, StabilityCalculator

    # Load hull geometry
    hull = Hull("ship.stl")
    vessel = Vessel(hull)

    # Calculate GZ curve
    stab = StabilityCalculator(vessel, water_density=1025.0)
    curve = stab.calculate_gz_curve(
        displacement_mass=8635000,  # kg
        cog=(71.67, 0.0, 7.555),    # LCG, TCG, VCG
        heels=[0, 10, 20, 30, 40, 50, 60]
    )

    for heel, gz in zip(curve.heels(), curve.values()):
        print(f"Heel {heel:5.1f}°: GZ = {gz:.3f}m")

.. toctree::
   :maxdepth: 2
   :hidden:

   installation
   quickstart
   api/python
   api/rust
   tutorials/hydrostatics
   tutorials/stability
