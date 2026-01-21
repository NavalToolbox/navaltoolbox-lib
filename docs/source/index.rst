NavalToolbox Documentation
===========================

.. raw:: html

   <div class="hero-section">
       <h1>NavalToolbox</h1>
       <p>High-performance naval architecture library for hydrostatics and stability calculations</p>
       <div class="hero-buttons">
           <a href="userguide/quickstart.html" class="hero-btn hero-btn-primary">Get Started</a>
           <a href="userguide/installation.html" class="hero-btn hero-btn-secondary">Install</a>
           <a href="api/index.html" class="hero-btn hero-btn-secondary">API Reference</a>
       </div>
   </div>

**NavalToolbox** is a Rust library with Python bindings designed for naval architects, 
marine engineers, and researchers who need accurate and fast hydrostatic and stability calculations.

Key Features
------------

.. raw:: html

   <div class="feature-grid">
       <div class="feature-card">
           <h3>⚡ High Performance</h3>
           <p>Rust core with parallel processing (Rayon) delivers blazing-fast calculations. GZ curves computed in seconds.</p>
       </div>
       <div class="feature-card">
           <h3>🎯 Accurate</h3>
           <p>Validated against DTMB 5415 reference hull with <3.5cm GZ error. Wall-sided formula verification included.</p>
       </div>
       <div class="feature-card">
           <h3>🐍 Python API</h3>
           <p>Easy-to-use Python bindings via PyO3. Ideal for scripts, Jupyter notebooks, and rapid prototyping.</p>
       </div>
       <div class="feature-card">
           <h3>📐 Complete</h3>
           <p>Hull, Vessel, Tanks, Hydrostatics, Stability modules. Supports multi-hull configurations.</p>
       </div>
   </div>

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
   :caption: Get started

   userguide/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Tutorials

   tutorials/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: API Reference

   api/index
