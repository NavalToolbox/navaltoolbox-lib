Get started
===========

This guide covers installation, getting started, and important information for using NavalToolbox.

.. grid:: 1 1 2 2
    :gutter: 3

    .. grid-item-card:: 📦 Installation
        :link: installation
        :link-type: doc

        Install NavalToolbox from PyPI or build from source.

    .. grid-item-card:: 🚀 Quick Start
        :link: quickstart
        :link-type: doc

        Get started quickly with basic examples.

    .. grid-item-card:: 📐 Coordinate System
        :link: coordinates
        :link-type: doc

        Reference for axes, rotations, and sign conventions.

Features
--------

- ⚓ **Hull Geometry**: Load and manipulate ship hulls from STL/VTK files
- 🚢 **Multi-hull Support**: Catamarans, trimarans, and arbitrary configurations
- 📊 **Hydrostatics**: Volume, COB, Waterplane (:math:`A_{wp}`, LCF, :math:`BM_t`, :math:`BM_l`), Free Surface Correction (:math:`GM_{dry}/GM_{wet}``)
- ⚖️ **Stability Analysis**: GZ curve calculation with automatic trim optimization
- 🌊 **Downflooding Detection**: Automatic detection of submerged openings
- 🛢️ **Tank Management**: Fill levels, free surface effects, sounding tables
- 💨 **Wind Heeling**: Silhouette-based wind calculations (DXF/VTK support)
- ⚡ **High Performance**: Rust backend with Python convenience

License
-------

This project is licensed under the **GNU Affero General Public License v3.0 or later (AGPL-3.0-or-later)**: https://github.com/NavalToolbox/navaltoolbox-lib/blob/main/LICENSE

Disclaimer
----------

.. warning::

   NavalToolbox has been developed with care to ensure that all models and methods 
   are correct, and that calculations reflect the most accurate results achievable 
   with the implemented algorithms. 
   
   However, **results must not be considered as a guarantee of performance**. 
   
   The author cannot be held responsible for any inaccuracies in the calculations 
   or for any consequences arising from the use of this software. Users are advised 
   to independently verify critical calculations and to use this software as a tool 
   to support, not replace, professional engineering judgment.

.. toctree::
   :maxdepth: 2
   :hidden:

   installation
   quickstart
   coordinates
   scripting
