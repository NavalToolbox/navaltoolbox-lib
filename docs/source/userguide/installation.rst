Installation
============

Python Package
--------------

Install from PyPI :

.. code-block:: bash

    pip install navaltoolbox

Development Installation
------------------------

For development, clone the repository and install with maturin:

.. code-block:: bash

    git clone https://github.com/NavalToolbox/navaltoolbox-lib.git
    cd navaltoolbox/navaltoobox-lib/python

    # Install maturin if needed
    pip install maturin

    # Build and install in development mode
    maturin develop --release

Rust Library
------------

Add to your ``Cargo.toml``:

.. code-block:: toml

    [dependencies]
    navaltoolbox = "0.1"

Requirements
------------

- Python >= 3.9
- Rust >= 1.75 (for building from source)
