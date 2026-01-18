# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Add Python package to path for autodoc
sys.path.insert(0, os.path.abspath("../../python"))

# -- Project information -----------------------------------------------------

project = "NavalToolbox"
copyright = "2026, Antoine ANCEAU"
author = "Antoine ANCEAU"
release = "0.1.0"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",  # Google/NumPy style docstrings
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",  # Copy button for code blocks
]

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

html_theme = "pydata_sphinx_theme"

html_theme_options = {
    "logo": {
        "image_light": "_static/logo.png",
        "image_dark": "_static/logo.png",
        "text": "NavalToolbox",
    },
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/NavalToolbox/navaltoolbox-lib",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/navaltoolbox/",
            "icon": "fa-solid fa-box",
        },
    ],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "show_toc_level": 2,
    "navigation_with_keys": True,
}

html_favicon = "_static/favicon.png"

html_static_path = ["_static"]
html_css_files = ["custom.css"]

# -- Extension configuration -------------------------------------------------

# Napoleon settings (NumPy style docstrings)
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_use_param = True
napoleon_use_rtype = True

# Autodoc settings
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

# Intersphinx
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}
