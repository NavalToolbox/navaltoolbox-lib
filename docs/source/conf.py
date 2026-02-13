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
release = "0.6.4"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",  # Google/NumPy style docstrings
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",  # Copy button for code blocks
    "sphinx_design",  # Grid and cards
]

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

html_theme = "pydata_sphinx_theme"

# NavalToolbox color scheme
# Light: primary=#125C8B, secondary=#0094D3
# Dark: primary=#E8F0F5, secondary=#6DCFF6, background=#1E3A4C
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
            "icon": "fa-brands fa-python",
        },
        {
            "name": "crates.io",
            "url": "https://crates.io/crates/navaltoolbox",
            "icon": "fa-solid fa-cube",
        },
    ],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "show_toc_level": 2,
    "navigation_with_keys": True,
    # Hide source links and prev/next buttons
    "show_prev_next": False,
    # Custom colors for NavalToolbox
    "primary_sidebar_end": [],
    "pygments_light_style": "default",
    "pygments_dark_style": "monokai",
}

html_favicon = "_static/favicon.png"

html_static_path = ["_static"]
html_css_files = ["custom.css"]

# Hide "Show Source" link
html_show_sourcelink = False

# Custom context for templates
html_context = {
    "default_mode": "auto",
}

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
