import os
import sys

sys.path.insert(0, os.path.abspath("../../src"))

# project info
project = "clibas"
copyright = "2026, Alexander Vinogradov"
author = "Alexander Vinogradov"
release = "0.4.2"

# general config
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_parser",
    "nbsphinx",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# MyST config
myst_enable_extensions = [
    "deflist",
    "colon_fence",
    "attrs_inline",
]

# napoleon settings for docstring parsing
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False

# autodoc config
autodoc_typehints = "description"
autodoc_member_order = "bysource"

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

# html styling
html_theme = "furo"
html_static_path = ["_static"]

html_css_files = ["custom.css"]
html_title = "clibas docs"
