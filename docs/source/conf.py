# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import re
import os
import sys
from importlib.metadata import version as _version, PackageNotFoundError

try:
    _ver = _version("GRiD")
except PackageNotFoundError:
    _ver = "1.0.0"

if _ver != "1.0.0":
    print(f"Building documentation for GRiD version {_ver}")
    _ver = "1.0.0"

project = 'GRiD'
copyright = '2026, Zachary Caterer'
author = 'Zachary Caterer'
release = _ver

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",
    "sphinx_click",
    "sphinx_copybutton"
]

myst_enable_extensions = [
    "dollarmath",   # enables $...$ and $$...$$ math syntax in .md files
    "amsmath",      # enables \begin{equation} ... \end{equation} environments
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown'
}

bibtex_bibfiles = ["../../paper.bib"]

# -- Paths ---------------------------------------------------------------
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# -- Autodoc options -----------------------------------------------------
autodoc_default_options = {
    'members': True,           # Document all members (functions, classes)
    'undoc-members': True,     # Include members without docstrings
    'private-members': True,   # Include _private members
    'special-members': True,   # Include __special__ methods
    'show-inheritance': True,  # Show class inheritance
}

autodoc_typehints = 'description'  # Show type hints in description

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
html_logo = "../../assets/grid_logo.png"
html_favicon = "_static/favicon.ico"


html_theme_options = {
    "light_css_variables": {
        "color-brand-primary": "#F47B7B",
        "color-brand-content": "#F47B7B",
        "color-link": "#F47B7B",
    },
    "dark_css_variables": {
        "color-brand-primary": "#F47B7B",
        "color-brand-content": "#F47B7B",
        "color-link": "#F47B7B",
    },
}
