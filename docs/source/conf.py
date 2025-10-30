# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'GRiD'
copyright = '2025, Zachary Caterer'
author = 'Zachary Caterer'
release = 'v0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_click",
    "sphinx_copybutton"
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown'
}

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
        "color-brand-primary": "#D91E36",
        "color-brand-content": "#D91E36",
        "color-link": "#D91E36",
    },
    "dark_css_variables": {
        "color-brand-primary": "#D91E36",
        "color-brand-content": "#D91E36",
        "color-link": "#D91E36",
    },
}