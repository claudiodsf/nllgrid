# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from unittest.mock import Mock
sys.modules['numpy'] = Mock()
sys.modules['scipy'] = Mock()
sys.modules['scipy.ndimage'] = Mock()
sys.modules['pyproj'] = Mock()
sys.path.insert(0, os.path.abspath('..'))
from nllgrid._version import get_versions #NOQA
__version__ = get_versions()['version']


# -- Project information -----------------------------------------------------

project = 'NLLGrid'
copyright = '2015-2024, Claudio Satriano'
author = 'Claudio Satriano'
release = __version__
version = __version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'autoclasstoc',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinx_autodoc_typehints',
]
autodoc_mock_imports = ['numpy', 'scipy', 'pyproj']
autodoc_member_order = 'bysource'
autodoc_default_options = {
    'members': True,
    'special-members': False,
    'private-members': False,
    'inherited-members': True,
    'undoc-members': True,
    'exclude-members': '__weakref__',
}
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('http://docs.scipy.org/doc/scipy/', None),
}
napoleon_use_param = True
napoleon_preprocess_types = True
napoleon_type_aliases = {
    'array_like': ':term:`array_like`',
    'ArrayLike': ':py:data:`~numpy.typing.ArrayLike`'
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
