# -*- coding: utf-8 -*-
"""Configuration for apexpy documentation."""
import json
import os
from pyproject_parser import PyProject

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'sphinx.ext.coverage',
              'sphinx.ext.viewcode',
              'sphinx.ext.githubpages',
              'sphinx.ext.napoleon',
              'sphinx.ext.extlinks',
              'autoapi.extension']

# General information about the project.
info = PyProject.load("../pyproject.toml")

# Define common elements
source_suffix = '.rst'
master_doc = 'index'
project = 'ApexPy'
year = '2024'
zenodo = json.loads(open('../.zenodo.json').read())
author = ' and '.join([zcreator['name'] for zcreator in zenodo['creators']])
copyright = ', '.join([year, author])
version = release = info.project['version'].base_version

# Configure autoapi
autoapi_type = 'python'
autoapi_dirs = ['../apexpy']
autoapi_keep_files = True
autoapi_root = 'autoapi/generated'

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = os.path.join(os.path.abspath('.'), 'apexpy.png')

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'sphinx_rtd_theme'
html_theme_path = ["_themes", ]

pygments_style = 'trac'
templates_path = ['.']
html_use_smartypants = True
html_last_updated_fmt = '%b %d, %Y'
html_split_index = True
html_sidebars = {'**': ['searchbox.html', 'globaltoc.html', 'sourcelink.html']}
html_short_title = '-'.join([project, version])
autodoc_member_order = 'bysource'
autodoc_mock_imports = ['apexpy']
napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = False

extlinks = {'doi': ('http://dx.doi.org/%s', 'doi:%s')}
