# -*- coding: utf-8 -*-

import json
import os
import re

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'sphinx.ext.coverage',
              'sphinx.ext.viewcode',
              'sphinx.ext.githubpages',
              'sphinx.ext.napoleon',
              'sphinx.ext.extlinks',
              'autoapi.extension']

# Define common elements

source_suffix = '.rst'
master_doc = 'index'
project = 'ApexPy'
year = '2022'
zenodo = json.loads(open('../.zenodo.json').read())
author = ' and '.join([zcreator['name'] for zcreator in zenodo['creators']])
copyright = ', '.join([year, author])

# Get version number from __init__.py
regex = r"(?<=__version__..\s)\S+"
with open('../apexpy/__init__.py', 'r') as fin:
    text = fin.read()
match = re.findall(regex, text)
version = release = match[0].strip("'")

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
