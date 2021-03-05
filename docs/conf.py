# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os
import re


extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.extlinks'
]

source_suffix = '.rst'
master_doc = 'index'
project = u'Apex Python library'
year = u'2015'
author = u'Christer van der Meeren, Angeline G. Burrell'
copyright = '{0}, {1}'.format(year, author)
# Get version number from __init__.py
here = os.path.abspath(os.path.dirname(__file__))
regex = "(?<=__version__..\s)\S+"
with open(os.path.join(here,'../src/apexpy/__init__.py'),'r', encoding='utf-8') as f:
    text = f.read()
match = re.findall(regex,text)
version = release = match[0].strip("'")
# on_rtd is whether we are on readthedocs.org
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:  # only import and set the theme if we're building docs locally
    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

pygments_style = 'trac'
templates_path = ['.']
html_use_smartypants = True
html_last_updated_fmt = '%b %d, %Y'
html_split_index = True
html_sidebars = {
   '**': ['searchbox.html', 'globaltoc.html', 'sourcelink.html'],
}
html_short_title = '%s-%s' % (project, version)
autodoc_member_order='bysource'
napoleon_use_ivar=True
napoleon_use_rtype=False
napoleon_use_param=False

extlinks = {'doi': ('http://dx.doi.org/%s', 'doi:')}
