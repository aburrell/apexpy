#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import

from glob import glob
from os import path, environ
from setuptools import setup, find_packages

# Include extensions only when not on readthedocs.org
if environ.get('READTHEDOCS', None) == 'True':
    extensions = []
else:
    from numpy.distutils.core import setup, Extension
    extensions = [
        Extension(name='apexpy.fortranapex',
                  sources=['src/fortranapex/magfld.f', 'src/fortranapex/apex.f',
                           'src/fortranapex/makeapexsh.f90',
                           'src/fortranapex/igrf.f90',
                           'src/fortranapex/apexsh.f90',
                           'src/fortranapex/checkapexsh.f90',
                           'src/fortranapex/fortranapex.pyf'])]

setup_kwargs = {'py_modules': [path.splitext(path.basename(pp))[0]
                               for pp in glob('src/*.py')],
                'ext_modules': extensions,
                'packages': find_packages(where='src')}

setup(**setup_kwargs)
