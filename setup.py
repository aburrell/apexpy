#!/usr/bin/env python
# -*- encoding: utf-8 -*-

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
                  sources=['fortranapex/igrf.f90',
                           'fortranapex/magfld.f90',
                           'fortranapex/apex.f90',
                           'fortranapex/makeapexsh.f90',
                           'fortranapex/apexsh.f90',
                           'fortranapex/checkapexsh.f90',
                           'fortranapex/fortranapex.pyf'])]

setup_kwargs = {'py_modules': [path.splitext(path.basename(pp))[0]
                               for pp in glob('*.py')],
                'ext_modules': extensions,
                'packages': find_packages()}

setup(**setup_kwargs)
