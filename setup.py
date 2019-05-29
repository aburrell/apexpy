#!/usr/bin/env python
# -*- encoding: utf-8 -*-
import os

# Include extensions only when not on readthedocs.org
if os.environ.get('READTHEDOCS', None) == 'True':
    from setuptools import setup
    from distutils.core import Extension
    extensions = []
else:
    from numpy.distutils.core import setup, Extension

    # add MinGW to Windows only, if needed
    if os.name == 'nt':
        sfn = os.path.join(os.path.dirname(__file__), 'setup.cfg')
        stxt = open(sfn).read()
        if '[build_ext]' not in stxt:
            with open(sfn, 'a') as f:
                f.write("\n[build_ext]\ncompiler = mingw32")

    extensions = [
        Extension(name='apexpy.fortranapex',
                  sources=['src/fortranapex/magfld.f', 'src/fortranapex/apex.f',
                           'src/fortranapex/makeapexsh.f90',
                           'src/fortranapex/apexsh.f90',
                           'src/fortranapex/checkapexsh.f90'])]

setup(ext_modules=extensions)
