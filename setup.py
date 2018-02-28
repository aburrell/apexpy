#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import, print_function

import io
import os
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import relpath
from os.path import splitext

from setuptools import find_packages

# on_rtd is whether we are on readthedocs.org
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if not on_rtd:
    from numpy.distutils.core import setup, Extension
    extensions = [
        Extension(name='apexpy.fortranapex',
                  sources=['src/fortranapex/magfld.f', 'src/fortranapex/apex.f',
                           'src/fortranapex/makeapexsh.f90',
                           'src/fortranapex/apexsh.f90',
                           'src/fortranapex/checkapexsh.f90'])]
else:
    from setuptools import setup
    from distutils.core import Extension
    extensions = []


def read(*names, **kwargs):
    return io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ).read()


if __name__ == "__main__":
    setup(
        name='apexpy',
        version='1.0.2',
        license='MIT',
        description='A Python wrapper for Apex coordinates',
        long_description='%s\n%s' % (read('README.rst'),
                                     re.sub(':[a-z]+:`~?(.*?)`', r'``\1``',
                                            read('CHANGELOG.rst'))),
        author='Christer van der Meeren; Angeline G. Burrell',
        author_email='agb073000@utdallas.edu',
        url='https://github.com/aburrell/apexpy',
        packages=find_packages('src'),
        package_dir={'': 'src'},
        py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
        package_data={'apexpy': ['apexsh.dat']},
        zip_safe=False,
        classifiers=[
            # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Operating System :: Unix',
            'Operating System :: POSIX',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: MacOS :: MacOS X',
            'Programming Language :: Python',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: Implementation :: CPython',
            'Topic :: Scientific/Engineering :: Physics',
            'Topic :: Utilities',
        ],
        keywords=[
            'apex',
            'modified apex',
            'quasi-dipole',
            'quasi dipole',
            'coordinates',
            'magnetic coordinates',
            'mlt',
            'magnetic local time',
            'conversion',
            'converting',
        ],
        install_requires=[
            'numpy',
        ],
        ext_modules=extensions,
        entry_points={
            'console_scripts': [
                'apexpy = apexpy.__main__:main',
            ]
        },
    )
