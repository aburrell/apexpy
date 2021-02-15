#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import

import io
import re
from glob import glob
from os import path, environ

from setuptools import find_packages

# Include extensions only when not on readthedocs.org
if environ.get('READTHEDOCS', None) == 'True':
    from setuptools import setup
    from distutils.core import Extension
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

def read(*names, **kwargs):
    return io.open(
        path.join(path.dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ).read()


if __name__ == "__main__":
    setup(
        name='apexpy',
        version='1.0.4',
        license='MIT',
        description='A Python wrapper for Apex coordinates',
        long_description='%s\n%s' % (read('README.rst'),
                                     re.sub(':[a-z]+:`~?(.*?)`', r'``\1``',
                                            read('CHANGELOG.rst'))),
        author='Christer van der Meeren; Angeline G. Burrell',
        author_email='angeline.burrell@nrl.navy.mil',
        url='https://github.com/aburrell/apexpy',
        packages=find_packages('src'),
        package_dir={'': 'src'},
        py_modules=[path.splitext(path.basename(pp))[0]
                    for pp in glob('src/*.py')],
        package_data={'apexpy': ['apexsh.dat', 'igrfcoeffs.txt']},
        zip_safe=False,
        classifiers=[
            # complete classifier list:
            # http://pypi.python.org/pypi?%3Aaction=list_classifiers
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
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
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
