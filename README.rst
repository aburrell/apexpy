========
Overview
========

|docs| |version| |doi|

This is a Python wrapper for the Apex fortran library by
Emmert et al. [2010] [1]_, which allows converting between geodetic, modified
apex, and quasi-dipole coordinates as well as getting modified apex and
quasi-dipole base vectors (Richmond [1995] [2]_). The geodetic system used here
is WGS84. MLT calculations are also included. The package is free software
(MIT license).

Quick start
===========

Install from PyPI using ``pip``::

    pip install apexpy

This assumes that the same version of libgfortran is installed in the same
location as when the pip wheel was built (if a wheel was used). If not, you may
have trouble importing apexpy.  If you run into trouble, try the command::

    pip install --no-binary :apexpy: apexpy

which requires both libgfortran and gfortran to be installed on your system.

Conversion is done by creating an ``Apex`` object and using its methods to
perform the desired calculations. Some simple examples::

    >>> from apexpy import Apex
    >>> from __future__ import print_function
    >>> A = Apex(date=2015.3)  # datetime objects are also supported
    >>> # geo to apex, scalar input
    >>> mlat, mlon = A.convert(60, 15, 'geo', 'apex', height=300)
    >>> print("{:.12f}, {:.12f}".format(mlat, mlon))
    57.477310180664, 93.590156555176
    >>> # apex to geo, array input
    >>> glat, glon = A.convert([90, -90], 0, 'apex', 'geo', height=0)
    >>> print(["{:.12f}, {:.12f}".format(ll, glon[i]) for i,ll in enumerate(glat)])
    ['83.103820800781, -84.526657104492', '-74.388252258301, 125.736274719238']
    >>> # geo to MLT
    >>> import datetime as dt
    >>> mlat, mlt = A.convert(60, 15, 'geo', 'mlt', datetime=dt.datetime(2015, 2, 10, 18, 0, 0))
    >>> print("{:.12f}, {:.12f}".format(mlat, mlt))
    56.598316192627, 19.107861709595
    >>> # can also convert magnetic longitude to mlt
    >>> mlt = A.mlon2mlt(120, dt.datetime(2015, 2, 10, 18, 0, 0))
    >>> print("{:.2f}".format(mlt))
    20.90

If you don't know or use Python, you can also use the command line. See details in the full documentation.

Documentation
=============

https://apexpy.readthedocs.org/

References
==========

.. [1] Emmert, J. T., A. D. Richmond, and D. P. Drob (2010),
       A computationally compact representation of Magnetic-Apex
       and Quasi-Dipole coordinates with smooth base vectors,
       J. Geophys. Res., 115(A8), A08322,
       `doi:10.1029/2010JA015326 <http://dx.doi.org/10.1029/2010JA015326>`_.

.. [2] Richmond, A. D. (1995), Ionospheric Electrodynamics Using
       Magnetic Apex Coordinates, Journal of geomagnetism and
       geoelectricity, 47(2), 191–212,
       `doi:10.5636/jgg.47.191 <http://dx.doi.org/10.5636/jgg.47.191>`_.

Badges
======

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |travis| |appveyor|
        | |coveralls| |requires|
        | |codeclimate| |scrutinizer| |codacy|
    * - package
      - | |version| |supported-versions|
        | |wheel| |supported-implementations|

.. |docs| image:: https://readthedocs.org/projects/apexpy/badge/?style=flat
    :target: https://readthedocs.org/projects/apexpy
    :alt: Documentation Status

.. |travis| image:: https://travis-ci.org/aburrell/apexpy.svg?branch=main
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/aburrell/apexpy

.. |appveyor| image:: https://ci.appveyor.com/api/projects/status/github/aburrell/apexpy?branch=main&svg=true
    :alt: AppVeyor Build Status
    :target: https://ci.appveyor.com/project/aburrell/apexpy

.. |requires| image:: https://requires.io/github/aburrell/apexpy/requirements.svg?branch=main
     :alt: Requirements Status
     :target: https://requires.io/github/aburrell/apexpy/requirements/?branch=main

.. |coveralls| image:: https://coveralls.io/repos/github/aburrell/apexpy/badge.svg?branch=main
    :alt: Coverage Status
    :target: https://coveralls.io/github/aburrell/apexpy?branch=main

.. |codacy| image:: https://api.codacy.com/project/badge/Grade/7d4c1a6c60e747ca95cdf97746c39cda
   :alt: Codacy Badge
   :target: https://app.codacy.com/gh/aburrell/apexpy?utm_source=github.com&utm_medium=referral&utm_content=aburrell/apexpy&utm_campaign=Badge_Grade

.. |codeclimate| image:: https://api.codeclimate.com/v1/badges/da1d972dee790da595f8/maintainability.svg
   :target: https://codeclimate.com/github/aburrell/apexpy
   :alt: CodeClimate Quality Status

.. |version| image:: https://img.shields.io/pypi/v/apexpy.svg?style=flat
    :alt: PyPI Package latest release
    :target: https://pypi.python.org/pypi/apexpy

.. |downloads| image:: https://img.shields.io/pypi/dm/apexpy.svg?style=flat
    :alt: PyPI Package monthly downloads
    :target: https://pypi.python.org/pypi/apexpy

.. |wheel| image:: https://img.shields.io/pypi/wheel/apexpy.svg?style=flat
    :alt: PyPI Wheel
    :target: https://pypi.python.org/pypi/apexpy

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/apexpy.svg?style=flat
    :alt: Supported versions
    :target: https://pypi.python.org/pypi/apexpy

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/apexpy.svg?style=flat
    :alt: Supported implementations
    :target: https://pypi.python.org/pypi/apexpy

.. |scrutinizer| image:: https://img.shields.io/scrutinizer/g/aburrell/apexpy/main.svg?style=flat
    :alt: Scrutinizer Status
    :target: https://scrutinizer-ci.com/g/aburrell/apexpy/

.. |doi| image:: https://www.zenodo.org/badge/46420037.svg
   :target: https://www.zenodo.org/badge/latestdoi/46420037
