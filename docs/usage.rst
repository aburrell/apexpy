==============
Usage examples
==============

Python library
==============

For full documentation of the functions, see the reference for :class:`apexpy.Apex`.

    >>> from apexpy import Apex
    >>> from __future__ import print_function
    >>> A = Apex(date=2015.3)  # datetime objects are also supported
    >>> # geo to apex, scalar input
    >>> mlat, mlon = A.convert(60, 15, 'geo', 'apex', height=300)
    >>> print("{:.12f}, {:.12f}".format(mlat, mlon))
    57.469573974609, 93.633583068848
    >>> # apex to geo, array input
    >>> glat, glon = A.convert([90, -90], 0, 'apex', 'geo', height=0)
    >>> print(["{:.12f}, {:.12f}".format(ll, glon[i]) for i,ll in enumerate(glat)])
    ['83.099594116211, -84.594589233398', '-74.388267517090, 125.714927673340']
    >>> # geo to MLT
    >>> import datetime as dt
    >>> mlat, mlt = A.convert(60, 15, 'geo', 'mlt', datetime=dt.datetime(2015, 2, 10, 18, 0, 0))
    >>> print("{:.12f}, {:.12f}".format(mlat, mlt))
    56.590423583984, 19.108103879293
    >>> # can also convert magnetic longitude to mlt
    >>> mlt = A.mlon2mlt(120, dt.datetime(2015, 2, 10, 18, 0, 0))
    >>> print("{:.2f}".format(mlt))
    20.89




Command-line interface
======================

.. highlight:: none

The Python package also installs a command called ``apexpy`` which allows using the :py:meth:`~apexpy.Apex.convert` method from the command line. The command-line interface allows you to make use of the Python library even if you don't know or use Python. See the reference for :doc:`Command-line interface <reference/cli>` for a list of arguments to the commands. Below are some simple usage examples.

Produce a file called e.g. ``input.txt`` with the input latitudes and longitudes on each row separated by whitespace::

    # lat lon
    # comment lines like these are ignored
    60 15
    61 15
    62 15

To convert this to apex using the magnetic field model for the date 2015-01-01 using a height of 300 km, run the command ``apexpy geo apex 20150101 --height 300 -i input.txt -o output.txt`` (in this specific example you could also just use ``2015`` or ``201501`` for the date). The output file will look like this::

    57.46954727 93.63981628
    58.52270126 94.04476166
    59.57146454 94.47725677

Alternatively, you can skip the files and just use command-line piping::

    $ echo 60 15 | apexpy geo apex 2015 --height 300
    57.46954727 93.63981628

MLT conversion works in much the same way, but requires both date and time (``YYYYMMDDHHMMSS``). For example, if you want to find the MLT (and the apex latitude) at geodetic coordinates (60, 15) for midnight on the day 2015-01-01, run ``echo 60 15 | apexpy geo mlt 20150101000000``. The output will look like this::

    56.59033585 1.03556417
