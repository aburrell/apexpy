==============
Usage examples
==============

Python library
==============

For full documentation of the functions, see the reference for :class:`apexpy.Apex`.

    >>> from apexpy import Apex
    >>> A = Apex(date=2015.3)  # datetime objects are also supported
    >>> # geo to apex, scalar input
    >>> mlat, mlon = A.convert(60, 15, 'geo', 'apex', height=300)
    >>> mlat
    57.469573974609375
    >>> mlon
    93.633583068847656
    >>> # apex to geo, array input
    >>> glat, glon = A.convert([90, -90], 0, 'apex', 'geo', height=0)
    >>> glat
    array([ 83.09959412, -74.38826752])
    >>> glon
    array([ -84.59458923,  125.71492767])
    >>> # geo to MLT
    >>> import datetime as dt
    >>> mlat, mlt = A.convert(60, 15, 'geo', 'mlt', datetime=dt.datetime(2015, 2, 10, 18, 0, 0))
    >>> mlat
    56.590423583984375
    >>> mlt
    19.107855542500815
    >>> # can also convert magnetic longitude to mlt
    >>> mlt = A.mlon2mlt(120, dt.datetime(2015, 2, 10, 18, 0, 0))
    >>> mlt
    20.893299166361491




Command-line interface
======================

.. highlight:: none

The Python package also installs a command called ``apexpy`` which allows using the :meth:`apexpy.Apex.convert <Apex.convert>` method from the command line. The command-line interface allows you to make use of the Python library even if you don't know or use Python. See the reference for :doc:`Command-line interface <reference/cli>` for a list of arguments to the commands. Below are some simple usage examples.

Produce a file called e.g. ``input.txt`` with the input latitudes, longitudes and altitudes on each row separated by whitespace::

    # lat lon
    # comment lines like these are ignored
    60 15
    61 15
    62 15

To convert this to apex using the magnetic field model for the date 2015-02-24 using a height of 300 km, run the command ``apexpy geo apex -i input.txt -o output.txt -d 20150224 -h 300``. The output file will look like this::

    57.47612194 93.55719875
    58.53323704 93.96069212
    59.58522105 94.38968625

Alternatively, you can skip the files and just use command-line piping::

    $ echo 60 15 | apexpy -d 20150224 -h 300
    57.47612194 93.55719875


MLT conversion works in much the same way, but requires both date and time in the ``-d`` argument. If the columns in the input file shown above are apex latitude, longitude and height and you want to convert to MLT, run e.g. ``apexpy apex mlt -d 20150224140015 -i input.txt -o output.txt`` (note that the date/time is a required parameter). The output file will then look like this::

    60 13.7
    61 13.7
    62 13.9
