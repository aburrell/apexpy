.. _ex-cmi:

Command-Line Interface
======================
The Python package also installs a command called ``apexpy`` which allows using
the :py:meth:`~apexpy.Apex.convert` method from the command line. The
command-line interface allows you to make use of the Python library even if you
don't know or use Python. See the reference for
:ref:`cli` for a list of arguments to the
commands. Below are some simple usage examples.


Running Many Locations
----------------------

You can convert many locations at a single time using an input file.  To follow
this example, create a file called ``input.txt`` with the input latitudes and
longitudes on each row separated by whitespace as shown below.
::

    # gdlat gdlon
    # comment lines like these are ignored
    60 15
    61 15
    62 15

To convert these geodetic coordinates to apex coordinates using the magnetic
field model for the date 2015-01-01 using a height of 300 km, run the command
``apexpy geo apex 20150101 --height 300 -i input.txt -o output.txt``
(in this specific example you could also just use ``2015`` or ``201501`` for
the date). This will create an output file named ``output.txt`` that looks like
this:
::

    57.46954727 93.63981628
    58.52270126 94.04476166
    59.57146454 94.47725677


Running One Location
--------------------
If you don't have a lot of data to process, you can skip the input and output
files using command-line piping.  The ``echo`` command will provide the
input geodetic latitude and longitude in this example and the output appears
on the screen:
::

    $ echo 60 15 | apexpy geo apex 2015 --height 300
    57.46954727 93.63981628


MLT Conversions
---------------
The previous examples all showed how to convert from geodetic to apex
coordinates, but you can convert to and from any of the coordinate systems
supported by :py:meth:`~apexpy.Apex.convert`.  In this example, we show how
to convert from geodetic latitude and longitude to apex latitude and magnetic
local time (MLT).
    
MLT conversion works in much the same way as any other coordinate conversion,
but requires both date and time (``YYYYMMDDHHMMSS``). For example, if you want
to find the MLT (and the apex latitude) at geodetic coordinates (60, 15) for
midnight on the day 2015-01-01, run
``echo 60 15 | apexpy geo mlt 20150101000000``. The output will look like this:
::

    56.59033585 1.03556417
