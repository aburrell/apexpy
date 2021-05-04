.. _ex-convert:

Using apexpy
============

The :py:class:`~apexpy.Apex` class contains most of the methods you'll want to
use when converting between geodetic and apex or quasi-dipole coordinates.
The :py:meth:`~apexpy.Apex.convert` method is designed to be the primary user
interface for coordinate conversion.  For full documentation of this and other
class methods, see the reference for :py:class:`~apexpy.Apex`.  Some simple
examples to get you going follow below.


Initialize the Apex class
-------------------------

The Apex class requires a date and time as well as a reference height when
initialized.  If you don't supply one, then the class will default to the
current time (in Universal Time) and a reference height of 0 km.  For this
example, we will specify the time.  Time can be supplied as either a decimal
year, :py:class:`~datetime.datetime` object, or a :py:class:`~datetime.date`
object.
::

  import apexpy
  apex_out = apexpy.Apex(date=2015.3)
  print(apex_out)


This yields:
::
  
  Apex class conversions performed with
  -------------------------------------
  Decimal year: 2015.30000000
  Reference height: 0.000 km
  Earth radius: 6371.009 km

  Coefficient file: '/path/to/programs/apexpy/src/apexpy/apexsh.dat'
  Cython Fortran library: '/path/to/programs/Git/apexpy/src/apexpy/fortranapex.cpython-37m-darwin.so'


Convert from Geodetic to Magnetic Coordinates 
----------------------------------------------

You can use the initialized :py:meth:`~apexpy.Apex` object to convert from
geodetic coordinates to magnetic coordinates and back again.  When converting to
and from apex coordinates, you need to supply the height of the observations.
::

  alat, alon = apex_out.convert(60, 15, 'geo', 'apex', height=300)
  print("{:.12f}, {:.12f}".format(alat, alon))
  57.477310180664, 93.590156555176

For quasi-dipole coordinates, this isn't necessary.
::
  
  qlat, qlon = apex_out.convert(60, 15, 'geo', 'qd')
  print("{:.12f}, {:.12f}".format(qlat, qlon))
  56.598316192627, 93.174751281738

You can calculate multiple locaitons at once using arrays, as long as the
inputs are broadcastable.  For example, you can provide a list or array of
different latitudes for a single longitude (and height).  It is also acceptable
to provide a list or array of the same shape that provides paired latitude,
longitude, and heights (if needed).  However, you can't provide mismatched array
or list inputs.  Here is an example where we convert from apex coordinates to
geodetic coordinates for two different latitudes at the same longitude and
height.
::

  glat, glon = apex_out.convert([90, -90], 0, 'apex', 'geo', height=0)
  print(["{:.12f}, {:.12f}".format(ll, glon[i]) for i,ll in enumerate(glat)])
  ['83.103820800781, -84.526657104492', '-74.388252258301, 125.736274719238']


Convert to Magnetic Local Time
------------------------------

When converting to magnetic local time (MLT), the convert function requires
a datetime input alongside a latitude and longitude.
::

  import datetime as dt
  utime = dt.datetime(2015, 2, 10, 18, 0, 0)
  mlat, mlt = apex_out.convert(60, 15, 'geo', 'mlt', datetime=utime)
  print("{:.12f}, {:.12f}".format(mlat, mlt))
  56.598316192627, 19.107861709595


If you already have magnetic longitude, you can also calculate MLT using
:py:meth:`~apexpy.Apex.mlon2mlt`.
::

  mlt = apex_out.mlon2mlt(120, utime)
  print("{:.2f}".format(mlt))
  20.90
