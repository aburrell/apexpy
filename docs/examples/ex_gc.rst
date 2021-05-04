.. _ex-gc:

Convert from various non-Magnetic Coordinate Systems
====================================================

Usually non-magnetic coordinates for space sciences are provided in geodetic
(geographic) coordinates, where the latitude is the angle between the
equatorial plane and the normal to the ellipsoid surface at the desired
location. There are several different reference ellipsoids that may be used,
but the most common (and the one used by :py:class:`~apexpy.Apex` is WGS84, the
World Geodetic System 1984 (as described in, for example, Snay and Soler [1]_).
However, it frequently makes sense for a particular instrument to use a
different reference ellipsoid or another type of coordinate system.

Different types of geodetic or geocentric coordinate conversions are not
supported within apexpy, because they are not a magnetic coordinate
transformations.  We recommend first converting your other non-magnetic
coordinates to geodetic WGS84 coordinates and then performing the desired
conversion to apex or quasi-dipole coordinates.

Geocentric Example
------------------

One commonly encountered alternative to geodetic latitude is geocentric
latitude.  Geocentric latitude is the angle between the equatorial plane and a
line joining the center of the Earth and the desired location. If you have data
in geocentric coordinates you can convert them to geodetic using a simple
equation layed out in equation 3-28 of Snyder [2]_.
::

   import numpy as np

   def gc_to_gd(gc_lat, e_sq):
       """Convert from geocentric to geodetic

       Parameters
       ----------
       gc_lat : array-like
          Geocentric latitude in degrees
       e_sq : float
           First eccentricity squared, unitless

       Returns
       -------
       gd_lat : array-like
           Geodetic latitude in degrees
       """
       gd_lat = np.arctan(np.tan(np.radians(gc_lat)) / (1.0 - e_sq))
       return np.degrees(gd_lat)

       
The function above requires the first eccentricity of the reference ellipsoid.
This example uses the ``pyproj`` library [3]_ to get the WGS84 ellipsoid data,
but the function shown will take any float.  This lets you decide the level of
precision you need in your coordinate transformation.
::

  import pyproj

   # Initialize the WGS84 reference ellipsoid
   wgs84_geod = pyproj.crs.GeographicCRS(name='WGS84').get_geod()
   print("{:.12f}".format(wgs84_geod.es))
   0.006694379990


Now, try converting from geocentric to quasi-dipole coordinates.  In this
example you need to supply height, because the coordinate transformation from
geodetic to quasi-dipole takes place by converting from geodetic to apex and
then from apex to quasi-dipole.
::

   import apexpy
   
   # Define the starting values
   year = 2015.3
   gc_lat = 45.0
   glon = 0.0
   height = 0.0

   # Get the quasi-dipole coordiantes
   apex_out = apexpy.Apex(date=year)
   qlat, qlon = apex_out.geo2qd(gc_to_gd(gc_lat, wgs84_geod.es), glon, height)
   print("{:.12f}, {:.12f}".format(qlat, qlon))
   39.852313995361, 76.711242675781


.. [1] Snay and Soler (1999) Modern Terrestrial Reference Systems (Part 1),
       `Professional Surveyor <https://www.ngs.noaa.gov/CORS/Articles/
       Reference-Systems-Part-1.pdf>`_.
.. [2] Snyder, J. P. Map projections — A working manual. Professional Paper
       1395, U.S. Geological Survey, 1987.
       `doi:10.3133/pp1395 <https://pubs.er.usgs.gov/publication/pp1395>`_.
.. [3] `pyproj GitHub page <https://github.com/pyproj4/pyproj>`_.      
