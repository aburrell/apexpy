.. _ex-apexh:

Calculate L-Shells
==================

L-shells are the apex heights of specified magnetic field lines in units of
Earth radii where L=0 corresponds to the center of the Earth.  You can calculate
the L-shells seen by a given instrument using the
:py:meth:`~apexpy.Apex.geo2apex` and :py:meth:`~apexpy.Apex.get_apex` methods.
An example of this is shown for a single point along the orbit of the
International Space Station (ISS) below.
::

   import apexpy
   import datetime as dt

   # Set the location of the ISS in geodetic coordinates at a single time
   stime = stime = dt.datetime(2021, 3, 15, 15, 6)
   iss_lat = 9.8
   iss_lon = 142.2
   iss_alt = 419.0

   # Get the apex lat
   alat, alon = apex_iss.geo2apex(iss_lat, iss_lon, iss_alt)

   # Get the apex height
   aalt = apex_iss.get_apex(alat, iss_alt)

   # Convert from apex height in km to L-shell
   L_iss = 1.0 + aalt / apex_iss.RE
   print("apex height={:.3f} km, L={:.2f}".format(aalt, L_iss))
   apex height=875.155 km, L=1.14
