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
   stime = dt.datetime(2021, 3, 15, 15, 6)
   iss_lat = 9.8
   iss_lon = 142.2
   iss_alt = 419.0

   # Get the apex lat
   apex_iss = apexpy.Apex(stime, refh=iss_alt)
   alat, alon = apex_iss.geo2apex(iss_lat, iss_lon, iss_alt)

   # Get the apex height
   aalt = apex_iss.get_apex(alat, iss_alt)

   # Convert from apex height in km to L-shell
   L_iss = 1.0 + aalt / apex_iss.RE
   print("apex height={:.3f} km, L={:.2f}".format(aalt, L_iss))
   apex height=428.007 km, L=1.07


Trace a Field Line
==================

It can be useful to trace a field line with a specified apex height across a
range of known latitudes.  This can then be useful when plotting field lines at
a particular magnetic meridian or used to gather data along the same field line,
but measured by different instruments.

::

   import numpy as np

   # Continue form the previous example, define the field line at all latitudes
   lats = np.linspace(-90, 90, 181)
   alts = apex_iss.get_height(lats, aalt)

   # Select the locations with positive altitudes (above the Earth's surface)
   iline = np.where(alts >= 0.0)[0]

   # Print the latitude limits for which the field line is above the surface
   # of the Earth
   print("lat={:.2f} deg, alt={:.2f} km; lat={:.2f} deg, alt={:.2f} km".format(
       lats[iline[0]], alts[iline[0]], lats[iline[-1]], alts[iline[-1]]))
   lat=-14.00 deg, alt=30.09 km; lat=14.00 deg, alt=30.09 km
