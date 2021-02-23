
Changelog
=========

1.1.0 (2021-03-05)
------------------
* Adapted Fortran to read IRGF coefficients from a file (updated to IGRF-13)
* Improved the subsol routine to allow array input
* Improved PEP8 compliance
* Added some missing docstrings to unit tests
* Fixed AppVeyor test environment
* Updated python test versions
* Updated community and package documentation
* Fixed bug where NaNs caused array input to crash
* Fixed bug in quasi-dipole to apex conversion at equator
* Removed duplicate CI services

1.0.4 (2019-04-05)
----------------------------------------
* Updated installation instructions
* Simplified warning tests
* Made some PEP8 changes

1.0.3 (2018-04-05)
-----------------------------------------
* Updated badges and added DOI
* Added tests for python 3.6
* Removed tests for python 3.3
* Made some PEP8 changes

1.0.2 (2018-02-27)
-----------------------------------------

* Extend character limit for allowable data file path, and update documentation
  to reflect a change in maintainers.  Also updated testing implimentation,
  reduced fortran compiler warnings, and improved PEP8 compliance.

1.0.1 (2016-03-10)
-----------------------------------------

* Remove geocentric to geodetic conversion of subsolar point based on feedback
  from Art Richmond. (The subsolar point is the same in geocentric and geodetic
  coordinates.) The helper function `gc2gdlat` have been kept to preserve
  backwards compatibility.


1.0.0 (2015-11-30)
-----------------------------------------

* Initial release
