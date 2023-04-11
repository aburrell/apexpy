
Changelog
=========

2.0.1 (2023-04-11)
------------------
* Expanded installation instructions in the documenation
* Added unit tests for todays date, ensuring that `apex.dat` is current
* Added cron unit test to GitHub Actions CI
* Added a logo
* Correct indexing bug in Fortran source that was causing array overflow and
  memory errors for extrapolated years beyond the latest formal IGRF fit

2.0.0 (2022-12-09)
------------------
* Update Fortran source code to Fortran 90 standards
* Removed Python 2 support
* Updated community and package documentation
* Updated unit test style to reduce duplication and better follow PEP8
* Updated code style using codacy suggestions and reduced complexity
* Added class representation strings to Apex
* Improved input testing for Apex methods
* Added more examples and installation help to the documentation
* Fixed missing microseconds bug in helpers.subsol
* Added function to calculate height along a field line
* Changed installation to use meson
* Added wheel creation to CI
* Updated flake8 ignore syntax

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
