Package Maintenance
===================

Updating IGRF
-------------

The `International Geomagnetic Reference Field <https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html>`_
is regularly updated to reflect the most recent changes to the Terrestrial
magnetic field. apexpy currently uses IRGF-13 coefficients, which are provided
in the ``apexpy/src/fortranapex/igrf13coeff.txt`` file. To change or update the
magnetic field coefficients used by apexpy, you need to update the fortran code.
Assuming your new coefficient file has the same format, the process is simple:

1. Clone the repository or your fork of the repository (see :ref:`contributing`)
2. Update ``apexpy/src/fortranapex/magfld.f`` variable ``FILENAME`` by setting
   it equal to the new target filename
3. Optional: Test the fortran compilation by running the test executable to
   ensure that the new coefficients are used and deviate an acceptable amount
   from prior coefficients::

     make
     apextest
     make clean

4. Install the python package from the command line
   (see :ref:`installation-cmd`)

Updating tests and style standards
-----------------------------------

apexpy is in the process of updating unit and integration tests to reduce code
duplication and implementing cleaner style standards. Additionally, some parts
of the fortran code adhere to older coding standards and raise warnings when
compiled with newer compilers. If you would like to assist in these efforts
(help would be appreciated), please discuss your potential contribution with
the current maintainer to ensure a minimal duplication of effort.
