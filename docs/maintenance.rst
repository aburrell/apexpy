Package Maintenance
===================

Providing Wheels with Releases
------------------------------

The Continuous Integration (CI) now saves wheels created for each tested Python
version and computer Operating System (OS) as artifacts. When preparing a new
PyPi release, these wheels may be downloaded from the release candidate.  We
currently don't include them, because the wheels only work when the installation
environment mirrors the CI environment.

Updating IGRF
-------------

The `International Geomagnetic Reference Field <https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html>`_
is regularly updated to reflect the most recent changes to the Terrestrial
magnetic field. apexpy currently uses IRGF-14 coefficients, which are provided
in the ``apexpy/apexpy/igrf14coeff.txt`` file. To change or update the
magnetic field coefficients used by apexpy, you need to update the python code,
then rerun the fortran program that builds ``apexpy/apexpy/apexsh.dat``. This 
is what makes apexpy performant. For more details, see Emmert et al. [2010] [1]_.

Assuming your new coefficient file has the same format, the process is simple:

1. Clone the repository or your fork of the repository (see :ref:`contributing`).
2. Update ``apexpy/apexpy/apex.py`` variable ``igrf_fn`` by setting
   it equal to the new IGRF coefficient filename (``igrf14coeff.txt``, for
   example).
3. In ``apexpy/fortranapex/checkapexsh.f90``, update the variable ``igrffilein``
   to the new IGRF coefficent filename.  Relative paths are allowed.
4. Modify ``checkapexsh.f90`` by adding the next 5 year epoch to the
   ``epochgrid`` variable and updating the ``nepochgrid`` variable as
   necessary. For example, if the newest IGRF coefficients are good up to 2030
   and ``epochgrid`` only has up to the year 2025, then add 2030 to
   ``epochgrid`` and then increment ``nepochgrid`` by 1.
5. Execute the ``apextest`` binary to generate the new ``apexsh.dat`` file.
6. Update the unit tests in the class ``TestApexMethodExtrapolateIGRF`` in 
   ``apexpy/apexpy/tests/test_Apex.py`` so that they check the methods are
   working correctly with dates after the latest IGRF epoch (i.e., if the
   latest epoch is 2025, set the test to initialize with the year 2030).  You
   will have to update the hard-coded confirmation values used by these tests.
7. Commit all changes and create a pull request on GitHub to integrate your 
   branch with updated IGRF into the main repository.

Modifying Fortran Source
------------------------
When modifying the fortran source code, it can be helpful to run a preliminary
validation of the fortran output independent of the python wrapper. This should
be done within the ``apexpy/fortranapex`` directory.

1. Remove any existing binaries by running the ``make clean`` command.
2. Build the ``apextest`` binary by running the ``make`` command.
3. Execute the ``apextext`` binary.
4. Confirm the output printed to the screen matches the test output shown in
   the comment block at the bottom of ``checkapexsh.f90``. The output may not
   match the test output exactly due to floating point errors and improvements
   in the precision of the calculation.
5. If the modifications involved adding or removing fortran source files, modify
   the list of extension sources in ``setup.cfg``.
6. Rebuild and install apexpy following the instructions in
   :ref:`installation-build`.

Updating tests and style standards
-----------------------------------

apexpy is in the process of updating unit and integration tests to reduce code
duplication and implementing cleaner style standards. Additionally, some parts
of the fortran code adhere to older coding standards and raise warnings when
compiled with newer compilers. If you would like to assist in these efforts
(help would be appreciated), please discuss your potential contribution with
the current maintainer to ensure a minimal duplication of effort.


.. [1] Emmert, J. T., A. D. Richmond, and D. P. Drob (2010),
       A computationally compact representation of Magnetic-Apex
       and Quasi-Dipole coordinates with smooth base vectors,
       J. Geophys. Res., 115(A8), A08322, doi:10.1029/2010JA015326.
