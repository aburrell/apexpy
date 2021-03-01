Package Maintenance
===================

Updating IGRF
-------------

The `International Geomagnetic Reference Field <https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html>`_
is regularly updated to reflect the most recent changes to the Terrestrial
magnetic field. apexpy currently uses IRGF-13 coefficients, which are provided
in the ``apexpy/src/apexpy/igrf13coeff.txt`` file. To change or update the
magnetic field coefficients used by apexpy, you need to update the python code.
Assuming your new coefficient file has the same format, the process is simple:

1. Clone the repository or your fork of the repository (see :ref:`contributing`)
2. Update ``apexpy/src/apexpy/apex.py`` variable ``igrf_fn`` by setting
   it equal to the new target filename
3. Install the python package from the command line
   (see :ref:`installation-cmd`)

Updating ``apexsh.dat``
-----------------------

After updating the IGRF coefficients file the ``apexsh.dat`` file needs to be
rebuilt. This file is what makes apexpy performant. For more details, see
Emmert et al. [2010] [1]_.

Updating ``apexsh.dat`` is done by modifying and building the fortran source code
in the ``apexpy/src/fortranapex`` directory. Working in that directory:

1. Make sure a copy of the latest IGRF coefficient file is present in the
   selfsame directory.
2. Modify the ``igrffilein`` in ``checkapexsh.f90`` to the name of the IGRF
   coefficient file (``igrf13coeff.txt``, for example).
3. Modify ``checkapexsh.f90`` by adding the next 5 year epoch to the 
   ``epochgrid`` variable and updating the ``nepochgrid`` variable as
   necessary. For example, if the newest IGRF coefficients are good up to 2025
   and ``epochgrid`` only has up to the year 2020, then add 2025 to
   ``epochgrid`` and then increment ``nepochgrid`` by 1.
4. Build the ``apextest`` binary by running the ``make`` command.
5. Execute the ``apextest`` binary to generate the new ``apexsh.dat`` file
4. Copy the new ``apexsh.dat`` file to the ``apexpy/src/apexpy`` directory.

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
       J. Geophys. Res., 115(A8), A08322, :doi:`10.1029/2010JA015326`.
