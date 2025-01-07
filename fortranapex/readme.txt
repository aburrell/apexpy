Auxiliary material for Paper 2010JA015326

A computationally compact representation of Magnetic-Apex and Quasi-Dipole
coordinates with smooth base vectors

J. T. Emmert and D. P. Drob
Space Science Division, U.S. Naval Research Laboratory, Washington, D.C., USA

A. D. Richmond
High Altitude Observatory, National Center for Atmospheric Research, Boulder,
Colorado, USA


Emmert, J. T., A. D. Richmond, and D. P. Drob (2010), A computationally compact
representation of magnetic Apex and Quasi-Dipole coordinates with smooth base
vectors, J. Geophys. Res., 115, A08322, doi:10.1029/2010JA015326.

Introduction
------------

Auxiliary material for this article contains Fortran code for converting from
geodetic to Quasi-Dipole (QD) and Modified Apex coordinates, for computing base
vectors in these coordinates, for converting from QD to geodetic coordinates,
and for generating expansion coefficients for these conversions, as described
in the paper. Updated versions of this code will be made available on the  CEDAR
Data System at http://cedarweb.hao.ucar.edu/.

1. 2010ja015326-txts01.txt 
checkapexsh.f90: Fortran-90 test driver and usage example for the QD 
coordinate conversion code package.

2. 2010ja015326-txts02.txt 
apexsh.f90: Fortran-90 code for converting among coordinate systems
and computing base vectors. This file contains all the routines necessary to do
the conversions. The last three files are only needed to generate a file
containing  the expansion coefficients.

3. 2010ja015326-txt03.txt 
makeapexsh.f90: Fortran-90 code for generating the expansion
coefficients used by apexsh.f90. Must be compiled together with apexsh.f90,
apex.f, and magfld.f.

4. 2010ja015326-txts04.txt 
apex.f: Fortran-77 code for computing the apex height of magnetic
field lines.

5. 2010ja015326-txts05.txt 
magfld.f: Fortran-77 code for evaluating IGRF.

##############################################

DEV NOTES:
----------

The above text is copy-pasted from the auxiliary materials for the Emmert et al.
(2010) paper.

In February 2021, the magfld.f source was modified so that it could read IGRF
coefficients from the text file instead of from hard-coded constants.
Specifically the COFRM subroutine was modified. This also required modifying
APEX subroutine in apex.f, the makeapxsh subroutine in makeapexsh.f90, and the
checkapexsh program in checkapexsh.f90.


MAINTENANCE NOTE:
-----------------

After updating the IGRF coefficients file, we need to rebuild the apexsh.dat
file. This is done by:

1) Adding the next 5 year epoch to "epochgrid" and updating the "nepochgrid"
   variable in checkapexsh.f90. For example, if epochgrid has up to the year
   2020 and the newest IGRF coefficients are good up to 2025, we should add
   2025 to epochgrid and then increment nepochgrid by 1.
2) checkapexsh.f90 also expects the IGRF coefficient file to be in the same
   directory with the name "igrf14coeffs.txt", so you may also need to update
   the "igrffilein" variable as well.
3) Building the "apextest" binary using the "make" command.
4) Copying the resulting "apexsh.dat" file to the apexpy/src/apexpy directory.

If the fortran code is updated, the .pyf file will also need to be recreated.
This can be done by:

$ f2py -m fortranapex apexsh.f90 igrf.f90 apex.f90 magfld.f90 makeapexsh.f90 -h fortranapex.pyf

Be sure you remove the existing `fortranapex.pyf` file before performing this
step.
