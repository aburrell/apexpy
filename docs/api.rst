Reference
=========

The :py:class:`apexpy.Apex` class is used for all the main functionality
(converting between coordinate systems, field line mapping, and calculating
base vectors). The :py:mod:`apexpy.helpers` sub-module includes additional
functions that may be useful, especially :py:func:`~apexpy.helpers.subsol`.
The :py:mod:`apexpy.fortranapex` module is the interface to the apex Fortran
library by Emmert et al. [2010] [1]_. The interface is not documented.
Use :class:`apexpy.Apex` for all conversions and calculations. You can find
some documentation of the actual Fortran library in the source file
`apexsh.f90 <https://github.com/aburrell/apexpy/blob/main/fortranapex/apexsh.f90>`_.  These functions may also be accessed through the command-line interface.

.. _api:

API
---

.. toctree::
   :maxdepth: 2

   autoapi/generated/apexpy/index.rst


.. _cli:

Command-line interface
----------------------

.. highlight:: none

When you install this package you will get a command called ``apexpy``, which
is an interface to the :py:meth:`~apexpy.Apex.convert` method. See the
documentation for this method for a more thorough explanation of arguments and
behaviour.

You can get help on the command by running ``apexpy -h``.

.. code::

    $ apexpy -h
    usage: apexpy [-h] [--height HEIGHT] [--refh REFH] [-i FILE_IN]
                  [-o FILE_OUT] SOURCE DEST DATE

    Converts between geodetic, modified apex, quasi-dipole and MLT

    positional arguments:
      SOURCE                Convert from {geo, apex, qd, mlt}
      DEST                  Convert to {geo, apex, qd, mlt}
      DATE                  YYYY[MM[DD[HHMMSS]]] date/time for IGRF
                            coefficients, time part required for MLT
			    calculations

    optional arguments:
      -h, --help            show this help message and exit
      --height HEIGHT       height for conversion
      --refh REFH           reference height for modified apex coordinates
      -i FILE_IN, --input FILE_IN
                            input file (stdin if none specified)
      -o FILE_OUT, --output FILE_OUT
                            output file (stdout if none specified)

.. [1] Emmert, J. T., A. D. Richmond, and D. P. Drob (2010),
       A computationally compact representation of Magnetic-Apex
       and Quasi-Dipole coordinates with smooth base vectors,
       J. Geophys. Res., 115(A8), A08322, doi:10.1029/2010JA015326.
