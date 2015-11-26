Command-line interface
======================

.. highlight:: none

When you install this package you will get a command called ``apexpy``, which is an interface to the :meth:`~apexpy.Apex.convert` method. See the documentation for this method for a more thorough explanation of arguments and behaviour.

You can get help on the command by running ``apexpy -h``.

.. code::

    $ apexpy -h
    usage: apexpy [-h] [--height HEIGHT] [--refh REFH] [-i FILE_IN] [-o FILE_OUT]
                  SOURCE DEST DATE

    Converts between geodetic, modified apex, quasi-dipole and MLT

    positional arguments:
      SOURCE                Convert from {geo, apex, qd, mlt}
      DEST                  Convert to {geo, apex, qd, mlt}
      DATE                  YYYY[MM[DD[HHMMSS]]] date/time for IGRF coefficients,
                            time part required for MLT calculations

    optional arguments:
      -h, --help            show this help message and exit
      --height HEIGHT       height for conversion
      --refh REFH           reference height for modified apex coordinates
      -i FILE_IN, --input FILE_IN
                            input file (stdin if none specified)
      -o FILE_OUT, --output FILE_OUT
                            output file (stdout if none specified)
