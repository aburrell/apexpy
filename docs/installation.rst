============
Installation
============

This package requires NumPy, which you can install alone or as a part of SciPy. `Some Python distributions <http://www.scipy.org/install.html#scientific-python-distributions>`_ come with NumPy/SciPy pre-installed. For Python distributions without NumPy/SciPy, Windows/Mac users should install `pre-compiled binaries of NumPy/SciPy <http://www.scipy.org/scipylib/download.html#official-source-and-binary-releases>`_, and Linux users may have NumPy/SciPy available in `their repositories <http://www.scipy.org/scipylib/download.html#third-party-vendor-package-managers>`_.

When you have NumPy, install this package at the command line using ``pip`` [1]_::

    pip install apexpy

On Linux this will build apexpy from source, which requires a fortran compiler such as gfortran.

The package has been tested with the following setups (others might work, too):

* Windows (32/64 bit Python) and Linux (64 bit)
* Python 2.7, 3.3, 3.4 (and 3.5 on Linux [2]_)

.. [1] pip is included with Python 2 from v2.7.9 and Python 3 from v3.4. If you don't have pip, `get it here <http://pip.readthedocs.org/en/stable/installing/>`_.
.. [2] I do not know of any way to compile the Fortran extension on Windows in a manner that is compatible with Python 3.5. If you get it working, let me know!
