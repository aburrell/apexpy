.. _installation:

Installation
============

This package requires NumPy, which you can install alone or as a part of SciPy.
`Some Python distributions <https://www.scipy.org/install.html>`_
come with NumPy/SciPy pre-installed. For Python distributions without
NumPy/SciPy, Windows/Mac users should install
`pre-compiled binaries of NumPy/SciPy <https://www.scipy.org/scipylib/download.html#official-source-and-binary-releases>`_, and Linux users may have
NumPy/SciPy available in
`their repositories <https://www.scipy.org/scipylib/download.html#third-party-vendor-package-managers>`_.

When you have NumPy, you may use either PyPi or GitHub to install this package.


.. _installation-pip:

PyPi
----
This is the most straigforward option!  From the command line use
``pip`` [1]_::

    pip install apexpy

This assumes that the same version of libgfortran is installed in the same
location as when the pip wheel was built (in the event that a wheel is
available). If not, you may have trouble importing apexpy.  In the event that
you run into issues, you can get around this problem by using
``pip`` [1]_::

    pip install --global-option='build_ext' apexpy

which requires both libgfortran and gfortran to be installed on your system.
This is the default option for Linux, and so should not be an issue there.

The package has been tested with the following setups (others might work, too):

* Windows (32/64 bit Python), Linux (64 bit), and Mac (64 bit)
* Python 2.7, 3.6, 3.7, 3.8, 3.9


.. _installation-cmd:

GitHub
------
If you are indending on modifying or contributing to apexpy, it's easier to
install apexpy by forking the repository and installing it locally or within
a virtual environment. After clonining the fork (see :ref:`contributing`),
you may install by::

  cd apexpy
  python setup.py --user


Another benefit of installing apexpy from the command line is specifying the
fortran compiler you would like to use.  By default, apexpy uses
`gfortran <https://gcc.gnu.org/wiki/GFortran>`_, but you can alter the
Makefile in ``src/fortranapex`` to use other compilers or specify different
compilation flags.  However, if using an Intel compiler, you will need to
uncomment a line at the top of ``src/fortranapex/igrf.f90`` to ensure all
necessary libraries are imported.

.. [1] pip is included with Python 2 from v2.7.9 and Python 3 from v3.4.
       If you don't have pip,
       `get it here <https://pip.pypa.io/en/stable/installing/>`_.
