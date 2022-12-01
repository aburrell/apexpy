.. _installation:

Installation
============

This package requires NumPy, which you can install alone or as a part of SciPy.
`Some Python distributions <https://scipy.org/install/>`_
come with NumPy/SciPy pre-installed. For Python distributions without
NumPy/SciPy, Windows/Mac users should install
`pre-compiled binaries of NumPy/SciPy <https://scipy.org/download/#official-source-and-binary-releases>`_, and Linux users may have
NumPy/SciPy available in `their repositories <https://scipy.org/download/>`_.

When you have NumPy, you may use either PyPI or GitHub to install this package.

The code behind this package is written in Fortran.  Because of this, you
**MUST** have a fortran compiler installed on your system before you attempt
the next step.  `Gfortran <https://gcc.gnu.org/wiki/GFortran>`_ is a free
compiler that can be installed, if one is not already available on your system.
If you are installling this or MinGW in Windows, make sure you install it
**after** installing the Windows Microsoft C++ Build tools.  The apexpy
installation has been tested successfully with gfortran 7 and some more recent
versions.  Earlier versions of gfortran may not work and are not recommended.


.. _installation-tested:

Tested environments
-------------------

The package has been tested with the following setups (others might work, too):

* Windows (64 bit Python), Linux (64 bit), and Mac (64 bit)
* Python 3.7, 3.8, 3.9, 3.10


.. _installation-pip:

PyPI
----
This is the most straightforward option!  From the command line use
``pip`` [1]_::

    pip install apexpy

Installation Issues
^^^^^^^^^^^^^^^^^^^

Because the base code is in Fortran, installation can be tricky and different
problems can arise even if you already have a compiler installed (but please
check that you do before trying these solutions).

ApexPy is not compatible with numpy version 1.19.X, older and newer versions
of numpy work without issue.

Many times, the following modification to ``pip`` [1]_ will solve the
installation problem, but it requires that both libgfortran and gfortran are
installed on your system.::

    pip install --no-binary :apexpy: apexpy

This is the default option for Linux, and so should not be an issue there. On
Windows with the Mingw32 compiler, you might find `this information <https://wiki.python.org/moin/WindowsCompilers#GCC_-_MinGW-w64_.28x86.2C_x64.29>`_
useful for helping build apexpy.

Problems have also been encountered when installing in a conda environment.
In this case, pip seems to ignore the installed numpy version when installing.
This appears to result in a successful installation that fails upon import.  In
this case, try::

  pip install apexpy --no-build-isolation

The specific fortran compiler to use can be specified with the `FC` flag.  This
can be useful your system has multiple versions of gfortran and the default
is not appropriate (ie., an older version).

::
   FC=/path/to/correct/gfortran pip install apexpy

Apple Silicon systems require certain compliation flags to deal with memory
problems.  In this case, the following command has worked (after first
installing both gfortran and numpy).

::
   CFLAGS="-falign-functions=8 ${CFLAGS}" pip install --no-binary :apexpy: apexpy

If you are on Apple and encounter a library error such as
``ld: library not found for -lm``, you will need to provide an additional
linking flag to the Mac OSX SDK library.  This example assumes you are building
locally from the cloned Git repository.  Issues on Mac OS have also been
encountered when using clang for ``CC`` alongside gfortran.  This resulted in a
seemly successful installation with apexpy reporting that fortranapex cannot be
imported.

::
   LDFLAGS="-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib ${LDFLAGS}" pip install apexpy

Windows systems are known to have issues with Fortran-based codes.  The Windows
testing we do uses miniconda, so we recommend using the Anaconda environment.
One problem that has been encountered is a lack of LAPACK/BLAS tools that
causes numpy to not behave as expected.  This can be fixed by installing
scipy before numpy and then installing apexpy.


.. _installation-cmd:

GitHub
------
If you are indending on modifying or contributing to apexpy, it's easier to
install apexpy by forking the repository and installing it locally or within
a virtual environment. After clonining the fork (see :ref:`contributing`),
you may install by::

  cd apexpy
  python -m build .
  pip install -e .


The ``-e`` flag uses pip to peform what used to be ``python setup.py develop``.

Another benefit of installing apexpy from the command line is specifying the
compilers you would like to use.  These can be changed by altering the ``CC``
and ``FC`` environment variables on your computer.  If using an Intel compiler,
you will need to uncomment a line at the top of ``src/fortranapex/igrf.f90`` to
ensure all necessary libraries are imported.

If the above command doesn't work for you (as may be the case for Windows), you
can try::

  cd apexpy
  meson setup build
  ninja -j 2 -C build
  cd build
  meson install

.. [1] pip is included with Python 2 from v2.7.9 and Python 3 from v3.4.
       If you don't have pip,
       `get it here <https://pip.pypa.io/en/stable/installing/>`_.
