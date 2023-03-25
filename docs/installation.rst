.. _installation:

Installation
============

.. _installation-cmd:

Standard Installation
---------------------

The recommended (and most straightforward) method for users to install
:py:mod:`apexpy` is through PyPI. From the command line use ``pip`` [1]_::

    pip install apexpy

You should be able to import :py:mod:`apexpy` and run basic conversions as shown in the
examples.  If you get errors or warnings upon importing, see below for more
advanced options and troubleshooting.


.. _installation-tested:

Tested environments
-------------------

The package has been tested with the following setups (others might work, too):

* Windows (64 bit Python), Linux (64 bit), and Mac (64 bit)
* Python 3.7, 3.8, 3.9, 3.10


Advanced Installation
---------------------

If you cannot install :py:mod:`apexpy` from the distribtuion available on PyPI,
you will have to use one of the following more advanced options. These are
generally only recommended if you are planning on developing or modifying the
:py:mod:`apexpy` source code.

The code behind this package is written in Fortran.  Because of this, you
**MUST** have a fortran compiler installed on your system before you attempt
the following steps.  `Gfortran <https://gcc.gnu.org/wiki/GFortran>`_ is a free
compiler that can be installed, if one is not already available on your system.
If you are installing this or MinGW in Windows, make sure you install it
**after** installing the Windows Microsoft C++ Build tools. You must also make
sure that the compilers and Python that are installed both use the same
processing standard (either 32-bit or 64-bit). The :py:mod:`apexpy`
installation has been tested successfully with gfortran 7 and some more recent
versions.  Earlier versions of gfortran may not work and are not recommended.

Installation also requires a C compiler of the same type as the fortran
compiler. `GCC <https://gcc.gnu.org/>`_ is a free compiler that works with
Gfortran and can be installed from a variety of sources and standard package
managers. It is recommended that you check to see if you have gcc available on
your system before installing as it is relatively common and multiple competing
versions may cause problems if paths are not managed carefully.

This package requires NumPy, which you can install alone or as a part of SciPy.
`Some Python distributions <https://scipy.org/install/>`_
come with NumPy/SciPy pre-installed. For Python distributions without
NumPy/SciPy, various package managers for different operating systems allow
for simple local installation (as directed on the SciPy installation page.
Pip should install NumPy automatically when installing :py:mod:`apexpy`, but if
not, install it manually before attempting to install :py:mod:`apexpy`.
:py:mod:`apexpy` is not compatible with NumPy version 1.19.X, older and newer
versions of NumPy work without issue.

**IMPORTANT:** If you are struggling with installing :py:mod:`apexpy` and trying
some of the following options, it is recommended that you user the
``--no-cache`` flag with pip to avoid repeatedly reinstalling the same
non-functional build.


Install from GitHub
^^^^^^^^^^^^^^^^^^^

:py:mod:`apexpy` can be installed from the source code on GitHub, so long as a
fortran compiler is available::

  pip install git+https://github.com/aburrell/apexpy.git

This is advantageous if you would like to install from a particular branch or
tag instead of the latest published stable release on PyPI.  Do this by
appending `@target-branch` to the end of the above command.  For instance, if
you would like to install from the develop branch instead of main, the
appropriate command would be::

  pip install git+https://github.com/aburrell/apexpy.git@develop


Install without Wheels
^^^^^^^^^^^^^^^^^^^^^^

Many times, skipping building wheels locally will solve installation problems,
but it requires that both libgfortran and gfortran are installed on your
system::

    pip install --no-binary :apexpy: apexpy

This is the default option for Linux, and so should not be an issue there. On
Windows with the Mingw32 compiler, you might find `this information <https://wiki.python.org/moin/WindowsCompilers#GCC_-_MinGW-w64_.28x86.2C_x64.29>`_
useful for helping build :py:mod:`apexpy`.

Install against an incompatible numpy version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  pip install apexpy --no-build-isolation --no-cache


Installation using CI Wheels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If your local set up is essentially identical to one of the CI test
environments, then you can use one of the wheel artifacts to install
:py:mod:`apexpy`. The list of artifacts may be found
`here <https://api.github.com/repos/aburrell/apexpy/actions/artifacts>`_.

To download an artifact:

1. If you don't have a GitHub Personal Access Token, follow
   `these instructions <https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token>`_
   to create one.
2. Run ``curl -v -H "Authorization: token <GITHUB-ACCESS-TOKEN>" https://api.github.com/repos/aburrell/apexpy/actions/artifacts/<ARTIFACT-ID>/zip``, where
   <ITEM> should be replaced with the appropriate item string.
3. Copy the URL from the ``Location`` output produced by the previous command
   into a browser, which will download a zip archive into your standard
   download location. Alternatively (or if this doesn't work) you can use
   `wget` to retrieve the archive.
4. Copy the zip archive into the ``apexpy/dist`` directory and unzip.
5. Check the archive for the expected matrix of ``*.whl`` objects

To install, use ``pip install .``

.. _installation-build:

Build from Source
^^^^^^^^^^^^^^^^^

If you intend to modify or contribute to :py:mod:`apexpy`, you should install
:py:mod:`apexpy` by forking the repository and installing it locally or within
a virtual environment. After cloning the fork (see :ref:`contributing`),
you may install by::

  cd apexpy
  python -m build .
  pip install .


Note that the ``-e`` flag for pip, which performs what used to be
``python setup.py develop``, isn't used here.  That's because meson currently
doesn't support develop style builds.

If the above command doesn't work for you (as may be the case for Windows), you
can try::

  cd apexpy
  meson setup build
  ninja -j 2 -C build
  cd build
  meson install


Specifying Compilers
^^^^^^^^^^^^^^^^^^^^

When you install :py:mod:`apexpy` from the command line you can specify the
compilers you would like to use.  These can be changed by altering the ``CC``
and ``FC`` environment variables on your computer::

  FC=/path/to/correct/gfortran CC=/path/to/correct/gcc python -m build
  pip install .

This can be useful your system has multiple versions of gfortran or gcc and the
default is not appropriate (ie., an older version). If using an Intel compiler,
you will need to clone the repository locally and uncomment a line at the top of
``src/fortranapex/igrf.f90`` to ensure all necessary libraries are imported.


When All Else Fails
^^^^^^^^^^^^^^^^^^^

Because the base code is in Fortran, installation can be tricky and different
problems can arise even if you already have a compiler installed.  The following
are a series of installation commands that users have reported working for
different system configurations.  We have not been able to reproduce some of
the issues users report and cannot fully explain why some of the options work,
none the less they are recorded here as they may be useful to other users.  If
you feel like you can provide more insight on the situations where these
commands are appropriate or discover a new installation process that works for
your system when none of the previously described standard approaches work,
please consider contributing to this documentation (see :ref:`contributing`).

Problems have been encountered when installing in a conda environment. In this
case, pip seems to ignore the installed NumPy version when installing. This
appears to result in a successful installation that fails upon import or causes
a RuntimeError. This happens when the version of NumPy used to build
:py:mod:`apexpy` is newer than the system version of NumPy (NumPy may not be
forwards compatible). In this case, try::

  pip install apexpy --no-build-isolation --no-cache


Apple Silicon systems require certain compilation flags to deal with memory
problems. :py:mod:`apexpy` may appear to install and import correctly, but then
fail with BUS errors when used. In this case, the following command has worked::

  CFLAGS="-falign-functions=8 ${CFLAGS}" pip install --no-binary :apexpy: apexpy


If you are on Apple and encounter a library error such as
``ld: library not found for -lm``, you will need to provide an additional
linking flag to the Mac OSX SDK library::

  LDFLAGS="-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib ${LDFLAGS}" pip install .

This example assumes you are building
locally from the cloned Git repository.  Issues on Mac OS have also been
encountered when using clang for ``CC`` alongside gfortran.  This resulted in a
seemly successful installation with :py:mod:`apexpy` reporting that fortranapex
cannot be imported.

Some users have reported unusual behavior when using Anaconda on Apple Silicon
systems.  Anaconda will attempt to build and install the Intel versions of
wheels instead of the M1 versions and run everything through Rosetta.  This
configuration has not been fully evaluated, but it results in a seemly
successful installation with :py:mod:`apexpy` reporting that fortranapex
cannot be imported.  Users should confirm that wheels created by conda (both for
apexpy and other packages) end in ``arm64.whl`` not ``osx-64.whl``.  If the
later is true, users should consider uninstalling anaconda completely, and
instead installing miniconda following
`these instructions <https://conda.io/projects/conda/en/stable/user-guide/install/macos.html>`_,
which has been confirmed to work. **WARNING:** This will remove any environments
you have set up and likely undo all IDE settings, so be cautious and consider
backing up your work first!


Windows systems are known to have issues with Fortran-based codes.  The Windows
testing we do uses miniconda, so we recommend using the Anaconda environment.
One problem that has been encountered is a lack of LAPACK/BLAS tools that
causes NumPy to not behave as expected.  This can be fixed by installing
scipy before NumPy and then installing :py:mod:`apexpy`.


.. [1] pip is included with Python 2 from v2.7.9 and Python 3 from v3.4.
       If you don't have pip,
       `get it here <https://pip.pypa.io/en/stable/installing/>`_.
