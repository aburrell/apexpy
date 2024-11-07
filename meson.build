project( 'apexpy', 'c',
  # Note that the git commit hash cannot be added dynamically here
  version: '2.0.1',
  license: 'MIT',
  meson_version: '>= 0.63.0',
  default_options: [
    'buildtype=debugoptimized',
    # TODO: the below -Wno flags are all needed to silence warnings in
    # f2py-generated code. This should be fixed in f2py itself.
    'c_args=-Wno-unused-function -Wno-conversion -Wno-misleading-indentation -Wno-incompatible-pointer-types',
    'fortran_args=-Wno-conversion',  # silence "conversion from REAL(8) to INTEGER(4)"
    'fortran_std=legacy',
  ],
)

# Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)


fc = meson.get_compiler('fortran')
cc = meson.get_compiler('c')

# Don't use the deprecated NumPy C API. Define this to a fixed version instead
# of NPY_API_VERSION in order not to break compilation for released versions
# when NumPy introduces a new deprecation. Use in a meson.build file::
#
#   py3.extension_module('_name',
#     'source_fname',
#     numpy_nodepr_api)
#
numpy_nodepr_api = '-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION'

# This argument is called -Wno-unused-but-set-variable by GCC, however Clang
# doesn't recognize that.
if cc.has_argument('-Wno-unused-but-set-variable')
  add_global_arguments('-Wno-unused-but-set-variable', language : 'c')
endif

# For Fortran code, Meson already adds `-lm`.
m_dep = cc.find_library('m', required : false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

# Add more link arguments
add_project_link_arguments('-lquadmath', language: ['c', 'fortran'])

# https://mesonbuild.com/Python-module.html
py_mod = import('python')
py3 = py_mod.find_installation()
py3_dep = py3.dependency()
message(py3.full_path())
message(py3.get_install_dir())

incdir_numpy = run_command(py3,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py3,
  ['-c', 'import os; os.chdir(".."); from numpy import f2py; print(f2py.get_include())'],
  check : true
).stdout().strip()

inc_dirs = include_directories(incdir_numpy, incdir_f2py)

# Don't use the deprecated NumPy C API. Define this to a fixed version instead of
# NPY_API_VERSION in order not to break compilation for released SciPy versions
# when NumPy introduces a new deprecation. Use in a meson.build file::
#
#   py3.extension_module('_name',
#     'source_fname',
#     numpy_nodepr_api)
#
c_flags = ['-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION']

# Platform detection to set more flags for Windows systems
is_windows = host_machine.system() == 'windows'
is_mac = host_machine.system() == 'darwin'
is_mingw = is_windows and cc.get_id() == 'gcc'

cython_c_args = []
if is_windows
  # For mingw-w64, link statically against the UCRT.
  gcc_link_args = ['-lucrt', '-static']
  if is_mingw
    add_global_link_arguments(gcc_link_args, language: ['c'])

    # Force gcc to float64 long doubles for compatibility with MSVC
    # builds, for C only.
    add_global_arguments('-mlong-double-64', language: 'c')

    # Make fprintf("%zd") work (see https://github.com/rgommers/scipy/issues/118)
    add_global_arguments('-D__USE_MINGW_ANSI_STDIO=1', language: ['c'])

    # Manual add of MS_WIN64 macro when not using MSVC.
    # https://bugs.python.org/issue28267
    bitness = run_command('apexpy/_gcc_build_bitness.py').stdout().strip()
    if bitness == '64'
      add_global_arguments('-DMS_WIN64', language: ['c', 'fortran'])
    endif

    # Silence warnings emitted by PyOS_snprintf for (%zd), see
    # https://github.com/rgommers/scipy/issues/118.
    # Use as c_args for extensions containing Cython code
    cython_c_args += ['-Wno-format-extra-args', '-Wno-format']

    # Windows may need the directory for Python.h added.  This location
    # does not have the same name for all installs.  This is the one
    # for Windows on GitHub actions
    incdir_py = run_command(py3,
      ['-c', 'import os; os.chdir(".."); import sys; print(os.path.join(sys.prefix, "include"))'],
      check : true).stdout().strip()
    inc_dirs = include_directories(incdir_numpy, incdir_f2py, incdir_py)
  endif

  if meson.get_compiler('fortran').get_id() == 'gcc'
    add_global_link_arguments(gcc_link_args, language: ['fortran'])

    # Flag needed to work around BLAS and LAPACK Gfortran dependence on
    # undocumented C feature when passing single character string
    # arguments.
    # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90329
    # https://github.com/wch/r-source/blob/838f9d5a7be08f2a8c08e47bcd28756f5d0aac90/src/gnuwin32/MkRules.rules#L121
    add_global_arguments('-fno-optimize-sibling-calls',
      language: ['fortran'])
  endif
elif is_mac
  c_flags += ['-FALIGN_FUNCTIONS=8']
endif

# Check the python headers
cc.check_header('Python.h', dependencies: [py3_dep], required: true)

# Unlike distutils, meson doesn't yet include some of the f2py stuff
fortranobject_c = incdir_f2py / 'fortranobject.c'

fortranobject_lib = static_library('_fortranobject',
  fortranobject_c,
  c_args: c_flags,
  dependencies: py3_dep,
  include_directories: [incdir_numpy, incdir_f2py])

fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  include_directories: [incdir_numpy, incdir_f2py])

# Make a custom target, output should not contain a path segment.
# The name 'fortranapex-f2pywrappers2.f90' is used because an empty
# 'fortranapex-f2pywrappers.f' is also created, making the F90 wrapper the
# second wrapper.  F2PY controls these names, so it may change at some point.
fortranapex_source = custom_target(
  'fortranapexmodule.c',
   input : ['fortranapex/magfld.f90', 'fortranapex/apex.f90',
            'fortranapex/makeapexsh.f90', 'fortranapex/igrf.f90',
	    'fortranapex/apexsh.f90', 'fortranapex/checkapexsh.f90'],
   output : ['fortranapexmodule.c', 'fortranapex-f2pywrappers2.f90'],
   command : [py3, '-m', 'numpy.f2py', '@INPUT@', '-m', 'fortranapex',
              '--lower'])

# Declare the fortran extension module
py3.extension_module('fortranapex',
  ['fortranapex/magfld.f90', 'fortranapex/apex.f90',
   'fortranapex/makeapexsh.f90', 'fortranapex/igrf.f90',
   'fortranapex/apexsh.f90', 'fortranapex/checkapexsh.f90',
   fortranapex_source, fortranobject_c],
  c_args: c_flags,
  include_directories: inc_dirs,
  link_with: fortranobject_lib,
  dependencies : [py3_dep, fortranobject_dep],
  subdir: 'apexpy',
  install : true)

# Declare the sources
py3.install_sources([
  'apexpy/__init__.py',
  'apexpy/__main__.py',
  'apexpy/apex.py',
  'apexpy/apexsh.dat',
  'apexpy/helpers.py',
  'apexpy/igrf13coeffs.txt'],
  pure: false,
  subdir: 'apexpy')
