project( 'apexpy', 'fortran', 'c',
  # Note that the git commit hash cannot be added dynamically here
  version: '1.1.0',
  license: 'MIT',
  meson_version: '>= 0.60',
  default_options: [
    'buildtype=debugoptimized',
    # TODO: the below -Wno flags are all needed to silence warnings in
    # f2py-generated code. This should be fixed in f2py itself.
    'c_args=-Wno-unused-function -Wno-conversion -Wno-misleading-indentation -Wno-incompatible-pointer-types',
    'fortran_args=-Wno-conversion',
    'fortran_std=legacy',
  ],
)

fc = meson.get_compiler('fortran')
cc = meson.get_compiler('c')

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
add_project_link_arguments('-lquadmath', language: 'c')

# Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)

copier = find_program(['cp', 'apexpy/_copyfiles.py'])

# https://mesonbuild.com/Python-module.html
py_mod = import('python')
py3 = py_mod.find_installation()
py3_dep = py3.dependency()

incdir_numpy = run_command(py3,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py3,
  ['-c', 'import os; os.chdir(".."); from numpy import f2py; print(f2py.get_include())'],
  check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

py3.extension_module('apexpy.fortranapex',
  'fortranapex/magfld.f', 'fortranapex/apex.f', 'fortranapex/makeapexsh.f90',
  'fortranapex/igrf.f90', 'fortranapex/apexsh.f90',
  'fortranapex/checkapexsh.f90', 'fortranapex/fortranapexmodule.c',
  include_directories: inc_np,
  dependencies : py3_dep,
  install : true)
