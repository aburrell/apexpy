project(
  'apexpy',
  'fortran',
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
    'blas=openblas',
    'lapack=openblas'
  ],
)

# For Fortran code, Meson already adds `-lm`.
fc = meson.get_compiler('fortran')

# Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)

cython = find_program('cython')
pythran = find_program('pythran')
copier = find_program(['cp', 'apexpy/_copyfiles.py'])

# https://mesonbuild.com/Python-module.html
py_mod = import('python')
py3 = py_mod.find_installation()
py3_dep = py3.dependency()

subdir('apexpy')