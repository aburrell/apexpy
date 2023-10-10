from sys import stderr

# Below try..catch required for autodoc to work on readthedocs
try:
    from apexpy import fortranapex  # noqa F401
except ImportError as ierr:
    stderr.write("".join(["fortranapex module could not be imported. ",
                          "apexpy probably won't work. Failed with error: ",
                          str(ierr)]))

from apexpy.apex import Apex, ApexHeightError  # noqa F401
from apexpy import helpers  # noqa F401

# Define the global variables
__version__ = "2.0.1"
__all__ = ['Apex', 'fortranapex', 'helpers', 'ApexHeightError']
