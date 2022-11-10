from sys import stderr

# Below try..catch required for autodoc to work on readthedocs
try:
    from apexpy import fortranapex
except ImportError as ierr:
    stderr.write("".join(["fortranapex module could not be imported. ",
                          "apexpy probably won't work. Failed with error: ",
                          str(ierr)]))

from apexpy.apex import Apex, ApexHeightError
from apexpy import helpers

# Define the global variables
__version__ = "1.1.0"
__all__ = ['Apex', 'fortranapex', 'helpers', 'ApexHeightError']
