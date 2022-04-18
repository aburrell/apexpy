from sys import stderr

from apexpy.apex import Apex, ApexHeightError
from apexpy import helpers

# Below try..catch required for autodoc to work on readthedocs
try:
    from apexpy import fortranapex
except ImportError:
    stderr.write("".join(["fortranapex module could not be imported. ",
                          "apexpy probably won't work"]))

# Define the global variables
__version__ = "1.1.0"
__all__ = ['Apex', 'fortranapex', 'helpers', 'ApexHeightError']
