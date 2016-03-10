from __future__ import absolute_import, division, print_function

from .apex import Apex, ApexHeightError
from . import helpers
# below try..catch required for autodoc to work on readthedocs
try:
    from . import fortranapex
except ImportError:
    print("ERROR: fortranapex module could not be imported. apexpy probably won't work")


__version__ = "1.0.1"

__all__ = ['Apex', 'fortranapex', 'helpers', 'ApexHeightError']
