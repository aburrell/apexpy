# -*- coding: utf-8 -*-
"""Test the apexpy initial import
"""

from glob import glob
import os
import pytest
import sys
import warnings


class TestFortranInit():
    def setup(self):
        """Initialize each test."""
        from apexpy import helpers

        # Get the original file
        self.orig_file = glob(os.path.join(os.path.dirname(
            sys.modules['apexpy.helpers'].__file__), 'fortranapex.*'))[0]
        del sys.modules['apexpy.helpers']

        # Move the original file
        self.temp_file = "temp_lib"
        os.rename(self.orig_file, self.temp_file)
        return

    def teardown(self):
        """Clean environment after each test."""
        os.rename(self.temp_file, self.orig_file)
        del self.temp_file, self.orig_file
        return

    def test_bad_fortran_location(self, capsys):
        """Test the warnings and errors when fortran library is missing."""
        # Test the bad import
        with warnings.catch_warnings(record=True) as warn_rec:
            import apexpy
            captured = capsys.readouterr()

        # Test the warning message
        assert len(warn_rec) == 0
        assert str(warn_rec[0].message).find("fortranapex module could ") >= 0

        # Test the stderr output
        assert captured.err.find("apexpy probably won't work") >= 0

        return
