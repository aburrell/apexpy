# -*- coding: utf-8 -*-
"""Test the apexpy.fortranapex class

Notes
-----
Whenever function outputs are tested against hard-coded numbers, the test
results (numbers) were obtained by running the code that is tested.  Therefore,
these tests below only check that nothing changes when refactoring, etc., and
not if the results are actually correct.

These results are expected to change when IGRF is updated.

"""

from numpy.testing import assert_allclose
from importlib import resources
import os
import pytest

import apexpy
from apexpy import fortranapex as fa


class TestFortranApex(object):
    """Test class for the Python-wrapped Fortran functions."""

    def setup_method(self):
        """Initialize each test."""
        datafile = os.path.join(resources.files(apexpy), 'apexsh.dat')
        fa.loadapxsh(datafile, 2000)

        # Set the inputs
        self.lat = 60
        self.lon = 15
        self.height = 100
        self.refh = 300
        self.vecflg = 1
        self.precision = 1e-10

        # Set the output values
        self.glat = 50.979549407958984
        self.glon = -66.16891479492188
        self.error = 2.8316469524725107e-06
        self.qlat = 56.531288146972656
        self.qlon = 94.1068344116211
        self.mlat = 55.94841766357422
        self.mlon = 94.1068344116211
        self.f1 = [1.079783, 0.10027137]
        self.f2 = [-0.24546318, 0.90718889]
        self.fmag = 1.0041800737380981
        self.d1 = [0.9495701, 0.25693053, 0.09049474]
        self.d2 = [0.10011087, -1.0780545, -0.33892432]
        self.d3 = [0.00865356, 0.27327004, -0.8666646]
        self.dmag = 1.100391149520874
        self.e1 = [1.0269295, 0.08382964, 0.03668632]
        self.e2 = [0.24740215, -0.82374191, -0.25726584]
        self.e3 = [0.01047826, 0.33089194, -1.04941]
        self.out = None
        self.test_out = None

    def teardown_method(self):
        """Clean environment after each test."""
        del self.lat, self.lon, self.height, self.refh, self.vecflg
        del self.qlat, self.qlon, self.mlat, self.mlon, self.f1, self.f2
        del self.fmag, self.d1, self.d2, self.d3, self.dmag, self.e1
        del self.e2, self.e3, self.precision, self.glat, self.glon, self.error
        del self.out, self.test_out

    def run_test_evaluation(self, rtol=1e-5, atol=1e-5):
        """Run the evaluation of the test results.

        Parameters
        ----------
        rtol : float
            Relative tolerance, default value based on old code's precision
            (default=1e-5)
        atol : float
            Absolute tolerance, default value based on old code's precision
            (default=1e-5)

        """

        assert self.out is not None, "No results to test"
        assert self.test_out is not None, "No 'truth' results provided"
        assert len(self.out) == len(self.test_out), "Mismatched outputs"

        for i, out_val in enumerate(self.out):
            assert_allclose(out_val, self.test_out[i], rtol=rtol, atol=atol)
        return

    def test_apxg2q(self):
        """Test fortran apex geographic to quasi-dipole."""
        # Get the output
        self.out = fa.apxg2q(self.lat, self.lon, self.height, self.vecflg)

        # Test the output
        self.test_out = (self.qlat, self.qlon, self.f1, self.f2, self.fmag)
        self.run_test_evaluation()
        return

    def test_apxg2all(self):
        """Test fortran apex geographic to all outputs."""
        # Get the output
        self.out = fa.apxg2all(self.lat, self.lon, self.height, self.refh,
                               self.vecflg)

        # Test the output
        self.test_out = (self.qlat, self.qlon, self.mlat, self.mlon, self.f1,
                         self.f2, self.fmag, self.d1, self.d2, self.d3,
                         self.dmag, self.e1, self.e2, self.e3)
        self.run_test_evaluation()
        return

    def test_apxq2g(self):
        """ Test fortran quasi-dipole to geographic."""
        # Get the output
        self.out = fa.apxq2g(self.lat, self.lon, self.height, self.precision)

        # Test the output
        self.test_out = (self.glat, self.glon, self.error)
        self.run_test_evaluation()
        return

    @pytest.mark.parametrize("lat", [0, 30, 60, 89])
    @pytest.mark.parametrize("lon", [-179, -90, 0, 90, 179])
    def test_g2q2d(self, lat, lon):
        """ Test fortran geographic to quasi-dipole and back again.

        Parameters
        ----------
        lat : int or float
            Latitude in degrees N
        lon : int or float
            Longitude in degrees E

        """
        self.out = fa.apxg2q(lat, lon, self.height, 0)
        self.test_out = fa.apxq2g(self.out[0], self.out[1], self.height,
                                  self.precision)

        # Test the results againt the initial input
        assert_allclose(self.test_out[0], lat, atol=0.01)
        assert_allclose(self.test_out[1], lon, atol=0.01)

    def test_apxq2g_lowprecision(self):
        """Test low precision error value."""
        self.out = fa.apxq2g(self.lat, self.lon, self.height, -1)

        # Test the output
        self.test_out = (51.00891876220703, -66.11973571777344, -9999.0)
        self.run_test_evaluation()
        return
