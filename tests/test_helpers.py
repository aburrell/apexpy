# -*- coding: utf-8 -*-
"""Test the apexpy.Apex class

Notes
-----
Whenever function outputs are tested against hard-coded numbers, the test
results (numbers) were obtained by running the code that is tested.  Therefore,
these tests below only check that nothing changes when refactoring, etc., and
not if the results are actually correct.

These results are expected to change when IGRF is updated.

"""

import datetime as dt
import numpy as np
import pytest

from apexpy import helpers


class TestHelpers():
    def setup(self):
        self.in_shape = None
        self.calc_val = None
        self.test_val = None

    def teardown(self):
        del self.in_shape, self.calc_val, self.test_val

    def eval_output(self, rtol=1e-7, atol=0.0):
        """Evaluate the values and shape of the calculated and expected output.
        """
        np.testing.assert_allclose(self.calc_val, self.test_val, rtol=rtol,
                                   atol=atol)
        assert np.asarray(self.calc_val).shape == self.in_shape
        return

    def datetime64_to_datetime(self, dt64):
        """Convert numpy datetime64 object to a datetime datetime object.

        Parameters
        ----------
        dt64 : np.datetime64
            Numpy datetime64 object

        Returns
        -------
        dt.datetime
            Equivalent datetime object with a resolution of days

        Notes
        -----
        Works outside 32 bit int second range of 1970

        """
        year_floor = dt64.astype('datetime64[Y]')
        month_floor = dt64.astype('datetime64[M]')
        day_floor = dt64.astype('datetime64[D]')
        year = year_floor.astype(int) + 1970
        month = (month_floor
                 - year_floor).astype('timedelta64[M]').astype(int) + 1
        day = (day_floor - month_floor).astype('timedelta64[D]').astype(int) + 1
        return dt.datetime(year, month, day)

    @pytest.mark.parametrize('lat', [90, 0, -90, np.nan])
    def test_checklat_scalar(self, lat):
        """Test good latitude check with scalars."""
        if np.isnan(lat):
            assert np.isnan(helpers.checklat(lat))
        else:
            assert helpers.checklat(lat) == lat
        return

    @pytest.mark.parametrize('lat', [(90 + 1e-5), (-90 - 1e-5)])
    def test_checklat_scalar_clip(self, lat):
        """Test good latitude check with scalars just beyond the lat limits."""
        assert helpers.checklat(lat) == np.sign(lat) * np.floor(abs(lat))
        return

    @pytest.mark.parametrize('in_args,msg',
                             [([90 + 1e-4], "lat must be in"),
                              ([-90 - 1e-4, 'glat'], "glat must be in"),
                              ([[-90 - 1e-5, -90, 0, 90, 90 + 1e-4], 'glat'],
                               "glat must be in"),
                              ([[-90 - 1e-4, -90, np.nan, np.nan, 90 + 1e-5]],
                               'lat must be in')])
    def test_checklat_error(self, in_args, msg):
        """Test bad latitude raises ValueError with appropriate message."""
        with pytest.raises(ValueError) as verr:
            helpers.checklat(*in_args)

        assert str(verr.value).startswith(msg)
        return

    @pytest.mark.parametrize('lat,test_lat',
                             [(np.linspace(-90 - 1e-5, 90 + 1e-5, 3),
                               [-90, 0, 90]),
                              (np.linspace(-90, 90, 3), [-90, 0, 90]),
                              ([-90 - 1e-5, 0, 90, np.nan],
                               [-90, 0, 90, np.nan])])
    def test_checklat_array(self, lat, test_lat):
        """Test good latitude with finite values."""
        self.calc_val = helpers.checklat(lat)
        self.in_shape = np.asarray(lat).shape
        self.test_val = test_lat
        self.eval_output(atol=1e-8)
        return

    @pytest.mark.parametrize('lat,test_sin', [
        (60, 0.96076892283052284), (10, 0.33257924500670238),
        ([60, 10], [0.96076892283052284, 0.33257924500670238]),
        ([[60, 10], [60, 10]], [[0.96076892283052284, 0.33257924500670238],
                                [0.96076892283052284, 0.33257924500670238]])])
    def test_getsinIm(self, lat, test_sin):
        """Test sin(Im) calculation for scalar and array inputs."""
        self.calc_val = helpers.getsinIm(lat)
        self.in_shape = np.asarray(lat).shape
        self.test_val = test_sin
        self.eval_output()
        return

    @pytest.mark.parametrize('lat,test_cos', [
        (60, 0.27735009811261463), (10, 0.94307531289434765),
        ([60, 10], [0.27735009811261463, 0.94307531289434765]),
        ([[60, 10], [60, 10]], [[0.27735009811261463, 0.94307531289434765],
                                [0.27735009811261463, 0.94307531289434765]])])
    def test_getcosIm(self, lat, test_cos):
        """Test cos(Im) calculation for scalar and array inputs."""
        self.calc_val = helpers.getcosIm(lat)
        self.in_shape = np.asarray(lat).shape
        self.test_val = test_cos
        self.eval_output()
        return

    @pytest.mark.parametrize('in_time,year', [
        (dt.datetime(2001, 1, 1), 2001), (dt.date(2001, 1, 1), 2001),
        (dt.datetime(2002, 1, 1), 2002),
        (dt.datetime(2005, 2, 3, 4, 5, 6), 2005.090877283105),
        (dt.datetime(2005, 12, 11, 10, 9, 8), 2005.943624682902)])
    def test_toYearFraction(self, in_time, year):
        """Test the datetime to fractional year calculation."""
        self.calc_val = helpers.toYearFraction(in_time)
        np.testing.assert_allclose(self.calc_val, year)
        return

    @pytest.mark.parametrize('gc_lat,gd_lat', [
        (0, 0), (90, 90), (30, 30.166923849507356), (60, 60.166364190170931),
        ([0, 90, 30], [0, 90, 30.166923849507356]),
        ([[0, 30], [90, 60]], [[0, 30.16692384950735],
                               [90, 60.166364190170931]])])
    def test_gc2gdlat(self, gc_lat, gd_lat):
        """Test geocentric to geodetic calculation."""
        self.calc_val = helpers.gc2gdlat(gc_lat)
        self.in_shape = np.asarray(gc_lat).shape
        self.test_val = gd_lat
        self.eval_output()
        return

    @pytest.mark.parametrize('in_time,test_loc', [
        (dt.datetime(2005, 2, 3, 4, 5, 6), (-16.505391672592904,
                                            122.17768157084515)),
        (dt.datetime(2010, 12, 11, 10, 9, 8), (-23.001554595838947,
                                               26.008999999955023)),
        (dt.datetime(1601, 1, 1, 0, 0, 0), (-23.06239721771427,
                                            -178.90131731228584)),
        (dt.datetime(2100, 12, 31, 23, 59, 59), (-23.021061422069053,
                                                 -179.23129780639425))])
    def test_subsol(self, in_time, test_loc):
        """Test the subsolar location calculation."""
        self.calc_val = helpers.subsol(in_time)
        np.testing.assert_allclose(self.calc_val, test_loc)
        return

    @pytest.mark.parametrize('in_time', [dt.datetime(1600, 12, 31, 23, 59, 59),
                                         dt.datetime(2101, 1, 1, 0, 0, 0)])
    def test_bad_subsol(self, in_time):
        """Test raises ValueError for bad time in subsolar calculation."""
        with pytest.raises(ValueError) as verr:
            helpers.subsol(in_time)

        assert str(verr.value).startswith('Year must be in')
        return

    def test_subsol_array(self):
        """Verify subsolar point calculation using an array of np.datetime64.

        Notes
        -----
        Tested by ensuring the array of np.datetime64 is equivalent to
        converting using single dt.datetime values

        """
        in_dates = np.arange(np.datetime64("1601"), np.datetime64("2100"),
                             np.timedelta64(100, 'D')).astype('datetime64[s]')
        sslat, sslon = helpers.subsol(in_dates)

        # Test the shape of the output
        assert sslat.shape == in_dates.shape
        assert sslon.shape == in_dates.shape

        # Test the values
        for i, in_date in enumerate(in_dates):
            dtime = self.datetime64_to_datetime(in_date)
            true_sslat, true_sslon = helpers.subsol(dtime)
            assert sslat[i] == true_sslat
            assert sslon[i] == true_sslon
        return
