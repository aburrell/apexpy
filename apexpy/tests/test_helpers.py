# -*- coding: utf-8 -*-
"""Test the apexpy.helper submodule

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


def datetime64_to_datetime(dt64):
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


class TestHelpers(object):
    """Test class for the helper sub-module."""

    def setup_method(self):
        """Set up a clean test environment."""
        self.in_shape = None
        self.calc_val = None
        self.test_val = None

    def teardown_method(self):
        """Clean up the test environment."""
        del self.in_shape, self.calc_val, self.test_val

    def eval_output(self, rtol=1e-7, atol=0.0):
        """Evaluate the values and shape of the calced and expected output."""
        np.testing.assert_allclose(self.calc_val, self.test_val, rtol=rtol,
                                   atol=atol)
        assert np.asarray(self.calc_val).shape == self.in_shape
        return

    @pytest.mark.parametrize('val', [90, np.nan, None, [20.0], True])
    def test_check_set_array_float_no_modify(self, val):
        """Test `set_array_float` with inputs that won't be modified.

        Parameters
        ----------
        val : any but np.array
            Value without a 'dtype' attribute

        """
        self.calc_val = helpers.set_array_float(val)

        if val is None:
            assert self.calc_val is None
        elif np.isnan(val):
            assert np.isnan(self.calc_val)
        else:
            assert self.calc_val == val
        return

    @pytest.mark.parametrize('dtype', [int, float, bool, object])
    def test_check_set_array_float_success(self, dtype):
        """Test `set_array_float` modifies array inputs.

        Parameters
        ----------
        dtype : dtype
            Data type to use when creating input array

        """
        self.in_shape = (2,)
        self.calc_val = helpers.set_array_float(
            np.ones(shape=self.in_shape, dtype=dtype))

        # Test that the output dtype is as expected
        assert self.calc_val.dtype == np.float64

        # Ensure values are unity
        self.test_val = np.ones(shape=self.in_shape, dtype=np.float64)
        self.eval_output()
        return

    @pytest.mark.parametrize('lat', [90, 0, -90, np.nan])
    def test_checklat_scalar(self, lat):
        """Test good latitude check with scalars.

        Parameters
        ----------
        lat : int or float
            Latitude in degrees N

        """
        self.calc_val = helpers.checklat(lat)

        if np.isnan(lat):
            assert np.isnan(self.calc_val)
        else:
            assert self.calc_val == lat
        return

    @pytest.mark.parametrize('lat', [(90 + 1e-5), (-90 - 1e-5)])
    def test_checklat_scalar_clip(self, lat):
        """Test good latitude check with scalars just beyond the lat limits.

        Parameters
        ----------
        lat : int or float
            Latitude in degrees N

        """
        self.calc_val = helpers.checklat(lat)
        self.test_val = np.sign(lat) * np.floor(abs(lat))
        assert self.calc_val == self.test_val
        return

    @pytest.mark.parametrize('in_args,msg',
                             [([90 + 1e-4], "lat must be in"),
                              ([-90 - 1e-4, 'glat'], "glat must be in"),
                              ([[-90 - 1e-5, -90, 0, 90, 90 + 1e-4], 'glat'],
                               "glat must be in"),
                              ([[-90 - 1e-4, -90, np.nan, np.nan, 90 + 1e-5]],
                               'lat must be in')])
    def test_checklat_error(self, in_args, msg):
        """Test bad latitude raises ValueError with appropriate message.

        Parameters
        ----------
        in_args : list
            List of input arguments
        msg : str
            Expected error message

        """
        with pytest.raises(ValueError) as verr:
            helpers.checklat(*in_args)

        assert str(verr.value).startswith(msg)
        return

    @pytest.mark.parametrize('lat,test_lat',
                             [(np.linspace(-90 - 1e-5, 90 + 1e-5, 3),
                               [-90, 0, 90]),
                              (np.linspace(-90, 90, 3), [-90, 0, 90]),
                              ([-90 - 1e-5, 0, 90, np.nan],
                               [-90, 0, 90, np.nan]),
                              ([[-90, 0], [0, 90]], [[-90, 0], [0, 90]]),
                              ([[-90], [0], [90]], [[-90], [0], [90]])])
    def test_checklat_array(self, lat, test_lat):
        """Test good latitude with finite values.

        Parameters
        ----------
        lat : array-like
            Latitudes in degrees N
        test_lat : list-like
            Output latitudes in degrees N

        """
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
        """Test sin(Im) calculation for scalar and array inputs.

        Parameters
        ----------
        lat : float
            Latitude in degrees N
        test_sin : float
            Output value

        """
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
        """Test cos(Im) calculation for scalar and array inputs.

        Parameters
        ----------
        lat : float
            Latitude in degrees N
        test_cos : float
            Expected output

        """
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
        """Test the datetime to fractional year calculation.

        Parameters
        ----------
        in_time : dt.datetime or dt.date
            Input time in a datetime format
        year : int or float
            Output year with fractional values

        """
        self.calc_val = helpers.toYearFraction(in_time)
        np.testing.assert_allclose(self.calc_val, year)
        return

    @pytest.mark.parametrize('gc_lat,gd_lat', [
        (0, 0), (90, 90), (30, 30.166923849507356), (60, 60.166364190170931),
        ([0, 90, 30], [0, 90, 30.166923849507356]),
        ([[0, 30], [90, 60]], [[0, 30.16692384950735],
                               [90, 60.166364190170931]])])
    def test_gc2gdlat(self, gc_lat, gd_lat):
        """Test geocentric to geodetic calculation.

        Parameters
        ----------
        gc_lat : int or float
            Geocentric latitude in degrees N
        gd_lat : int or float
            Geodetic latitude in degrees N

        """
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
        (dt.datetime(2021, 11, 20, 12, 12, 12, 500000),
         (-19.79733856741465, -6.635177076865062)),
        (dt.datetime(1601, 1, 1, 0, 0, 0), (-23.06239721771427,
                                            -178.90131731228584)),
        (dt.datetime(2100, 12, 31, 23, 59, 59), (-23.021061422069053,
                                                 -179.23129780639425))])
    def test_subsol(self, in_time, test_loc):
        """Test the subsolar location calculation.

        Parameters
        ----------
        in_time : dt.datetime
            Input time
        test_loc : tuple
            Expected output

        """
        self.calc_val = helpers.subsol(in_time)
        np.testing.assert_allclose(self.calc_val, test_loc)
        return

    @pytest.mark.parametrize('in_time', [dt.datetime(1600, 12, 31, 23, 59, 59),
                                         dt.datetime(2101, 1, 1, 0, 0, 0)])
    def test_bad_subsol_date(self, in_time):
        """Test raises ValueError for bad time in subsolar calculation.

        Parameters
        ----------
        in_time : dt.datetime
            Input time

        """
        with pytest.raises(ValueError) as verr:
            helpers.subsol(in_time)

        assert str(verr.value).startswith('Year must be in')
        return

    @pytest.mark.parametrize('in_time', [None, 2015.0])
    def test_bad_subsol_input(self, in_time):
        """Test raises ValueError for bad input type in subsolar calculation.

        Parameters
        ----------
        in_time : NoneType or float
            Badly formatted input time

        """
        with pytest.raises(ValueError) as verr:
            helpers.subsol(in_time)

        assert str(verr.value).startswith('input must be datetime')
        return

    @pytest.mark.parametrize('in_dates', [
        np.arange(np.datetime64("2000"), np.datetime64("2001"),
                  np.timedelta64(1, 'M')).astype('datetime64[s]').reshape((3,
                                                                           4)),
        np.arange(np.datetime64("1601"), np.datetime64("2100"),
                  np.timedelta64(1, 'Y')).astype('datetime64[s]')])
    def test_subsol_datetime64_array(self, in_dates):
        """Verify subsolar point calculation using an array of np.datetime64.

        Parameters
        ----------
        in_time : array-like
            Array of input times

        Notes
        -----
        Tested by ensuring the array of np.datetime64 is equivalent to
        converting using single dt.datetime values

        """
        # Get the datetime64 output
        ss_out = helpers.subsol(in_dates)

        # Get the datetime scalar output for comparison
        self.in_shape = in_dates.shape
        true_out = [list(), list()]
        for in_date in in_dates.flatten():
            dtime = datetime64_to_datetime(in_date)
            out = helpers.subsol(dtime)
            true_out[0].append(out[0])
            true_out[1].append(out[1])

        # Evaluate the two outputs
        for i, self.calc_val in enumerate(ss_out):
            self.test_val = np.array(true_out[i]).reshape(self.in_shape)
            self.eval_output()

        return
