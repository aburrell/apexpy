# -*- coding: utf-8 -*-

from __future__ import division, absolute_import, unicode_literals

import datetime as dt

import numpy as np
import pytest
from numpy.testing import assert_allclose

from apexpy import helpers


# ----------------------------------------------------------------------------
# NOTE: whenever function outputs are tested against hard-coded numbers, the
# test results (numbers) were obtained by running the code that is tested.
# Therefore these tests below only check that nothing changes when refactoring,
# etc., and not if the results are actually correct
# ----------------------------------------------------------------------------


# ============================================================================
#  Test checklat
# ============================================================================

@pytest.mark.parametrize('lat', [(90), (0), (-90), (np.nan)])
def test_checklat_scalar(lat):
    """Test good latitude check with scalars."""
    if np.isnan(lat):
        assert np.isnan(helpers.checklat(lat))
    else:
        assert helpers.checklat(lat) == lat


@pytest.mark.parametrize('lat', [(90 + 1e-5), (-90 - 1e-5)])
def test_checklat_scalar_clip(lat):
    """Test good latitude check with scalars just beyond the lat limits."""
    assert helpers.checklat(lat) == np.sign(lat) * np.floor(abs(lat))


@pytest.mark.parametrize('in_args,msg',
                         [([90 + 1e-4], "lat must be in"),
                          ([-90 - 1e-4, 'glat'], "glat must be in"),
                          ([[-90 - 1e-5, -90, 0, 90, 90 + 1e-4], 'glat'],
                           "glat must be in"),
                          ([[-90 - 1e-4, -90, np.nan, np.nan, 90 + 1e-5]],
                           'lat must be in')])
def test_checklat_error(in_args, msg):
    """Test bad latitude raises ValueError with appropriate message."""
    with pytest.raises(ValueError) as verr:
        helpers.checklat(*in_args)

    assert str(verr.value).startswith(msg)


def test_checklat_array():
    """Test good latitude with finite values."""
    assert_allclose(helpers.checklat([-90 - 1e-5, -90, 0, 90, 90 + 1e-5]),
                    np.array([-90, -90, 0, 90, 90]), rtol=0, atol=1e-8)

    return


def test_checklat_array_withnan():
    """Test good latitude input mixed with NaNs."""

    in_lat = np.array([-90 - 1e-5, -90, 0, 90, 90 + 1e-5, np.nan, np.nan])
    fin_mask = np.isfinite(in_lat)
    out_lat = helpers.checklat(in_lat)
    assert_allclose(np.array([-90, -90, 0, 90, 90]), out_lat[fin_mask],
                    rtol=0, atol=1e-8)

    assert np.all(np.isnan(out_lat[~fin_mask]))


# ============================================================================
#  Test getsinIm
# ============================================================================

def test_getsinIm_scalar():
    assert_allclose(helpers.getsinIm(60), 0.96076892283052284)
    assert_allclose(helpers.getsinIm(10), 0.33257924500670238)
    assert type(helpers.getsinIm(60)) != np.ndarray


def test_getsinIm_1Darray():
    assert_allclose(helpers.getsinIm([60, 10]),
                    [0.96076892283052284, 0.33257924500670238])


def test_getsinIm_2Darray():
    assert_allclose(helpers.getsinIm([[60, 10], [60, 10]]),
                    [[0.96076892283052284, 0.33257924500670238],
                     [0.96076892283052284, 0.33257924500670238]])


# ============================================================================
#  Test getcosIm
# ============================================================================


def test_getcosIm_scalar():
    assert_allclose(helpers.getcosIm(60), 0.27735009811261463)
    assert_allclose(helpers.getcosIm(10), 0.94307531289434765)
    assert type(helpers.getcosIm(60)) != np.ndarray


def test_getcosIm_1Darray():
    assert_allclose(helpers.getcosIm([60, 10]),
                    [0.27735009811261463, 0.94307531289434765])


def test_getcosIm_2Darray():
    assert_allclose(helpers.getcosIm([[60, 10], [60, 10]]),
                    [[0.27735009811261463, 0.94307531289434765],
                     [0.27735009811261463, 0.94307531289434765]])


# ============================================================================
#  Test toYearFraction
# ============================================================================


def test_toYearFraction():
    assert_allclose(helpers.toYearFraction(dt.datetime(2001, 1, 1, 0, 0, 0)),
                    2001)
    assert_allclose(helpers.toYearFraction(dt.date(2001, 1, 1)), 2001)
    assert_allclose(helpers.toYearFraction(dt.datetime(2002, 1, 1, 0, 0, 0)),
                    2002)
    assert_allclose(helpers.toYearFraction(dt.datetime(2005, 2, 3, 4, 5, 6)),
                    2005.090877283105)
    assert_allclose(helpers.toYearFraction(dt.datetime(2005, 12, 11, 10, 9, 8)),
                    2005.943624682902)


# ============================================================================
#  Test gc2gdlat
# ============================================================================


def test_gc2gdlat():
    assert_allclose(helpers.gc2gdlat(0), 0)
    assert_allclose(helpers.gc2gdlat(90), 90)
    assert_allclose(helpers.gc2gdlat(30), 30.166923849507356)
    assert_allclose(helpers.gc2gdlat(60), 60.166364190170931)


# ============================================================================
#  Test subsol
# ============================================================================

def test_subsol():
    assert_allclose(helpers.subsol(dt.datetime(2005, 2, 3, 4, 5, 6)),
                    (-16.505391672592904, 122.17768157084515))
    assert_allclose(helpers.subsol(dt.datetime(2010, 12, 11, 10, 9, 8)),
                    (-23.001554595838947, 26.008999999955023))

    with pytest.raises(ValueError):
        helpers.subsol(dt.datetime(1600, 12, 31, 23, 59, 59))
    assert_allclose(helpers.subsol(dt.datetime(1601, 1, 1, 0, 0, 0)),
                    (-23.06239721771427, -178.90131731228584))
    with pytest.raises(ValueError):
        helpers.subsol(dt.datetime(2101, 1, 1, 0, 0, 0))
    assert_allclose(helpers.subsol(dt.datetime(2100, 12, 31, 23, 59, 59)),
                    (-23.021061422069053, -179.23129780639425))


def datetime64_to_datetime(dt64):
    """Convert numpy datetime64 object to a datetime datetime object.

    Notes
    -----
    Works outside 32 bit int second range of 1970

    """
    year_floor = dt64.astype('datetime64[Y]')
    month_floor = dt64.astype('datetime64[M]')
    day_floor = dt64.astype('datetime64[D]')
    year = year_floor.astype(int) + 1970
    month = (month_floor - year_floor).astype('timedelta64[M]').astype(int) + 1
    day = (day_floor - month_floor).astype('timedelta64[D]').astype(int) + 1
    return dt.datetime(year, month, day)


def test_subsol_array():
    """Verify subsolar point calculation using an array of np.datetime64.

    Notes
    -----
    Tested by ensuring the array of np.datetime64 is equivalent to converting
    using single dt.datetime values

    """
    dates = np.arange(np.datetime64("1601"), np.datetime64("2100"),
                      np.timedelta64(100, 'D')).astype('datetime64[s]')
    sslat, sslon = helpers.subsol(dates)
    for i, date in enumerate(dates):
        datetime = datetime64_to_datetime(date)
        true_sslat, true_sslon = helpers.subsol(datetime)
        assert sslat[i] == true_sslat
        assert sslon[i] == true_sslon


if __name__ == '__main__':
    pytest.main()
