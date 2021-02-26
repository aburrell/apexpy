# -*- coding: utf-8 -*-

from __future__ import division, absolute_import, unicode_literals

import datetime as dt
import itertools
import numpy as np
from numpy.testing import assert_allclose
import os
import pytest

from apexpy import fortranapex as fa
from apexpy import Apex, ApexHeightError, helpers


# ----------------------------------------------------------------------------
# NOTE: whenever function outputs are tested against hard-coded numbers, the
# test results (numbers) were obtained by running the code that is tested.
# Therefore these tests below only check that nothing changes when refactoring,
# etc., and not if the results are actually correct.
# -----------------------------------------------------------------------------


# ============================================================================
#  Test initiating the Apex class
# ============================================================================


def test_init_defaults():
    Apex()


def test_init_date_int():
    apex_out = Apex(date=2015)
    assert apex_out.year == 2015


def test_init_date_float():
    apex_out = Apex(date=2015.5)
    assert apex_out.year == 2015.5


def test_init_date():
    date = dt.date(2015, 1, 1)
    apex_out = Apex(date=date)
    assert apex_out.year == helpers.toYearFraction(date)


def test_init_datetime():
    datetime = dt.datetime(2015, 6, 1, 18, 23, 45)
    apex_out = Apex(date=datetime)
    assert apex_out.year == helpers.toYearFraction(datetime)


def test_init_datafile_IOError():
    with pytest.raises(IOError):
        Apex(date=2015, datafile='foo/path/to/datafile.blah')


# ============================================================================
#  Test the low-level interfaces to the fortran wrappers
# ============================================================================

def test__geo2qd_scalar():
    apex_out = Apex(date=2000, refh=300)
    for lat in [0, 30, 60, 89]:
        for lon in [-179, -90, 0, 90, 180]:
            assert_allclose(apex_out._geo2qd(lat, lon, 100),
                            fa.apxg2q(lat, lon, 100, 0)[:2])


def test__geo2qd_array():
    apex_out = Apex(date=2000, refh=300)
    lats, lons = apex_out._geo2qd([[0, 30], [60, 90]], 15,
                                  [[100, 200], [300, 400]])
    lat1, lon1 = fa.apxg2q(0, 15, 100, 0)[:2]
    lat2, lon2 = fa.apxg2q(30, 15, 200, 0)[:2]
    lat3, lon3 = fa.apxg2q(60, 15, 300, 0)[:2]
    lat4, lon4 = fa.apxg2q(90, 15, 400, 0)[:2]
    assert_allclose(lats.astype(float), np.array([[lat1, lat2], [lat3, lat4]],
                                                 dtype=float))
    assert_allclose(lons.astype(float), np.array([[lon1, lon2], [lon3, lon4]],
                                                 dtype=float))


def test__geo2qd_longitude():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out._geo2qd(60, 180, 100),
                    fa.apxg2q(60, 180, 100, 0)[:2])
    assert_allclose(apex_out._geo2qd(60, -180, 100),
                    fa.apxg2q(60, -180, 100, 0)[:2])
    assert_allclose(apex_out._geo2qd(60, -180, 100),
                    apex_out._geo2qd(60, 180, 100))
    for i in range(-5, 5):
        for lat in [0, 30, 60, 90]:
            assert_allclose(apex_out._geo2qd(lat, 15 + i * 360, 100),
                            fa.apxg2q(lat, 15, 100, 0)[:2])


def test__geo2apex_scalar():
    apex_out = Apex(date=2000, refh=300)
    for lat in [0, 30, 60, 89]:
        for lon in [-179, -90, 0, 90, 180]:
            assert_allclose(apex_out._geo2apex(lat, lon, 100),
                            fa.apxg2all(lat, lon, 100, 300, 0)[2:4])


def test__geo2apex_array():
    apex_out = Apex(date=2000, refh=300)
    lats, lons = apex_out._geo2apex([[0, 30], [60, 90]], 15,
                                    [[100, 200], [300, 400]])
    lat1, lon1 = fa.apxg2all(0, 15, 100, 300, 0)[2:4]
    lat2, lon2 = fa.apxg2all(30, 15, 200, 300, 0)[2:4]
    lat3, lon3 = fa.apxg2all(60, 15, 300, 300, 0)[2:4]
    lat4, lon4 = fa.apxg2all(90, 15, 400, 300, 0)[2:4]
    assert_allclose(lats.astype(float), np.array([[lat1, lat2], [lat3, lat4]],
                                                 dtype=float))
    assert_allclose(lons.astype(float), np.array([[lon1, lon2], [lon3, lon4]],
                                                 dtype=float))


def test__geo2apex_longitude():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out._geo2apex(60, 180, 100),
                    fa.apxg2all(60, 180, 100, 300, 0)[2:4])
    assert_allclose(apex_out._geo2apex(60, -180, 100),
                    fa.apxg2all(60, -180, 100, 300, 0)[2:4])
    assert_allclose(apex_out._geo2apex(60, -180, 100),
                    apex_out._geo2apex(60, 180, 100))
    for i in range(-5, 5):
        for lat in [0, 30, 60, 90]:
            assert_allclose(apex_out._geo2apex(lat, 15 + i * 360, 100),
                            fa.apxg2all(lat, 15, 100, 300, 0)[2:4])


def test__geo2apexall_scalar():
    apex_out = Apex(date=2000, refh=300)
    for lat in [0, 30, 60, 89]:
        for lon in [-179, -90, 0, 90, 180]:
            ret1 = apex_out._geo2apexall(lat, lon, 100)
            ret2 = fa.apxg2all(lat, lon, 100, 300, 1)
            for r1, r2 in zip(ret1, ret2):
                assert_allclose(r1, r2)


def test__geo2apexall_array():
    apex_out = Apex(date=2000, refh=300)
    ret = apex_out._geo2apexall([[0, 30], [60, 90]], 15,
                                [[100, 200], [300, 400]])
    ret1 = fa.apxg2all(0, 15, 100, 300, 1)
    ret2 = fa.apxg2all(30, 15, 200, 300, 1)
    ret3 = fa.apxg2all(60, 15, 300, 300, 1)
    ret4 = fa.apxg2all(90, 15, 400, 300, 1)
    for i in range(len(ret)):
        try:
            # ret[i] is array of floats
            assert_allclose(ret[i].astype(float),
                            np.array([[ret1[i], ret2[i]], [ret3[i], ret4[i]]],
                                     dtype=float))
        except ValueError:
            # ret[i] is array of arrays
            assert_allclose(ret[i][0, 0], ret1[i])
            assert_allclose(ret[i][0, 1], ret2[i])
            assert_allclose(ret[i][1, 0], ret3[i])
            assert_allclose(ret[i][1, 1], ret4[i])


def test__qd2geo_scalar():
    apex_out = Apex(date=2000, refh=300)
    for lat in [0, 30, 60, 89]:
        for lon in [-179, -90, 0, 90, 180]:
            for prec in [-1, 1e-2, 1e-10]:
                assert_allclose(apex_out._qd2geo(lat, lon, 100, prec),
                                fa.apxq2g(lat, lon, 100, prec))


def test__qd2geo_array():
    apex_out = Apex(date=2000, refh=300)
    lats, lons, errs = apex_out._qd2geo([[0, 30], [60, 90]], 15,
                                        [[100, 200], [300, 400]], 1e-2)
    lat1, lon1, err1 = fa.apxq2g(0, 15, 100, 1e-2)
    lat2, lon2, err2 = fa.apxq2g(30, 15, 200, 1e-2)
    lat3, lon3, err3 = fa.apxq2g(60, 15, 300, 1e-2)
    lat4, lon4, err4 = fa.apxq2g(90, 15, 400, 1e-2)
    assert_allclose(lats.astype(float), np.array([[lat1, lat2], [lat3, lat4]],
                                                 dtype=float))
    assert_allclose(lons.astype(float), np.array([[lon1, lon2], [lon3, lon4]],
                                                 dtype=float))
    assert_allclose(errs.astype(float), np.array([[err1, err2], [err3, err4]],
                                                 dtype=float))


def test__qd2geo_longitude():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out._qd2geo(60, 180, 100, 1e-2),
                    fa.apxq2g(60, 180, 100, 1e-2))
    assert_allclose(apex_out._qd2geo(60, -180, 100, 1e-2),
                    fa.apxq2g(60, -180, 100, 1e-2))
    assert_allclose(apex_out._qd2geo(60, -180, 100, 1e-2),
                    apex_out._qd2geo(60, 180, 100, 1e-2))
    for i in range(-5, 5):
        for lat in [0, 30, 60, 90]:
            assert_allclose(apex_out._qd2geo(lat, 15 + i * 360, 100, 1e-2),
                            fa.apxq2g(lat, 15, 100, 1e-2))


def test__basevec_scalar():
    apex_out = Apex(date=2000, refh=300)
    for lat in [0, 30, 60, 89]:
        for lon in [-179, -90, 0, 90, 180]:
            assert_allclose(apex_out._basevec(lat, lon, 100),
                            fa.apxg2q(lat, lon, 100, 1)[2:4])


def test__basevec_array():
    apex_out = Apex(date=2000, refh=300)
    f1s, f2s = apex_out._basevec([[0, 30], [60, 90]], 15,
                                 [[100, 200], [300, 400]])
    f11, f21 = fa.apxg2q(0, 15, 100, 1)[2:4]
    f12, f22 = fa.apxg2q(30, 15, 200, 1)[2:4]
    f13, f23 = fa.apxg2q(60, 15, 300, 1)[2:4]
    f14, f24 = fa.apxg2q(90, 15, 400, 1)[2:4]
    assert_allclose(f1s[0, 0], f11)
    assert_allclose(f1s[0, 1], f12)
    assert_allclose(f1s[1, 0], f13)
    assert_allclose(f1s[1, 1], f14)
    assert_allclose(f2s[0, 0], f21)
    assert_allclose(f2s[0, 1], f22)
    assert_allclose(f2s[1, 0], f23)
    assert_allclose(f2s[1, 1], f24)


def test__basevec_longitude():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out._basevec(60, 180, 100),
                    fa.apxg2q(60, 180, 100, 1)[2:4])
    assert_allclose(apex_out._basevec(60, -180, 100),
                    fa.apxg2q(60, -180, 100, 1)[2:4])
    assert_allclose(apex_out._basevec(60, -180, 100),
                    apex_out._basevec(60, 180, 100))
    for i in range(-5, 5):
        for lat in [0, 30, 60, 90]:
            assert_allclose(apex_out._basevec(lat, 15 + i * 360, 100),
                            fa.apxg2q(lat, 15, 100, 1)[2:4])


# ============================================================================
#  Test the convert() method
# ============================================================================


def test_convert_geo2apex():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.convert(60, 15, 'geo', 'apex', height=100),
                    apex_out.geo2apex(60, 15, 100))


def test_convert_geo2qd():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.convert(60, 15, 'geo', 'qd', height=100),
                    apex_out.geo2qd(60, 15, 100))


def test_convert_geo2mlt_nodate():
    apex_out = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        apex_out.convert(60, 15, 'geo', 'mlt')


def test_convert_geo2mlt():
    datetime = dt.datetime(2000, 3, 9, 14, 25, 58)
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.convert(60, 15, 'geo', 'mlt', height=100,
                                     ssheight=2e5, datetime=datetime)[1],
                    apex_out.mlon2mlt(apex_out.geo2apex(60, 15, 100)[1],
                                      datetime, ssheight=2e5))


def test_convert_apex2geo():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.convert(60, 15, 'apex', 'geo', height=100,
                                     precision=1e-2),
                    apex_out.apex2geo(60, 15, 100, precision=1e-2)[:-1])


def test_convert_apex2qd():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.convert(60, 15, 'apex', 'qd', height=100),
                    apex_out.apex2qd(60, 15, height=100))


def test_convert_apex2mlt():
    datetime = dt.datetime(2000, 3, 9, 14, 25, 58)
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.convert(60, 15, 'apex', 'mlt', height=100,
                                     datetime=datetime, ssheight=2e5)[1],
                    apex_out.mlon2mlt(15, datetime, ssheight=2e5))


def test_convert_qd2geo():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.convert(60, 15, 'qd', 'geo', height=100,
                                     precision=1e-2),
                    apex_out.qd2geo(60, 15, 100, precision=1e-2)[:-1])


def test_convert_qd2apex():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.convert(60, 15, 'qd', 'apex', height=100),
                    apex_out.qd2apex(60, 15, height=100))


def test_convert_qd2apex_at_equator():
    """Test the quasi-dipole to apex conversion at the magnetic equator."""
    apex_out = Apex(date=2000, refh=80)
    elat, elon = apex_out.convert(lat=0.0, lon=0, source='qd', dest='apex',
                                  height=120.0)
    clat, clon = apex_out.convert(lat=0.001, lon=0, source='qd', dest='apex',
                                  height=120.0)
    assert_allclose([elat, elon], [clat, clon], atol=1e-4)


def test_convert_qd2mlt():
    datetime = dt.datetime(2000, 3, 9, 14, 25, 58)
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.convert(60, 15, 'qd', 'mlt', height=100,
                                     datetime=datetime, ssheight=2e5)[1],
                    apex_out.mlon2mlt(15, datetime, ssheight=2e5))


def test_convert_mlt2geo():
    datetime = dt.datetime(2000, 3, 9, 14, 25, 58)
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.convert(60, 15, 'mlt', 'geo', height=100,
                                     datetime=datetime, precision=1e-2,
                                     ssheight=2e5),
                    apex_out.apex2geo(60, apex_out.mlt2mlon(15, datetime,
                                                            ssheight=2e5), 100,
                                      precision=1e-2)[:-1])


def test_convert_mlt2apex():
    datetime = dt.datetime(2000, 3, 9, 14, 25, 58)
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.convert(60, 15, 'mlt', 'apex', height=100,
                                     datetime=datetime, ssheight=2e5),
                    (60, apex_out.mlt2mlon(15, datetime, ssheight=2e5)))


def test_convert_mlt2qd():
    datetime = dt.datetime(2000, 3, 9, 14, 25, 58)
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.convert(60, 15, 'mlt', 'qd', height=100,
                                     datetime=datetime, ssheight=2e5),
                    apex_out.apex2qd(60, apex_out.mlt2mlon(15, datetime,
                                                           ssheight=2e5),
                                     height=100))


def test_convert_invalid_lat():
    apex_out = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        apex_out.convert(91, 0, 'geo', 'geo')
    with pytest.raises(ValueError):
        apex_out.convert(-91, 0, 'geo', 'geo')
    apex_out.convert(90, 0, 'geo', 'geo')
    apex_out.convert(-90, 0, 'geo', 'geo')

    assert_allclose(apex_out.convert(90 + 1e-5, 0, 'geo', 'apex'),
                    apex_out.convert(90, 0, 'geo', 'apex'), rtol=0, atol=1e-8)


def test_convert_invalid_transformation():
    apex_out = Apex(date=2000, refh=300)
    with pytest.raises(NotImplementedError):
        apex_out.convert(0, 0, 'foobar', 'geo')
    with pytest.raises(NotImplementedError):
        apex_out.convert(0, 0, 'geo', 'foobar')


coord_names = ['geo', 'apex', 'qd']


@pytest.mark.parametrize('transform', itertools.product(coord_names,
                                                        coord_names))
def test_convert_withnan(transform):
    """Test Apex.convert success with NaN input."""
    num_nans = 5
    in_lat = np.arange(0, 10, dtype=float)
    in_lat[:num_nans] = np.nan
    in_lon = np.arange(0, 10, dtype=float)
    in_lon[:num_nans] = np.nan
    src, dest = transform
    apex_out = Apex(date=2000, refh=80)
    out_lat, out_lon = apex_out.convert(in_lat, in_lon, src, dest, height=120)
    assert np.all(np.isnan(out_lat[:num_nans]))
    assert np.all(np.isnan(out_lon[:num_nans]))
    assert np.all(np.isfinite(out_lat[num_nans:]))
    assert np.all(np.isfinite(out_lat[num_nans:]))


# ============================================================================
#  Test the geo2apex() method
# ============================================================================


def test_geo2apex():
    apex_out = Apex(date=2000, refh=300)
    lat, lon = apex_out.geo2apex(60, 15, 100)
    assert_allclose((lat, lon), apex_out._geo2apex(60, 15, 100))
    assert type(lat) != np.ndarray
    assert type(lon) != np.ndarray


def test_geo2apex_vectorization():
    apex_out = Apex(date=2000, refh=300)
    assert apex_out.geo2apex([60, 60], 15, 100)[0].shape == (2,)
    assert apex_out.geo2apex(60, [15, 15], 100)[0].shape == (2,)
    assert apex_out.geo2apex(60, 15, [100, 100])[0].shape == (2,)


def test_geo2apex_invalid_lat():
    apex_out = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        apex_out.geo2apex(91, 0, 0)
    with pytest.raises(ValueError):
        apex_out.geo2apex(-91, 0, 0)
    apex_out.geo2apex(90, 0, 0)
    apex_out.geo2apex(-90, 0, 0)

    assert_allclose(apex_out.geo2apex(90 + 1e-5, 0, 0),
                    apex_out.geo2apex(90, 0, 0), rtol=0, atol=1e-8)


def test_geo2apex_undefined_warning(recwarn):
    """Test warning and fill values for an undefined location."""
    apex_out = Apex(date=2000, refh=10000)
    ret = apex_out.geo2apex(0, 0, 0)

    assert np.isnan(ret[0])
    assert len(recwarn) == 1
    assert issubclass(recwarn[-1].category, UserWarning)
    assert 'set to NaN where' in str(recwarn[-1].message)


# ============================================================================
#  Test the apex2geo() method
# ============================================================================


def test_apex2geo():
    apex_out = Apex(date=2000, refh=300)
    lat, lon, error = apex_out.apex2geo(60, 15, 100, precision=1e-2)
    assert_allclose((lat, lon, error),
                    apex_out.qd2geo(*apex_out.apex2qd(60, 15, 100), height=100,
                                    precision=1e-2))

    assert type(lat) != np.ndarray
    assert type(lon) != np.ndarray
    assert type(error) != np.ndarray


def test_apex2geo_vectorization():
    apex_out = Apex(date=2000, refh=300)
    assert apex_out.apex2geo([60, 60], 15, 100)[0].shape == (2,)
    assert apex_out.apex2geo(60, [15, 15], 100)[0].shape == (2,)
    assert apex_out.apex2geo(60, 15, [100, 100])[0].shape == (2,)


def test_apex2geo_invalid_lat():
    apex_out = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        apex_out.apex2geo(91, 0, 0, 1e-2)
    with pytest.raises(ValueError):
        apex_out.apex2geo(-91, 0, 0, 1e-2)
    apex_out.apex2geo(90, 0, 0, 1e-2)
    apex_out.apex2geo(-90, 0, 0, 1e-2)

    assert_allclose(apex_out.apex2geo(90 + 1e-5, 0, 0, 1e-2),
                    apex_out.apex2geo(90, 0, 0, 1e-2), rtol=0, atol=1e-8)


# ============================================================================
#  Test the geo2qd() method
# ============================================================================


def test_geo2qd():
    apex_out = Apex(date=2000, refh=300)
    lat, lon = apex_out.geo2qd(60, 15, 100)
    assert_allclose((lat, lon), apex_out._geo2qd(60, 15, 100))
    assert type(lat) != np.ndarray
    assert type(lon) != np.ndarray


def test_geo2qd_vectorization():
    apex_out = Apex(date=2000, refh=300)
    assert apex_out.geo2qd([60, 60], 15, 100)[0].shape == (2,)
    assert apex_out.geo2qd(60, [15, 15], 100)[0].shape == (2,)
    assert apex_out.geo2qd(60, 15, [100, 100])[0].shape == (2,)


def test_geo2qd_invalid_lat():
    apex_out = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        apex_out.geo2qd(91, 0, 0)
    with pytest.raises(ValueError):
        apex_out.geo2qd(-91, 0, 0)
    apex_out.geo2qd(90, 0, 0)
    apex_out.geo2qd(-90, 0, 0)

    assert_allclose(apex_out.geo2qd(90 + 1e-5, 0, 0), apex_out.geo2qd(90, 0, 0),
                    rtol=0, atol=1e-8)


# ============================================================================
#  Test the qd2geo() method
# ============================================================================


def test_qd2geo():
    apex_out = Apex(date=2000, refh=300)
    lat, lon, error = apex_out.qd2geo(60, 15, 100, precision=1e-2)
    assert_allclose((lat, lon, error), apex_out._qd2geo(60, 15, 100, 1e-2))
    assert type(lat) != np.ndarray
    assert type(lon) != np.ndarray
    assert type(error) != np.ndarray


def test_qd2geo_vectorization():
    apex_out = Apex(date=2000, refh=300)
    assert apex_out.qd2geo([60, 60], 15, 100)[0].shape == (2,)
    assert apex_out.qd2geo(60, [15, 15], 100)[0].shape == (2,)
    assert apex_out.qd2geo(60, 15, [100, 100])[0].shape == (2,)


def test_qd2geo_invalid_lat():
    apex_out = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        apex_out.qd2geo(91, 0, 0, precision=1e-2)
    with pytest.raises(ValueError):
        apex_out.qd2geo(-91, 0, 0, precision=1e-2)
    apex_out.qd2geo(90, 0, 0, precision=1e-2)
    apex_out.qd2geo(-90, 0, 0, precision=1e-2)

    assert_allclose(apex_out.qd2geo(90 + 1e-5, 0, 0, 1e-2),
                    apex_out.qd2geo(90, 0, 0, 1e-2), rtol=0, atol=1e-8)


# ============================================================================
#  Test the apex2qd() method
# ============================================================================


def test_apex2qd():
    apex_out = Apex(date=2000, refh=300)
    lat, lon = apex_out.apex2qd(60, 15, 100)
    assert_allclose((lat, lon),
                    [60.498401, 15])
    assert type(lat) != np.ndarray
    assert type(lon) != np.ndarray


def test_apex2qd_vectorization():
    apex_out = Apex(date=2000, refh=300)
    assert apex_out.apex2qd([60, 60], 15, 100)[0].shape == (2,)
    assert apex_out.apex2qd(60, [15, 15], 100)[0].shape == (2,)
    assert apex_out.apex2qd(60, 15, [100, 100])[0].shape == (2,)


def test_apex2qd_invalid_lat():
    apex_out = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        apex_out.apex2qd(91, 0, 0)
    with pytest.raises(ValueError):
        apex_out.apex2qd(-91, 0, 0)
    apex_out.apex2qd(90, 0, 0)
    apex_out.apex2qd(-90, 0, 0)

    assert_allclose(apex_out.apex2qd(90 + 1e-5, 0, 0),
                    apex_out.apex2qd(90, 0, 0), rtol=0, atol=1e-8)


def test_apex2qd_apexheight_close():
    apex_out = Apex(date=2000, refh=300)
    apex_out.apex2qd(0, 15, 300 + 1e-6)


def test_apex2qd_apexheight_over():
    apex_out = Apex(date=2000, refh=300)
    with pytest.raises(ApexHeightError):
        apex_out.apex2qd(0, 15, 301)


# ============================================================================
#  Test the qd2apex() method
# ============================================================================


def test_qd2apex():
    apex_out = Apex(date=2000, refh=300)
    lat, lon = apex_out.qd2apex(60, 15, 100)
    assert_allclose((lat, lon),
                    [59.491381, 15])
    assert type(lat) != np.ndarray
    assert type(lon) != np.ndarray


def test_qd2apex_vectorization():
    apex_out = Apex(date=2000, refh=300)
    assert apex_out.qd2apex([60, 60], 15, 100)[0].shape == (2,)
    assert apex_out.qd2apex(60, [15, 15], 100)[0].shape == (2,)
    assert apex_out.qd2apex(60, 15, [100, 100])[0].shape == (2,)


def test_qd2apex_invalid_lat():
    apex_out = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        apex_out.qd2apex(91, 0, 0)
    with pytest.raises(ValueError):
        apex_out.qd2apex(-91, 0, 0)
    apex_out.qd2apex(90, 0, 0)
    apex_out.qd2apex(-90, 0, 0)

    assert_allclose(apex_out.qd2apex(90 + 1e-5, 0, 0),
                    apex_out.qd2apex(90, 0, 0), rtol=0, atol=1e-8)


def test_qd2apex_apexheight_close():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.qd2apex(0, 15, 300 - 1e-5),
                    apex_out.qd2apex(0, 15, 300))


def test_qd2apex_apexheight_over():
    apex_out = Apex(date=2000, refh=300)
    with pytest.raises(ApexHeightError):
        apex_out.qd2apex(0, 15, 299)


# ============================================================================
#  Test mlon2mlt()
# ============================================================================


def test_mlon2mlt_scalar():
    apex_out = Apex(date=2000, refh=300)
    mlon = apex_out.mlon2mlt(0, dt.datetime(2000, 2, 3, 4, 5, 6))
    assert_allclose(mlon, 23.019629923502603)
    assert type(mlon) != np.ndarray


def test_mlon2mlt_ssheight():
    apex_out = Apex(date=2000, refh=300)
    mlt = apex_out.mlon2mlt(0, dt.datetime(2000, 2, 3, 4, 5, 6),
                            ssheight=50 * 2000)
    assert_allclose(mlt, 23.026712036132814)


def test_mlon2mlt_1Darray():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.mlon2mlt([0, 180],
                                      dt.datetime(2000, 2, 3, 4, 5, 6)),
                    [23.019261, 11.019261], rtol=1e-4)


def test_mlon2mlt_2Darray():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.mlon2mlt([[0, 180], [0, 180]],
                                      dt.datetime(2000, 2, 3, 4, 5, 6)),
                    [[23.019261, 11.019261], [23.019261, 11.019261]], rtol=1e-4)


def test_mlon2mlt_diffdates():
    apex_out = Apex(date=2000, refh=300)
    dtime1 = dt.datetime(2000, 2, 3, 4, 5, 6)
    dtime2 = dt.datetime(2000, 2, 3, 5, 5, 6)
    assert apex_out.mlon2mlt(0, dtime1) != apex_out.mlon2mlt(0, dtime2)


def test_mlon2mlt_offset():
    apex_out = Apex(date=2000, refh=300)
    date = dt.datetime(2000, 2, 3, 4, 5, 6)
    assert_allclose(apex_out.mlon2mlt(0, date),
                    apex_out.mlon2mlt(-15, date) + 1)
    assert_allclose(apex_out.mlon2mlt(0, date),
                    apex_out.mlon2mlt(-10 * 15, date) + 10)


def test_mlon2mlt_range():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.mlon2mlt(range(0, 361, 30),
                                      dt.datetime(2000, 2, 3, 4, 5, 6)),
                    [23.01963, 1.01963, 3.01963, 5.01963, 7.01963,
                     9.01963, 11.01963, 13.01963, 15.01963, 17.01963,
                     19.01963, 21.01963, 23.01963],
                    rtol=1e-4)


# ============================================================================
#  Test mlt2mlon()
# ============================================================================


def test_mlt2mlon_scalar():
    apex_out = Apex(date=2000, refh=300)
    mlt = apex_out.mlt2mlon(0, dt.datetime(2000, 2, 3, 4, 5, 6))
    assert_allclose(mlt, 14.705535888671875)
    assert type(mlt) != np.ndarray


def test_mlt2mlon_ssheight():
    apex_out = Apex(date=2000, refh=300)
    mlt = apex_out.mlt2mlon(0, dt.datetime(2000, 2, 3, 4, 5, 6),
                            ssheight=50 * 2000)
    assert_allclose(mlt, 14.599319458007812)


def test_mlt2mlon_1Darray():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.mlt2mlon([0, 12],
                                      dt.datetime(2000, 2, 3, 4, 5, 6)),
                    [14.705551, 194.705551], rtol=1e-4)


def test_mlt2mlon_2Darray():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.mlt2mlon([[0, 12], [0, 12]],
                                      dt.datetime(2000, 2, 3, 4, 5, 6)),
                    [[14.705551, 194.705551], [14.705551, 194.705551]],
                    rtol=1e-4)


def test_mlt2mlon_diffdates():
    apex_out = Apex(date=2000, refh=300)
    dtime1 = dt.datetime(2000, 2, 3, 4, 5, 6)
    dtime2 = dt.datetime(2000, 2, 3, 5, 5, 6)
    assert apex_out.mlt2mlon(0, dtime1) != apex_out.mlt2mlon(0, dtime2)


def test_mlt2mlon_offset():
    apex_out = Apex(date=2000, refh=300)
    date = dt.datetime(2000, 2, 3, 4, 5, 6)
    assert_allclose(apex_out.mlt2mlon(0, date), apex_out.mlt2mlon(1, date) - 15)
    assert_allclose(apex_out.mlt2mlon(0, date),
                    apex_out.mlt2mlon(10, date) - 150)


def test_mlt2mlon_range():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.mlt2mlon(range(0, 25, 2),
                                      dt.datetime(2000, 2, 3, 4, 5, 6)),
                    [14.705551, 44.705551, 74.705551, 104.705551, 134.705551,
                     164.705551, 194.705551, 224.705551, 254.705551, 284.705551,
                     314.705551, 344.705551, 14.705551],
                    rtol=1e-4)


# ============================================================================
#  Test mlt/mlon back and forth
# ============================================================================


def test_mlon2mlt2mlon():
    apex_out = Apex(date=2000, refh=300)
    date = dt.datetime(2000, 2, 3, 4, 5, 6)
    assert_allclose(apex_out.mlon2mlt(apex_out.mlt2mlon(0, date), date), 0)
    assert_allclose(apex_out.mlon2mlt(apex_out.mlt2mlon(6, date), date), 6)
    assert_allclose(apex_out.mlon2mlt(apex_out.mlt2mlon(12, date), date), 12)
    assert_allclose(apex_out.mlon2mlt(apex_out.mlt2mlon(18, date), date), 18)
    assert_allclose(apex_out.mlon2mlt(apex_out.mlt2mlon(24, date), date), 0)


def test_mlt2mlon2mlt():
    apex_out = Apex(date=2000, refh=300)
    date = dt.datetime(2000, 2, 3, 4, 5, 6)
    assert_allclose(apex_out.mlt2mlon(apex_out.mlon2mlt(0, date), date), 0)
    assert_allclose(apex_out.mlt2mlon(apex_out.mlon2mlt(90, date), date), 90)
    assert_allclose(apex_out.mlt2mlon(apex_out.mlon2mlt(180, date), date), 180)
    assert_allclose(apex_out.mlt2mlon(apex_out.mlon2mlt(270, date), date), 270)
    assert_allclose(apex_out.mlt2mlon(apex_out.mlon2mlt(360, date), date), 0)


# ============================================================================
#  Test the map_to_height() method
# ============================================================================


def test_map_to_height():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.map_to_height(60, 15, 100, 10000, conjugate=False,
                                           precision=1e-10),
                    (31.841466903686523, 17.916635513305664,
                     1.7075473124350538e-6))
    assert_allclose(apex_out.map_to_height(30, 170, 100, 500, conjugate=False,
                                           precision=1e-2),
                    (25.727270126342773, 169.60546875, 0.00017573432705830783))


def test_map_to_height_same_height():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.map_to_height(60, 15, 100, 100, conjugate=False,
                                           precision=1e-10),
                    (60.0, 15.000003814697266, 0.0), rtol=1e-5)


def test_map_to_height_conjugate():
    """Test results of map_to_height using conjugacy."""
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.map_to_height(60, 15, 100, 10000, conjugate=True,
                                           precision=1e-10),
                    (-25.424888610839844, 27.310426712036133,
                     1.2074182222931995e-6), atol=1e-6)
    assert_allclose(apex_out.map_to_height(30, 170, 100, 500, conjugate=True,
                                           precision=1e-2),
                    (-13.76642894744873, 164.24259948730469,
                     0.00056820799363777041), atol=1e-6)


def test_map_to_height_vectorization():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.map_to_height([60, 60], 15, 100, 100),
                    ([60] * 2, [15.00000381] * 2, [0] * 2), rtol=1e-5)
    assert_allclose(apex_out.map_to_height(60, [15, 15], 100, 100),
                    ([60] * 2, [15.00000381] * 2, [0] * 2), rtol=1e-5)
    assert_allclose(apex_out.map_to_height(60, 15, [100, 100], 100),
                    ([60] * 2, [15.00000381] * 2, [0] * 2), rtol=1e-5)
    assert_allclose(apex_out.map_to_height(60, 15, 100, [100, 100]),
                    ([60] * 2, [15.00000381] * 2, [0] * 2), rtol=1e-5)


def test_map_to_height_ApexHeightError():
    apex_out = Apex(date=2000, refh=300)
    with pytest.raises(ApexHeightError):
        apex_out.map_to_height(0, 15, 100, 10000)


# ============================================================================
#  Test the map_E_to_height() method
# ============================================================================


def test_map_E_to_height():
    apex_out = Apex(date=2000, refh=300)
    out_60_15_100_500 = [0.71152183, 2.35624876, 0.57260784]
    out_60_15_100_500_234 = [1.56028502, 3.43916636, 0.78235384]
    out_60_15_100_1000 = [0.67796492, 2.08982134, 0.55860785]
    out_60_15_200_500 = [0.72377397, 2.42737471, 0.59083726]
    out_60_30_100_500 = [0.68626344, 2.37530133, 0.60060124]
    out_70_15_100_500 = [0.72760378, 2.18082305, 0.29141979]

    # scalar
    assert_allclose(apex_out.map_E_to_height(60, 15, 100, 500, [1, 2, 3]),
                    out_60_15_100_500, rtol=1e-5)
    assert_allclose(apex_out.map_E_to_height(60, 15, 100, 500, [2, 3, 4]),
                    out_60_15_100_500_234, rtol=1e-5)
    assert_allclose(apex_out.map_E_to_height(60, 15, 100, 1000, [1, 2, 3]),
                    out_60_15_100_1000, rtol=1e-5)
    assert_allclose(apex_out.map_E_to_height(60, 15, 200, 500, [1, 2, 3]),
                    out_60_15_200_500, rtol=1e-5)
    assert_allclose(apex_out.map_E_to_height(60, 30, 100, 500, [1, 2, 3]),
                    out_60_30_100_500, rtol=1e-5)
    assert_allclose(apex_out.map_E_to_height(70, 15, 100, 500, [1, 2, 3]),
                    out_70_15_100_500, rtol=1e-5)

    # vectorize lat
    assert_allclose(apex_out.map_E_to_height([60, 70], 15, 100, 500,
                                             np.array([[1, 2, 3]] * 2).T),
                    np.array([out_60_15_100_500, out_70_15_100_500]).T,
                    rtol=1e-5)

    # vectorize lon
    assert_allclose(apex_out.map_E_to_height(60, [15, 30], 100, 500,
                                             np.array([[1, 2, 3]] * 2).T),
                    np.array([out_60_15_100_500, out_60_30_100_500]).T,
                    rtol=1e-5)

    # vectorize height
    assert_allclose(apex_out.map_E_to_height(60, 15, [100, 200], 500,
                                             np.array([[1, 2, 3]] * 2).T),
                    np.array([out_60_15_100_500, out_60_15_200_500]).T,
                    rtol=1e-5)

    # vectorize newheight
    assert_allclose(apex_out.map_E_to_height(60, 15, 100, [500, 1000],
                                             np.array([[1, 2, 3]] * 2).T),
                    np.array([out_60_15_100_500, out_60_15_100_1000]).T,
                    rtol=1e-5)

    # vectorize E
    assert_allclose(apex_out.map_E_to_height(60, 15, 100, 500,
                                             np.array([[1, 2, 3],
                                                       [2, 3, 4]]).T),
                    np.array([out_60_15_100_500, out_60_15_100_500_234]).T,
                    rtol=1e-5)


# ============================================================================
#  Test the map_V_to_height() method
# ============================================================================


def test_map_V_to_height():
    apex_out = Apex(date=2000, refh=300)
    out_60_15_100_500 = [0.81971957, 2.84512495, 0.69545001]
    out_60_15_100_500_234 = [1.83027746, 4.14346436, 0.94764179]
    out_60_15_100_1000 = [0.92457698, 3.14997661, 0.85135187]
    out_60_15_200_500 = [0.80388262, 2.79321504, 0.68285158]
    out_60_30_100_500 = [0.76141245, 2.87884673, 0.73655941]
    out_70_15_100_500 = [0.84681866, 2.5925821,  0.34792655]

    # scalar
    assert_allclose(apex_out.map_V_to_height(60, 15, 100, 500, [1, 2, 3]),
                    out_60_15_100_500, rtol=1e-5)
    assert_allclose(apex_out.map_V_to_height(60, 15, 100, 500, [2, 3, 4]),
                    out_60_15_100_500_234, rtol=1e-5)
    assert_allclose(apex_out.map_V_to_height(60, 15, 100, 1000, [1, 2, 3]),
                    out_60_15_100_1000, rtol=1e-5)
    assert_allclose(apex_out.map_V_to_height(60, 15, 200, 500, [1, 2, 3]),
                    out_60_15_200_500, rtol=1e-5)
    assert_allclose(apex_out.map_V_to_height(60, 30, 100, 500, [1, 2, 3]),
                    out_60_30_100_500, rtol=1e-5)
    assert_allclose(apex_out.map_V_to_height(70, 15, 100, 500, [1, 2, 3]),
                    out_70_15_100_500, rtol=1e-5)

    # vectorize lat
    assert_allclose(apex_out.map_V_to_height([60, 70], 15, 100, 500,
                                             np.array([[1, 2, 3]] * 2).T),
                    np.array([out_60_15_100_500, out_70_15_100_500]).T,
                    rtol=1e-5)

    # vectorize lon
    assert_allclose(apex_out.map_V_to_height(60, [15, 30], 100, 500,
                                             np.array([[1, 2, 3]] * 2).T),
                    np.array([out_60_15_100_500, out_60_30_100_500]).T,
                    rtol=1e-5)

    # vectorize height
    assert_allclose(apex_out.map_V_to_height(60, 15, [100, 200], 500,
                                             np.array([[1, 2, 3]] * 2).T),
                    np.array([out_60_15_100_500, out_60_15_200_500]).T,
                    rtol=1e-5)

    # vectorize newheight
    assert_allclose(apex_out.map_V_to_height(60, 15, 100, [500, 1000],
                                             np.array([[1, 2, 3]] * 2).T),
                    np.array([out_60_15_100_500, out_60_15_100_1000]).T,
                    rtol=1e-5)

    # vectorize E
    assert_allclose(apex_out.map_V_to_height(60, 15, 100, 500,
                                             np.array([[1, 2, 3],
                                                       [2, 3, 4]]).T),
                    np.array([out_60_15_100_500, out_60_15_100_500_234]).T,
                    rtol=1e-5)


# ============================================================================
#  Test basevectors_qd()
# ============================================================================


# test coords

def test_basevectors_qd_scalar_geo():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.basevectors_qd(60, 15, 100, coords='geo'),
                    apex_out._basevec(60, 15, 100))


def test_basevectors_qd_scalar_apex():
    apex_out = Apex(date=2000, refh=300)
    glat, glon, _ = apex_out.apex2geo(60, 15, 100, precision=1e-2)
    assert_allclose(apex_out.basevectors_qd(60, 15, 100, coords='apex',
                                            precision=1e-2),
                    apex_out._basevec(glat, glon, 100))


def test_basevectors_qd_scalar_qd():
    apex_out = Apex(date=2000, refh=300)
    glat, glon, _ = apex_out.qd2geo(60, 15, 100, precision=1e-2)
    assert_allclose(apex_out.basevectors_qd(60, 15, 100, coords='qd',
                                            precision=1e-2),
                    apex_out._basevec(glat, glon, 100))


# test shapes and vectorization of arguments
def test_basevectors_qd_scalar_shape():
    apex_out = Apex(date=2000, refh=300)
    ret = apex_out.basevectors_qd(60, 15, 100)
    for r in ret:
        assert r.shape == (2,)


def test_basevectors_qd_vectorization():
    apex_out = Apex(date=2000, refh=300)
    ret = apex_out.basevectors_qd([60, 60, 60, 60], 15, 100, coords='geo')
    for r in ret:
        assert r.shape == (2, 4)
    ret = apex_out.basevectors_qd(60, [15, 15, 15, 15], 100, coords='geo')
    for r in ret:
        assert r.shape == (2, 4)
    ret = apex_out.basevectors_qd(60, 15, [100, 100, 100, 100], coords='geo')
    for r in ret:
        assert r.shape == (2, 4)


# test array return values

def test_basevectors_qd_array():
    apex_out = Apex(date=2000, refh=300)
    f1, f2 = apex_out.basevectors_qd([0, 30], 15, 100, coords='geo')
    f1_lat0, f2_lat0 = apex_out._basevec(0, 15, 100)
    f1_lat30, f2_lat30 = apex_out._basevec(30, 15, 100)
    assert_allclose(f1[:, 0], f1_lat0)
    assert_allclose(f2[:, 0], f2_lat0)
    assert_allclose(f1[:, 1], f1_lat30)
    assert_allclose(f2[:, 1], f2_lat30)


# ============================================================================
#  Test basevectors_apex()
# ============================================================================


# test against return from _geo2apexall for different coords

def test_basevectors_apex_scalar_geo():
    apex_out = Apex(date=2000, refh=300)

    (f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2,
     e3) = apex_out.basevectors_apex(60, 15, 100, coords='geo')

    (_, _, _, _, f1_, f2_, _, d1_, d2_, d3_, _, e1_, e2_,
     e3_) = apex_out._geo2apexall(60, 15, 100)

    assert_allclose(f1, f1_)
    assert_allclose(f2, f2_)
    assert_allclose(d1, d1_)
    assert_allclose(d2, d2_)
    assert_allclose(d3, d3_)
    assert_allclose(e1, e1_)
    assert_allclose(e2, e2_)
    assert_allclose(e3, e3_)


def test_basevectors_apex_scalar_apex():
    apex_out = Apex(date=2000, refh=300)

    (f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2,
     e3) = apex_out.basevectors_apex(60, 15, 100, coords='apex', precision=1e-2)

    glat, glon, _ = apex_out.apex2geo(60, 15, 100, precision=1e-2)
    (_, _, _, _, f1_, f2_, _, d1_, d2_, d3_, _, e1_, e2_,
     e3_) = apex_out._geo2apexall(glat, glon, 100)

    assert_allclose(f1, f1_)
    assert_allclose(f2, f2_)
    assert_allclose(d1, d1_)
    assert_allclose(d2, d2_)
    assert_allclose(d3, d3_)
    assert_allclose(e1, e1_)
    assert_allclose(e2, e2_)
    assert_allclose(e3, e3_)


def test_basevectors_apex_scalar_qd():
    apex_out = Apex(date=2000, refh=300)

    (f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2,
     e3) = apex_out.basevectors_apex(60, 15, 100, coords='qd', precision=1e-2)

    glat, glon, _ = apex_out.qd2geo(60, 15, 100, precision=1e-2)
    (_, _, _, _, f1_, f2_, _, d1_, d2_, d3_, _, e1_, e2_,
     e3_) = apex_out._geo2apexall(glat, glon, 100)

    assert_allclose(f1, f1_)
    assert_allclose(f2, f2_)
    assert_allclose(d1, d1_)
    assert_allclose(d2, d2_)
    assert_allclose(d3, d3_)
    assert_allclose(e1, e1_)
    assert_allclose(e2, e2_)
    assert_allclose(e3, e3_)


# test shapes and vectorization of arguments

def test_basevectors_apex_scalar_shape():
    apex_out = Apex(date=2000, refh=300)
    ret = apex_out.basevectors_apex(60, 15, 100, precision=1e-2)
    for r in ret[:2]:
        assert r.shape == (2,)
    for r in ret[2:]:
        assert r.shape == (3,)


def test_basevectors_apex_vectorization():
    apex_out = Apex(date=2000, refh=300)
    ret = apex_out.basevectors_apex([60, 60, 60, 60], 15, 100)
    for r in ret[:2]:
        assert r.shape == (2, 4)
    for r in ret[2:]:
        assert r.shape == (3, 4)
    ret = apex_out.basevectors_apex(60, [15, 15, 15, 15], 100)
    for r in ret[:2]:
        assert r.shape == (2, 4)
    for r in ret[2:]:
        assert r.shape == (3, 4)
    ret = apex_out.basevectors_apex(60, 15, [100, 100, 100, 100])
    for r in ret[:2]:
        assert r.shape == (2, 4)
    for r in ret[2:]:
        assert r.shape == (3, 4)


# test correct vectorization of height
def test_basevectors_apex_vectorization_height():
    apex_out = Apex(date=2000, refh=0)
    (f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2,
     e3) = apex_out.basevectors_apex(60, 15, [200, 400], coords='geo')
    (_, _, _, _, f1_1, f2_1, _, d1_1, d2_1, d3_1, _, e1_1, e2_1,
     e3_1) = apex_out._geo2apexall(60, 15, 200)
    (_, _, _, _, f1_2, f2_2, _, d1_2, d2_2, d3_2, _, e1_2, e2_2,
     e3_2) = apex_out._geo2apexall(60, 15, 400)

    assert_allclose(f1[:, 0], f1_1)
    assert_allclose(f2[:, 0], f2_1)
    assert_allclose(d1[:, 0], d1_1)
    assert_allclose(d2[:, 0], d2_1)
    assert_allclose(d3[:, 0], d3_1)
    assert_allclose(e1[:, 0], e1_1)
    assert_allclose(e2[:, 0], e2_1)
    assert_allclose(e3[:, 0], e3_1)

    assert_allclose(f3[:, 0], np.array([-0.088671, -0.018272, 0.993576]),
                    rtol=1e-4)
    assert_allclose(g1[:, 0], np.array([0.903098, 0.245273, 0.085107]),
                    rtol=1e-4)
    assert_allclose(g2[:, 0], np.array([-0.103495, 1.072078, 0.01048]),
                    rtol=1e-4)
    assert_allclose(g3[:, 0], np.array([0, 0, 1.006465]), rtol=1e-4)

    assert_allclose(f1[:, 1], f1_2)
    assert_allclose(f2[:, 1], f2_2)
    assert_allclose(d1[:, 1], d1_2)
    assert_allclose(d2[:, 1], d2_2)
    assert_allclose(d3[:, 1], d3_2)
    assert_allclose(e1[:, 1], e1_2)
    assert_allclose(e2[:, 1], e2_2)
    assert_allclose(e3[:, 1], e3_2)

    assert_allclose(f3[:, 1], np.array([-0.085415, -0.021176, 0.989645]),
                    rtol=1e-4)
    assert_allclose(g1[:, 1], np.array([0.902695, 0.246919, 0.083194]),
                    rtol=1e-4)
    assert_allclose(g2[:, 1], np.array([-0.11051, 1.066094, 0.013274]),
                    rtol=1e-4)
    assert_allclose(g3[:, 1], np.array([0, 0, 1.010463]), rtol=1e-4)


# test scalar return values

def test_basevectors_apex_scalar():
    apex_out = Apex(date=2000, refh=300)

    (f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2,
     e3) = apex_out.basevectors_apex(0, 15, 100, coords='geo')
    (_, _, _, _, f1_1, f2_1, _, d1_1, d2_1, d3_1, _, e1_1, e2_1,
     e3_1) = apex_out._geo2apexall(0, 15, 100)

    assert_allclose(f1, f1_1)
    assert_allclose(f2, f2_1)
    assert_allclose(d1, d1_1)
    assert_allclose(d2, d2_1)
    assert_allclose(d3, d3_1)
    assert_allclose(e1, e1_1)
    assert_allclose(e2, e2_1)
    assert_allclose(e3, e3_1)

    assert_allclose(f3, np.array([0.092637, -0.245951, 0.938848]), rtol=1e-4)
    assert_allclose(g1, np.array([0.939012, 0.073416, -0.07342]), rtol=1e-4)
    assert_allclose(g2, np.array([0.055389, 1.004155, 0.257594]), rtol=1e-4)
    assert_allclose(g3, np.array([0, 0, 1.065135]), rtol=1e-4)


# test 1D array return values

def test_basevectors_apex_array():
    apex_out = Apex(date=2000, refh=300)
    (f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2,
     e3) = apex_out.basevectors_apex([0, 30], 15, 100, coords='geo')
    (_, _, _, _, f1_1, f2_1, _, d1_1, d2_1, d3_1, _, e1_1, e2_1,
     e3_1) = apex_out._geo2apexall(0, 15, 100)
    (_, _, _, _, f1_2, f2_2, _, d1_2, d2_2, d3_2, _, e1_2, e2_2,
     e3_2) = apex_out._geo2apexall(30, 15, 100)

    assert_allclose(f1[:, 0], f1_1)
    assert_allclose(f2[:, 0], f2_1)
    assert_allclose(d1[:, 0], d1_1)
    assert_allclose(d2[:, 0], d2_1)
    assert_allclose(d3[:, 0], d3_1)
    assert_allclose(e1[:, 0], e1_1)
    assert_allclose(e2[:, 0], e2_1)
    assert_allclose(e3[:, 0], e3_1)

    assert_allclose(f3[:, 0], np.array([0.092637, -0.245951, 0.938848]),
                    rtol=1e-4)
    assert_allclose(g1[:, 0], np.array([0.939012, 0.073416, -0.07342]),
                    rtol=1e-4)
    assert_allclose(g2[:, 0], np.array([0.055389, 1.004155, 0.257594]),
                    rtol=1e-4)
    assert_allclose(g3[:, 0], np.array([0, 0, 1.065135]), rtol=1e-4)

    assert_allclose(f1[:, 1], f1_2)
    assert_allclose(f2[:, 1], f2_2)
    assert_allclose(d1[:, 1], d1_2)
    assert_allclose(d2[:, 1], d2_2)
    assert_allclose(d3[:, 1], d3_2)
    assert_allclose(e1[:, 1], e1_2)
    assert_allclose(e2[:, 1], e2_2)
    assert_allclose(e3[:, 1], e3_2)

    assert_allclose(f3[:, 1], np.array([-0.036618, -0.071019, 0.861604]),
                    rtol=1e-4)
    assert_allclose(g1[:, 1], np.array([0.844391, 0.015353, 0.037152]),
                    rtol=1e-4)
    assert_allclose(g2[:, 1], np.array([0.050808, 1.02131, 0.086342]),
                    rtol=1e-4)
    assert_allclose(g3[:, 1], np.array([0, 0, 1.160625]), rtol=1e-4)


# test that vectors are calculated correctly

def test_basevectors_apex_delta():
    apex_out = Apex(date=2000, refh=300)
    for lat in range(0, 90, 10):
        for lon in range(0, 360, 15):
            (f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2,
             e3) = apex_out.basevectors_apex(lat, lon, 500)
            f = [np.append(f1, 0), np.append(f2, 0), f3]
            g = [g1, g2, g3]
            d = [d1, d2, d3]
            e = [e1, e2, e3]
            for i, j in [(i, j) for i in range(3) for j in range(3)]:
                delta = 1 if i == j else 0
                assert_allclose(np.sum(f[i] * g[j]), delta, rtol=0, atol=1e-5)
                assert_allclose(np.sum(d[i] * e[j]), delta, rtol=0, atol=1e-5)


def test_basevectors_apex_invalid_scalar(recwarn):
    """Test warning and fill values for calculating base vectors with bad value.
    """
    apex_out = Apex(date=2000, refh=10000)
    base_vecs = apex_out.basevectors_apex(0, 0, 0)

    assert issubclass(recwarn[-1].category, UserWarning)
    assert 'set to NaN where' in str(recwarn[-1].message)

    invalid = np.ones(3) * np.nan
    for i, bvec in enumerate(base_vecs):
        if i < 2:
            assert not np.allclose(bvec, invalid[:2])
        else:
            assert_allclose(bvec, invalid)


# ============================================================================
#  Test the get_apex() method
# ============================================================================


def test_get_apex():
    apex_out = Apex(date=2000, refh=300)
    assert_allclose(apex_out.get_apex(10), 507.409702543805)
    assert_allclose(apex_out.get_apex(60), 20313.026999999987)


def test_get_apex_invalid_lat():
    apex_out = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        apex_out.get_apex(91)
    with pytest.raises(ValueError):
        apex_out.get_apex(-91)
    apex_out.get_apex(90)
    apex_out.get_apex(-90)

    assert_allclose(apex_out.get_apex(90 + 1e-5), apex_out.get_apex(90),
                    rtol=0, atol=1e-8)


# ============================================================================
#  Test the set_epoch() method
# ============================================================================


def test_set_epoch():
    """Test successful setting of Apex epoch."""
    apex_out = Apex(date=2000.2, refh=300)
    assert_allclose(apex_out.year, 2000.2)
    ret_2000_2_py = apex_out._geo2apex(60, 15, 100)
    apex_out.set_epoch(2000.8)
    assert_allclose(apex_out.year, 2000.8)
    ret_2000_8_py = apex_out._geo2apex(60, 15, 100)

    assert ret_2000_2_py != ret_2000_8_py

    fa.loadapxsh(apex_out.datafile, 2000.2)
    ret_2000_2_apex = fa.apxg2all(60, 15, 100, 300, 0)[2:4]
    fa.loadapxsh(apex_out.datafile, 2000.8)
    ret_2000_8_apex = fa.apxg2all(60, 15, 100, 300, 0)[2:4]

    assert ret_2000_2_apex != ret_2000_8_apex

    assert_allclose(ret_2000_2_py, ret_2000_2_apex)
    assert_allclose(ret_2000_8_py, ret_2000_8_apex)


@pytest.fixture()
def igrf_file():
    # Ensure the coefficient file exists
    original_file = os.path.join(os.path.dirname(helpers.__file__),
                                 'igrf13coeffs.txt')
    tmp_file = "temp_coeff.txt"
    assert os.path.isfile(original_file)
    # Move the coefficient file
    os.rename(original_file, tmp_file)
    yield original_file
    # Move the coefficient file back
    os.rename(tmp_file, original_file)


def test_set_epoch_file_error(igrf_file):
    """Test raises OSError when IGRF coefficient file is missing."""
    # Test missing coefficient file failure
    with pytest.raises(OSError) as oerr:
        Apex(date=2000.2, refh=300)
    error_string = "File {:} does not exist".format(igrf_file)
    assert str(oerr.value).startswith(error_string)


# ============================================================================
#  Test the set_refh() method
# ============================================================================


def test_set_refh():
    apex_out = Apex(date=2000, refh=300)
    assert apex_out.refh, 300
    ret_300 = apex_out._geo2apex(60, 15, 100)
    apex_out.set_refh(500)
    assert apex_out.refh == 500
    ret_500 = apex_out._geo2apex(60, 15, 100)

    assert_allclose(ret_300, fa.apxg2all(60, 15, 100, 300, 0)[2:4])
    assert_allclose(ret_500, fa.apxg2all(60, 15, 100, 500, 0)[2:4])


# ============================================================================
#  Test the get_babs() method
# ============================================================================


def test_get_babs():
    inputs = [[[80], [100], [300]], [range(50, 90, 8), range(0, 360, 80),
                                     [300] * 5], [90.0, 0, 1000]]
    temp1 = np.array([4.22045410e-05, 5.15672743e-05, 4.98150200e-05,
                     5.06769359e-05, 4.91028428e-05])
    expected = [[5.1303124427795412e-05], temp1, [3.793962299823761e-05]]

    apex_out = Apex(date=2018.1, refh=0)
    for i in range(len(inputs)):
        outputs = apex_out.get_babs(*inputs[i])
        if isinstance(outputs, np.float64):
            outputs = [outputs]
        for j, output in enumerate(outputs):
            assert_allclose(output, expected[i][j], rtol=0, atol=1e-5)


# ============================================================================
#  Test the bvectors_apex() method
# ============================================================================


def test_bvectors_apex():
    inputs = [[80, 81], [100, 120], [100, 200]]

    expected = (np.array([5.95166171e-05, 5.95958974e-05]),
                np.array([[0.0191583, 0.0020023],
                          [0.03547136, 0.03392595],
                          [-0.9412518, -0.8991005]]),
                np.array([5.28257734e-05, 4.82450628e-05]),
                np.array([[0.02158486, 0.00247339],
                          [0.03996412, 0.04190787],
                          [-1.0604696, -1.110636]]))

    apex_out = Apex(date=2018.1, refh=0)

    outputs = apex_out.bvectors_apex(*inputs, coords='geo', precision=1e-10)
    for i, output in enumerate(outputs):
        for j in range(output.size):
            assert_allclose(output.ravel()[j], expected[i].ravel()[j], rtol=0,
                            atol=1e-5)


if __name__ == '__main__':
    pytest.main()
