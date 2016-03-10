# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals

import datetime as dt
import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

from apexpy import fortranapex as fa
from apexpy import Apex, ApexHeightError, helpers


##############################################################################
# NOTE: whenever function outputs are tested against hard-coded numbers,     #
# the test results (numbers) were obtained by running the code that is       #
# tested. Therefore these tests below only check that nothing changes when   #
# refactoring etc., and not if the results are actually correct              #
##############################################################################


###============================================================================
### Test initiating the Apex class
###============================================================================


def test_init_defaults():
    Apex()


def test_init_date_int():
    A = Apex(date=2015)
    assert A.year == 2015


def test_init_date_float():
    A = Apex(date=2015.5)
    assert A.year == 2015.5


def test_init_date():
    date = dt.date(2015, 1, 1)
    A = Apex(date=date)
    assert A.year == helpers.toYearFraction(date)


def test_init_datetime():
    datetime = dt.datetime(2015, 6, 1, 18, 23, 45)
    A = Apex(date=datetime)
    assert A.year == helpers.toYearFraction(datetime)


def test_init_datafile_IOError():
    with pytest.raises(IOError):
        Apex(date=2015, datafile='foo/path/to/datafile.blah')


###============================================================================
### Test the low-level interfaces to the fortran wrappers
###============================================================================

def test__geo2qd_scalar():
    A = Apex(date=2000, refh=300)
    for lat in [0, 30, 60, 89]:
        for lon in [-179, -90, 0, 90, 180]:
            assert_allclose(A._geo2qd(lat, lon, 100), fa.apxg2q(lat, lon, 100, 0)[:2])


def test__geo2qd_array():
    A = Apex(date=2000, refh=300)
    lats, lons = A._geo2qd([[0, 30], [60, 90]], 15, [[100, 200], [300, 400]])
    lat1, lon1 = fa.apxg2q(0, 15, 100, 0)[:2]
    lat2, lon2 = fa.apxg2q(30, 15, 200, 0)[:2]
    lat3, lon3 = fa.apxg2q(60, 15, 300, 0)[:2]
    lat4, lon4 = fa.apxg2q(90, 15, 400, 0)[:2]
    assert_allclose(lats.astype(float), np.array([[lat1, lat2], [lat3, lat4]], dtype=float))
    assert_allclose(lons.astype(float), np.array([[lon1, lon2], [lon3, lon4]], dtype=float))


def test__geo2qd_longitude():
    A = Apex(date=2000, refh=300)
    assert_allclose(A._geo2qd(60, 180, 100), fa.apxg2q(60, 180, 100, 0)[:2])
    assert_allclose(A._geo2qd(60, -180, 100), fa.apxg2q(60, -180, 100, 0)[:2])
    assert_allclose(A._geo2qd(60, -180, 100), A._geo2qd(60, 180, 100))
    for i in range(-5, 5):
        for lat in [0, 30, 60, 90]:
            assert_allclose(A._geo2qd(lat, 15+i*360, 100), fa.apxg2q(lat, 15, 100, 0)[:2])


def test__geo2apex_scalar():
    A = Apex(date=2000, refh=300)
    for lat in [0, 30, 60, 89]:
        for lon in [-179, -90, 0, 90, 180]:
            assert_allclose(A._geo2apex(lat, lon, 100), fa.apxg2all(lat, lon, 100, 300, 0)[2:4])


def test__geo2apex_array():
    A = Apex(date=2000, refh=300)
    lats, lons = A._geo2apex([[0, 30], [60, 90]], 15, [[100, 200], [300, 400]])
    lat1, lon1 = fa.apxg2all(0, 15, 100, 300, 0)[2:4]
    lat2, lon2 = fa.apxg2all(30, 15, 200, 300, 0)[2:4]
    lat3, lon3 = fa.apxg2all(60, 15, 300, 300, 0)[2:4]
    lat4, lon4 = fa.apxg2all(90, 15, 400, 300, 0)[2:4]
    assert_allclose(lats.astype(float), np.array([[lat1, lat2], [lat3, lat4]], dtype=float))
    assert_allclose(lons.astype(float), np.array([[lon1, lon2], [lon3, lon4]], dtype=float))


def test__geo2apex_longitude():
    A = Apex(date=2000, refh=300)
    assert_allclose(A._geo2apex(60, 180, 100), fa.apxg2all(60, 180, 100, 300, 0)[2:4])
    assert_allclose(A._geo2apex(60, -180, 100), fa.apxg2all(60, -180, 100, 300, 0)[2:4])
    assert_allclose(A._geo2apex(60, -180, 100), A._geo2apex(60, 180, 100))
    for i in range(-5, 5):
        for lat in [0, 30, 60, 90]:
            assert_allclose(A._geo2apex(lat, 15+i*360, 100), fa.apxg2all(lat, 15, 100, 300, 0)[2:4])


def test__geo2apexall_scalar():
    A = Apex(date=2000, refh=300)
    for lat in [0, 30, 60, 89]:
        for lon in [-179, -90, 0, 90, 180]:
            ret1 = A._geo2apexall(lat, lon, 100)
            ret2 = fa.apxg2all(lat, lon, 100, 300, 1)
            for r1, r2 in zip(ret1, ret2):
                assert_allclose(r1, r2)


def test__geo2apexall_array():
    A = Apex(date=2000, refh=300)
    ret = A._geo2apexall([[0, 30], [60, 90]], 15, [[100, 200], [300, 400]])
    ret1 = fa.apxg2all(0, 15, 100, 300, 1)
    ret2 = fa.apxg2all(30, 15, 200, 300, 1)
    ret3 = fa.apxg2all(60, 15, 300, 300, 1)
    ret4 = fa.apxg2all(90, 15, 400, 300, 1)
    for i in range(len(ret)):
        try:
            # ret[i] is array of floats
            assert_allclose(ret[i].astype(float), np.array([[ret1[i], ret2[i]], [ret3[i], ret4[i]]], dtype=float))
        except:
            # ret[i] is array of arrays
            assert_allclose(ret[i][0, 0], ret1[i])
            assert_allclose(ret[i][0, 1], ret2[i])
            assert_allclose(ret[i][1, 0], ret3[i])
            assert_allclose(ret[i][1, 1], ret4[i])


def test__qd2geo_scalar():
    A = Apex(date=2000, refh=300)
    for lat in [0, 30, 60, 89]:
        for lon in [-179, -90, 0, 90, 180]:
            for prec in [-1, 1e-2, 1e-10]:
                assert_allclose(A._qd2geo(lat, lon, 100, prec), fa.apxq2g(lat, lon, 100, prec))


def test__qd2geo_array():
    A = Apex(date=2000, refh=300)
    lats, lons, errs = A._qd2geo([[0, 30], [60, 90]], 15, [[100, 200], [300, 400]], 1e-2)
    lat1, lon1, err1 = fa.apxq2g(0, 15, 100, 1e-2)
    lat2, lon2, err2 = fa.apxq2g(30, 15, 200, 1e-2)
    lat3, lon3, err3 = fa.apxq2g(60, 15, 300, 1e-2)
    lat4, lon4, err4 = fa.apxq2g(90, 15, 400, 1e-2)
    assert_allclose(lats.astype(float), np.array([[lat1, lat2], [lat3, lat4]], dtype=float))
    assert_allclose(lons.astype(float), np.array([[lon1, lon2], [lon3, lon4]], dtype=float))
    assert_allclose(errs.astype(float), np.array([[err1, err2], [err3, err4]], dtype=float))


def test__qd2geo_longitude():
    A = Apex(date=2000, refh=300)
    assert_allclose(A._qd2geo(60, 180, 100, 1e-2), fa.apxq2g(60, 180, 100, 1e-2))
    assert_allclose(A._qd2geo(60, -180, 100, 1e-2), fa.apxq2g(60, -180, 100, 1e-2))
    assert_allclose(A._qd2geo(60, -180, 100, 1e-2), A._qd2geo(60, 180, 100, 1e-2))
    for i in range(-5, 5):
        for lat in [0, 30, 60, 90]:
            assert_allclose(A._qd2geo(lat, 15+i*360, 100, 1e-2), fa.apxq2g(lat, 15, 100, 1e-2))


def test__basevec_scalar():
    A = Apex(date=2000, refh=300)
    for lat in [0, 30, 60, 89]:
        for lon in [-179, -90, 0, 90, 180]:
            assert_allclose(A._basevec(lat, lon, 100), fa.apxg2q(lat, lon, 100, 1)[2:4])


def test__basevec_array():
    A = Apex(date=2000, refh=300)
    f1s, f2s = A._basevec([[0, 30], [60, 90]], 15, [[100, 200], [300, 400]])
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
    A = Apex(date=2000, refh=300)
    assert_allclose(A._basevec(60, 180, 100), fa.apxg2q(60, 180, 100, 1)[2:4])
    assert_allclose(A._basevec(60, -180, 100), fa.apxg2q(60, -180, 100, 1)[2:4])
    assert_allclose(A._basevec(60, -180, 100), A._basevec(60, 180, 100))
    for i in range(-5, 5):
        for lat in [0, 30, 60, 90]:
            assert_allclose(A._basevec(lat, 15+i*360, 100), fa.apxg2q(lat, 15, 100, 1)[2:4])


###============================================================================
### Test the convert() method
###============================================================================


def test_convert_geo2apex():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.convert(60, 15, 'geo', 'apex', height=100), A.geo2apex(60, 15, 100))


def test_convert_geo2qd():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.convert(60, 15, 'geo', 'qd', height=100), A.geo2qd(60, 15, 100))


def test_convert_geo2mlt_nodate():
    A = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        A.convert(60, 15, 'geo', 'mlt')


def test_convert_geo2mlt():
    datetime = dt.datetime(2000, 3, 9, 14, 25, 58)
    A = Apex(date=2000, refh=300)
    assert_allclose(A.convert(60, 15, 'geo', 'mlt', height=100, ssheight=2e5, datetime=datetime)[1], A.mlon2mlt(A.geo2apex(60, 15, 100)[1], datetime, ssheight=2e5))


def test_convert_apex2geo():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.convert(60, 15, 'apex', 'geo', height=100, precision=1e-2), A.apex2geo(60, 15, 100, precision=1e-2)[:-1])


def test_convert_apex2qd():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.convert(60, 15, 'apex', 'qd', height=100), A.apex2qd(60, 15, height=100))


def test_convert_apex2mlt():
    datetime = dt.datetime(2000, 3, 9, 14, 25, 58)
    A = Apex(date=2000, refh=300)
    assert_allclose(A.convert(60, 15, 'apex', 'mlt', height=100, datetime=datetime, ssheight=2e5)[1], A.mlon2mlt(15, datetime, ssheight=2e5))


def test_convert_qd2geo():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.convert(60, 15, 'qd', 'geo', height=100, precision=1e-2), A.qd2geo(60, 15, 100, precision=1e-2)[:-1])


def test_convert_qd2apex():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.convert(60, 15, 'qd', 'apex', height=100), A.qd2apex(60, 15, height=100))


def test_convert_qd2mlt():
    datetime = dt.datetime(2000, 3, 9, 14, 25, 58)
    A = Apex(date=2000, refh=300)
    assert_allclose(A.convert(60, 15, 'qd', 'mlt', height=100, datetime=datetime, ssheight=2e5)[1], A.mlon2mlt(15, datetime, ssheight=2e5))


def test_convert_mlt2geo():
    datetime = dt.datetime(2000, 3, 9, 14, 25, 58)
    A = Apex(date=2000, refh=300)
    assert_allclose(A.convert(60, 15, 'mlt', 'geo', height=100, datetime=datetime, precision=1e-2, ssheight=2e5), A.apex2geo(60, A.mlt2mlon(15, datetime, ssheight=2e5), 100, precision=1e-2)[:-1])


def test_convert_mlt2apex():
    datetime = dt.datetime(2000, 3, 9, 14, 25, 58)
    A = Apex(date=2000, refh=300)
    assert_allclose(A.convert(60, 15, 'mlt', 'apex', height=100, datetime=datetime, ssheight=2e5), (60, A.mlt2mlon(15, datetime, ssheight=2e5)))


def test_convert_mlt2qd():
    datetime = dt.datetime(2000, 3, 9, 14, 25, 58)
    A = Apex(date=2000, refh=300)
    assert_allclose(A.convert(60, 15, 'mlt', 'qd', height=100, datetime=datetime, ssheight=2e5), A.apex2qd(60, A.mlt2mlon(15, datetime, ssheight=2e5), height=100))


def test_convert_invalid_lat():
    A = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        A.convert(91, 0, 'geo', 'geo')
    with pytest.raises(ValueError):
        A.convert(-91, 0, 'geo', 'geo')
    A.convert(90, 0, 'geo', 'geo')
    A.convert(-90, 0, 'geo', 'geo')

    assert_allclose(A.convert(90+1e-5, 0, 'geo', 'apex'), A.convert(90, 0, 'geo', 'apex'), rtol=0, atol=1e-8)


def test_convert_invalid_transformation():
    A = Apex(date=2000, refh=300)
    with pytest.raises(NotImplementedError):
        A.convert(0, 0, 'foobar', 'geo')
    with pytest.raises(NotImplementedError):
        A.convert(0, 0, 'geo', 'foobar')


###============================================================================
### Test the geo2apex() method
###============================================================================


def test_geo2apex():
    A = Apex(date=2000, refh=300)
    lat, lon = A.geo2apex(60, 15, 100)
    assert_allclose((lat, lon), A._geo2apex(60, 15, 100))
    assert type(lat) != np.ndarray
    assert type(lon) != np.ndarray


def test_geo2apex_vectorization():
    A = Apex(date=2000, refh=300)
    assert A.geo2apex([60, 60], 15, 100)[0].shape == (2,)
    assert A.geo2apex(60, [15, 15], 100)[0].shape == (2,)
    assert A.geo2apex(60, 15, [100, 100])[0].shape == (2,)


def test_geo2apex_invalid_lat():
    A = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        A.geo2apex(91, 0, 0)
    with pytest.raises(ValueError):
        A.geo2apex(-91, 0, 0)
    A.geo2apex(90, 0, 0)
    A.geo2apex(-90, 0, 0)

    assert_allclose(A.geo2apex(90+1e-5, 0, 0), A.geo2apex(90, 0, 0), rtol=0, atol=1e-8)


def test_geo2apex_undefined_warning():
    A = Apex(date=2000, refh=10000)
    with warnings.catch_warnings(record=True) as w:
        ret = A.geo2apex(0, 0, 0)
        A.geo2apex(0, 0, 0)
        assert ret[0] == -9999
        assert len(w) == 2
        assert issubclass(w[-1].category, UserWarning)
        assert 'set to -9999 where' in str(w[-1].message)


###============================================================================
### Test the apex2geo() method
###============================================================================


def test_apex2geo():
    A = Apex(date=2000, refh=300)
    lat, lon, error = A.apex2geo(60, 15, 100, precision=1e-2)
    assert_allclose((lat, lon, error),
                    A.qd2geo(*A.apex2qd(60, 15, 100), height=100, precision=1e-2))
    assert type(lat) != np.ndarray
    assert type(lon) != np.ndarray
    assert type(error) != np.ndarray


def test_apex2geo_vectorization():
    A = Apex(date=2000, refh=300)
    assert A.apex2geo([60, 60], 15, 100)[0].shape == (2,)
    assert A.apex2geo(60, [15, 15], 100)[0].shape == (2,)
    assert A.apex2geo(60, 15, [100, 100])[0].shape == (2,)


def test_apex2geo_invalid_lat():
    A = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        A.apex2geo(91, 0, 0, 1e-2)
    with pytest.raises(ValueError):
        A.apex2geo(-91, 0, 0, 1e-2)
    A.apex2geo(90, 0, 0, 1e-2)
    A.apex2geo(-90, 0, 0, 1e-2)

    assert_allclose(A.apex2geo(90+1e-5, 0, 0, 1e-2), A.apex2geo(90, 0, 0, 1e-2), rtol=0, atol=1e-8)


###============================================================================
### Test the geo2qd() method
###============================================================================


def test_geo2qd():
    A = Apex(date=2000, refh=300)
    lat, lon = A.geo2qd(60, 15, 100)
    assert_allclose((lat, lon), A._geo2qd(60, 15, 100))
    assert type(lat) != np.ndarray
    assert type(lon) != np.ndarray


def test_geo2qd_vectorization():
    A = Apex(date=2000, refh=300)
    assert A.geo2qd([60, 60], 15, 100)[0].shape == (2,)
    assert A.geo2qd(60, [15, 15], 100)[0].shape == (2,)
    assert A.geo2qd(60, 15, [100, 100])[0].shape == (2,)


def test_geo2qd_invalid_lat():
    A = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        A.geo2qd(91, 0, 0)
    with pytest.raises(ValueError):
        A.geo2qd(-91, 0, 0)
    A.geo2qd(90, 0, 0)
    A.geo2qd(-90, 0, 0)

    assert_allclose(A.geo2qd(90+1e-5, 0, 0), A.geo2qd(90, 0, 0), rtol=0, atol=1e-8)


###============================================================================
### Test the qd2geo() method
###============================================================================


def test_qd2geo():
    A = Apex(date=2000, refh=300)
    lat, lon, error = A.qd2geo(60, 15, 100, precision=1e-2)
    assert_allclose((lat, lon, error), A._qd2geo(60, 15, 100, 1e-2))
    assert type(lat) != np.ndarray
    assert type(lon) != np.ndarray
    assert type(error) != np.ndarray


def test_qd2geo_vectorization():
    A = Apex(date=2000, refh=300)
    assert A.qd2geo([60, 60], 15, 100)[0].shape == (2,)
    assert A.qd2geo(60, [15, 15], 100)[0].shape == (2,)
    assert A.qd2geo(60, 15, [100, 100])[0].shape == (2,)


def test_qd2geo_invalid_lat():
    A = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        A.qd2geo(91, 0, 0, precision=1e-2)
    with pytest.raises(ValueError):
        A.qd2geo(-91, 0, 0, precision=1e-2)
    A.qd2geo(90, 0, 0, precision=1e-2)
    A.qd2geo(-90, 0, 0, precision=1e-2)

    assert_allclose(A.qd2geo(90+1e-5, 0, 0, 1e-2), A.qd2geo(90, 0, 0, 1e-2), rtol=0, atol=1e-8)


###============================================================================
### Test the apex2qd() method
###============================================================================


def test_apex2qd():
    A = Apex(date=2000, refh=300)
    lat, lon = A.apex2qd(60, 15, 100)
    assert_allclose((lat, lon),
                    [60.498401, 15])
    assert type(lat) != np.ndarray
    assert type(lon) != np.ndarray


def test_apex2qd_vectorization():
    A = Apex(date=2000, refh=300)
    assert A.apex2qd([60, 60], 15, 100)[0].shape == (2,)
    assert A.apex2qd(60, [15, 15], 100)[0].shape == (2,)
    assert A.apex2qd(60, 15, [100, 100])[0].shape == (2,)


def test_apex2qd_invalid_lat():
    A = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        A.apex2qd(91, 0, 0)
    with pytest.raises(ValueError):
        A.apex2qd(-91, 0, 0)
    A.apex2qd(90, 0, 0)
    A.apex2qd(-90, 0, 0)

    assert_allclose(A.apex2qd(90+1e-5, 0, 0), A.apex2qd(90, 0, 0), rtol=0, atol=1e-8)


def test_apex2qd_apexheight_close():
    A = Apex(date=2000, refh=300)
    A.apex2qd(0, 15, 300+1e-6)


def test_apex2qd_apexheight_over():
    A = Apex(date=2000, refh=300)
    with pytest.raises(ApexHeightError):
        A.apex2qd(0, 15, 301)


###============================================================================
### Test the qd2apex() method
###============================================================================


def test_qd2apex():
    A = Apex(date=2000, refh=300)
    lat, lon = A.qd2apex(60, 15, 100)
    assert_allclose((lat, lon),
                    [59.491381, 15])
    assert type(lat) != np.ndarray
    assert type(lon) != np.ndarray


def test_qd2apex_vectorization():
    A = Apex(date=2000, refh=300)
    assert A.qd2apex([60, 60], 15, 100)[0].shape == (2,)
    assert A.qd2apex(60, [15, 15], 100)[0].shape == (2,)
    assert A.qd2apex(60, 15, [100, 100])[0].shape == (2,)


def test_qd2apex_invalid_lat():
    A = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        A.qd2apex(91, 0, 0)
    with pytest.raises(ValueError):
        A.qd2apex(-91, 0, 0)
    A.qd2apex(90, 0, 0)
    A.qd2apex(-90, 0, 0)

    assert_allclose(A.qd2apex(90+1e-5, 0, 0), A.qd2apex(90, 0, 0), rtol=0, atol=1e-8)


def test_qd2apex_apexheight_close():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.qd2apex(0, 15, 300-1e-5), A.qd2apex(0, 15, 300))


def test_qd2apex_apexheight_over():
    A = Apex(date=2000, refh=300)
    with pytest.raises(ApexHeightError):
        A.qd2apex(0, 15, 299)


###============================================================================
### Test mlon2mlt()
###============================================================================


def test_mlon2mlt_scalar():
    A = Apex(date=2000, refh=300)
    mlon = A.mlon2mlt(0, dt.datetime(2000, 2, 3, 4, 5, 6))
    assert_allclose(mlon, 23.019629923502603)
    assert type(mlon) != np.ndarray


def test_mlon2mlt_ssheight():
    A = Apex(date=2000, refh=300)
    mlt = A.mlon2mlt(0, dt.datetime(2000, 2, 3, 4, 5, 6), ssheight=50*2000)
    assert_allclose(mlt, 23.026712036132814)


def test_mlon2mlt_1Darray():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.mlon2mlt([0, 180], dt.datetime(2000, 2, 3, 4, 5, 6)), [23.019261, 11.019261], rtol=1e-4)


def test_mlon2mlt_2Darray():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.mlon2mlt([[0, 180], [0, 180]], dt.datetime(2000, 2, 3, 4, 5, 6)), [[23.019261, 11.019261], [23.019261, 11.019261]], rtol=1e-4)


def test_mlon2mlt_diffdates():
    A = Apex(date=2000, refh=300)
    assert A.mlon2mlt(0, dt.datetime(2000, 2, 3, 4, 5, 6)) != A.mlon2mlt(0, dt.datetime(2000, 2, 3, 5, 5, 6))


def test_mlon2mlt_offset():
    A = Apex(date=2000, refh=300)
    date = dt.datetime(2000, 2, 3, 4, 5, 6)
    assert_allclose(A.mlon2mlt(0, date), A.mlon2mlt(-15, date) + 1)
    assert_allclose(A.mlon2mlt(0, date), A.mlon2mlt(-10*15, date) + 10)


def test_mlon2mlt_range():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.mlon2mlt(range(0, 361, 30), dt.datetime(2000, 2, 3, 4, 5, 6)),
                    [23.01963, 1.01963, 3.01963, 5.01963, 7.01963,
                     9.01963, 11.01963, 13.01963, 15.01963, 17.01963,
                     19.01963, 21.01963, 23.01963],
                    rtol=1e-4)


###============================================================================
### Test mlt2mlon()
###============================================================================


def test_mlt2mlon_scalar():
    A = Apex(date=2000, refh=300)
    mlt = A.mlt2mlon(0, dt.datetime(2000, 2, 3, 4, 5, 6))
    assert_allclose(mlt, 14.705551147460938)
    assert type(mlt) != np.ndarray


def test_mlt2mlon_ssheight():
    A = Apex(date=2000, refh=300)
    mlt = A.mlt2mlon(0, dt.datetime(2000, 2, 3, 4, 5, 6), ssheight=50*2000)
    assert_allclose(mlt, 14.599319458007812)


def test_mlt2mlon_1Darray():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.mlt2mlon([0, 12], dt.datetime(2000, 2, 3, 4, 5, 6)), [14.705551, 194.705551], rtol=1e-4)


def test_mlt2mlon_2Darray():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.mlt2mlon([[0, 12], [0, 12]], dt.datetime(2000, 2, 3, 4, 5, 6)), [[14.705551, 194.705551], [14.705551, 194.705551]], rtol=1e-4)


def test_mlt2mlon_diffdates():
    A = Apex(date=2000, refh=300)
    assert A.mlt2mlon(0, dt.datetime(2000, 2, 3, 4, 5, 6)) != A.mlt2mlon(0, dt.datetime(2000, 2, 3, 5, 5, 6))


def test_mlt2mlon_offset():
    A = Apex(date=2000, refh=300)
    date = dt.datetime(2000, 2, 3, 4, 5, 6)
    assert_allclose(A.mlt2mlon(0, date), A.mlt2mlon(1, date) - 15)
    assert_allclose(A.mlt2mlon(0, date), A.mlt2mlon(10, date) - 150)


def test_mlt2mlon_range():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.mlt2mlon(range(0, 25, 2), dt.datetime(2000, 2, 3, 4, 5, 6)),
                    [14.705551, 44.705551, 74.705551, 104.705551, 134.705551,
                     164.705551, 194.705551, 224.705551, 254.705551, 284.705551,
                     314.705551, 344.705551, 14.705551],
                    rtol=1e-4)


###============================================================================
### Test mlt/mlon back and forth
###============================================================================


def test_mlon2mlt2mlon():
    A = Apex(date=2000, refh=300)
    date = dt.datetime(2000, 2, 3, 4, 5, 6)
    assert_allclose(A.mlon2mlt(A.mlt2mlon(0, date), date), 0)
    assert_allclose(A.mlon2mlt(A.mlt2mlon(6, date), date), 6)
    assert_allclose(A.mlon2mlt(A.mlt2mlon(12, date), date), 12)
    assert_allclose(A.mlon2mlt(A.mlt2mlon(18, date), date), 18)
    assert_allclose(A.mlon2mlt(A.mlt2mlon(24, date), date), 0)


def test_mlt2mlon2mlt():
    A = Apex(date=2000, refh=300)
    date = dt.datetime(2000, 2, 3, 4, 5, 6)
    assert_allclose(A.mlt2mlon(A.mlon2mlt(0, date), date), 0)
    assert_allclose(A.mlt2mlon(A.mlon2mlt(90, date), date), 90)
    assert_allclose(A.mlt2mlon(A.mlon2mlt(180, date), date), 180)
    assert_allclose(A.mlt2mlon(A.mlon2mlt(270, date), date), 270)
    assert_allclose(A.mlt2mlon(A.mlon2mlt(360, date), date), 0)


###============================================================================
### Test the map_to_height() method
###============================================================================


def test_map_to_height():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.map_to_height(60, 15, 100, 10000, conjugate=False, precision=1e-10), (31.841459274291992, 17.916629791259766, 0))
    assert_allclose(A.map_to_height(30, 170, 100, 500, conjugate=False, precision=1e-2), (25.727252960205078, 169.60546875, 0.00017655163537710905))


def test_map_to_height_same_height():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.map_to_height(60, 15, 100, 100, conjugate=False, precision=1e-10), (60, 15, 3.4150946248701075e-6), rtol=1e-5)


def test_map_to_height_conjugate():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.map_to_height(60, 15, 100, 10000, conjugate=True, precision=1e-10), (-25.424892425537109, 27.310417175292969, 1.2074182222931995e-6))
    assert_allclose(A.map_to_height(30, 170, 100, 500, conjugate=True, precision=1e-2), (-13.76642894744873, 164.24259948730469, 0.00056820799363777041))


def test_map_to_height_vectorization():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.map_to_height([60, 60], 15, 100, 100), ([60]*2, [15]*2, [3.4150946248701075e-6]*2), rtol=1e-5)
    assert_allclose(A.map_to_height(60, [15, 15], 100, 100), ([60]*2, [15]*2, [3.4150946248701075e-6]*2), rtol=1e-5)
    assert_allclose(A.map_to_height(60, 15, [100, 100], 100), ([60]*2, [15]*2, [3.4150946248701075e-6]*2), rtol=1e-5)
    assert_allclose(A.map_to_height(60, 15, 100, [100, 100]), ([60]*2, [15]*2, [3.4150946248701075e-6]*2), rtol=1e-5)


def test_map_to_height_ApexHeightError():
    A = Apex(date=2000, refh=300)
    with pytest.raises(ApexHeightError):
        A.map_to_height(0, 15, 100, 10000)


###============================================================================
### Test the map_E_to_height() method
###============================================================================


def test_map_E_to_height():
    A = Apex(date=2000, refh=300)
    out_60_15_100_500 = [0.7115211, 2.3562392, 0.57259707]
    out_60_15_100_500_234 = [1.560284, 3.439154, 0.782339]
    out_60_15_100_1000 = [0.677964, 2.089811, 0.558601]
    out_60_15_200_500 = [0.723773, 2.427366, 0.590826]
    out_60_30_100_500 = [0.686265, 2.375296, 0.600594]
    out_70_15_100_500 = [0.727605, 2.180817, 0.291414]

    # scalar
    assert_allclose(A.map_E_to_height(60, 15, 100, 500, [1, 2, 3]), out_60_15_100_500, rtol=1e-5)
    assert_allclose(A.map_E_to_height(60, 15, 100, 500, [2, 3, 4]), out_60_15_100_500_234, rtol=1e-5)
    assert_allclose(A.map_E_to_height(60, 15, 100, 1000, [1, 2, 3]), out_60_15_100_1000, rtol=1e-5)
    assert_allclose(A.map_E_to_height(60, 15, 200, 500, [1, 2, 3]), out_60_15_200_500, rtol=1e-5)
    assert_allclose(A.map_E_to_height(60, 30, 100, 500, [1, 2, 3]), out_60_30_100_500, rtol=1e-5)
    assert_allclose(A.map_E_to_height(70, 15, 100, 500, [1, 2, 3]), out_70_15_100_500, rtol=1e-5)

    # vectorize lat
    assert_allclose(A.map_E_to_height([60, 70], 15, 100, 500, np.array([[1, 2, 3]]*2).T),
                    np.array([out_60_15_100_500, out_70_15_100_500]).T, rtol=1e-5)

    # vectorize lon
    assert_allclose(A.map_E_to_height(60, [15, 30], 100, 500, np.array([[1, 2, 3]]*2).T),
                    np.array([out_60_15_100_500, out_60_30_100_500]).T, rtol=1e-5)

    # vectorize height
    assert_allclose(A.map_E_to_height(60, 15, [100, 200], 500, np.array([[1, 2, 3]]*2).T),
                    np.array([out_60_15_100_500, out_60_15_200_500]).T, rtol=1e-5)

    # vectorize newheight
    assert_allclose(A.map_E_to_height(60, 15, 100, [500, 1000], np.array([[1, 2, 3]]*2).T),
                    np.array([out_60_15_100_500, out_60_15_100_1000]).T, rtol=1e-5)

    # vectorize E
    assert_allclose(A.map_E_to_height(60, 15, 100, 500, np.array([[1, 2, 3], [2, 3, 4]]).T),
                    np.array([out_60_15_100_500, out_60_15_100_500_234]).T, rtol=1e-5)


###============================================================================
### Test the map_V_to_height() method
###============================================================================


def test_map_V_to_height():
    A = Apex(date=2000, refh=300)
    out_60_15_100_500 = [0.819719, 2.845114, 0.695437]
    out_60_15_100_500_234 = [1.830277, 4.14345, 0.947624]
    out_60_15_100_1000 = [0.924577, 3.149964, 0.851343]
    out_60_15_200_500 = [0.803882, 2.793206, 0.682839]
    out_60_30_100_500 = [0.761412, 2.878837, 0.736549]
    out_70_15_100_500 = [0.846819, 2.592572, 0.347919]

    # scalar
    assert_allclose(A.map_V_to_height(60, 15, 100, 500, [1, 2, 3]), out_60_15_100_500, rtol=1e-5)
    assert_allclose(A.map_V_to_height(60, 15, 100, 500, [2, 3, 4]), out_60_15_100_500_234, rtol=1e-5)
    assert_allclose(A.map_V_to_height(60, 15, 100, 1000, [1, 2, 3]), out_60_15_100_1000, rtol=1e-5)
    assert_allclose(A.map_V_to_height(60, 15, 200, 500, [1, 2, 3]), out_60_15_200_500, rtol=1e-5)
    assert_allclose(A.map_V_to_height(60, 30, 100, 500, [1, 2, 3]), out_60_30_100_500, rtol=1e-5)
    assert_allclose(A.map_V_to_height(70, 15, 100, 500, [1, 2, 3]), out_70_15_100_500, rtol=1e-5)

    # vectorize lat
    assert_allclose(A.map_V_to_height([60, 70], 15, 100, 500, np.array([[1, 2, 3]]*2).T),
                    np.array([out_60_15_100_500, out_70_15_100_500]).T, rtol=1e-5)

    # vectorize lon
    assert_allclose(A.map_V_to_height(60, [15, 30], 100, 500, np.array([[1, 2, 3]]*2).T),
                    np.array([out_60_15_100_500, out_60_30_100_500]).T, rtol=1e-5)

    # vectorize height
    assert_allclose(A.map_V_to_height(60, 15, [100, 200], 500, np.array([[1, 2, 3]]*2).T),
                    np.array([out_60_15_100_500, out_60_15_200_500]).T, rtol=1e-5)

    # vectorize newheight
    assert_allclose(A.map_V_to_height(60, 15, 100, [500, 1000], np.array([[1, 2, 3]]*2).T),
                    np.array([out_60_15_100_500, out_60_15_100_1000]).T, rtol=1e-5)

    # vectorize E
    assert_allclose(A.map_V_to_height(60, 15, 100, 500, np.array([[1, 2, 3], [2, 3, 4]]).T),
                    np.array([out_60_15_100_500, out_60_15_100_500_234]).T, rtol=1e-5)


###============================================================================
### Test basevectors_qd()
###============================================================================


# test coords

def test_basevectors_qd_scalar_geo():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.basevectors_qd(60, 15, 100, coords='geo'), A._basevec(60, 15, 100))


def test_basevectors_qd_scalar_apex():
    A = Apex(date=2000, refh=300)
    glat, glon, _ = A.apex2geo(60, 15, 100, precision=1e-2)
    assert_allclose(A.basevectors_qd(60, 15, 100, coords='apex', precision=1e-2), A._basevec(glat, glon, 100))


def test_basevectors_qd_scalar_qd():
    A = Apex(date=2000, refh=300)
    glat, glon, _ = A.qd2geo(60, 15, 100, precision=1e-2)
    assert_allclose(A.basevectors_qd(60, 15, 100, coords='qd', precision=1e-2), A._basevec(glat, glon, 100))


# test shapes and vectorization of arguments

def test_basevectors_qd_scalar_shape():
    A = Apex(date=2000, refh=300)
    ret = A.basevectors_qd(60, 15, 100)
    for r in ret:
        assert r.shape == (2,)


def test_basevectors_qd_vectorization():
    A = Apex(date=2000, refh=300)
    ret = A.basevectors_qd([60, 60, 60, 60], 15, 100, coords='geo')
    for r in ret:
        assert r.shape == (2, 4)
    ret = A.basevectors_qd(60, [15, 15, 15, 15], 100, coords='geo')
    for r in ret:
        assert r.shape == (2, 4)
    ret = A.basevectors_qd(60, 15, [100, 100, 100, 100], coords='geo')
    for r in ret:
        assert r.shape == (2, 4)


# test array return values

def test_basevectors_qd_array():
    A = Apex(date=2000, refh=300)
    f1, f2 = A.basevectors_qd([0, 30], 15, 100, coords='geo')
    f1_lat0, f2_lat0 = A._basevec(0, 15, 100)
    f1_lat30, f2_lat30 = A._basevec(30, 15, 100)
    assert_allclose(f1[:, 0], f1_lat0)
    assert_allclose(f2[:, 0], f2_lat0)
    assert_allclose(f1[:, 1], f1_lat30)
    assert_allclose(f2[:, 1], f2_lat30)


###============================================================================
### Test basevectors_apex()
###============================================================================


# test against return from _geo2apexall for different coords

def test_basevectors_apex_scalar_geo():
    A = Apex(date=2000, refh=300)

    f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = A.basevectors_apex(60, 15, 100, coords='geo')

    _, _, _, _, f1_, f2_, _, d1_, d2_, d3_, _, e1_, e2_, e3_ = A._geo2apexall(60, 15, 100)

    assert_allclose(f1, f1_)
    assert_allclose(f2, f2_)
    assert_allclose(d1, d1_)
    assert_allclose(d2, d2_)
    assert_allclose(d3, d3_)
    assert_allclose(e1, e1_)
    assert_allclose(e2, e2_)
    assert_allclose(e3, e3_)


def test_basevectors_apex_scalar_apex():
    A = Apex(date=2000, refh=300)

    f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = A.basevectors_apex(60, 15, 100, coords='apex', precision=1e-2)

    glat, glon, _ = A.apex2geo(60, 15, 100, precision=1e-2)
    _, _, _, _, f1_, f2_, _, d1_, d2_, d3_, _, e1_, e2_, e3_ = A._geo2apexall(glat, glon, 100)

    assert_allclose(f1, f1_)
    assert_allclose(f2, f2_)
    assert_allclose(d1, d1_)
    assert_allclose(d2, d2_)
    assert_allclose(d3, d3_)
    assert_allclose(e1, e1_)
    assert_allclose(e2, e2_)
    assert_allclose(e3, e3_)


def test_basevectors_apex_scalar_qd():
    A = Apex(date=2000, refh=300)

    f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = A.basevectors_apex(60, 15, 100, coords='qd', precision=1e-2)

    glat, glon, _ = A.qd2geo(60, 15, 100, precision=1e-2)
    _, _, _, _, f1_, f2_, _, d1_, d2_, d3_, _, e1_, e2_, e3_ = A._geo2apexall(glat, glon, 100)

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
    A = Apex(date=2000, refh=300)
    ret = A.basevectors_apex(60, 15, 100, precision=1e-2)
    for r in ret[:2]:
        assert r.shape == (2,)
    for r in ret[2:]:
        assert r.shape == (3,)


def test_basevectors_apex_vectorization():
    A = Apex(date=2000, refh=300)
    ret = A.basevectors_apex([60, 60, 60, 60], 15, 100)
    for r in ret[:2]:
        assert r.shape == (2, 4)
    for r in ret[2:]:
        assert r.shape == (3, 4)
    ret = A.basevectors_apex(60, [15, 15, 15, 15], 100)
    for r in ret[:2]:
        assert r.shape == (2, 4)
    for r in ret[2:]:
        assert r.shape == (3, 4)
    ret = A.basevectors_apex(60, 15, [100, 100, 100, 100])
    for r in ret[:2]:
        assert r.shape == (2, 4)
    for r in ret[2:]:
        assert r.shape == (3, 4)


# test correct vectorization of height
def test_basevectors_apex_vectorization_height():
    A = Apex(date=2000, refh=0)
    f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = A.basevectors_apex(60, 15, [200, 400], coords='geo')
    _, _, _, _, f1_1, f2_1, _, d1_1, d2_1, d3_1, _, e1_1, e2_1, e3_1 = A._geo2apexall(60, 15, 200)
    _, _, _, _, f1_2, f2_2, _, d1_2, d2_2, d3_2, _, e1_2, e2_2, e3_2 = A._geo2apexall(60, 15, 400)

    assert_allclose(f1[:, 0], f1_1)
    assert_allclose(f2[:, 0], f2_1)
    assert_allclose(d1[:, 0], d1_1)
    assert_allclose(d2[:, 0], d2_1)
    assert_allclose(d3[:, 0], d3_1)
    assert_allclose(e1[:, 0], e1_1)
    assert_allclose(e2[:, 0], e2_1)
    assert_allclose(e3[:, 0], e3_1)

    assert_allclose(f3[:, 0], np.array([-0.088671, -0.018272, 0.993576]), rtol=1e-4)
    assert_allclose(g1[:, 0], np.array([0.903098, 0.245273, 0.085107]), rtol=1e-4)
    assert_allclose(g2[:, 0], np.array([-0.103495, 1.072078, 0.01048]), rtol=1e-4)
    assert_allclose(g3[:, 0], np.array([0, 0, 1.006465]), rtol=1e-4)

    assert_allclose(f1[:, 1], f1_2)
    assert_allclose(f2[:, 1], f2_2)
    assert_allclose(d1[:, 1], d1_2)
    assert_allclose(d2[:, 1], d2_2)
    assert_allclose(d3[:, 1], d3_2)
    assert_allclose(e1[:, 1], e1_2)
    assert_allclose(e2[:, 1], e2_2)
    assert_allclose(e3[:, 1], e3_2)

    assert_allclose(f3[:, 1], np.array([-0.085415, -0.021176, 0.989645]), rtol=1e-4)
    assert_allclose(g1[:, 1], np.array([0.902695, 0.246919, 0.083194]), rtol=1e-4)
    assert_allclose(g2[:, 1], np.array([-0.11051, 1.066094, 0.013274]), rtol=1e-4)
    assert_allclose(g3[:, 1], np.array([0, 0, 1.010463]), rtol=1e-4)


# test scalar return values

def test_basevectors_apex_scalar():
    A = Apex(date=2000, refh=300)

    f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = A.basevectors_apex(0, 15, 100, coords='geo')
    _, _, _, _, f1_1, f2_1, _, d1_1, d2_1, d3_1, _, e1_1, e2_1, e3_1 = A._geo2apexall(0, 15, 100)

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
    A = Apex(date=2000, refh=300)
    f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = A.basevectors_apex([0, 30], 15, 100, coords='geo')
    _, _, _, _, f1_1, f2_1, _, d1_1, d2_1, d3_1, _, e1_1, e2_1, e3_1 = A._geo2apexall(0, 15, 100)
    _, _, _, _, f1_2, f2_2, _, d1_2, d2_2, d3_2, _, e1_2, e2_2, e3_2 = A._geo2apexall(30, 15, 100)

    assert_allclose(f1[:, 0], f1_1)
    assert_allclose(f2[:, 0], f2_1)
    assert_allclose(d1[:, 0], d1_1)
    assert_allclose(d2[:, 0], d2_1)
    assert_allclose(d3[:, 0], d3_1)
    assert_allclose(e1[:, 0], e1_1)
    assert_allclose(e2[:, 0], e2_1)
    assert_allclose(e3[:, 0], e3_1)

    assert_allclose(f3[:, 0], np.array([0.092637, -0.245951, 0.938848]), rtol=1e-4)
    assert_allclose(g1[:, 0], np.array([0.939012, 0.073416, -0.07342]), rtol=1e-4)
    assert_allclose(g2[:, 0], np.array([0.055389, 1.004155, 0.257594]), rtol=1e-4)
    assert_allclose(g3[:, 0], np.array([0, 0, 1.065135]), rtol=1e-4)

    assert_allclose(f1[:, 1], f1_2)
    assert_allclose(f2[:, 1], f2_2)
    assert_allclose(d1[:, 1], d1_2)
    assert_allclose(d2[:, 1], d2_2)
    assert_allclose(d3[:, 1], d3_2)
    assert_allclose(e1[:, 1], e1_2)
    assert_allclose(e2[:, 1], e2_2)
    assert_allclose(e3[:, 1], e3_2)

    assert_allclose(f3[:, 1], np.array([-0.036618, -0.071019, 0.861604]), rtol=1e-4)
    assert_allclose(g1[:, 1], np.array([0.844391, 0.015353, 0.037152]), rtol=1e-4)
    assert_allclose(g2[:, 1], np.array([0.050808, 1.02131, 0.086342]), rtol=1e-4)
    assert_allclose(g3[:, 1], np.array([0, 0, 1.160625]), rtol=1e-4)


# test that vectors are calculated correctly

def test_basevectors_apex_delta():
    A = Apex(date=2000, refh=300)
    for lat in range(0, 90, 10):
        for lon in range(0, 360, 15):
            f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = A.basevectors_apex(lat, lon, 500)
            f = [np.append(f1, 0), np.append(f2, 0), f3]
            g = [g1, g2, g3]
            d = [d1, d2, d3]
            e = [e1, e2, e3]
            for i, j in [(i, j) for i in range(3) for j in range(3)]:
                delta = 1 if i == j else 0
                assert_allclose(np.sum(f[i]*g[j]), delta, rtol=0, atol=1e-5)
                assert_allclose(np.sum(d[i]*e[j]), delta, rtol=0, atol=1e-5)


def test_basevectors_apex_invalid_scalar():
    A = Apex(date=2000, refh=10000)
    with warnings.catch_warnings(record=True) as w:
        f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = A.basevectors_apex(0, 0, 0)
        A.basevectors_apex(0, 0, 0)
        assert len(w) == 2
        assert issubclass(w[-1].category, UserWarning)
        assert 'set to -9999 where' in str(w[-1].message)

    invalid = [-9999, -9999, -9999]
    assert not np.allclose(f1, invalid[:2])
    assert not np.allclose(f2, invalid[:2])
    assert_allclose(f3, invalid)
    assert_allclose(g1, invalid)
    assert_allclose(g2, invalid)
    assert_allclose(g3, invalid)
    assert_allclose(d1, invalid)
    assert_allclose(d2, invalid)
    assert_allclose(d3, invalid)
    assert_allclose(e1, invalid)
    assert_allclose(e2, invalid)
    assert_allclose(e3, invalid)


###============================================================================
### Test the get_apex() method
###============================================================================


def test_get_apex():
    A = Apex(date=2000, refh=300)
    assert_allclose(A.get_apex(10), 507.409702543805)
    assert_allclose(A.get_apex(60), 20313.026999999987)


def test_get_apex_invalid_lat():
    A = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        A.get_apex(91)
    with pytest.raises(ValueError):
        A.get_apex(-91)
    A.get_apex(90)
    A.get_apex(-90)

    assert_allclose(A.get_apex(90+1e-5), A.get_apex(90), rtol=0, atol=1e-8)


###============================================================================
### Test the set_epoch() method
###============================================================================


def test_set_epoch():
    A = Apex(date=2000.2, refh=300)
    assert_allclose(A.year, 2000.2)
    ret_2000_2_py = A._geo2apex(60, 15, 100)
    A.set_epoch(2000.8)
    assert_allclose(A.year, 2000.8)
    ret_2000_8_py = A._geo2apex(60, 15, 100)

    assert ret_2000_2_py != ret_2000_8_py

    fa.loadapxsh(A.datafile, 2000.2)
    ret_2000_2_apex = fa.apxg2all(60, 15, 100, 300, 0)[2:4]
    fa.loadapxsh(A.datafile, 2000.8)
    ret_2000_8_apex = fa.apxg2all(60, 15, 100, 300, 0)[2:4]

    assert ret_2000_2_apex != ret_2000_8_apex

    assert_allclose(ret_2000_2_py, ret_2000_2_apex)
    assert_allclose(ret_2000_8_py, ret_2000_8_apex)


###============================================================================
### Test the set_refh() method
###============================================================================


def test_set_refh():
    A = Apex(date=2000, refh=300)
    assert A.refh, 300
    ret_300 = A._geo2apex(60, 15, 100)
    A.set_refh(500)
    assert A.refh == 500
    ret_500 = A._geo2apex(60, 15, 100)

    assert_allclose(ret_300, fa.apxg2all(60, 15, 100, 300, 0)[2:4])
    assert_allclose(ret_500, fa.apxg2all(60, 15, 100, 500, 0)[2:4])


if __name__ == '__main__':
    pytest.main()
