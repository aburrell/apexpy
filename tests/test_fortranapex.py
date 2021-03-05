# -*- coding: utf-8 -*-

from numpy.testing import assert_allclose
import os
import pytest

import apexpy
from apexpy import fortranapex as fa

# ----------------------------------------------------------------------------
# NOTE: the test results (numbers) were obtained by running the code that is
# tested, therefore the tests below only check that nothing changes when
# refactoring etc., and not if the results are actually correct
# ----------------------------------------------------------------------------


def test_apxg2q():
    """Test fortran apex geographic to quasi-dipole
    """
    fa.loadapxsh(os.path.join(os.path.dirname(apexpy.__file__), 'apexsh.dat'),
                 2000)
    qlat, qlon, f1, f2, f = fa.apxg2q(60, 15, 100, 1)
    assert_allclose(qlat, 56.531288146972656)
    assert_allclose(qlon, 94.1068344116211)
    assert_allclose(f1, [1.079783, 0.10027137], rtol=1e-6)
    assert_allclose(f2, [-0.24546318, 0.90718889], rtol=1e-6)
    assert_allclose(f, 1.0041800737380981)


def test_apxg2all():
    fa.loadapxsh(os.path.join(os.path.dirname(apexpy.__file__), 'apexsh.dat'),
                 2000)
    qlat, qlon, mlat, mlon, f1, f2, f, d1, d2, d3, d, e1, e2, e3 = fa.apxg2all(
        60, 15, 100, 300, 1)
    assert_allclose(qlat, 56.531288146972656, rtol=1e-6)
    assert_allclose(qlon, 94.1068344116211, rtol=1e-6)
    assert_allclose(mlat, 55.94841766357422, rtol=1e-6)
    assert_allclose(mlon, 94.1068344116211, rtol=1e-6)
    assert_allclose(f1, [1.079783, 0.10027137], rtol=1e-6)
    assert_allclose(f2, [-0.24546318, 0.90718889], rtol=1e-6)
    assert_allclose(f, 1.0041800737380981, rtol=1e-6)
    assert_allclose(d1, [0.9495701, 0.25693053, 0.09049474], rtol=1e-6)
    assert_allclose(d2, [0.10011087, -1.0780545, -0.33892432], rtol=1e-6)
    assert_allclose(d3, [0.00865356, 0.27327004, -0.8666646], rtol=1e-6)
    assert_allclose(d, 1.100391149520874, rtol=1e-6)
    assert_allclose(e1, [1.0269295, 0.08382964, 0.03668632], rtol=1e-6)
    assert_allclose(e2, [0.24740215, -0.82374191, -0.25726584], rtol=1e-6)
    assert_allclose(e3, [0.01047826,  0.33089194, -1.04941], rtol=1e-6)


def test_apxq2g():
    """ Test fortran quasi-dipole to geographic
    """
    fa.loadapxsh(os.path.join(os.path.dirname(apexpy.__file__), 'apexsh.dat'),
                 2000)
    glat, glon, error = fa.apxq2g(60, 15, 100, 1e-2)
    assert_allclose([glat, glon, error],
                    [50.97946548461914, -66.16902923583984,
                     0.00010020843910751864], atol=1e-6)


def test_g2q2d():
    fa.loadapxsh(os.path.join(os.path.dirname(apexpy.__file__), 'apexsh.dat'),
                 2000)
    for lat in [0, 30, 60, 89]:
        for lon in [-179, -90, 0, 90, 179]:
            qlat, qlon, _, _, _ = fa.apxg2q(lat, lon, 100, 0)
            glat, glon, _ = fa.apxq2g(qlat, qlon, 100, 1e-10)
            assert_allclose(glat, lat, atol=0.01)
            assert_allclose(glon, lon, atol=0.01)


def test_apxq2g_lowprecision():
    fa.loadapxsh(os.path.join(os.path.dirname(apexpy.__file__), 'apexsh.dat'),
                 2000)
    glat, glon, error = fa.apxq2g(60, 15, 100, -1)
    assert_allclose(glat, 51.00891876220703)
    assert_allclose(glon, -66.11973571777344)
    assert_allclose(error, -9999.0)


if __name__ == '__main__':
    pytest.main()
