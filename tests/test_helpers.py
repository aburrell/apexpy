# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals

import datetime as dt

import numpy as np
import pytest
from numpy.testing import assert_allclose

from apexpy import helpers


##############################################################################
# NOTE: whenever function outputs are tested against hard-coded numbers,     #
# the test results (numbers) were obtained by running the code that is       #
# tested. Therefore these tests below only check that nothing changes when   #
# refactoring etc., and not if the results are actually correct              #
##############################################################################


###============================================================================
### Test constants
###============================================================================


def test_constants():
    assert_allclose(helpers.d2r, np.pi/180)
    assert_allclose(helpers.r2d, 180/np.pi)


###============================================================================
### Test checklat
###============================================================================

def test_checklat_scalar():
    assert helpers.checklat(90) == 90
    assert helpers.checklat(0) == 0
    assert helpers.checklat(-90) == -90

    assert helpers.checklat(90+1e-5) == 90
    assert helpers.checklat(-90-1e-5) == -90

    assert type(helpers.checklat(0.)) == float
    assert type(helpers.checklat(0)) == int
    assert type(helpers.checklat(90+1e-5)) == int

    with pytest.raises(ValueError):
        helpers.checklat(90+1e-4)
    with pytest.raises(ValueError):
        helpers.checklat(-90-1e-4)


def test_checklat_message():
    with pytest.raises(ValueError) as excinfo:
        helpers.checklat(100)
    assert str(excinfo.value).startswith('lat must be in')
    with pytest.raises(ValueError) as excinfo:
        helpers.checklat(100, name='glat')
    assert str(excinfo.value).startswith('glat')


def test_checklat_array():
    assert_allclose(helpers.checklat([-90-1e-5, -90, 0, 90, 90+1e-5]), np.array([-90, -90, 0, 90, 90]), rtol=0, atol=1e-8)

    assert type(helpers.checklat([0])) == list
    assert type(helpers.checklat(np.array([0]))) == np.ndarray

    with pytest.raises(ValueError):
        helpers.checklat([-90-1e-4, -90, 0, 90, 90+1e-5])

    with pytest.raises(ValueError):
        helpers.checklat([-90-1e-5, -90, 0, 90, 90+1e-4])


###============================================================================
### Test getsinIm
###============================================================================

def test_getsinIm_scalar():
    assert_allclose(helpers.getsinIm(60), 0.96076892283052284)
    assert_allclose(helpers.getsinIm(10), 0.33257924500670238)
    assert type(helpers.getsinIm(60)) != np.ndarray


def test_getsinIm_1Darray():
    assert_allclose(helpers.getsinIm([60, 10]), [0.96076892283052284, 0.33257924500670238])


def test_getsinIm_2Darray():
    assert_allclose(helpers.getsinIm([[60, 10], [60, 10]]), [[0.96076892283052284, 0.33257924500670238],
                                                             [0.96076892283052284, 0.33257924500670238]])


###============================================================================
### Test getcosIm
###============================================================================


def test_getcosIm_scalar():
    assert_allclose(helpers.getcosIm(60), 0.27735009811261463)
    assert_allclose(helpers.getcosIm(10), 0.94307531289434765)
    assert type(helpers.getcosIm(60)) != np.ndarray


def test_getcosIm_1Darray():
    assert_allclose(helpers.getcosIm([60, 10]), [0.27735009811261463, 0.94307531289434765])


def test_getcosIm_2Darray():
    assert_allclose(helpers.getcosIm([[60, 10], [60, 10]]), [[0.27735009811261463, 0.94307531289434765],
                                                             [0.27735009811261463, 0.94307531289434765]])


###============================================================================
### Test toYearFraction
###============================================================================


def test_toYearFraction():
    assert_allclose(helpers.toYearFraction(dt.datetime(2001, 1, 1, 0, 0, 0)), 2001)
    assert_allclose(helpers.toYearFraction(dt.date(2001, 1, 1)), 2001)
    assert_allclose(helpers.toYearFraction(dt.datetime(2002, 1, 1, 0, 0, 0)), 2002)
    assert_allclose(helpers.toYearFraction(dt.datetime(2005, 2, 3, 4, 5, 6)), 2005.090877283105)
    assert_allclose(helpers.toYearFraction(dt.datetime(2005, 12, 11, 10, 9, 8)), 2005.943624682902)


###============================================================================
### Test gc2gdlat
###============================================================================


def test_gc2gdlat():
    assert_allclose(helpers.gc2gdlat(0), 0)
    assert_allclose(helpers.gc2gdlat(90), 90)
    assert_allclose(helpers.gc2gdlat(30), 30.166923849507356)
    assert_allclose(helpers.gc2gdlat(60), 60.166364190170931)


###============================================================================
### Test subsol
###============================================================================

def test_subsol():
    assert_allclose(helpers.subsol(dt.datetime(2005, 2, 3, 4, 5, 6)), (-16.505391672592904, 122.17768157084515))
    assert_allclose(helpers.subsol(dt.datetime(2010, 12, 11, 10, 9, 8)), (-23.001554595838947, 26.008999999955023))

    with pytest.raises(ValueError):
        helpers.subsol(dt.datetime(1600, 12, 31, 23, 59, 59))
    assert_allclose(helpers.subsol(dt.datetime(1601, 1, 1, 0, 0, 0)), (-23.06239721771427, -178.90131731228584))
    with pytest.raises(ValueError):
        helpers.subsol(dt.datetime(2101, 1, 1, 0, 0, 0))
    assert_allclose(helpers.subsol(dt.datetime(2100, 12, 31, 23, 59, 59)), (-23.021061422069053, -179.23129780639425))

if __name__ == '__main__':
    pytest.main()
