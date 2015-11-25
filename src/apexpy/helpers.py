# -*- coding: utf-8 -*-

'''This module contains helper functions used by :class:`~apexpy.Apex`.'''

from __future__ import division, print_function, absolute_import

import time
import datetime as dt

import numpy as np

# constants
d2r = np.pi/180
r2d = 180/np.pi


def checklat(lat, name='lat'):
    '''Makes sure the latitude is inside [-90, 90], clipping close values (tolerance 1e-4).

    Parameters
    ==========
    lat : array_like
        latitude
    name : str, optional
        parameter name to use in the exception message

    Returns
    =======
    lat : ndarray or float
        Same as input where values just outside the range have been
        clipped to [-90, 90]

    Raises
    ======
    ValueError
        if any values are too far outside the range [-90, 90]
    '''

    if np.all(np.float64(lat) >= -90) and np.all(np.float64(lat) <= 90):
        return lat

    if np.isscalar(lat):
        if lat > 90 and np.isclose(lat, 90, rtol=0, atol=1e-4):
            lat = 90
            return lat
        elif lat < -90 and np.isclose(lat, -90, rtol=0, atol=1e-4):
            lat = -90
            return lat
    else:
        lat = np.float64(lat)  # make sure we have an array, not list
        lat[(lat > 90) & (np.isclose(lat, 90, rtol=0, atol=1e-4))] = 90
        lat[(lat < -90) & (np.isclose(lat, -90, rtol=0, atol=1e-4))] = -90
        if np.all(lat >= -90) and np.all(lat <= 90):
            return lat

    # we haven't returned yet, so raise exception
    raise ValueError(name + ' must be in [-90, 90]')


def getsinIm(alat):
    '''Computes sinIm from modified apex latitude.

    Parameters
    ==========
    alat : array_like
        Modified apex latitude

    Returns
    =======
    sinIm : ndarray or float

    '''

    alat = np.float64(alat)

    return 2*np.sin(alat*d2r)/np.sqrt(4 - 3*np.cos(alat*d2r)**2)


def getcosIm(alat):
    '''Computes cosIm from modified apex latitude.

    Parameters
    ==========
    alat : array_like
        Modified apex latitude

    Returns
    =======
    cosIm : ndarray or float

    '''

    alat = np.float64(alat)

    return np.cos(alat*d2r)/np.sqrt(4 - 3*np.cos(alat*d2r)**2)


def toYearFraction(date):
    '''Converts :class:`datetime.date` or :class:`datetime.datetime` to decimal year.

    Parameters
    ==========
    date : :class:`datetime.date` or :class:`datetime.datetime`

    Returns
    =======
    year : float
        Decimal year

    Notes
    =====
    The algorithm is taken from http://stackoverflow.com/a/6451892/2978652

    '''

    def sinceEpoch(date):
        '''returns seconds since epoch'''
        return time.mktime(date.timetuple())
    year = date.year
    startOfThisYear = dt.datetime(year=year, month=1, day=1)
    startOfNextYear = dt.datetime(year=year+1, month=1, day=1)

    yearElapsed = sinceEpoch(date) - sinceEpoch(startOfThisYear)
    yearDuration = sinceEpoch(startOfNextYear) - sinceEpoch(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction


def gc2gdlat(gclat):
    '''Converts geocentric latitude to geodetic latitude using WGS84.

    Parameters
    ==========
    gclat : array_like
        Geocentric latitude

    Returns
    =======
    gdlat : ndarray or float
        Geodetic latitude

    '''
    WGS84_e2 = 0.006694379990141317  # WGS84 first eccentricity squared
    return np.rad2deg(-np.arctan(np.tan(np.deg2rad(gclat))/(WGS84_e2 - 1)))


def subsol(datetime):
    '''Finds subsolar geocentric latitude and longitude.

    Parameters
    ==========
    datetime : :class:`datetime.datetime`

    Returns
    =======
    sbsllat : float
        Latitude of subsolar point
    sbsllon : float
        Longitude of subsolar point

    Notes
    =====
    Based on formulas in Astronomical Almanac for the year 1996, p. C24.
    (U.S. Government Printing Office, 1994). Usable for years 1601-2100,
    inclusive. According to the Almanac, results are good to at least 0.01
    degree latitude and 0.025 degrees longitude between years 1950 and 2050.
    Accuracy for other years has not been tested. Every day is assumed to have
    exactly 86400 seconds; thus leap seconds that sometimes occur on December
    31 are ignored (their effect is below the accuracy threshold of the
    algorithm).

    After Fortran code by A. D. Richmond, NCAR. Translated from IDL
    by K. Laundal.

    '''

    from numpy import sin, cos, pi, arctan2, arcsin

    # convert to year, day of year and seconds since midnight
    year = datetime.year
    doy = datetime.timetuple().tm_yday
    ut = datetime.hour*3600 + datetime.minute*60 + datetime.second

    if not 1601 <= year <= 2100:
        raise ValueError('Year must be in [1601, 2100]')

    yr = year - 2000

    nleap = np.floor((year-1601)/4)
    nleap = nleap - 99
    if year <= 1900:
        ncent = np.floor((year-1601)/100)
        ncent = 3 - ncent
        nleap = nleap + ncent

    l0 = -79.549 + (-0.238699*(yr-4*nleap) + 3.08514e-2*nleap)

    g0 = -2.472 + (-0.2558905*(yr-4*nleap) - 3.79617e-2*nleap)

    # Days (including fraction) since 12 UT on January 1 of IYR:
    df = (ut/86400 - 1.5) + doy

    # Addition to Mean longitude of Sun since January 1 of IYR:
    lf = 0.9856474*df

    # Addition to Mean anomaly since January 1 of IYR:
    gf = 0.9856003*df

    # Mean longitude of Sun:
    l = l0 + lf

    # Mean anomaly:
    g = g0 + gf
    grad = g*pi/180

    # Ecliptic longitude:
    lmbda = l + 1.915*sin(grad) + 0.020*sin(2*grad)
    lmrad = lmbda*pi/180
    sinlm = sin(lmrad)

    # Days (including fraction) since 12 UT on January 1 of 2000:
    n = df + 365*yr + nleap

    # Obliquity of ecliptic:
    epsilon = 23.439 - 4e-7*n
    epsrad = epsilon*pi/180

    # Right ascension:
    alpha = arctan2(cos(epsrad)*sinlm, cos(lmrad)) * 180/pi

    # Declination:
    delta = arcsin(sin(epsrad)*sinlm) * 180/pi

    # Subsolar latitude:
    sslat = delta

    # Equation of time (degrees):
    etdeg = l - alpha
    nrot = round(etdeg/360)
    etdeg = etdeg - 360*nrot

    # Apparent time (degrees):
    aptime = ut/240 + etdeg    # Earth rotates one degree every 240 s.

    # Subsolar longitude:
    sslon = 180 - aptime
    nrot = round(sslon/360)
    sslon = sslon - 360*nrot

    return sslat, sslon
