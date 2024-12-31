# -*- coding: utf-8 -*-
"""This module contains helper functions used by :class:`~apexpy.Apex`."""

import datetime as dt
import numpy as np
import time


def set_array_float(in_val):
    """Set array data type to float.

    Parameters
    ----------
    in_val : any
        Input value, only modified if it is a np.ndarray

    Returns
    -------
    out_val : any
        Output value, if `in_val` was an array, `out_val` will be an array of
        type `np.float64`.

    """

    if isinstance(in_val, np.ndarray):
        out_val = in_val.astype(np.float64)
    else:
        out_val = in_val

    return out_val


def checklat(lat, name='lat'):
    """Makes sure the latitude is inside [-90, 90], clipping close values
    (tolerance 1e-4).

    Parameters
    ----------
    lat : array-like
        latitude
    name : str, optional
        parameter name to use in the exception message

    Returns
    -------
    lat : ndarray or float
        Same as input where values just outside the range have been
        clipped to [-90, 90]

    Raises
    ------
    ValueError
        if any values are too far outside the range [-90, 90]

    """
    if np.any(np.abs(lat) > 90 + 1e-5):
        raise ValueError('{:} must be in [-90, 90]'.format(name))

    return np.clip(lat, -90.0, 90.0)


def getsinIm(alat):
    """Computes sinIm from modified apex latitude.

    Parameters
    ----------
    alat : array-like
        Modified apex latitude

    Returns
    -------
    sinIm : ndarray or float

    """

    alat = np.float64(alat)

    return 2 * np.sin(np.radians(alat)) / np.sqrt(4 - 3
                                                  * np.cos(np.radians(alat))**2)


def getcosIm(alat):
    """Computes cosIm from modified apex latitude.

    Parameters
    ----------
    alat : array-like
        Modified apex latitude

    Returns
    -------
    cosIm : ndarray or float

    """

    alat = np.float64(alat)

    return np.cos(np.radians(alat)) / np.sqrt(4 - 3
                                              * np.cos(np.radians(alat))**2)


def toYearFraction(date):
    """Converts :class:`datetime.date` or :class:`datetime.datetime` to decimal
    year.

    Parameters
    ----------
    date : :class:`datetime.date` or :class:`datetime.datetime`
        Input date or datetime object

    Returns
    -------
    year : float
        Decimal year

    Notes
    -----
    The algorithm is taken from http://stackoverflow.com/a/6451892/2978652

    """

    def sinceEpoch(date):
        """returns seconds since epoch"""
        return time.mktime(date.timetuple())

    year = date.year
    startOfThisYear = dt.datetime(year=year, month=1, day=1)
    startOfNextYear = dt.datetime(year=year + 1, month=1, day=1)

    yearElapsed = sinceEpoch(date) - sinceEpoch(startOfThisYear)
    yearDuration = sinceEpoch(startOfNextYear) - sinceEpoch(startOfThisYear)
    fraction = yearElapsed / yearDuration

    year += fraction

    return year


def gc2gdlat(gclat):
    """Converts geocentric latitude to geodetic latitude using WGS84.

    Parameters
    ---------
    gclat : array-like
        Geocentric latitude

    Returns
    -------
    gdlat : ndarray or float
        Geodetic latitude

    """
    WGS84_e2 = 0.006694379990141317  # WGS84 first eccentricity squared
    return np.rad2deg(-np.arctan(np.tan(np.deg2rad(gclat)) / (WGS84_e2 - 1)))


def subsol(datetime):
    """Finds subsolar geocentric latitude and longitude.

    Parameters
    ----------
    datetime : :class:`datetime.datetime` or :class:`numpy.ndarray[datetime64]`
        Date and time in UTC (naive objects are treated as UTC)

    Returns
    -------
    sbsllat : float
        Latitude of subsolar point
    sbsllon : float
        Longitude of subsolar point

    Notes
    -----
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

    """
    # Convert to year, day of year and seconds since midnight
    if isinstance(datetime, dt.datetime):
        year = np.asanyarray([datetime.year])
        doy = np.asanyarray([datetime.timetuple().tm_yday])
        ut = np.asanyarray([datetime.hour * 3600 + datetime.minute * 60
                            + datetime.second + datetime.microsecond / 1.0e6])
    elif isinstance(datetime, np.ndarray):
        # This conversion works for datetime of wrong precision or unit epoch
        times = datetime.astype('datetime64[us]')
        year_floor = times.astype('datetime64[Y]')
        day_floor = times.astype('datetime64[D]')
        year = year_floor.astype(int) + 1970
        doy = (day_floor - year_floor).astype(int) + 1
        ut = (times.astype('datetime64[us]') - day_floor).astype(float)
        ut /= 1e6
    else:
        raise ValueError("input must be datetime.datetime or numpy array")

    if not (np.all(1601 <= year) and np.all(year <= 2100)):
        raise ValueError('Year must be in [1601, 2100]')

    yr = year - 2000

    nleap = np.floor((year - 1601.0) / 4.0).astype(int)
    nleap -= 99
    mask_1900 = year <= 1900
    if np.any(mask_1900):
        ncent = np.floor((year[mask_1900] - 1601.0) / 100.0).astype(int)
        ncent = 3 - ncent
        nleap[mask_1900] = nleap[mask_1900] + ncent

    l0 = -79.549 + (-0.238699 * (yr - 4.0 * nleap) + 3.08514e-2 * nleap)
    g0 = -2.472 + (-0.2558905 * (yr - 4.0 * nleap) - 3.79617e-2 * nleap)

    # Days (including fraction) since 12 UT on January 1 of IYR:
    df = (ut / 86400.0 - 1.5) + doy

    # Mean longitude of Sun:
    lmean = l0 + 0.9856474 * df

    # Mean anomaly in radians:
    grad = np.radians(g0 + 0.9856003 * df)

    # Ecliptic longitude:
    lmrad = np.radians(lmean + 1.915 * np.sin(grad)
                       + 0.020 * np.sin(2.0 * grad))
    sinlm = np.sin(lmrad)

    # Obliquity of ecliptic in radians:
    epsrad = np.radians(23.439 - 4e-7 * (df + 365 * yr + nleap))

    # Right ascension:
    alpha = np.degrees(np.arctan2(np.cos(epsrad) * sinlm, np.cos(lmrad)))

    # Declination, which is also the subsolar latitude:
    sslat = np.degrees(np.arcsin(np.sin(epsrad) * sinlm))

    # Equation of time (degrees):
    etdeg = lmean - alpha
    nrot = np.round(etdeg / 360.0)
    etdeg = etdeg - 360.0 * nrot

    # Subsolar longitude calculation. Earth rotates one degree every 240 s.
    sslon = 180.0 - (ut / 240.0 + etdeg)
    nrot = np.round(sslon / 360.0)
    sslon = sslon - 360.0 * nrot

    # Return a single value from the output if the input was a single value
    if isinstance(datetime, dt.datetime):
        return sslat[0], sslon[0]
    return sslat, sslon
