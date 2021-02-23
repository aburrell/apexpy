# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import

import datetime as dt
import numpy as np
import os
import warnings

from . import helpers

# Below try..catch required for autodoc to work on readthedocs
try:
    from . import fortranapex as fa
except ImportError as err:
    warnings.warn("".join(["fortranapex module could not be imported, so ",
                           "apexpy probably won't work.  Make sure you have ",
                           "a gfortran compiler. Wheels installation ",
                           "assumes your compiler lives in /opt/local/bin"]))
    raise err

# make sure invalid warnings are always shown
warnings.filterwarnings('always', message='.*set to NaN where*',
                        module='apexpy.apex')


class ApexHeightError(ValueError):
    """Specialized error type definition
    """
    pass


class Apex(object):
    """Performs coordinate conversions, field-line mapping and base vector
    calculations.

    Parameters
    ----------
    date : float, :class:`dt.date`, or :class:`dt.datetime`, optional
        Determines which IGRF coefficients are used in conversions. Uses
        current date as default.  If float, use decimal year.
    refh : float, optional
        Reference height in km for apex coordinates (the field lines are mapped
        to this height)
    datafile : str, optional
        Path to custom coefficient file

    Attributes
    ----------
    year : float
        Decimal year used for the IGRF model
    refh : float
        Reference height in km for apex coordinates
    datafile : str
        Path to coefficient file

    Notes
    -----
    The calculations use IGRF-13 with coefficients from 1900 to 2025 [1]_.

    The geodetic reference ellipsoid is WGS84.

    References
    ----------

    .. [1] Thébault, E. et al. (2015), International Geomagnetic Reference
           Field: the 12th generation, Earth, Planets and Space, 67(1), 79,
           :doi:`10.1186/s40623-015-0228-9`.

    """

    def __init__(self, date=None, refh=0, datafile=None, fortranlib=None):

        if datafile is None:
            datafile = os.path.join(os.path.dirname(__file__), 'apexsh.dat')

        if fortranlib is None:
            fortranlib = fa.__file__

        self.RE = 6371.009  # mean Earth radius
        self.set_refh(refh)  # reference height

        if date is None:
            self.year = helpers.toYearFraction(dt.datetime.utcnow())
        else:
            try:
                # Convert date/datetime object to decimal year
                self.year = helpers.toYearFraction(date)
            except AttributeError:
                # Failed while finding datetime attribute, so
                # date is probably an int or float; use directly
                self.year = date

        if not os.path.isfile(datafile):
            raise IOError('Data file does not exist: {}'.format(datafile))

        if not os.path.isfile(fortranlib):
            raise IOError('Fortran library does not exist: {}'.format(
                fortranlib))

        self.datafile = datafile
        self.fortranlib = fortranlib

        self.set_epoch(self.year)

        # vectorize fortran functions
        self._geo2qd = np.frompyfunc(
            lambda glat, glon, height: fa.apxg2q(glat, (glon + 180) % 360 - 180,
                                                 height, 0)[:2], 3, 2)
        self._geo2apex = np.frompyfunc(
            lambda glat, glon, height: fa.apxg2all(glat, (glon + 180) % 360
                                                   - 180, height, self.refh,
                                                   0)[2:4], 3, 2)
        self._geo2apexall = np.frompyfunc(
            lambda glat, glon, height: fa.apxg2all(glat, (glon + 180) % 360
                                                   - 180, height, self.refh,
                                                   1), 3, 14)
        self._qd2geo = np.frompyfunc(
            lambda qlat, qlon, height, precision: fa.apxq2g(qlat, (qlon + 180)
                                                            % 360 - 180, height,
                                                            precision), 4, 3)
        self._basevec = np.frompyfunc(
            lambda glat, glon, height: fa.apxg2q(glat, (glon + 180) % 360 - 180,
                                                 height, 1)[2:4], 3, 2)

        # vectorize other nonvectorized functions
        self._apex2qd = np.frompyfunc(self._apex2qd_nonvectorized, 3, 2)
        self._qd2apex = np.frompyfunc(self._qd2apex_nonvectorized, 3, 2)
        self._get_babs = np.frompyfunc(self._get_babs_nonvectorized, 3, 1)

    def convert(self, lat, lon, source, dest, height=0, datetime=None,
                precision=1e-10, ssheight=50 * 6371):
        """Converts between geodetic, modified apex, quasi-dipole and MLT.

        Parameters
        ----------
        lat : array_like
            Latitude
        lon : array_like
            Longitude/MLT
        source : {'geo', 'apex', 'qd', 'mlt'}
            Input coordinate system
        dest : {'geo', 'apex', 'qd', 'mlt'}
            Output coordinate system
        height : array_like, optional
            Altitude in km
        datetime : :class:`datetime.datetime`
            Date and time for MLT conversions (required for MLT conversions)
        precision : float, optional
            Precision of output (degrees) when converting to geo. A negative
            value of this argument produces a low-precision calculation of
            geodetic lat/lon based only on their spherical harmonic
            representation.
            A positive value causes the underlying Fortran routine to iterate
            until feeding the output geo lat/lon into geo2qd (APXG2Q) reproduces
            the input QD lat/lon to within the specified precision (all
            coordinates being converted to geo are converted to QD first and
            passed through APXG2Q).
        ssheight : float, optional
            Altitude in km to use for converting the subsolar point from
            geographic to magnetic coordinates. A high altitude is used
            to ensure the subsolar point is mapped to high latitudes, which
            prevents the South-Atlantic Anomaly (SAA) from influencing the MLT.

        Returns
        -------
        lat : ndarray or float
            Converted latitude (if converting to MLT, output latitude is apex)
        lat : ndarray or float
            Converted longitude/MLT

        """

        if datetime is None and ('mlt' in [source, dest]):
            raise ValueError('datetime must be given for MLT calculations')

        lat = helpers.checklat(lat)

        if source == dest:
            return lat, lon
        # from geo
        elif source == 'geo' and dest == 'apex':
            lat, lon = self.geo2apex(lat, lon, height)
        elif source == 'geo' and dest == 'qd':
            lat, lon = self.geo2qd(lat, lon, height)
        elif source == 'geo' and dest == 'mlt':
            lat, lon = self.geo2apex(lat, lon, height)
            lon = self.mlon2mlt(lon, datetime, ssheight=ssheight)
        # from apex
        elif source == 'apex' and dest == 'geo':
            lat, lon, _ = self.apex2geo(lat, lon, height, precision=precision)
        elif source == 'apex' and dest == 'qd':
            lat, lon = self.apex2qd(lat, lon, height=height)
        elif source == 'apex' and dest == 'mlt':
            lon = self.mlon2mlt(lon, datetime, ssheight=ssheight)
        # from qd
        elif source == 'qd' and dest == 'geo':
            lat, lon, _ = self.qd2geo(lat, lon, height, precision=precision)
        elif source == 'qd' and dest == 'apex':
            lat, lon = self.qd2apex(lat, lon, height=height)
        elif source == 'qd' and dest == 'mlt':
            lat, lon = self.qd2apex(lat, lon, height=height)
            lon = self.mlon2mlt(lon, datetime, ssheight=ssheight)
        # from mlt (input latitude assumed apex)
        elif source == 'mlt' and dest == 'geo':
            lon = self.mlt2mlon(lon, datetime, ssheight=ssheight)
            lat, lon, _ = self.apex2geo(lat, lon, height, precision=precision)
        elif source == 'mlt' and dest == 'apex':
            lon = self.mlt2mlon(lon, datetime, ssheight=ssheight)
        elif source == 'mlt' and dest == 'qd':
            lon = self.mlt2mlon(lon, datetime, ssheight=ssheight)
            lat, lon = self.apex2qd(lat, lon, height=height)
        # no other transformations are implemented
        else:
            estr = 'Unknown coordinate transformation: '
            estr += '{} -> {}'.format(source, dest)
            raise NotImplementedError(estr)

        return lat, lon

    def geo2apex(self, glat, glon, height):
        """Converts geodetic to modified apex coordinates.

        Parameters
        ----------
        glat : array_like
            Geodetic latitude
        glon : array_like
            Geodetic longitude
        height : array_like
            Altitude in km

        Returns
        -------
        alat : ndarray or float
            Modified apex latitude
        alon : ndarray or float
            Modified apex longitude

        """

        glat = helpers.checklat(glat, name='glat')

        alat, alon = self._geo2apex(glat, glon, height)

        if np.any(alat == -9999):
            warnings.warn('Apex latitude set to NaN where undefined '
                          '(apex height may be < reference height)')
            if np.isscalar(alat):
                alat = np.nan
            else:
                alat[alat == -9999] = np.nan

        # if array is returned, dtype is object, so convert to float
        return np.float64(alat), np.float64(alon)

    def apex2geo(self, alat, alon, height, precision=1e-10):
        """Converts modified apex to geodetic coordinates.

        Parameters
        ----------
        alat : array_like
            Modified apex latitude
        alon : array_like
            Modified apex longitude
        height : array_like
            Altitude in km
        precision : float, optional
            Precision of output (degrees). A negative value of this argument
            produces a low-precision calculation of geodetic lat/lon based only
            on their spherical harmonic representation. A positive value causes
            the underlying Fortran routine to iterate until feeding the output
            geo lat/lon into geo2qd (APXG2Q) reproduces the input QD lat/lon to
            within the specified precision.

        Returns
        -------
        glat : ndarray or float
            Geodetic latitude
        glon : ndarray or float
            Geodetic longitude
        error : ndarray or float
            The angular difference (degrees) between the input QD coordinates
            and the qlat/qlon produced by feeding the output glat and glon
            into geo2qd (APXG2Q)

        """

        alat = helpers.checklat(alat, name='alat')

        qlat, qlon = self.apex2qd(alat, alon, height=height)
        glat, glon, error = self.qd2geo(qlat, qlon, height, precision=precision)

        return glat, glon, error

    def geo2qd(self, glat, glon, height):
        """Converts geodetic to quasi-dipole coordinates.

        Parameters
        ----------
        glat : array_like
            Geodetic latitude
        glon : array_like
            Geodetic longitude
        height : array_like
            Altitude in km

        Returns
        -------
        qlat : ndarray or float
            Quasi-dipole latitude
        qlon : ndarray or float
            Quasi-dipole longitude

        """

        glat = helpers.checklat(glat, name='glat')

        qlat, qlon = self._geo2qd(glat, glon, height)

        # if array is returned, dtype is object, so convert to float
        return np.float64(qlat), np.float64(qlon)

    def qd2geo(self, qlat, qlon, height, precision=1e-10):
        """Converts quasi-dipole to geodetic coordinates.

        Parameters
        ----------
        qlat : array_like
            Quasi-dipole latitude
        qlon : array_like
            Quasi-dipole longitude
        height : array_like
            Altitude in km
        precision : float, optional
            Precision of output (degrees). A negative value of this argument
            produces a low-precision calculation of geodetic lat/lon based only
            on their spherical harmonic representation. A positive value causes
            the underlying Fortran routine to iterate until feeding the output
            geo lat/lon into geo2qd (APXG2Q) reproduces the input QD lat/lon to
            within the specified precision.

        Returns
        -------
        glat : ndarray or float
            Geodetic latitude
        glon : ndarray or float
            Geodetic longitude
        error : ndarray or float
            The angular difference (degrees) between the input QD coordinates
            and the qlat/qlon produced by feeding the output glat and glon
            into geo2qd (APXG2Q)

        """

        qlat = helpers.checklat(qlat, name='qlat')

        glat, glon, error = self._qd2geo(qlat, qlon, height, precision)

        # if array is returned, dtype is object, so convert to float
        return np.float64(glat), np.float64(glon), np.float64(error)

    def _apex2qd_nonvectorized(self, alat, alon, height):
        """Convert from apex to quasi-dipole (not-vectorised)

        Parameters
        -----------
        alat : (float)
            Apex latitude in degrees
        alon : (float)
            Apex longitude in degrees
        height : (float)
            Height in km

        Returns
        ---------
        qlat : (float)
            Quasi-dipole latitude in degrees
        qlon : (float)
            Quasi-diplole longitude in degrees
        """

        alat = helpers.checklat(alat, name='alat')

        # convert modified apex to quasi-dipole:
        qlon = alon

        # apex height
        hA = self.get_apex(alat)

        if hA < height:
            if np.isclose(hA, height, rtol=0, atol=1e-5):
                # allow for values that are close
                hA = height
            else:
                estr = 'height {:.3g} is > apex height'.format(np.max(height))\
                       + ' {:.3g} for alat {:.3g}'.format(hA, alat)
                raise ApexHeightError(estr)

        salat = np.sign(alat) if alat != 0 else 1
        qlat = salat * np.degrees(np.arccos(np.sqrt((self.RE + height) /
                                                    (self.RE + hA))))

        return qlat, qlon

    def apex2qd(self, alat, alon, height):
        """Converts modified apex to quasi-dipole coordinates.

        Parameters
        ----------
        alat : array_like
            Modified apex latitude
        alon : array_like
            Modified apex longitude
        height : array_like
            Altitude in km

        Returns
        -------
        qlat : ndarray or float
            Quasi-dipole latitude
        qlon : ndarray or float
            Quasi-dipole longitude

        Raises
        ------
        ApexHeightError
            if `height` > apex height

        """

        qlat, qlon = self._apex2qd(alat, alon, height)

        # if array is returned, the dtype is object, so convert to float
        return np.float64(qlat), np.float64(qlon)

    def _qd2apex_nonvectorized(self, qlat, qlon, height):

        qlat = helpers.checklat(qlat, name='qlat')

        alon = qlon
        hA = self.get_apex(qlat, height)  # apex height

        if hA < self.refh:
            if np.isclose(hA, self.refh, rtol=0, atol=1e-5):
                # allow for values that are close
                hA = self.refh
            else:
                estr = 'apex height ({:.3g}) is < reference height '.format(hA)
                estr += '({:.3g}) for qlat {:.3g}'.format(self.refh, qlat)
                raise ApexHeightError(estr)

        sqlat = np.sign(qlat) if qlat != 0 else 1
        alat = sqlat * np.degrees(np.arccos(np.sqrt((self.RE + self.refh) /
                                                    (self.RE + hA))))

        return alat, alon

    def qd2apex(self, qlat, qlon, height):
        """Converts quasi-dipole to modified apex coordinates.

        Parameters
        ----------
        qlat : array_like
            Quasi-dipole latitude
        qlon : array_like
            Quasi-dipole longitude
        height : array_like
            Altitude in km

        Returns
        -------
        alat : ndarray or float
            Modified apex latitude
        alon : ndarray or float
            Modified apex longitude

        Raises
        ------
        ApexHeightError
            if apex height < reference height

        """

        alat, alon = self._qd2apex(qlat, qlon, height)

        # if array is returned, the dtype is object, so convert to float
        return np.float64(alat), np.float64(alon)

    def mlon2mlt(self, mlon, datetime, ssheight=50 * 6371):
        """Computes the magnetic local time at the specified magnetic longitude
        and UT.

        Parameters
        ----------
        mlon : array_like
            Magnetic longitude (apex and quasi-dipole longitude are always
            equal)
        datetime : :class:`datetime.datetime`
            Date and time
        ssheight : float, optional
            Altitude in km to use for converting the subsolar point from
            geographic to magnetic coordinates. A high altitude is used
            to ensure the subsolar point is mapped to high latitudes, which
            prevents the South-Atlantic Anomaly (SAA) from influencing the MLT.

        Returns
        -------
        mlt : ndarray or float
            Magnetic local time [0, 24)

        Notes
        -----
        To compute the MLT, we find the apex longitude of the subsolar point at
        the given time. Then the MLT of the given point will be computed from
        the separation in magnetic longitude from this point (1 hour = 15
        degrees).

        """
        ssglat, ssglon = helpers.subsol(datetime)
        ssalat, ssalon = self.geo2apex(ssglat, ssglon, ssheight)

        # np.float64 will ensure lists are converted to arrays
        return (180 + np.float64(mlon) - ssalon) / 15 % 24

    def mlt2mlon(self, mlt, datetime, ssheight=50 * 6371):
        """Computes the magnetic longitude at the specified magnetic local time
        and UT.

        Parameters
        ----------
        mlt : array_like
            Magnetic local time
        datetime : :class:`datetime.datetime`
            Date and time
        ssheight : float, optional
            Altitude in km to use for converting the subsolar point from
            geographic to magnetic coordinates. A high altitude is used
            to ensure the subsolar point is mapped to high latitudes, which
            prevents the South-Atlantic Anomaly (SAA) from influencing the MLT.

        Returns
        -------
        mlon : ndarray or float
            Magnetic longitude [0, 360) (apex and quasi-dipole longitude are
            always equal)

        Notes
        -----
        To compute the magnetic longitude, we find the apex longitude of the
        subsolar point at the given time. Then the magnetic longitude of the
        given point will be computed from the separation in magnetic local time
        from this point (1 hour = 15 degrees).
        """

        ssglat, ssglon = helpers.subsol(datetime)
        ssalat, ssalon = self.geo2apex(ssglat, ssglon, ssheight)

        # np.float64 will ensure lists are converted to arrays
        return (15 * np.float64(mlt) - 180 + ssalon + 360) % 360

    def map_to_height(self, glat, glon, height, newheight, conjugate=False,
                      precision=1e-10):
        """Performs mapping of points along the magnetic field to the closest
        or conjugate hemisphere.

        Parameters
        ----------
        glat : array_like
            Geodetic latitude
        glon : array_like
            Geodetic longitude
        height : array_like
            Source altitude in km
        newheight : array_like
            Destination altitude in km
        conjugate : bool, optional
            Map to `newheight` in the conjugate hemisphere instead of the
            closest hemisphere
        precision : float, optional
            Precision of output (degrees). A negative value of this argument
            produces a low-precision calculation of geodetic lat/lon based only
            on their spherical harmonic representation. A positive value causes
            the underlying Fortran routine to iterate until feeding the output
            geo lat/lon into geo2qd (APXG2Q) reproduces the input QD lat/lon to
            within the specified precision.

        Returns
        -------
        newglat : ndarray or float
            Geodetic latitude of mapped point
        newglon : ndarray or float
            Geodetic longitude of mapped point
        error : ndarray or float
            The angular difference (degrees) between the input QD coordinates
            and the qlat/qlon produced by feeding the output glat and glon
            into geo2qd (APXG2Q)

        Notes
        -----
        The mapping is done by converting glat/glon/height to modified apex
        lat/lon, and converting back to geographic using newheight (if
        conjugate, use negative apex latitude when converting back)

        """

        alat, alon = self.geo2apex(glat, glon, height)
        if conjugate:
            alat = -alat
        try:
            newglat, newglon, error = self.apex2geo(alat, alon, newheight,
                                                    precision=precision)
        except ApexHeightError:
            raise ApexHeightError("newheight is > apex height")

        return newglat, newglon, error

    def _map_EV_to_height(self, alat, alon, height, newheight, X, EV):

        # make sure X is array of correct shape
        if (not (np.ndim(X) == 1 and np.size(X) == 3) and not (
                np.ndim(X) == 2 and np.shape(X)[0] == 3)):
            # raise ValueError because if passing e.g. a (6,) ndarray the
            # reshape below will work even though the input is invalid
            raise ValueError(EV + ' must be (3, N) or (3,) ndarray')
        X = np.reshape(X, (3, np.size(X) // 3))

        _, _, _, _, _, _, d1, d2, _, e1, e2, _ = self.basevectors_apex(
            alat, alon, height, coords='apex')

        if EV == 'E':
            v1 = e1
            v2 = e2
        else:
            v1 = d1
            v2 = d2

        # make sure v1 and v2 have shape (3, N)
        v1 = np.reshape(v1, (3, v1.size // 3))
        v2 = np.reshape(v2, (3, v2.size // 3))

        X1 = np.sum(X * v1, axis=0)  # E dot e1 or V dot d1
        X2 = np.sum(X * v2, axis=0)  # E dot e2 or V dot d2

        _, _, _, _, _, _, d1, d2, _, e1, e2, _ = self.basevectors_apex(
            alat, alon, newheight, coords='apex')

        if EV == 'E':
            v1 = d1
            v2 = d2
        else:
            v1 = e1
            v2 = e2

        # make sure v1 and v2 have shape (3, N)
        v1 = np.reshape(v1, (3, v1.size // 3))
        v2 = np.reshape(v2, (3, v2.size // 3))

        X_mapped = X1[np.newaxis, :] * v1 + X2[np.newaxis, :] * v2

        return np.squeeze(X_mapped)

    def map_E_to_height(self, alat, alon, height, newheight, E):
        """Performs mapping of electric field along the magnetic field.

        It is assumed that the electric field is perpendicular to B.

        Parameters
        ----------
        alat : (N,) array_like or float
            Modified apex latitude
        alon : (N,) array_like or float
            Modified apex longitude
        height : (N,) array_like or float
            Source altitude in km
        newheight : (N,) array_like or float
            Destination altitude in km
        E : (3,) or (3, N) array_like
            Electric field (at `alat`, `alon`, `height`) in geodetic east,
            north, and up components

        Returns
        -------
        E : (3, N) or (3,) ndarray
            The electric field at `newheight` (geodetic east, north, and up
            components)

        """
        return self._map_EV_to_height(alat, alon, height, newheight, E, 'E')

    def map_V_to_height(self, alat, alon, height, newheight, V):
        """Performs mapping of electric drift velocity along the magnetic field.

        It is assumed that the electric field is perpendicular to B.

        Parameters
        ----------
        alat : (N,) array_like or float
            Modified apex latitude
        alon : (N,) array_like or float
            Modified apex longitude
        height : (N,) array_like or float
            Source altitude in km
        newheight : (N,) array_like or float
            Destination altitude in km
        V : (3,) or (3, N) array_like
            Electric drift velocity (at `alat`, `alon`, `height`) in geodetic
            east, north, and up components

        Returns
        -------
        V : (3, N) or (3,) ndarray
            The electric drift velocity at `newheight` (geodetic east, north,
            and up components)

        """

        return self._map_EV_to_height(alat, alon, height, newheight, V, 'V')

    def basevectors_qd(self, lat, lon, height, coords='geo', precision=1e-10):
        """Returns quasi-dipole base vectors f1 and f2 at the specified
        coordinates.

        The vectors are described by Richmond [1995] [2]_ and
        Emmert et al. [2010] [3]_.  The vector components are geodetic east and
        north.

        Parameters
        ----------
        lat : (N,) array_like or float
            Latitude
        lon : (N,) array_like or float
            Longitude
        height : (N,) array_like or float
            Altitude in km
        coords : {'geo', 'apex', 'qd'}, optional
            Input coordinate system
        precision : float, optional
            Precision of output (degrees) when converting to geo. A negative
            value of this argument produces a low-precision calculation of
            geodetic lat/lon based only on their spherical harmonic
            representation.
            A positive value causes the underlying Fortran routine to iterate
            until feeding the output geo lat/lon into geo2qd (APXG2Q) reproduces
            the input QD lat/lon to within the specified precision (all
            coordinates being converted to geo are converted to QD first and
            passed through APXG2Q).

        Returns
        -------
        f1 : (2, N) or (2,) ndarray
        f2 : (2, N) or (2,) ndarray

        References
        ----------
        .. [2] Richmond, A. D. (1995), Ionospheric Electrodynamics Using
               Magnetic Apex Coordinates, Journal of geomagnetism and
               geoelectricity, 47(2), 191–212, :doi:`10.5636/jgg.47.191`.

        .. [3] Emmert, J. T., A. D. Richmond, and D. P. Drob (2010),
               A computationally compact representation of Magnetic-Apex
               and Quasi-Dipole coordinates with smooth base vectors,
               J. Geophys. Res., 115(A8), A08322, :doi:`10.1029/2010JA015326`.

        """

        glat, glon = self.convert(lat, lon, coords, 'geo', height=height,
                                  precision=precision)

        f1, f2 = self._basevec(glat, glon, height)

        # if inputs are not scalar, each vector is an array of arrays,
        # so reshape to a single array
        if f1.dtype == object:
            f1 = np.vstack(f1).T
            f2 = np.vstack(f2).T

        return f1, f2

    def basevectors_apex(self, lat, lon, height, coords='geo', precision=1e-10):
        """Returns base vectors in quasi-dipole and apex coordinates.

        The vectors are described by Richmond [1995] [4]_ and
        Emmert et al. [2010] [5]_.  The vector components are geodetic east,
        north, and up (only east and north for `f1` and `f2`).

        Parameters
        ----------
        lat : (N,) array_like or float
            Latitude
        lon : (N,) array_like or float
            Longitude
        height : (N,) array_like or float
            Altitude in km
        coords : {'geo', 'apex', 'qd'}, optional
            Input coordinate system
        precision : float, optional
            Precision of output (degrees) when converting to geo. A negative
            value of this argument produces a low-precision calculation of
            geodetic lat/lon based only on their spherical harmonic
            representation.
            A positive value causes the underlying Fortran routine to iterate
            until feeding the output geo lat/lon into geo2qd (APXG2Q) reproduces
            the input QD lat/lon to within the specified precision (all
            coordinates being converted to geo are converted to QD first and
            passed through APXG2Q).

        Returns
        -------
        f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 : (3, N) or (3,) ndarray

        Notes
        -----
        `f3`, `g1`, `g2`, and `g3` are not part of the Fortran code
        by Emmert et al. [2010] [5]_. They are calculated by this
        Python library according to the following equations in
        Richmond [1995] [4]_:

        * `g1`: Eqn. 6.3
        * `g2`: Eqn. 6.4
        * `g3`: Eqn. 6.5
        * `f3`: Eqn. 6.8

        References
        ----------

        .. [4] Richmond, A. D. (1995), Ionospheric Electrodynamics Using
               Magnetic Apex Coordinates, Journal of geomagnetism and
               geoelectricity, 47(2), 191–212, :doi:`10.5636/jgg.47.191`.

        .. [5] Emmert, J. T., A. D. Richmond, and D. P. Drob (2010),
               A computationally compact representation of Magnetic-Apex
               and Quasi-Dipole coordinates with smooth base vectors,
               J. Geophys. Res., 115(A8), A08322, :doi:`10.1029/2010JA015326`.

        """

        glat, glon = self.convert(lat, lon, coords, 'geo', height=height,
                                  precision=precision)

        returnvals = self._geo2apexall(glat, glon, height)
        qlat = np.float64(returnvals[0])
        alat = np.float64(returnvals[2])
        f1, f2 = returnvals[4:6]
        d1, d2, d3 = returnvals[7:10]
        e1, e2, e3 = returnvals[11:14]

        # if inputs are not scalar, each vector is an array of arrays,
        # so reshape to a single array
        if f1.dtype == object:
            f1 = np.vstack(f1).T
            f2 = np.vstack(f2).T
            d1 = np.vstack(d1).T
            d2 = np.vstack(d2).T
            d3 = np.vstack(d3).T
            e1 = np.vstack(e1).T
            e2 = np.vstack(e2).T
            e3 = np.vstack(e3).T

        # make sure arrays are 2D
        f1 = f1.reshape((2, f1.size // 2))
        f2 = f2.reshape((2, f2.size // 2))
        d1 = d1.reshape((3, d1.size // 3))
        d2 = d2.reshape((3, d2.size // 3))
        d3 = d3.reshape((3, d3.size // 3))
        e1 = e1.reshape((3, e1.size // 3))
        e2 = e2.reshape((3, e2.size // 3))
        e3 = e3.reshape((3, e3.size // 3))

        # compute f3, g1, g2, g3
        F1 = np.vstack((f1, np.zeros_like(f1[0])))
        F2 = np.vstack((f2, np.zeros_like(f2[0])))
        F = np.cross(F1.T, F2.T).T[-1]
        cosI = helpers.getcosIm(alat)
        k = np.array([0, 0, 1], dtype=np.float64).reshape((3, 1))
        g1 = ((self.RE + np.float64(height))
              / (self.RE + self.refh)) ** (3 / 2) * d1 / F
        g2 = -1.0 / (2.0 * F * np.tan(np.radians(qlat))) * (
            k + ((self.RE + np.float64(height))
                 / (self.RE + self.refh)) * d2 / cosI)
        g3 = k * F
        f3 = np.cross(g1.T, g2.T).T

        if np.any(alat == -9999):
            warnings.warn(''.join(['Base vectors g, d, e, and f3 set to NaN ',
                                   'where apex latitude is undefined (apex ',
                                   'height may be < reference height)']))
            mask = alat == -9999
            f3 = np.where(mask, np.nan, f3)
            g1 = np.where(mask, np.nan, g1)
            g2 = np.where(mask, np.nan, g2)
            g3 = np.where(mask, np.nan, g3)
            d1 = np.where(mask, np.nan, d1)
            d2 = np.where(mask, np.nan, d2)
            d3 = np.where(mask, np.nan, d3)
            e1 = np.where(mask, np.nan, e1)
            e2 = np.where(mask, np.nan, e2)
            e3 = np.where(mask, np.nan, e3)

        return tuple(np.squeeze(x) for x in
                     [f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3])

    def get_apex(self, lat, height=None):
        """ Calculate apex height

        Parameters
        -----------
        lat : (float)
            Latitude in degrees
        height : (float or NoneType)
            Height above the surface of the earth in km or NoneType to use
            reference height (default=None)

        Returns
        ----------
        apex_height : (float)
            Height of the field line apex in km
        """
        lat = helpers.checklat(lat, name='alat')
        if height is None:
            height = self.refh

        cos_lat_squared = np.cos(np.radians(lat)) ** 2
        apex_height = (self.RE + height) / cos_lat_squared - self.RE

        return apex_height

    def set_epoch(self, year):
        """Updates the epoch for all subsequent conversions.

        Parameters
        ----------
        year : float
            Decimal year

        """
        # f2py
        self.year = np.float64(year)
        fa.loadapxsh(self.datafile, self.year)
        igrf_fn = os.path.join(os.path.dirname(__file__), 'igrf13coeffs.txt')
        if not os.path.exists(igrf_fn):
            raise OSError("File {} does not exist".format(igrf_fn))
        fa.cofrm(self.year, igrf_fn)

    def set_refh(self, refh):
        """Updates the apex reference height for all subsequent conversions.

        Parameters
        ----------
        refh : float
            Apex reference height in km

        Notes
        -----
        The reference height is the height to which field lines will be mapped,
        and is only relevant for conversions involving apex (not quasi-dipole).

        """
        self.refh = refh

    def _get_babs_nonvectorized(self, glat, glon, height):
        bnorth, beast, bdown, babs = fa.feldg(1, glat, glon, height)
        # BABS is in guass, so convert to tesla
        return babs / 10000.0

    def get_babs(self, glat, glon, height):
        """Returns the magnitude of the IGRF magnetic field in tesla.

        Parameters
        ----------
        glat : array_like
            Geodetic latitude
        glon : array_like
            Geodetic longitude
        height : array_like
            Altitude in km

        Returns
        -------
        babs : ndarray or float
            Magnitude of the IGRF magnetic field

        """

        babs = self._get_babs(glat, glon, height)

        # if array is returned, the dtype is object, so convert to float
        return np.float64(babs)

    def bvectors_apex(self, lat, lon, height, coords='geo', precision=1e-10):
        """Returns the magnetic field vectors in apex coordinates.

        The apex magnetic field vectors described by Richmond [1995] [4]_ and
        Emmert et al. [2010] [5]_, specfically the Be3 and Bd3 components. The
        vector components are geodetic east, north, and up.

        Parameters
        ----------
        lat : (N,) array_like or float
            Latitude
        lon : (N,) array_like or float
            Longitude
        height : (N,) array_like or float
            Altitude in km
        coords : {'geo', 'apex', 'qd'}, optional
            Input coordinate system
        precision : float, optional
            Precision of output (degrees) when converting to geo. A negative
            value of this argument produces a low-precision calculation of
            geodetic lat/lon based only on their spherical harmonic
            representation.
            A positive value causes the underlying Fortran routine to iterate
            until feeding the output geo lat/lon into geo2qd (APXG2Q) reproduces
            the input QD lat/lon to within the specified precision (all
            coordinates being converted to geo are converted to QD first and
            passed through APXG2Q).

        Returns
        -------
        Be3: (1, N) or (1,) ndarray
        e3 : (3, N) or (3,) ndarray
        Bd3: (1, N) or (1,) ndarray
        d3 : (3, N) or (3,) ndarray

        Notes
        -----
        Be3 is not equivalent to the magnitude of the IGRF magnitude, but is
        instead equal to the IGRF magnitude divided by a scaling factor, D.
        Similarly, Bd3 is the IGRF magnitude multiplied by D.

        See Richmond, A. D. (1995) [4]_ equations 3.13 and 3.14

        References
        ----------
        Richmond, A. D. (1995) [4]_
        Emmert, J. T. et al. (2010) [5]_

        """
        glat, glon = self.convert(lat, lon, coords, 'geo', height=height,
                                  precision=precision)

        babs = self.get_babs(glat, glon, height)

        _, _, _, _, _, _, d1, d2, d3, _, _, e3 = self.basevectors_apex(
            glat, glon, height, coords='geo')
        d1_cross_d2 = np.cross(d1.T, d2.T).T
        D = np.sqrt(np.sum(d1_cross_d2 ** 2, axis=0))

        Be3 = babs / D
        Bd3 = babs * D

        return Be3, e3, Bd3, d3
