# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import

import os
import datetime as dt

import numpy as np

from . import helpers
from .helpers import d2r, r2d

# below try..catch required for autodoc to work on readthedocs
try:
    from . import fortranapex as fa
except:
    print("ERROR: fortranapex module could not be imported. apexpy probably won't work")


class ApexHeightError(ValueError):
    pass


class Apex(object):
    '''Performs coordinate conversions, field-line mapping and base vector calculations.

    Parameters
    ==========
    date : float (decimal year) or instance of :class:`datetime.date` or :class:`datetime.datetime`, optional
        Determines which IGRF coefficients are used in conversions. Uses current date as default.
    refh : float, optional
        Reference height in km for apex coordinates (the field lines are mapped to this height)
    datafile : str, optional
        Path to custom coefficient file

    Attributes
    ==========
    year : float
        Decimal year used for the IGRF model
    refh : float
        Reference height in km for apex coordinates
    datafile : str
        Path to coefficient file

    Notes
    =====
    The calculations use IGRF-12 with coefficients from 1900 to 2020 [1]_.

    References
    ==========

    .. [1] Thébault, E. et al. (2015), International Geomagnetic Reference
           Field: the 12th generation, Earth, Planets and Space, 67(1), 79,
           doi:10.1186/s40623-015-0228-9.

    '''

    def __init__(self, date=None, refh=0, datafile=None):

        if datafile is None:
            datafile = os.path.join(os.path.dirname(__file__), 'apexsh.dat')

        self.RE = 6371.009  # mean Earth radius
        self.set_refh(refh)  # reference height

        if date is None:
            self.year = helpers.toYearFraction(dt.datetime.now())
        else:
            try:
                # convert date/datetime object to decimal year
                self.year = helpers.toYearFraction(date)
            except:
                # failed so date is probably int/float, use directly
                self.year = date

        if not os.path.isfile(datafile):
            raise IOError('Datafile does not exist: {}'.format(datafile))

        self.datafile = datafile
        self.set_epoch(self.year)

        # vectorize fortran functions
        self._geo2qd = np.frompyfunc(
            lambda glat, glon, height: fa.apxg2q(glat, (glon + 180) % 360 - 180, height, 0)[:2], 3, 2)
        self._geo2apex = np.frompyfunc(
            lambda glat, glon, height: fa.apxg2all(glat, (glon + 180) % 360 - 180, height, self.refh, 0)[2:4], 3, 2)
        self._geo2apexall = np.frompyfunc(
            lambda glat, glon, height: fa.apxg2all(glat, (glon + 180) % 360 - 180, height, self.refh, 1), 3, 14)
        self._qd2geo = np.frompyfunc(
            lambda qlat, qlon, height, precision: fa.apxq2g(qlat, (qlon + 180) % 360 - 180, height, precision), 4, 3)
        self._basevec = np.frompyfunc(
            lambda glat, glon, height: fa.apxg2q(glat, (glon + 180) % 360 - 180, height, 1)[2:4], 3, 2)

        # vectorize other nonvectorized functions
        self._apex2qd = np.frompyfunc(self._apex2qd_nonvectorized, 3, 2)
        self._qd2apex = np.frompyfunc(self._qd2apex_nonvectorized, 3, 2)

    def convert(self, lat, lon, source, dest, height=0, datetime=None, precision=1e-10, ssheight=50*6371):
        '''Converts between geodetic, modified apex, quasi-dipole and MLT.

        Parameters
        ==========
        lat, lon : array_like
            Latitude and longitude/MLT
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
            geodetic lat/lon based only on their spherical harmonic representation.
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
        =======
        lat, lon : ndarray or float
            Converted latitude and longitude/MLT

        .. note::
            If converting to MLT, output latitude is apex.

        '''

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
            raise NotImplementedError('Unknown coordinate transformation: {} -> {}'.format(source, dest))

        return lat, lon

    def geo2apex(self, glat, glon, height):
        '''Converts geodetic to modified apex coordinates.

        Parameters
        ==========
        glat, glon : array_like
            Geodetic latitude and longitude
        height : array_like
            Altitude in km

        Returns
        =======
        alat, alon : ndarray or float
            Modified apex latitude and longitude

        '''

        glat = helpers.checklat(glat, name='glat')

        alat, alon = self._geo2apex(glat, glon, height)

        # if array is returned, dtype is object, so convert to float
        return np.float64(alat), np.float64(alon)

    def apex2geo(self, alat, alon, height, precision=1e-10):
        '''Converts modified apex to geodetic coordinates.

        Parameters
        ==========
        alat, alon : array_like
            Modified apex latitude and longitude
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
        =======
        glat, glon : ndarray or float
            Geodetic latitude and longitude
        error : ndarray or float
            The angular difference (degrees) between the input QD coordinates
            and the qlat/qlon produced by feeding the output glat and glon
            into geo2qd (APXG2Q)

        '''

        alat = helpers.checklat(alat, name='alat')

        qlat, qlon = self.apex2qd(alat, alon, height=height)
        glat, glon, error = self.qd2geo(qlat, qlon, height, precision=precision)

        return glat, glon, error

    def geo2qd(self, glat, glon, height):
        '''Converts geodetic to quasi-dipole coordinates.

        Parameters
        ==========
        glat, glon : array_like
            Geodetic latitude and longitude
        height : array_like
            Altitude in km

        Returns
        =======
        qlat, qlon : ndarray or float
            Quasi-dipole latitude and longitude

        '''

        glat = helpers.checklat(glat, name='glat')

        qlat, qlon = self._geo2qd(glat, glon, height)

        # if array is returned, dtype is object, so convert to float
        return np.float64(qlat), np.float64(qlon)

    def qd2geo(self, qlat, qlon, height, precision=1e-10):
        '''Converts quasi-dipole to geodetic coordinates.

        Parameters
        ==========
        qlat, qlon : array_like
            Quasi-dipole latitude and longitude
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
        =======
        glat, glon : ndarray or float
            Geodetic latitude and longitude
        error : ndarray or float
            The angular difference (degrees) between the input QD coordinates
            and the qlat/qlon produced by feeding the output glat and glon
            into geo2qd (APXG2Q)

        '''

        qlat = helpers.checklat(qlat, name='qlat')

        glat, glon, error = self._qd2geo(qlat, qlon, height, precision)

        # if array is returned, dtype is object, so convert to float
        return np.float64(glat), np.float64(glon), np.float64(error)

    def _apex2qd_nonvectorized(self, alat, alon, height):

        alat = helpers.checklat(alat, name='alat')

        # convert modified apex to quasi-dipole:
        qlon = alon
        hA = (self.RE + self.refh)/(np.cos(alat*d2r)**2) - self.RE  # apex height

        if hA < height:
            if np.isclose(hA, height, rtol=0, atol=1e-5):  # allow for values that are close
                hA = height
            else:
                raise ApexHeightError(
                    'height {:.3g} is > apex height {:.3g} for alat {:.3g}'.format(np.max(height), hA, alat))

        qlat = np.sign(alat) * np.arccos(np.sqrt((self.RE + height)/(self.RE + hA)))*r2d

        return qlat, qlon

    def apex2qd(self, alat, alon, height):
        '''Converts modified apex to quasi-dipole coordinates.

        Parameters
        ==========
        alat, alon : array_like
            Modified apex latitude and longitude
        height : array_like
            Altitude in km

        Returns
        =======
        qlat, qlon : ndarray or float
            Quasi-dipole latitude and longitude

        Raises
        ======
        ApexHeightError
            if `height` > apex height

        '''

        qlat, qlon = self._apex2qd(alat, alon, height)

        # if array is returned, the dtype is object, so convert to float
        return np.float64(qlat), np.float64(qlon)

    def _qd2apex_nonvectorized(self, qlat, qlon, height):

        qlat = helpers.checklat(qlat, name='qlat')

        alon = qlon
        hA = (self.RE + height)/(np.cos(qlat*d2r)**2) - self.RE  # apex height

        if hA < self.refh:
            if np.isclose(hA, self.refh, rtol=0, atol=1e-5):  # allow for values that are close
                hA = self.refh
            else:
                raise ApexHeightError(
                    'apex height ({:.3g}) is < reference height ({:.3g}) for qlat {:.3g}'.format(hA, self.refh, qlat))

        alat = np.sign(qlat) * np.arccos(np.sqrt((self.RE + self.refh)/(self.RE + hA)))*r2d

        return alat, alon

    def qd2apex(self, qlat, qlon, height):
        '''Converts quasi-dipole to modified apex coordinates.

        Parameters
        ==========
        qlat, qlon : array_like
            Quasi-dipole latitude and longitude
        height : array_like
            Altitude in km

        Returns
        =======
        alat, alon : ndarray or float
            Modified apex latitude and longitude

        Raises
        ======
        ApexHeightError
            if apex height < reference height

        '''

        alat, alon = self._qd2apex(qlat, qlon, height)

        # if array is returned, the dtype is object, so convert to float
        return np.float64(alat), np.float64(alon)

    def mlon2mlt(self, mlon, datetime, ssheight=50*6371):
        '''Computes the magnetic local time at the specified magnetic longitude and UT.

        To compute the MLT, we find the apex longitude of the subsolar point at
        the given time. Then the MLT of the given point will be computed from
        the separation in magnetic longitude from this point (1 hour = 15 degrees).

        Parameters
        ==========
        mlon : array_like
            Magnetic longitude (apex and quasi-dipole longitude are always equal)
        datetime : :class:`datetime.datetime`
            Date and time
        ssheight : float, optional
            Altitude in km to use for converting the subsolar point from
            geographic to magnetic coordinates. A high altitude is used
            to ensure the subsolar point is mapped to high latitudes, which
            prevents the South-Atlantic Anomaly (SAA) from influencing the MLT.

        Returns
        =======
        mlt : ndarray or float
            Magnetic local time [0, 24)

        '''
        ssglat, ssglon = helpers.subsol(datetime)
        ssglat = helpers.gc2gdlat(ssglat)
        ssalat, ssalon = self.geo2apex(ssglat, ssglon, ssheight)

        # np.float64 will ensure lists are converted to arrays
        return (180 + np.float64(mlon) - ssalon)/15 % 24

    def mlt2mlon(self, mlt, datetime, ssheight=50*6371):
        '''Computes the magnetic longitude at the specified magnetic local time and UT.

        To compute the magnetic longitude, we find the apex longitude of the subsolar point at
        the given time. Then the magnetic longitude of the given point will be computed from
        the separation in magnetic local time from this point (1 hour = 15 degrees).

        Parameters
        ==========
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
        =======
        mlon : ndarray or float
            Magnetic longitude [0, 360) (apex and quasi-dipole longitude are always equal)
        '''

        ssglat, ssglon = helpers.subsol(datetime)
        ssglat = helpers.gc2gdlat(ssglat)
        ssalat, ssalon = self.geo2apex(ssglat, ssglon, ssheight)

        # np.float64 will ensure lists are converted to arrays
        return (15*np.float64(mlt) - 180 + ssalon + 360) % 360

    def map_to_height(self, glat, glon, height, newheight, conjugate=False, precision=1e-10):
        '''Performs mapping of points along the magnetic field to the closest or conjugate hemisphere.

        Parameters
        ==========
        glat, glon : array_like
            Geodetic latitude and longitude
        height : array_like
            Source altitude in km
        newheight : array_like
            Destination altitude in km
        conjugate : bool, optional
            Map to `newheight` in the conjugate hemisphere instead of the closest hemisphere
        precision : float, optional
            Precision of output (degrees). A negative value of this argument
            produces a low-precision calculation of geodetic lat/lon based only
            on their spherical harmonic representation. A positive value causes
            the underlying Fortran routine to iterate until feeding the output
            geo lat/lon into geo2qd (APXG2Q) reproduces the input QD lat/lon to
            within the specified precision.

        Returns
        =======
        newglat : ndarray or float
            Geodetic latitude of mapped point
        newglon : ndarray or float
            Geodetic longitude of mapped point
        error : ndarray or float
            The angular difference (degrees) between the input QD coordinates
            and the qlat/qlon produced by feeding the output glat and glon
            into geo2qd (APXG2Q)

        Notes
        =====
        The mapping is done by converting glat/glon/height to modified apex lat/lon, and converting back
        to geographic using newheight (if conjugate, use negative apex latitude when converting back)

        '''

        alat, alon = self.geo2apex(glat, glon, height)
        if conjugate:
            alat = -alat
        try:
            newglat, newglon, error = self.apex2geo(alat, alon, newheight, precision=precision)
        except ApexHeightError:
            raise ApexHeightError("newheight is > apex height")

        return newglat, newglon, error

    def _map_EV_to_height(self, alat, alon, height, newheight, X, EV):

        # make sure X is array of correct shape
        if not (np.ndim(X) == 1 and np.size(X) == 3) and not (np.ndim(X) == 2 and np.shape(X)[0] == 3):
            # raise ValueError because if passing e.g. a (6,) ndarray the reshape below will work
            # even though the input is invalid
            raise ValueError(EV + ' must be (3, N) or (3,) ndarray')
        X = np.reshape(X, (3, np.size(X)/3))

        _, _, d1, d2, _, e1, e2, _ = self.basevectors_apex(alat, alon, height, coords='apex')

        if EV == 'E':
            v1 = e1
            v2 = e2
        else:
            v1 = d1
            v2 = d2

        # make sure v1 and v2 have shape (3, N)
        v1 = np.reshape(v1, (3, v1.size/3))
        v2 = np.reshape(v2, (3, v2.size/3))

        X1 = np.sum(X*v1, axis=0)  # E dot e1 or V dot d1
        X2 = np.sum(X*v2, axis=0)  # E dot e2 or V dot d2

        _, _, d1, d2, _, e1, e2, _ = self.basevectors_apex(alat, alon, newheight, coords='apex')

        if EV == 'E':
            v1 = d1
            v2 = d2
        else:
            v1 = e1
            v2 = e2

        # make sure v1 and v2 have shape (3, N)
        v1 = np.reshape(v1, (3, v1.size/3))
        v2 = np.reshape(v2, (3, v2.size/3))

        X_mapped = X1[np.newaxis, :]*v1 + X2[np.newaxis, :]*v2

        return np.squeeze(X_mapped)

    def map_E_to_height(self, alat, alon, height, newheight, E):
        '''Performs mapping of electric field along the magnetic field.

        It is assumed that the electric field is perpendicular to B

        Parameters
        ==========
        alat, alon : (N,) array_like or float
            Apex latitude and longitude
        height : (N,) array_like or float
            Source altitude in km
        newheight : (N,) array_like or float
            Destination altitude in km
        E : (3,) or (3, N) array_like
            Electric field (at `alat`, `alon`, `height`) in geodetic east, north, and up components

        Returns
        =======
        E : (3, N) or (3,) ndarray
            The electric field at `newheight` (geodetic east, north, and up components)

        '''

        return self._map_EV_to_height(alat, alon, height, newheight, E, 'E')

    def map_V_to_height(self, alat, alon, height, newheight, V):
        '''Performs mapping of electric drivt velocity along the magnetic field.

        It is assumed that the electric field is perpendicular to B

        Parameters
        ==========
        alat, alon : (N,) array_like or float
            Apex latitude and longitude
        height : (N,) array_like or float
            Source altitude in km
        newheight : (N,) array_like or float
            Destination altitude in km
        V : (3,) or (3, N) array_like
            Electric drift velocity (at `alat`, `alon`, `height`) in geodetic east, north, and up components

        Returns
        =======
        V : (3, N) or (3,) ndarray
            The electric drivt velocity at `newheight` (geodetic east, north, and up components)

        '''

        return self._map_EV_to_height(alat, alon, height, newheight, V, 'V')

    def basevectors_qd(self, lat, lon, height, coords='geo', precision=1e-10):
        '''Returns quasi-dipole base vectors f1 and f2 at the specified coordinates.

        The vectors are described by Richmond [2005] [1]_ and Emmert et al. [2010] [2]_.

        Parameters
        ==========
        lat, lon : array_like
            Latitude and longitude
        height : array_like
            Altitude in km
        coords : {'geo', 'apex', 'qd'}, optional
            Input coordinate system
        precision : float, optional
            Precision of output (degrees) when converting to geo. A negative
            value of this argument produces a low-precision calculation of
            geodetic lat/lon based only on their spherical harmonic representation.
            A positive value causes the underlying Fortran routine to iterate
            until feeding the output geo lat/lon into geo2qd (APXG2Q) reproduces
            the input QD lat/lon to within the specified precision (all
            coordinates being converted to geo are converted to QD first and
            passed through APXG2Q).

        Returns
        =======

        f1, f2 : ndarray (see note below on output shapes)

        .. note::

            * If the inputs are scalar, the outputs are vectors with 2 components.

            * If the inputs broadcast to 1D with length N, the outputs are 2xN arrays
              where the columns are the vectors, i.e. ``f1[:, 0]`` is the f1 vector
              corresponding to the first index in the broadcasted input.

            * If the inputs broadcast to 2D with shape NxM, the outputs are 2xNxM arrays
              where ``f1[:, 0, 0]`` is the f1 vector corresponding to the index [0, 0] in
              the broadcasted input.

            * Higher dimensions are untested.

        References
        ==========

        .. [1] Richmond, A. D. (1995), Ionospheric Electrodynamics Using
               Magnetic Apex Coordinates, Journal of geomagnetism and
               geoelectricity, 47(2), 191–212, doi:10.5636/jgg.47.191.

        .. [2] Emmert, J. T., A. D. Richmond, and D. P. Drob (2010),
               A computationally compact representation of Magnetic-Apex
               and Quasi-Dipole coordinates with smooth base vectors,
               J. Geophys. Res., 115(A8), A08322, doi:10.1029/2010JA015326.

        '''

        glat, glon = self.convert(lat, lon, coords, 'geo', height=height, precision=precision)

        f1, f2 = self._basevec(glat, glon, height)

        # if inputs are not scalar, each vector is an array of arrays,
        # so reshape to a single array
        if f1.dtype == object:
            f1 = np.dstack(f1.flat).reshape([2] + list(f1.shape))
            f2 = np.dstack(f2.flat).reshape([2] + list(f2.shape))

        return f1, f2

    def basevectors_apex(self, lat, lon, height, coords='geo', return_all=False, precision=1e-10):
        '''Returns base vectors in quasi-dipole and apex coordinates.

        The vectors are described by Richmond [2005] [1]_ and Emmert et al. [2010] [2]_.

        Parameters
        ==========
        lat, lon : array_like
            Latitude and longitude
        height : array_like
            Altitude in km
        coords : {'geo', 'apex', 'qd'}, optional
            Input coordinate system
        return_all : bool, optional
            Will also return f3, g1, g2, and g3, and f1 and f2 have 3 components
            (the last component is zero). Requires `lat`, `lon`, and `height` to be
            broadcast to 1D (at least one of the parameters must be 1D and the
            other two parameters must be 1D or 0D).
        precision : float, optional
            Precision of output (degrees) when converting to geo. A negative
            value of this argument produces a low-precision calculation of
            geodetic lat/lon based only on their spherical harmonic representation.
            A positive value causes the underlying Fortran routine to iterate
            until feeding the output geo lat/lon into geo2qd (APXG2Q) reproduces
            the input QD lat/lon to within the specified precision (all
            coordinates being converted to geo are converted to QD first and
            passed through APXG2Q).

        Returns
        =======

        f1, f2, d1, d2, d3, e1, e2, e3 : ndarray
            if `return_all` is False
        f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 : ndarray
            if `return_all` is True

        .. note::

            * If the inputs are scalar (only if `return_all` is False), the outputs
              are vectors with 2 (f) and 3 (d and e) components.

            * If the inputs broadcast to 1D with length N, the outputs are 2xN (f) and
              3xN (d and e) arrays where the columns are the vectors, i.e. ``f1[:, 0]`` is the
              f1 vector corresponding to the first index in the broadcasted input. If
              if `return_all` is True, f is 3xN instead of 2xN.

            * If the inputs broadcast to 2D with shape NxM (only if `return_all` is False),
              the outputs are 2xNxM and 3xNxM arrays where ``f1[:, 0, 0]`` is the f1
              vector corresponding to the index [0, 0] in the broadcasted input.

            * Higher dimensions are untested.

        References
        ==========

        .. [1] Richmond, A. D. (1995), Ionospheric Electrodynamics Using
               Magnetic Apex Coordinates, Journal of geomagnetism and
               geoelectricity, 47(2), 191–212, doi:10.5636/jgg.47.191.

        .. [2] Emmert, J. T., A. D. Richmond, and D. P. Drob (2010),
               A computationally compact representation of Magnetic-Apex
               and Quasi-Dipole coordinates with smooth base vectors,
               J. Geophys. Res., 115(A8), A08322, doi:10.1029/2010JA015326.

        '''

        glat, glon = self.convert(lat, lon, coords, 'geo', height=height, precision=precision)

        returnvals = self._geo2apexall(glat, glon, height)
        f1, f2 = returnvals[4:6]
        d1, d2, d3 = returnvals[7:10]
        e1, e2, e3 = returnvals[11:14]

        # if inputs are not scalar, each vector is an array of arrays,
        # so reshape to a single array
        if f1.dtype == object:
            f1 = np.dstack(f1.flat).reshape([2] + list(f1.shape))
            f2 = np.dstack(f2.flat).reshape([2] + list(f2.shape))
            d1 = np.dstack(d1.flat).reshape([3] + list(d1.shape))
            d2 = np.dstack(d2.flat).reshape([3] + list(d2.shape))
            d3 = np.dstack(d3.flat).reshape([3] + list(d3.shape))
            e1 = np.dstack(e1.flat).reshape([3] + list(e1.shape))
            e2 = np.dstack(e2.flat).reshape([3] + list(e2.shape))
            e3 = np.dstack(e3.flat).reshape([3] + list(e3.shape))

        if not return_all:
            return f1, f2, d1, d2, d3, e1, e2, e3
        else:

            # check that inputs are of the correct shape
            if np.broadcast(lat, lon, height).nd != 1:
                raise ValueError(
                    ('return_all=True requires lat, lon, height to be broadcast '
                     'to 1D (at least one parameter must be 1D and the other '
                     'parameters must be 1D or 0D'))

            # compute f3, g1, g2, g3
            qdlat, qdlon = self.geo2qd(lat, lon, height)
            alat, alon = self.geo2apex(lat, lon, height)
            F1 = np.vstack((f1, np.zeros_like(f1[0])))
            F2 = np.vstack((f2, np.zeros_like(f2[0])))
            F = np.cross(F1.T, F2.T).T[-1]
            cosI = helpers.getcosIm(alat)
            k = np.array([0, 0, 1], dtype=np.float64).reshape((3, 1))
            g1 = ((self.RE + height)/(self.RE + self.refh))**(3/2)*d1/F
            g2 = -1/(2*F*np.tan(qdlat*d2r))*(k + ((self.RE + height)/(self.RE + self.refh))*d2/cosI)
            g3 = k*F
            f3 = np.cross(g1.T, g2.T).T

            return F1, F2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3

    def get_apex(self, alat):
        '''Computes the field line apex for a given modified apex latitude.

        Parameters
        ==========
        alat : array_like
            Modified apex latitude

        Returns
        =======
        apex : ndarray or float
            Field line apex in km
        '''

        alat = helpers.checklat(alat, name='alat')

        apex = (self.RE + self.refh)/np.cos(d2r*alat)**2 - self.RE

        return apex

    def set_epoch(self, year):
        '''Updates the epoch for all subsequent conversions.

        Parameters
        ==========
        year : float
            Decimal year

        '''

        fa.loadapxsh(self.datafile, np.float(year))
        self.year = year

    def set_refh(self, refh):
        '''Updates the apex reference height for all subsequent conversions.

        Parameters
        ==========
        refh : float
            Apex reference height in km

        Notes
        =====
        The reference height is the height to which field lines will be mapped,
        and is only relevant for conversions involving apex (not quasi-dipole).

        '''

        self.refh = refh
