# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals
import os
import numpy as np

from . import helpers
import datetime as dt

from .helpers import d2r, r2d

# below try..catch required for autodoc to work on readthedocs
try:
    from . import fortranapex as fa
except:
    print("ERROR: fortranapex module could not be imported. apexpy probably won't work")


class ApexHeightError(ValueError):
    pass


class Apex(object):
    '''Performs coordinate conversions

    Parameters
    ==========
    date : float (decimal year) or instance of :class:`~datetime.date` or :class:`~datetime.datetime`
        IGRF coefficients are used in conversions. Uses current date as default.
    refh : float
        Reference height in km for apex coordinates (the field lines are mapped to this height)
    datafile : str
        Path to custom coefficient file

    Methods
    =======
    :meth:`convert`
        High-level, general-purpose conversion between geodetic, modified apex, quasi-dipole and MLT
    :meth:`geo2apex`, :meth:`apex2geo`, :meth:`geo2qd`, :meth:`qd2geo`, \
    :meth:`apex2qd`, :meth:`qd2apex`, :meth:`mlon2mlt`, :meth:`mlt2mlon`
        conversion functions for specific coordinate systems (called by :meth:`convert`)
    :meth:`map_to_height`
        Maps geodetic coordinates along the magnetic field to a new height in the closest or conjugate hemisphere
        (for finding footprints, conjugate points, etc.)
    :meth:`basevectors_qd`, :meth:`basevectors_apex`
        Calculate base vectors
    :meth:`get_apex`
        Compute field line apex from apex latitude
    :meth:`set_epoch`, :meth:`set_refh`
        Change epoch and reference height for subsequent conversions

    Attributes
    ==========
    year : float
        Decimal year used for the IGRF model
    refh : float
        Reference height in km for apex coordinates
    datafile : str
        Path to coefficient file

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
        '''Converts latitude and longitude/MLT between two coordinate systems

        `lat`, `lon`, `height` must be broadcastable to the same shape.

        Parameters
        ==========
        lat : array_like
            Latitude to convert
        lon : array_like
            Longitude to convert (or MLT if converting from MLT)
        source : {'geo', 'apex', 'qd', 'mlt'}
            Input coordinate system
        dest : {'geo', 'apex', 'qd', 'mlt'}
            Output coordinate system
        height : array_like
            Altitude in km
        datetime : :class:`~datetime.datetime`
            Date and time for MLT conversions (required for MLT conversions)
        precision : float
            Precision of output (degrees) when converting to geo. A negative
            value of this argument produces a low-precision calculation of
            geodetic lat/lon based only on their spherical harmonic representation.
            A positive value causes the underlying Fortran routine to iterate
            until feeding the output geo lat/lon into geo2qd (APXG2Q) reproduces
            the input QD lat/lon to within the specified precision (all
            coordinates being converted to geo are converted to QD first and
            passed through APXG2Q).
        ssheight : float
            Altitude in km to use for converting the subsolar point from
            geographic to magnetic coordinates. A high altitude is used
            to ensure the subsolar point is mapped to high latitudes, which
            prevents the South-Atlantic Anomaly (SAA) from influencing the MLT.

        Returns
        =======
        If input `lat`, `lon` or `height` have dimension > 0, outputs are arrays.

        lat : float or array
            Converted latitude
        lon : float or array
            Converted longitude (or MLT if converting to MLT)

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
        '''Converts geodetic to modified apex coordinates

        Parameters
        ==========
        glat : array_like
            Geodetic latitude
        glon : array_like
            Geodetic longitude
        height : array_like
            Altitude in km

        Returns
        =======
        If any inputs have dimension > 0, outputs are arrays.

        alat : float or array
            Modified apex latitude
        alon : float or array
            Modified apex longitude

        '''

        glat = helpers.checklat(glat, name='glat')

        alat, alon = self._geo2apex(glat, glon, height)

        # if array is returned, dtype is object, so convert to float
        return np.float64(alat), np.float64(alon)

    def apex2geo(self, alat, alon, height, precision=1e-10):
        '''Converts modified apex to geodetic coordinates

        Parameters
        ==========
        alat : array_like
            Modified apex latitude
        alon : array_like
            Modified apex longitude
        height : array_like
            Altitude in km
        precision : float
            Precision of output (degrees). A negative value of this argument
            produces a low-precision calculation of geodetic lat/lon based only
            on their spherical harmonic representation. A positive value causes
            the underlying Fortran routine to iterate until feeding the output
            geo lat/lon into geo2qd (APXG2Q) reproduces the input QD lat/lon to
            within the specified precision.

        Returns
        =======
        If `alat`, `alon` or `height` have dimension > 0, outputs are arrays.

        glat : float or array
            Geodetic latitude
        glon : float or array
            Geodetic longitude
        error : float or array
            The angular difference (degrees) between the input QD coordinates
            and the qlat/qlon produced by feeding the output glat and glon
            into geo2qd (APXG2Q)

        '''

        alat = helpers.checklat(alat, name='alat')

        qlat, qlon = self.apex2qd(alat, alon, height=height)
        glat, glon, error = self.qd2geo(qlat, qlon, height, precision=precision)

        return glat, glon, error

    def geo2qd(self, glat, glon, height):
        '''Converts geodetic to quasi-dipole coordinates

        Parameters
        ==========
        glat : array_like
            Geodetic latitude
        glon : array_like
            Geodetic longitude
        height : array_like
            Altitude in km

        Returns
        =======
        If any inputs have dimension > 0, outputs are arrays.

        qlat : float or array
            Quasi-dipole latitude
        qlon : float or array
            Quasi-dipole longitude

        '''

        glat = helpers.checklat(glat, name='glat')

        qlat, qlon = self._geo2qd(glat, glon, height)

        # if array is returned, dtype is object, so convert to float
        return np.float64(qlat), np.float64(qlon)

    def qd2geo(self, qlat, qlon, height, precision=1e-10):
        '''Converts quasi-dipole to geodetic coordinates

        Parameters
        ==========
        qlat : array_like
            Quasi-dipole latitude
        qlon : array_like
            Quasi-dipole longitude
        height : array_like
            Altitude in km
        precision : float
            Precision of output (degrees). A negative value of this argument
            produces a low-precision calculation of geodetic lat/lon based only
            on their spherical harmonic representation. A positive value causes
            the underlying Fortran routine to iterate until feeding the output
            geo lat/lon into geo2qd (APXG2Q) reproduces the input QD lat/lon to
            within the specified precision.

        Returns
        =======
        If `qlat`, `qlon` or `height` have dimension > 0, outputs are arrays.

        glat : float or array
            Geodetic latitude
        glon : float or array
            Geodetic longitude
        error : float or array
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
        '''Converts modified apex to quasi-dipole coordinates

        Parameters
        ==========
        alat : array_like
            Modified apex latitude
        alon : array_like
            Modified apex longitude
        height : array_like
            Altitude in km

        Returns
        =======
        If any inputs have dimension > 0, outputs are arrays.

        qlat : float or array
            Quasi-dipole latitude
        qlon : float or array
            Quasi-dipole longitude

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
        '''Converts quasi-dipole to modified apex coordinates

        Parameters
        ==========
        qlat : array_like
            Quasi-dipole latitude
        qlon : array_like
            Quasi-dipole longitude
        height : array_like
            Altitude in km

        Returns
        =======
        If any inputs have dimension > 0, outputs are arrays.

        alat : float or array
            Modified apex latitude
        alon : float or array
            Modified apex longitude

        Raises
        ======
        ApexHeightError
            if apex height is < reference height

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
            magnetic longitude (apex and quasi-dipole longitude are always equal)
        datetime : :class:`~datetime.datetime`
            Date and time
        ssheight : float
            Altitude in km to use for converting the subsolar point from
            geographic to magnetic coordinates. A high altitude is used
            to ensure the subsolar point is mapped to high latitudes, which
            prevents the South-Atlantic Anomaly (SAA) from influencing the MLT.

        Returns
        =======
        If `mlon` has dimension > 0, output is array.

        mlt : float or array
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
            magnetic local time
        datetime : :class:`~datetime.datetime`
            Date and time
        ssheight : float
            Altitude in km to use for converting the subsolar point from
            geographic to magnetic coordinates. A high altitude is used
            to ensure the subsolar point is mapped to high latitudes, which
            prevents the South-Atlantic Anomaly (SAA) from influencing the MLT.

        Returns
        =======
        If `mlt` has dimension > 0, output is array.

        mlon : float or array
            Magnetic longitude [0, 360) (apex and quasi-dipole longitude are always equal)
        '''

        ssglat, ssglon = helpers.subsol(datetime)
        ssglat = helpers.gc2gdlat(ssglat)
        ssalat, ssalon = self.geo2apex(ssglat, ssglon, ssheight)

        # np.float64 will ensure lists are converted to arrays
        return (15*np.float64(mlt) - 180 + ssalon + 360) % 360

    def map_to_height(self, glat, glon, height, newheight, conjugate=False, precision=1e-10):
        '''Maps geodetic coordinates along the magnetic field to `newheight` in the closest or conjugate hemisphere

        This is done by converting glat/glon/height to modified apex lat/lon, and converting back
        to geographic using newheight (if conjugate, use negative apex latitude when converting back)

        Parameters
        ==========
        glat : array_like
            Geodetic latitude
        glon : array_like
            Geodetic longitude
        height : array_like
            Source altitude in km
        newheight : array_like
            Destination altitude in km
        conjugate : bool
            Map to `newheight` in the conjugate hemisphere instead of the closest hemisphere
        precision : float
            Precision of output (degrees). A negative value of this argument
            produces a low-precision calculation of geodetic lat/lon based only
            on their spherical harmonic representation. A positive value causes
            the underlying Fortran routine to iterate until feeding the output
            geo lat/lon into geo2qd (APXG2Q) reproduces the input QD lat/lon to
            within the specified precision.

        Returns
        =======
        If `glat`, `glon` or `height` have dimension > 0, outputs are arrays.

        newglat : float or array
            Geodetic latitude of mapped point
        newglon : float or array
            Geodetic longitude of mapped point
        error : float or array
            The angular difference (degrees) between the input QD coordinates
            and the qlat/qlon produced by feeding the output glat and glon
            into geo2qd (APXG2Q)

        '''

        alat, alon = self.geo2apex(glat, glon, height)
        if conjugate:
            alat = -alat
        try:
            newglat, newglon, error = self.apex2geo(alat, alon, newheight, precision=precision)
        except ApexHeightError:
            raise ApexHeightError("newheight is > apex height")

        return newglat, newglon, error

    def basevectors_qd(self, lat, lon, height, coords='geo', precision=1e-10):
        '''Returns quasi-dipole base vectors f1 and f2 at the specified coordinates

        The vectors are described by Richmond [2005] [1]_.

        Parameters
        ==========
        lat : array_like
            Latitude
        lon : array_like
            Longitude
        height : array_like
            Altitude in km
        coords : {'geo', 'apex', 'qd'}
            Input coordinate system
        precision : float
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

        f1 : array
            ..
        f2 : array
            ..

        Output shapes:

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
        '''Returns base vectors in quasi-dipole (f) and apex (d and e) coordinates at the specified coordinates

        The vectors are described by Richmond [2005] [1]_.

        Parameters
        ==========
        lat : array_like
            Latitude
        lon : array_like
            Longitude
        height : array_like
            Altitude in km
        coords : {'geo', 'apex', 'qd'}
            Input coordinate system
        return_all : bool
            Will also return f3, g1, g2, and g3, and f1 and f2 have 3 components
            (the last one is zero). Requires `lat`, `lon`, `height` to be
            broadcast to 1D (at least one of the parameters must be 1D and the
            other two parameters must be 1D or 0D).
        precision : float
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

        Brackets indicate calling with ``return_all=True``.

        f1 : array
            ..
        f2 : array
            ..
        [f3] : array
             ..
        [g1] : array
            ..
        [g2] : array
            ..
        [g3] : array
            ..
        d1 : array
            ..
        d2 : array
            ..
        d3 : array
            ..
        e1 : array
            ..
        e2 : array
            ..
        e3 : array
            ..

        Output shapes:

        * If the inputs are scalar (only if ``return_all=False``), the outputs
          are vectors with 2 (f) and 3 (d and e) components.

        * If the inputs broadcast to 1D with length N, the outputs are 2xN (f) and
          3xN (d and e) arrays where the columns are the vectors, i.e. ``f1[:, 0]`` is the
          f1 vector corresponding to the first index in the broadcasted input. If
          ``return_all=True``, f is 3xN instead of 2xN.

        * If the inputs broadcast to 2D with shape NxM (only if ``return_all=False``),
          the outputs are 2xNxM and 3xNxM arrays where ``f1[:, 0, 0]`` is the f1
          vector corresponding to the index [0, 0] in the broadcasted input.

        * Higher dimensions are untested.

        References
        ==========

        .. [1] Richmond, A. D. (1995), Ionospheric Electrodynamics Using
               Magnetic Apex Coordinates, Journal of geomagnetism and
               geoelectricity, 47(2), 191–212, doi:10.5636/jgg.47.191

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
        '''Computes the field line apex for given modified apex latitude

        Parameters
        ==========
        alat : array_like
            Modified apex latitude

        Returns
        =======
        If `alat` has dimension > 0, output is array.

        apex : float or array
            Field line apex
        '''

        alat = helpers.checklat(alat, name='alat')

        apex = (self.RE + self.refh)/np.cos(d2r*alat)**2 - self.RE

        return apex

    def set_epoch(self, year):
        '''Updates the epoch for all subsequent conversions

        Parameters
        ==========
        year : float
            Decimal year

        '''

        fa.loadapxsh(self.datafile, np.float(year))
        self.year = year

    def set_refh(self, refh):
        '''Updates the apex reference height for all subsequent conversions

        The reference height is the height to which field lines will be mapped,
        and is only relevant for conversions involving apex (not quasi-dipole).

        Parameters
        ==========
        refh : float
            Apex reference height in km

        '''

        self.refh = refh
