# -*- coding: utf-8 -*-
"""Classes that make up the core of apexpy."""

import datetime as dt
import numpy as np
import os
import warnings

from apexpy import helpers

# Below try..catch required for autodoc to work on readthedocs
try:
    from apexpy import fortranapex as fa
except ImportError as ierr:
    warnings.warn("".join(["fortranapex module could not be imported, so ",
                           "apexpy probably won't work.  Make sure you have ",
                           "a gfortran compiler. Failed with error: ",
                           str(ierr)]))

# Make sure invalid warnings are always shown
warnings.filterwarnings('always', message='.*set to NaN where*',
                        module='apexpy.apex')


class ApexHeightError(ValueError):
    """Specialized ValueError definition, to be used when apex height is wrong.
    """
    pass


class Apex(object):
    """Calculates coordinate conversions, field-line mapping, and base vectors.

    Parameters
    ----------
    date : NoneType, float, :class:`dt.date`, or :class:`dt.datetime`, optional
        Determines which IGRF coefficients are used in conversions. Uses
        current date as default.  If float, use decimal year.  If None, uses
        current UTC. (default=None)
    refh : float, optional
        Reference height in km for apex coordinates, the field lines are mapped
        to this height. (default=0)
    datafile : str or NoneType, optional
        Path to custom coefficient file, if None uses `apexsh.dat` file
        (default=None)
    fortranlib : str or NoneType, optional
        Path to Fortran Apex CPython library, if None uses linked library file
        (default=None)

    Attributes
    ----------
    year : float
        Decimal year used for the IGRF model
    RE : float
        Earth radius in km, defaults to mean Earth radius
    refh : float
        Reference height in km for apex coordinates
    datafile : str
        Path to coefficient file
    fortranlib : str
        Path to Fortran Apex CPython library
    igrf_fn : str
        IGRF coefficient filename

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

    # ------------------------
    # Define the magic methods

    def __init__(self, date=None, refh=0, datafile=None, fortranlib=None):

        if datafile is None:
            datafile = os.path.join(os.path.dirname(__file__), 'apexsh.dat')

        if fortranlib is None:
            fortranlib = fa.__file__

        self.RE = 6371.009  # Mean Earth radius in km
        self.set_refh(refh)  # Reference height in km

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

        # Set the IGRF coefficient text file name
        self.igrf_fn = os.path.join(os.path.dirname(__file__),
                                    'igrf13coeffs.txt')
        if not os.path.exists(self.igrf_fn):
            raise OSError("File {} does not exist".format(self.igrf_fn))

        # Update the Fortran epoch using the year defined above
        self.set_epoch(self.year)

        # Vectorize the fortran functions
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

        # Vectorize other nonvectorized functions
        self._apex2qd = np.frompyfunc(self._apex2qd_nonvectorized, 3, 2)
        self._qd2apex = np.frompyfunc(self._qd2apex_nonvectorized, 3, 2)
        self._get_babs = np.frompyfunc(self._get_babs_nonvectorized, 3, 1)

        return

    def __repr__(self):
        """Produce an evaluatable representation of the Apex class."""
        rstr = "".join(["apexpy.Apex(date=", self.year.__repr__(), ", refh=",
                        self.refh.__repr__(), ", datafile=",
                        self.datafile.__repr__(), ", fortranlib=",
                        self.fortranlib.__repr__(), ")"])

        return rstr

    def __str__(self):
        """Produce a user-friendly representation of the Apex class."""

        out_str = '\n'.join(["Apex class conversions performed with",
                             "-------------------------------------",
                             "Decimal year: {:.8f}".format(self.year),
                             "Reference height: {:.3f} km".format(self.refh),
                             "Earth radius: {:.3f} km".format(self.RE), '',
                             "Coefficient file: '{:s}'".format(self.datafile),
                             "Cython Fortran library: '{:s}'".format(
                                 self.fortranlib)])
        return out_str

    def __eq__(self, comp_obj):
        """Performs equivalency evaluation.

        Parameters
        ----------
        comp_obj
            Object of any time to be compared to the class object

        Returns
        -------
        bool or NotImplemented
            True if self and comp_obj are identical, False if they are not,
            and NotImplemented if the classes are not the same

        """

        if isinstance(comp_obj, self.__class__):
            # Objects are the same if all the attributes are the same
            for apex_attr in ['year', 'refh', 'RE', 'datafile', 'fortranlib',
                              'igrf_fn']:
                bad_attr = False
                if hasattr(self, apex_attr):
                    aval = getattr(self, apex_attr)

                    if hasattr(comp_obj, apex_attr):
                        cval = getattr(comp_obj, apex_attr)

                        if cval != aval:
                            # Not equal, as the attribute values differ
                            bad_attr = True
                    else:
                        # The comparison object is missing an attribute
                        bad_attr = True
                else:
                    if hasattr(comp_obj, apex_attr):
                        # The current object is missing an attribute
                        bad_attr = True

                if bad_attr:
                    return False

            # If we arrive here, all expected attributes exist in both classes
            # and have the same values
            return True

        return NotImplemented

    # -------------------------
    # Define the hidden methods

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
        # Evaluate the latitude
        alat = helpers.checklat(alat, name='alat')

        # Convert modified apex to quasi-dipole, longitude is the same
        qlon = alon

        # Get the apex height
        h_apex = self.get_apex(alat)

        if h_apex < height:
            if np.isclose(h_apex, height, rtol=0, atol=1e-5):
                # Allow for values that are close
                h_apex = height
            else:
                estr = ''.join(['height {:.3g} is > '.format(np.max(height)),
                                'apex height {:.3g} for alat {:.3g}'.format(
                                    h_apex, alat)])
                raise ApexHeightError(estr)

        # Convert the latitude
        salat = np.sign(alat) if alat != 0 else 1
        qlat = salat * np.degrees(np.arccos(np.sqrt((self.RE + height)
                                                    / (self.RE + h_apex))))

        return qlat, qlon

    def _qd2apex_nonvectorized(self, qlat, qlon, height):
        """Converts quasi-dipole to modified apex coordinates.

        Parameters
        ----------
        qlat : float
            Quasi-dipole latitude
        qlon : float
            Quasi-dipole longitude
        height : float
            Altitude in km

        Returns
        -------
        alat : float
            Modified apex latitude
        alon : float
            Modified apex longitude

        Raises
        ------
        ApexHeightError
            if apex height < reference height

        """
        # Evaluate the latitude
        qlat = helpers.checklat(qlat, name='qlat')

        # Get the longitude and apex height
        alon = qlon
        h_apex = self.get_apex(qlat, height)

        if h_apex < self.refh:
            if np.isclose(h_apex, self.refh, rtol=0, atol=1e-5):
                # Allow for values that are close
                h_apex = self.refh
            else:
                estr = ''.join(['apex height ({:.3g}) is < '.format(h_apex),
                                'reference height ({:.3g})'.format(self.refh),
                                ' for qlat {:.3g}'.format(qlat)])
                raise ApexHeightError(estr)

        # Convert the latitude
        sqlat = np.sign(qlat) if qlat != 0 else 1
        alat = sqlat * np.degrees(np.arccos(np.sqrt((self.RE + self.refh)
                                                    / (self.RE + h_apex))))

        return alat, alon

    def _map_EV_to_height(self, alat, alon, height, newheight, data, ev_flag):
        """Maps electric field related values to a desired height

        Parameters
        ----------
        alat : array-like
            Apex latitude in degrees.
        alon : array-like
            Apex longitude in degrees.
        height : array-like
            Current altitude in km.
        new_height : array-like
            Desired altitude to which EV values will be mapped in km.
        data : array-like
            3D value(s) for the electric field or electric drift
        ev_flag : str
           Specify if value is an electric field ('E') or electric drift ('V')

        Returns
        -------
        data_mapped : array-like
            Data mapped along the magnetic field from the old height to the
            new height.

        Raises
        ------
        ValueError
            If an unknown `ev_flag` or badly shaped data input is supplied.

        """
        # Ensure the E-V flag is correct
        ev_flag = ev_flag.upper()
        if ev_flag not in ['E', 'V']:
            raise ValueError('unknown electric field/drift flag')

        # Make sure data is array of correct shape
        if not (np.ndim(data) == 1
                and np.size(data) == 3) and not (np.ndim(data) == 2
                                                 and np.shape(data)[0] == 3):
            # Raise ValueError because if passing e.g. a (6,) ndarray the
            # reshape below will work even though the input is invalid
            raise ValueError('{:} must be (3, N) or (3,) ndarray'.format(
                ev_flag))

        data = np.reshape(data, (3, np.size(data) // 3))

        # Get the necessary base vectors at the current and new height
        v1 = list()
        v2 = list()
        for i, alt in enumerate([height, newheight]):
            _, _, _, _, _, _, d1, d2, _, e1, e2, _ = self.basevectors_apex(
                alat, alon, alt, coords='apex')

            if ev_flag == 'E' and i == 0 or ev_flag == 'V' and i == 1:
                v1.append(e1)
                v2.append(e2)
            else:
                v1.append(d1)
                v2.append(d2)

            # Make sure v1 and v2 have shape (3, N)
            v1[-1] = np.reshape(v1[-1], (3, v1[-1].size // 3))
            v2[-1] = np.reshape(v2[-1], (3, v2[-1].size // 3))

        # Take the dot product between the data value and each base vector
        # at the current height
        data1 = np.sum(data * v1[0], axis=0)
        data2 = np.sum(data * v2[0], axis=0)

        # Map the data to the new height, removing any axes of length 1
        # after the calculation
        data_mapped = np.squeeze(data1[np.newaxis, :] * v1[1]
                                 + data2[np.newaxis, :] * v2[1])

        return data_mapped

    def _get_babs_nonvectorized(self, glat, glon, height):
        """Get the absolute value of the B-field in Tesla

        Parameters
        ----------
        glat : float
            Geodetic latitude in degrees
        glon : float
            Geodetic longitude in degrees
        height : float
            Altitude in km

        Returns
        -------
        babs : float
            Absolute value of the magnetic field in Tesla

        """
        # Evaluate the latitude
        glat = helpers.checklat(glat, name='qlat')

        # Get the magnetic field output: North, East, Down, Absolute value
        bout = fa.feldg(1, glat, glon, height)

        # Convert the absolute value from Gauss to Tesla
        babs = bout[-1] / 10000.0

        return babs

    # -----------------------
    # Define the user methods

    def convert(self, lat, lon, source, dest, height=0, datetime=None,
                precision=1e-10, ssheight=50 * 6371):
        """Converts between geodetic, modified apex, quasi-dipole and MLT.

        Parameters
        ----------
        lat : array_like
            Latitude in degrees
        lon : array_like
            Longitude in degrees or MLT in hours
        source : str
            Input coordinate system, accepts 'geo', 'apex', 'qd', or 'mlt'
        dest : str
            Output coordinate system, accepts 'geo', 'apex', 'qd', or 'mlt'
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
        lon : ndarray or float
            Converted longitude or MLT

        Raises
        ------
        ValueError
            For unknown source or destination coordinate system or a missing
            or badly formed latitude or datetime input

        """
        # Test the input values
        source = source.lower()
        dest = dest.lower()

        if source not in ['geo', 'apex', 'qd', 'mlt'] \
           or dest not in ['geo', 'apex', 'qd', 'mlt']:
            estr = 'Unknown coordinate transformation: {} -> {}'.format(source,
                                                                        dest)
            raise ValueError(estr)

        if datetime is None and ('mlt' in [source, dest]):
            raise ValueError('datetime must be given for MLT calculations')

        lat = helpers.checklat(lat)

        if source == dest:
            # The source and destination are the same, no conversion necessary
            return lat, lon
        else:
            # Select the correct conversion method
            if source == 'geo':
                if dest == 'qd':
                    lat, lon = self.geo2qd(lat, lon, height)
                else:
                    lat, lon = self.geo2apex(lat, lon, height)
                    if dest == 'mlt':
                        lon = self.mlon2mlt(lon, datetime, ssheight=ssheight)
            elif source == 'apex':
                if dest == 'geo':
                    lat, lon, _ = self.apex2geo(lat, lon, height,
                                                precision=precision)
                elif dest == 'qd':
                    lat, lon = self.apex2qd(lat, lon, height=height)
                elif dest == 'mlt':
                    lon = self.mlon2mlt(lon, datetime, ssheight=ssheight)
            elif source == 'qd':
                if dest == 'geo':
                    lat, lon, _ = self.qd2geo(lat, lon, height,
                                              precision=precision)
                else:
                    lat, lon = self.qd2apex(lat, lon, height=height)
                    if dest == 'mlt':
                        lon = self.mlon2mlt(lon, datetime, ssheight=ssheight)
            elif source == 'mlt':
                # From MLT means that the input latitude is assumed to be apex,
                # so we don't need to update latitude for apex conversions.
                lon = self.mlt2mlon(lon, datetime, ssheight=ssheight)
                if dest == 'geo':
                    lat, lon, _ = self.apex2geo(lat, lon, height,
                                                precision=precision)
                elif dest == 'qd':
                    lat, lon = self.apex2qd(lat, lon, height=height)

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
            warnings.warn(''.join(['Apex latitude set to NaN where undefined ',
                                   '(apex height may be < reference height)']))
            if np.isscalar(alat):
                alat = np.nan
            else:
                alat[alat == -9999] = np.nan

        # If array is returned, dtype is object, so convert to float
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
        # Evaluate the latitude
        alat = helpers.checklat(alat, name='alat')

        # Perform the two-step convertion to geodetic coordinates
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
        # Evaluate the latitude
        glat = helpers.checklat(glat, name='glat')

        # Convert to quasi-dipole coordinates
        qlat, qlon = self._geo2qd(glat, glon, height)

        # If array is returned, dtype is object, so convert to float
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
        # Evaluate the latitude
        qlat = helpers.checklat(qlat, name='qlat')

        # Convert to geodetic coordinates
        glat, glon, error = self._qd2geo(qlat, qlon, height, precision)

        # If array is returned, dtype is object, so convert to float
        return np.float64(glat), np.float64(glon), np.float64(error)

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
        # Convert the latitude to apex, latitude check performed in the
        # hidden method
        qlat, qlon = self._apex2qd(alat, alon, height)

        # If array is returned, the dtype is object, so convert to float
        return np.float64(qlat), np.float64(qlon)

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
        # Perform the conversion from quasi-dipole to apex coordinates
        alat, alon = self._qd2apex(qlat, qlon, height)

        # If array is returned, the dtype is object, so convert to float
        return np.float64(alat), np.float64(alon)

    def mlon2mlt(self, mlon, dtime, ssheight=318550):
        """Computes the magnetic local time at the specified magnetic longitude
        and UT.

        Parameters
        ----------
        mlon : array_like
            Magnetic longitude (apex and quasi-dipole longitude are always
            equal)
        dtime : :class:`datetime.datetime`
            Date and time
        ssheight : float, optional
            Altitude in km to use for converting the subsolar point from
            geographic to magnetic coordinates. A high altitude is used
            to ensure the subsolar point is mapped to high latitudes, which
            prevents the South-Atlantic Anomaly (SAA) from influencing the MLT.
            The current default is  50 * 6371, roughly 50 RE. (default=318550)

        Returns
        -------
        mlt : ndarray or float
            Magnetic local time in hours [0, 24)

        Notes
        -----
        To compute the MLT, we find the apex longitude of the subsolar point at
        the given time. Then the MLT of the given point will be computed from
        the separation in magnetic longitude from this point (1 hour = 15
        degrees).

        """
        # Get the subsolar location
        ssglat, ssglon = helpers.subsol(dtime)

        # Convert the subsolar location to apex coordinates
        _, ssalon = self.geo2apex(ssglat, ssglon, ssheight)

        # Calculate the magnetic local time (0-24 h range) from apex longitude.
        # np.float64 will ensure lists are converted to arrays
        mlt = (180 + np.float64(mlon) - ssalon) / 15 % 24

        return mlt

    def mlt2mlon(self, mlt, dtime, ssheight=318550):
        """Computes the magnetic longitude at the specified MLT and UT.

        Parameters
        ----------
        mlt : array_like
            Magnetic local time
        dtime : :class:`datetime.datetime`
            Date and time
        ssheight : float, optional
            Altitude in km to use for converting the subsolar point from
            geographic to magnetic coordinates. A high altitude is used
            to ensure the subsolar point is mapped to high latitudes, which
            prevents the South-Atlantic Anomaly (SAA) from influencing the MLT.
            The current default is  50 * 6371, roughly 50 RE. (default=318550)

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
        # Get the location of the subsolar point at this time
        ssglat, ssglon = helpers.subsol(dtime)

        # Convert the location of the subsolar point to apex coordinates
        _, ssalon = self.geo2apex(ssglat, ssglon, ssheight)

        # Calculate the magnetic longitude (0-360 h range) from MLT.
        # np.float64 will ensure lists are converted to arrays
        mlon = (15 * np.float64(mlt) - 180 + ssalon + 360) % 360

        return mlon

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
            raise ApexHeightError("input height is > apex height")

        return newglat, newglon, error

    def map_E_to_height(self, alat, alon, height, newheight, edata):
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
        edata : (3,) or (3, N) array_like
            Electric field (at `alat`, `alon`, `height`) in geodetic east,
            north, and up components

        Returns
        -------
        out : (3, N) or (3,) ndarray
            The electric field at `newheight` (geodetic east, north, and up
            components)

        """
        # Call hidden mapping method with flag for the electric field
        out = self._map_EV_to_height(alat, alon, height, newheight, edata, 'E')

        return out

    def map_V_to_height(self, alat, alon, height, newheight, vdata):
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
        vdata : (3,) or (3, N) array_like
            Electric drift velocity (at `alat`, `alon`, `height`) in geodetic
            east, north, and up components

        Returns
        -------
        out : (3, N) or (3,) ndarray
            The electric drift velocity at `newheight` (geodetic east, north,
            and up components)

        """
        # Call hidden mapping method with flag for the electric drift velocities
        out = self._map_EV_to_height(alat, alon, height, newheight, vdata, 'V')

        return out

    def basevectors_qd(self, lat, lon, height, coords='geo', precision=1e-10):
        """Get quasi-dipole base vectors f1 and f2 at the specified coordinates.

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

        Notes
        -----
        The vectors are described by Richmond [1995] [2]_ and
        Emmert et al. [2010] [3]_.  The vector components are geodetic east and
        north.

        References
        ----------
        .. [2] Richmond, A. D. (1995), Ionospheric Electrodynamics Using
               Magnetic Apex Coordinates, Journal of geomagnetism and
               geoelectricity, 47(2), 191–212, :doi:`10.5636/jgg.47.191`.

        .. [3] Emmert, J. T., A. D. Richmond, and D. P. Drob (2010),
               A computationally compact representation of Magnetic-Apex
               and Quasi-Dipole coordinates with smooth base vectors,
               J. Geophys. Res., 115(A8), A08322, doi:10.1029/2010JA015326.

        """
        # Convert from current coordinates to geodetic coordinates
        glat, glon = self.convert(lat, lon, coords, 'geo', height=height,
                                  precision=precision)

        # Get the geodetic base vectors
        f1, f2 = self._basevec(glat, glon, height)

        # If inputs are not scalar, each vector is an array of arrays,
        # so reshape to a single array
        if f1.dtype == object:
            f1 = np.vstack(f1).T
            f2 = np.vstack(f2).T

        return f1, f2

    def basevectors_apex(self, lat, lon, height, coords='geo', precision=1e-10):
        """Returns base vectors in quasi-dipole and apex coordinates.

        Parameters
        ----------
        lat : array_like or float
            Latitude in degrees; input must be broadcastable with `lon` and
            `height`.
        lon : array_like or float
            Longitude in degrees; input must be broadcastable with `lat` and
            `height`.
        height : array_like or float
            Altitude in km; input must be broadcastable with `lon` and `lat`.
        coords : str, optional
            Input coordinate system, expects one of 'geo', 'apex', or 'qd'
            (default='geo')
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
        f1 : (3, N) or (3,) ndarray
            Quasi-dipole base vector equivalent to e1, tangent to contours of
            constant lambdaA
        f2 : (3, N) or (3,) ndarray
            Quasi-dipole base vector equivalent to e2
        f3 : (3, N) or (3,) ndarray
            Quasi-dipole base vector equivalent to e3, tangent to contours of
            PhiA
        g1 : (3, N) or (3,) ndarray
            Quasi-dipole base vector equivalent to d1
        g2 : (3, N) or (3,) ndarray
            Quasi-dipole base vector equivalent to d2
        g3 : (3, N) or (3,) ndarray
            Quasi-dipole base vector equivalent to d3
        d1 : (3, N) or (3,) ndarray
            Apex base vector normal to contours of constant PhiA
        d2 : (3, N) or (3,) ndarray
            Apex base vector that completes the right-handed system
        d3 : (3, N) or (3,) ndarray
            Apex base vector normal to contours of constant lambdaA
        e1 : (3, N) or (3,) ndarray
            Apex base vector tangent to contours of constant V0
        e2 : (3, N) or (3,) ndarray
            Apex base vector that completes the right-handed system
        e3 : (3, N) or (3,) ndarray
            Apex base vector tangent to contours of constant PhiA

        Notes
        -----
        The vectors are described by Richmond [1995] [4]_ and
        Emmert et al. [2010] [5]_.  The vector components are geodetic east,
        north, and up (only east and north for `f1` and `f2`).

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
               J. Geophys. Res., 115(A8), A08322, doi:10.1029/2010JA015326.

        """
        # Convert to geodetic coordinates from current coordinate system
        glat, glon = self.convert(lat, lon, coords, 'geo', height=height,
                                  precision=precision)

        # Retrieve the desired magnetic locations and base vectors
        returnvals = list(self._geo2apexall(glat, glon, height))
        qlat = np.float64(returnvals[0])
        alat = np.float64(returnvals[2])
        bvec_ind = [4, 5, 7, 8, 9, 11, 12, 13]

        # If inputs are not scalar, each vector is an array of arrays,
        # so reshape to a single array
        if returnvals[4].dtype == object:
            for i in bvec_ind:
                returnvals[i] = np.vstack(returnvals[i]).T

        # Make sure the base vector arrays are 2D
        for i in bvec_ind:
            if i in [4, 5]:
                rsize = (2, returnvals[i].size // 2)
            else:
                rsize = (3, returnvals[i].size // 3)
            returnvals[i] = returnvals[i].reshape(rsize)

        # Assign the reshaped base vectors
        f1, f2 = returnvals[4:6]
        d1, d2, d3 = returnvals[7:10]
        e1, e2, e3 = returnvals[11:14]

        # Compute f3, g1, g2, g3 (outstanding quasi-dipole base vectors)
        #
        # Start by calculating the D equivalent, F (F_scalar)
        f1_stack = np.vstack((f1, np.zeros_like(f1[0])))
        f2_stack = np.vstack((f2, np.zeros_like(f2[0])))
        F_scalar = np.cross(f1_stack.T, f2_stack.T).T[-1]

        # Get the cosine of the magnetic inclination
        cos_mag_inc = helpers.getcosIm(alat)

        # Define the k base vectors
        k_unit = np.array([0, 0, 1], dtype=np.float64).reshape((3, 1))

        # Calculate the remaining quasi-dipole base vectors
        g1 = ((self.RE + np.float64(height))
              / (self.RE + self.refh)) ** (3 / 2) * d1 / F_scalar
        g2 = -1.0 / (2.0 * F_scalar * np.tan(np.radians(qlat))) * (
            k_unit + ((self.RE + np.float64(height))
                      / (self.RE + self.refh)) * d2 / cos_mag_inc)
        g3 = k_unit * F_scalar
        f3 = np.cross(g1.T, g2.T).T

        # Reshape the output
        out = [f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3]

        if np.any(alat == -9999):
            warnings.warn(''.join(['Base vectors g, d, e, and f3 set to NaN ',
                                   'where apex latitude is undefined (apex ',
                                   'height may be < reference height)']))
            mask = alat == -9999
            for i in range(len(out) - 2):
                out[i + 2] = np.where(mask, np.nan, out[i + 2])

        out = tuple(np.squeeze(bvec) for bvec in out)

        return out

    def get_apex(self, lat, height=None):
        """ Calculate apex height

        Parameters
        ----------
        lat : float
            Apex latitude in degrees
        height : float or NoneType
            Height above the surface of the Earth in km or NoneType to use
            reference height (default=None)

        Returns
        -------
        apex_height : float
            Height of the field line apex in km

        """
        # Check the latitude
        lat = helpers.checklat(lat, name='alat')

        # Ensure the height is set
        if height is None:
            height = self.refh

        # Calculate the apex height
        cos_lat_squared = np.cos(np.radians(lat)) ** 2
        apex_height = (self.RE + height) / cos_lat_squared - self.RE

        return apex_height

    def get_height(self, lat, apex_height):
        """Calculate the height given an apex latitude and apex height.

        Parameters
        ----------
        lat : float
            Apex latitude in degrees
        apex_height : float
            Maximum height of the apex field line above the surface of the
            Earth in km

        Returns
        -------
        height : float
            Height above the surface of the Earth in km

        """
        # Check the latitude
        lat = helpers.checklat(lat, name='alat')

        # Calculate the height from the apex height
        cos_lat_squared = np.cos(np.radians(lat)) ** 2
        height = cos_lat_squared * (apex_height + self.RE) - self.RE

        return height

    def set_epoch(self, year):
        """Updates the epoch for all subsequent conversions.

        Parameters
        ----------
        year : float
            Decimal year

        """
        # Set the year and load the data file
        self.year = np.float64(year)
        fa.loadapxsh(self.datafile, self.year)

        # Call the Fortran routine to set time
        fa.cofrm(self.year, self.igrf_fn)
        return

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

    def get_babs(self, glat, glon, height):
        """Returns the magnitude of the IGRF magnetic field in tesla.

        Parameters
        ----------
        glat : array_like
            Geodetic latitude in degrees
        glon : array_like
            Geodetic longitude in degrees
        height : array_like
            Altitude in km

        Returns
        -------
        babs : ndarray or float
            Magnitude of the IGRF magnetic field in Tesla

        """
        # Get the absolute value of the magnetic field at the desired location
        babs = self._get_babs(glat, glon, height)

        # If array is returned, the dtype is object, so convert to float
        return np.float64(babs)

    def bvectors_apex(self, lat, lon, height, coords='geo', precision=1e-10):
        """Returns the magnetic field vectors in apex coordinates.

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
        main_mag_e3: (1, N) or (1,) ndarray
           IGRF magnitude divided by a scaling factor, D (d_scale) to give
           the main B field magnitude along the e3 base vector
        e3 : (3, N) or (3,) ndarray
           Base vector tangent to the contours of constant V_0 and Phi_A
        main_mag_d3: (1, N) or (1,) ndarray
           IGRF magnitude multiplied by a scaling factor, D (d_scale) to give
           the main B field magnitudee along the d3 base vector
        d3 : (3, N) or (3,) ndarray
           Base vector equivalent to the scaled main field unit vector

        Notes
        -----
        See Richmond, A. D. (1995) [4]_ equations 3.8-3.14

        The apex magnetic field vectors described by Richmond [1995] [4]_ and
        Emmert et al. [2010] [5]_, specfically the Be3 (main_mag_e3) and Bd3
        (main_mag_d3) components. The vector components are geodetic east,
        north, and up.

        References
        ----------
        Richmond, A. D. (1995) [4]_
        Emmert, J. T. et al. (2010) [5]_

        """
        # Convert the current coordinates to geodetic coordinates
        glat, glon = self.convert(lat, lon, coords, 'geo', height=height,
                                  precision=precision)

        # Get the magnitude of the magnetic field at the desired location
        babs = self.get_babs(glat, glon, height)

        # Retrieve the necessary base vectors
        _, _, _, _, _, _, d1, d2, d3, _, _, e3 = self.basevectors_apex(
            glat, glon, height, coords='geo')

        # Perform the calculations described in [4]
        d1_cross_d2 = np.cross(d1.T, d2.T).T
        d_scale = np.sqrt(np.sum(d1_cross_d2 ** 2, axis=0))  # D in [4] 3.13

        main_mag_e3 = babs / d_scale  # Solve for b0 in [4] 3.13
        main_mag_d3 = babs * d_scale  # Solve for b0 in [4] 3.10

        return main_mag_e3, e3, main_mag_d3, d3
