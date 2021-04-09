# -*- coding: utf-8 -*-
"""Test the apexpy.Apex class

Notes
-----
Whenever function outputs are tested against hard-coded numbers, the test
results (numbers) were obtained by running the code that is tested.  Therefore,
these tests below only check that nothing changes when refactoring, etc., and
not if the results are actually correct.

These results are expected to change when IGRF is updated.

"""

import datetime as dt
import numpy as np
import os
import pytest
import warnings

from apexpy import fortranapex as fa
from apexpy import Apex, ApexHeightError, helpers


class TestApexInit():
    def setup(self):
        self.apex_out = None
        self.test_date = dt.datetime.utcnow()
        self.test_refh = 0

    def teardown(self):
        del self.apex_out, self.test_date, self.test_refh

    def eval_date(self):
        """Evaluate the times in self.test_date and self.apex_out."""
        if isinstance(self.test_date, dt.datetime) \
           or isinstance(self.test_date, dt.date):
            self.test_date = helpers.toYearFraction(self.test_date)

        # Assert the times are the same on the order of tens of seconds.
        # Necessary to evaluate the current UTC
        np.testing.assert_almost_equal(self.test_date, self.apex_out.year, 6)
        return

    def eval_refh(self):
        """Evaluate the reference height in self.refh and self.apex_out."""
        eval_str = "".join(["expected reference height [",
                            "{:}] not equal to Apex ".format(self.test_refh),
                            "reference height ",
                            "[{:}]".format(self.apex_out.refh)])
        assert self.test_refh == self.apex_out.refh, eval_str
        return

    def test_init_defaults(self):
        """Test Apex class default initialization."""
        self.apex_out = Apex()
        self.eval_date()
        self.eval_refh()
        return

    @pytest.mark.parametrize("in_date",
                             [(2015), (2015.5), (dt.date(2015, 1, 1)),
                              (dt.datetime(2015, 6, 1, 18, 23, 45))])
    def test_init_date(self, in_date):
        """Test Apex class with date initialization."""
        self.test_date = in_date
        self.apex_out = Apex(date=self.test_date)
        self.eval_date()
        self.eval_refh()
        return

    @pytest.mark.parametrize("in_refh", [(0.0), (300.0), (30000.0), (-1.0)])
    def test_init_refh(self, in_refh):
        """Test Apex class with reference height initialization."""
        self.test_refh = in_refh
        self.apex_out = Apex(refh=self.test_refh)
        self.eval_date()
        self.eval_refh()
        return

    def test_init_with_bad_datafile(self):
        """Test raises IOError with non-existent datafile input."""
        with pytest.raises(IOError):
            Apex(datafile='foo/path/to/datafile.blah')

        return


class TestApexMethod():
    """Test the Apex methods."""
    def setup(self):
        """Initialize all tests."""
        self.apex_out = Apex(date=2000, refh=300)

    def teardown(self):
        """Clean up after each test."""
        del self.apex_out

    def get_input_args(self, method_name, lat, lon, alt, precision=0.0):
        """Set the input arguments for the different Apex methods.

        Parameters
        ----------
        method_name : str
            Name of the Apex class method
        lat : float or array-like
            Value for the latitude
        lon : float or array-like
            Value for the longitude
        alt : float or array-like
            Value for the altitude
        precision : float
            Value for the precision (default=0.0)

        Returns
        -------
        in_args : list
            List of the appropriate input arguments

        """
        in_args = [lat, lon, alt]

        # Add precision, if needed
        if method_name in ["_qd2geo", "apxq2g", "apex2geo", "qd2geo",
                           "_apex2geo"]:
            in_args.append(precision)

        # Add a reference height, if needed
        if method_name in ["apxg2all"]:
            in_args.append(300)

        # Add a vector flag, if needed
        if method_name in ["apxg2all", "apxg2q"]:
            in_args.append(1)

        return in_args

    @pytest.mark.parametrize("apex_method,fortran_method,fslice",
                             [("_geo2qd", "apxg2q", slice(0, 2, 1)),
                              ("_geo2apex", "apxg2all", slice(2, 4, 1)),
                              ("_qd2geo", "apxq2g", slice(None)),
                              ("_basevec", "apxg2q", slice(2, 4, 1))])
    @pytest.mark.parametrize("lat", [(0), (30), (60), (89)])
    @pytest.mark.parametrize("lon", [(-179), (-90), (0), (90), (180)])
    def test_fortran_scalar_input(self, apex_method, fortran_method, fslice,
                                  lat, lon):
        """Tests Apex/fortran interface consistency for scalars."""
        # Get the Apex class method and the fortran function call
        apex_func = getattr(self.apex_out, apex_method)
        fortran_func = getattr(fa, fortran_method)

        # Get the appropriate input arguments
        apex_args = self.get_input_args(apex_method, lat, lon, 100)
        fortran_args = self.get_input_args(fortran_method, lat, lon, 100)

        # Evaluate the equivalent function calls
        np.testing.assert_allclose(apex_func(*apex_args),
                                   fortran_func(*fortran_args)[fslice])
        return

    @pytest.mark.parametrize("apex_method,fortran_method,fslice",
                             [("_geo2qd", "apxg2q", slice(0, 2, 1)),
                              ("_geo2apex", "apxg2all", slice(2, 4, 1)),
                              ("_qd2geo", "apxq2g", slice(None)),
                              ("_basevec", "apxg2q", slice(2, 4, 1))])
    @pytest.mark.parametrize("lat", [(0), (30), (60), (89)])
    @pytest.mark.parametrize("lon1,lon2", [(180, 180), (-180, -180),
                                           (180, -180), (-180, 180),
                                           (-345, 15), (375, 15)])
    def test_fortran_longitude_rollover(self, apex_method, fortran_method,
                                        fslice, lat, lon1, lon2):
        """Tests Apex/fortran interface consistency for longitude rollover."""
        # Get the Apex class method and the fortran function call
        apex_func = getattr(self.apex_out, apex_method)
        fortran_func = getattr(fa, fortran_method)

        # Get the appropriate input arguments
        apex_args = self.get_input_args(apex_method, lat, lon1, 100)
        fortran_args = self.get_input_args(fortran_method, lat, lon2, 100)

        # Evaluate the equivalent function calls
        np.testing.assert_allclose(apex_func(*apex_args),
                                   fortran_func(*fortran_args)[fslice])
        return

    @pytest.mark.parametrize("apex_method,fortran_method,fslice",
                             [("_geo2qd", "apxg2q", slice(0, 2, 1)),
                              ("_geo2apex", "apxg2all", slice(2, 4, 1)),
                              ("_qd2geo", "apxq2g", slice(None)),
                              ("_basevec", "apxg2q", slice(2, 4, 1))])
    def test_fortran_array_input(self, apex_method, fortran_method, fslice):
        """Tests Apex/fortran interface consistency for array input."""
        # Get the Apex class method and the fortran function call
        apex_func = getattr(self.apex_out, apex_method)
        fortran_func = getattr(fa, fortran_method)

        # Set up the input arrays
        in_lats = np.array([0, 30, 60, 90])
        in_lon = 15
        in_alts = np.array([100, 200, 300, 400])
        apex_args = self.get_input_args(apex_method, in_lats.reshape((2, 2)),
                                        in_lon, in_alts.reshape((2, 2)))

        # Get the Apex class results
        aret = apex_func(*apex_args)

        # Get the fortran function results
        flats = list()
        flons = list()

        for i, lat in enumerate(in_lats):
            fortran_args = self.get_input_args(fortran_method, lat, in_lon,
                                               in_alts[i])
            fret = fortran_func(*fortran_args)[fslice]
            flats.append(fret[0])
            flons.append(fret[1])

        flats = np.array(flats)
        flons = np.array(flons)

        # Evaluate results
        try:
            # This returned value is array of floats
            np.testing.assert_allclose(aret[0].astype(float),
                                       flats.reshape((2, 2)).astype(float))
            np.testing.assert_allclose(aret[1].astype(float),
                                       flons.reshape((2, 2)).astype(float))
        except ValueError:
            # This returned value is array of arrays
            alats = aret[0].reshape((4,))
            alons = aret[1].reshape((4,))
            for i, flat in enumerate(flats):
                np.testing.assert_array_almost_equal(alats[i], flat, 2)
                np.testing.assert_array_almost_equal(alons[i], flons[i], 2)

        return

    @pytest.mark.parametrize("lat", [(0), (30), (60), (89)])
    @pytest.mark.parametrize("lon", [(-179), (-90), (0), (90), (180)])
    def test_geo2apexall_scalar(self, lat, lon):
        """Test Apex/fortran geo2apexall interface consistency for scalars."""
        # Get the Apex and Fortran results
        aret = self.apex_out._geo2apexall(lat, lon, 100)
        fret = fa.apxg2all(lat, lon, 100, 300, 1)

        # Evaluate each element in the results
        for aval, fval in zip(aret, fret):
            np.testing.assert_allclose(aval, fval)

    def test_geo2apexall_array(self):
        """Test Apex/fortran geo2apexall interface consistency for arrays."""
        # Set the input
        in_lats = np.array([0, 30, 60, 90])
        in_lon = 15
        in_alts = np.array([100, 200, 300, 400])

        # Get the Apex class results
        aret = self.apex_out._geo2apexall(in_lats.reshape((2, 2)), in_lon,
                                          in_alts.reshape((2, 2)))

        # For each lat/alt pair, get the Fortran results
        fret = list()
        for i, lat in enumerate(in_lats):
            fret.append(fa.apxg2all(in_lats[i], in_lon, in_alts[i], 300, 1))

        # Cycle through all returned values
        for i, ret in enumerate(aret):
            try:
                # This returned value is array of floats
                np.testing.assert_allclose(ret.astype(float),
                                           np.array([[fret[0][i], fret[1][i]],
                                                     [fret[2][i], fret[3][i]]],
                                                    dtype=float))
            except ValueError:
                # This returned value is array of arrays
                ret = ret.reshape((4,))
                for j, single_fret in enumerate(fret):
                    np.testing.assert_allclose(ret[j], single_fret[i])
        return

    @pytest.mark.parametrize("in_coord", ["geo", "apex", "qd"])
    @pytest.mark.parametrize("out_coord", ["geo", "apex", "qd"])
    def test_convert_consistency(self, in_coord, out_coord):
        """Test the self-consistency of the Apex convert method."""
        if in_coord == out_coord:
            pytest.skip("Test not needed for same src and dest coordinates")

        # Define the method name
        method_name = "2".join([in_coord, out_coord])

        # Get the method and method inputs
        convert_kwargs = {'height': 100, 'precision': 0.0}
        apex_args = self.get_input_args(method_name, 60, 15, 100)
        apex_method = getattr(self.apex_out, method_name)

        # Define the slice needed to get equivalent output from the named method
        mslice = slice(0, -1, 1) if out_coord == "geo" else slice(None)

        # Get output using convert and named method
        convert_out = self.apex_out.convert(60, 15, in_coord, out_coord,
                                            **convert_kwargs)
        method_out = apex_method(*apex_args)[mslice]

        # Compare both outputs, should be identical
        np.testing.assert_allclose(convert_out, method_out)
        return

    @pytest.mark.parametrize("bound_lat", [(90), (-90)])
    @pytest.mark.parametrize("in_coord", ["geo", "apex", "qd"])
    @pytest.mark.parametrize("out_coord", ["geo", "apex", "qd"])
    def test_convert_at_lat_boundary(self, bound_lat, in_coord, out_coord):
        """Test the conversion at the latitude boundary, with allowed excess."""
        excess_lat = np.sign(bound_lat) * (abs(bound_lat) + 1.0e-5)

        # Get the two outputs, slight tolerance outside of boundary allowed
        bound_out = self.apex_out.convert(bound_lat, 0, in_coord, out_coord)
        excess_out = self.apex_out.convert(excess_lat, 0, in_coord, out_coord)

        # Test the outputs
        np.testing.assert_allclose(excess_out, bound_out, rtol=0, atol=1e-8)
        return

    def test_convert_qd2apex_at_equator(self):
        """Test the quasi-dipole to apex conversion at the magnetic equator."""
        eq_out = self.apex_out.convert(lat=0.0, lon=0, source='qd', dest='apex',
                                       height=320.0)
        close_out = self.apex_out.convert(lat=0.001, lon=0, source='qd',
                                          dest='apex', height=320.0)
        np.testing.assert_allclose(eq_out, close_out, atol=1e-4)

    @pytest.mark.parametrize("src", ["geo", "apex", "qd"])
    @pytest.mark.parametrize("dest", ["geo", "apex", "qd"])
    def test_convert_withnan(self, src, dest):
        """Test Apex.convert success with NaN input."""
        if src == dest:
            pytest.skip("Test not needed for same src and dest coordinates")

        num_nans = 5
        in_loc = np.arange(0, 10, dtype=float)
        in_loc[:num_nans] = np.nan

        out_loc = self.apex_out.convert(in_loc, in_loc, src, dest, height=320)

        for out in out_loc:
            assert np.all(np.isnan(out[:num_nans])), "NaN output expected"
            assert np.all(np.isfinite(out[num_nans:])), "Finite output expected"

        return

    @pytest.mark.parametrize("bad_lat", [(91), (-91)])
    def test_convert_invalid_lat(self, bad_lat):
        """Test convert raises ValueError for invalid latitudes."""

        with pytest.raises(ValueError):
            self.apex_out.convert(bad_lat, 0, 'geo', 'geo')
        return

    @pytest.mark.parametrize("coords", [("foobar", "geo"), ("geo", "foobar")])
    def test_convert_invalid_transformation(self, coords):
        """Test raises NotImplementedError for bad coordinates."""
        with pytest.raises(NotImplementedError):
            self.apex_out.convert(0, 0, *coords)
        return

    @pytest.mark.parametrize("method_name, out_comp",
                             [("geo2apex",
                               (55.94841766357422, 94.10684204101562)),
                              ("apex2geo",
                               (51.476322174072266, -66.22817993164062,
                                5.727287771151168e-06)),
                              ("geo2qd",
                               (56.531288146972656, 94.10684204101562)),
                              ("apex2qd", (60.498401178276744, 15.0)),
                              ("qd2apex", (59.49138097045895, 15.0))])
    def test_method_scalar_input(self, method_name, out_comp):
        """Test the user method against set values with scalars."""
        # Get the desired methods
        user_method = getattr(self.apex_out, method_name)

        # Get the user output
        user_out = user_method(60, 15, 100)

        # Evaluate the user output
        np.testing.assert_allclose(user_out, out_comp)

        for out_val in user_out:
            assert np.asarray(out_val).shape == (), "output is not a scalar"
        return

    @pytest.mark.parametrize("in_coord", ["geo", "apex", "qd"])
    @pytest.mark.parametrize("out_coord", ["geo", "apex", "qd"])
    @pytest.mark.parametrize("method_args, out_shape",
                             [([[60, 60], 15, 100], (2,)),
                              ([60, [15, 15], 100], (2,)),
                              ([60, 15, [100, 100]], (2,)),
                              ([[50, 60], [15, 16], [100, 200]], (2,))])
    def test_method_broadcast_input(self, in_coord, out_coord, method_args,
                                    out_shape):
        """Test the user method with inputs that require some broadcasting."""
        if in_coord == out_coord:
            pytest.skip("Test not needed for same src and dest coordinates")

        # Get the desired methods
        method_name = "2".join([in_coord, out_coord])
        user_method = getattr(self.apex_out, method_name)

        # Get the user output
        user_out = user_method(*method_args)

        # Evaluate the user output
        for out_val in user_out:
            assert hasattr(out_val, 'shape'), "output coordinate isn't np.array"
            assert out_val.shape == out_shape
        return

    @pytest.mark.parametrize("in_coord", ["geo", "apex", "qd"])
    @pytest.mark.parametrize("out_coord", ["geo", "apex", "qd"])
    @pytest.mark.parametrize("bad_lat", [(91), (-91)])
    def test_method_invalid_lat(self, in_coord, out_coord, bad_lat):
        """Test convert raises ValueError for invalid latitudes."""
        if in_coord == out_coord:
            pytest.skip("Test not needed for same src and dest coordinates")

        # Get the desired methods
        method_name = "2".join([in_coord, out_coord])
        user_method = getattr(self.apex_out, method_name)

        with pytest.raises(ValueError):
            user_method(bad_lat, 15, 100)
        return

    @pytest.mark.parametrize("in_coord", ["geo", "apex", "qd"])
    @pytest.mark.parametrize("out_coord", ["geo", "apex", "qd"])
    @pytest.mark.parametrize("bound_lat", [(90), (-90)])
    def test_method_at_lat_boundary(self, in_coord, out_coord, bound_lat):
        """Test user methods at the latitude boundary, with allowed excess."""
        if in_coord == out_coord:
            pytest.skip("Test not needed for same src and dest coordinates")

        # Get the desired methods
        method_name = "2".join([in_coord, out_coord])
        user_method = getattr(self.apex_out, method_name)

        # Get a latitude just beyond the limit
        excess_lat = np.sign(bound_lat) * (abs(bound_lat) + 1.0e-5)

        # Get the two outputs, slight tolerance outside of boundary allowed
        bound_out = user_method(bound_lat, 0, 100)
        excess_out = user_method(excess_lat, 0, 100)

        # Test the outputs
        np.testing.assert_allclose(excess_out, bound_out, rtol=0, atol=1e-8)
        return

    def test_geo2apex_undefined_warning(self):
        """Test geo2apex warning and fill values for an undefined location."""

        # Update the apex object
        self.apex_out = Apex(date=2000, refh=10000)

        # Get the output and the warnings
        with warnings.catch_warnings(record=True) as warn_rec:
            user_lat, user_lon = self.apex_out.geo2apex(0, 0, 0)

        assert np.isnan(user_lat)
        assert np.isfinite(user_lon)
        assert len(warn_rec) == 1
        assert issubclass(warn_rec[-1].category, UserWarning)
        assert 'latitude set to NaN where' in str(warn_rec[-1].message)
        return

    @pytest.mark.parametrize("method_name", ["apex2qd", "qd2apex"])
    def test_quasidipole_apexheight_close(self, method_name):
        """Test quasi-dipole success with a height close to the reference."""
        qd_method = getattr(self.apex_out, method_name)
        in_args = [0, 15, self.apex_out.refh + 1e-6]
        out_coords = qd_method(*in_args)

        for i, out_val in enumerate(out_coords):
            np.testing.assert_almost_equal(out_val, in_args[i], decimal=3)
        return

    @pytest.mark.parametrize("method_name, hinc, msg",
                             [("apex2qd", 1.0, "is > apex height"),
                              ("qd2apex", -1.0, "is < reference height")])
    def test_quasidipole_raises_apexheight(self, method_name, hinc, msg):
        """Quasi-dipole raises ApexHeightError when height above reference."""
        qd_method = getattr(self.apex_out, method_name)

        with pytest.raises(ApexHeightError) as aerr:
            qd_method(0, 15, self.apex_out.refh + hinc)

        assert str(aerr).find(msg) > 0
        return


class TestApexMLTMethods():
    """Test the Apex Magnetic Local Time (MLT) methods."""
    def setup(self):
        """Initialize all tests."""
        self.apex_out = Apex(date=2000, refh=300)
        self.in_time = dt.datetime(2000, 2, 3, 4, 5, 6)

    def teardown(self):
        """Clean up after each test."""
        del self.apex_out, self.in_time

    @pytest.mark.parametrize("in_coord", ["geo", "apex", "qd"])
    def test_convert_to_mlt(self, in_coord):
        """Test the conversions to MLT using Apex convert."""

        # Get the magnetic longitude from the appropriate method
        if in_coord == "geo":
            apex_method = getattr(self.apex_out, "{:s}2apex".format(in_coord))
            mlon = apex_method(60, 15, 100)[1]
        else:
            mlon = 15

        # Get the output MLT values
        convert_mlt = self.apex_out.convert(60, 15, in_coord, 'mlt',
                                            height=100, ssheight=2e5,
                                            datetime=self.in_time)[1]
        method_mlt = self.apex_out.mlon2mlt(mlon, self.in_time, ssheight=2e5)

        # Test the outputs
        np.testing.assert_allclose(convert_mlt, method_mlt)
        return

    @pytest.mark.parametrize("out_coord", ["geo", "apex", "qd"])
    def test_convert_mlt_to_lon(self, out_coord):
        """Test the conversions from MLT using Apex convert."""
        # Get the output longitudes
        convert_out = self.apex_out.convert(60, 15, 'mlt', out_coord,
                                            height=100, ssheight=2e5,
                                            datetime=self.in_time,
                                            precision=1e-2)
        mlon = self.apex_out.mlt2mlon(15, self.in_time, ssheight=2e5)

        if out_coord == "geo":
            method_out = self.apex_out.apex2geo(60, mlon, 100,
                                                precision=1e-2)[:-1]
        elif out_coord == "qd":
            method_out = self.apex_out.apex2qd(60, mlon, 100)
        else:
            method_out = (60, mlon)

        # Evaluate the outputs
        np.testing.assert_allclose(convert_out, method_out)
        return

    def test_convert_geo2mlt_nodate(self):
        """Test convert from geo to MLT raises ValueError with no datetime."""
        with pytest.raises(ValueError):
            self.apex_out.convert(60, 15, 'geo', 'mlt')
        return

    @pytest.mark.parametrize("mlon_kwargs,test_mlt",
                             [({}, 23.019629923502603),
                              ({"ssheight": 100000}, 23.026712036132814)])
    def test_mlon2mlt_scalar_inputs(self, mlon_kwargs, test_mlt):
        """Test mlon2mlt with scalar inputs."""
        mlt = self.apex_out.mlon2mlt(0, self.in_time, **mlon_kwargs)

        np.testing.assert_allclose(mlt, test_mlt)
        assert np.asarray(mlt).shape == ()
        return

    @pytest.mark.parametrize("mlt_kwargs,test_mlon",
                             [({}, 14.705535888671875),
                              ({"ssheight": 100000}, 14.599319458007812)])
    def test_mlt2mlon_scalar_inputs(self, mlt_kwargs, test_mlon):
        """Test mlt2mlon with scalar inputs."""
        mlon = self.apex_out.mlt2mlon(0, self.in_time, **mlt_kwargs)

        np.testing.assert_allclose(mlon, test_mlon)
        assert np.asarray(mlon).shape == ()
        return

    @pytest.mark.parametrize("mlon,test_mlt",
                             [([0, 180], [23.019261, 11.019261]),
                              (np.array([0, 180]), [23.019261, 11.019261]),
                              ([[0, 180], [0, 180]], [[23.019261, 11.019261],
                                                      [23.019261, 11.019261]]),
                              (range(0, 361, 30),
                               [23.01963, 1.01963, 3.01963, 5.01963, 7.01963,
                                9.01963, 11.01963, 13.01963, 15.01963, 17.01963,
                                19.01963, 21.01963, 23.01963])])
    def test_mlon2mlt_array(self, mlon, test_mlt):
        """Test mlon2mlt with array inputs."""
        mlt = self.apex_out.mlon2mlt(mlon, self.in_time)

        assert mlt.shape == np.asarray(test_mlt).shape
        np.testing.assert_allclose(mlt, test_mlt, rtol=1e-4)
        return

    @pytest.mark.parametrize("mlt,test_mlon",
                             [([0, 12], [14.705551, 194.705551]),
                              (np.array([0, 12]), [14.705551, 194.705551]),
                              ([[0, 12], [0, 12]], [[14.705551, 194.705551],
                                                    [14.705551, 194.705551]]),
                              (range(0, 25, 2),
                               [14.705551, 44.705551, 74.705551, 104.705551,
                                134.705551, 164.705551, 194.705551, 224.705551,
                                254.705551, 284.705551, 314.705551, 344.705551,
                                14.705551])])
    def test_mlt2mlon_array(self, mlt, test_mlon):
        """Test mlt2mlon with array inputs."""
        mlon = self.apex_out.mlt2mlon(mlt, self.in_time)

        assert mlon.shape == np.asarray(test_mlon).shape
        np.testing.assert_allclose(mlon, test_mlon, rtol=1e-4)
        return

    @pytest.mark.parametrize("method_name", ["mlon2mlt", "mlt2mlon"])
    def test_mlon2mlt_diffdates(self, method_name):
        """Test that MLT varies with universal time."""
        apex_method = getattr(self.apex_out, method_name)
        mlt1 = apex_method(0, self.in_time)
        mlt2 = apex_method(0, self.in_time + dt.timedelta(hours=1))

        assert mlt1 != mlt2
        return

    @pytest.mark.parametrize("mlt_offset", [1.0, 10.0])
    def test_mlon2mlt_offset(self, mlt_offset):
        """Test the time wrapping logic for the MLT."""
        mlt1 = self.apex_out.mlon2mlt(0.0, self.in_time)
        mlt2 = self.apex_out.mlon2mlt(-15.0 * mlt_offset,
                                      self.in_time) + mlt_offset

        np.testing.assert_allclose(mlt1, mlt2)
        return

    @pytest.mark.parametrize("mlon_offset", [15.0, 150.0])
    def test_mlt2mlon_offset(self, mlon_offset):
        """Test the time wrapping logic for the magnetic longitude."""
        mlon1 = self.apex_out.mlt2mlon(0, self.in_time)
        mlon2 = self.apex_out.mlt2mlon(mlon_offset / 15.0,
                                       self.in_time) - mlon_offset

        np.testing.assert_allclose(mlon1, mlon2)
        return

    @pytest.mark.parametrize("order", [["mlt", "mlon"], ["mlon", "mlt"]])
    @pytest.mark.parametrize("start_val", [0, 6, 12, 18, 22])
    def test_convert_and_return(self, order, start_val):
        """Test the conversion to magnetic longitude or MLT and back again."""
        first_method = getattr(self.apex_out, "2".join(order))
        second_method = getattr(self.apex_out, "2".join([order[1], order[0]]))

        middle_val = first_method(start_val, self.in_time)
        end_val = second_method(middle_val, self.in_time)

        np.testing.assert_allclose(start_val, end_val)
        return


class TestApexMapMethods():
    """Test the Apex height mapping methods."""
    def setup(self):
        """Initialize all tests."""
        self.apex_out = Apex(date=2000, refh=300)

    def teardown(self):
        """Clean up after each test."""
        del self.apex_out

    @pytest.mark.parametrize("in_args,test_mapped",
                             [([60, 15, 100, 10000],
                               [31.841466903686523, 17.916635513305664,
                                1.7075473124350538e-6]),
                              ([30, 170, 100, 500, False, 1e-2],
                               [25.727270126342773, 169.60546875,
                                0.00017573432705830783]),
                              ([60, 15, 100, 10000, True],
                               [-25.424888610839844, 27.310426712036133,
                                1.2074182222931995e-6]),
                              ([30, 170, 100, 500, True, 1e-2],
                               [-13.76642894744873, 164.24259948730469,
                                0.00056820799363777041])])
    def test_map_to_height(self, in_args, test_mapped):
        """Test the map_to_height function."""
        mapped = self.apex_out.map_to_height(*in_args)
        np.testing.assert_allclose(mapped, test_mapped, atol=1e-6)
        return

    def test_map_to_height_same_height(self):
        """Test the map_to_height function when mapping to same height."""
        mapped = self.apex_out.map_to_height(60, 15, 100, 100, conjugate=False,
                                             precision=1e-10)
        np.testing.assert_allclose(mapped, (60.0, 15.000003814697266, 0.0),
                                   rtol=1e-5)
        return

    @pytest.mark.parametrize('ivec', range(0, 4))
    def test_map_to_height_vectorization(self, ivec):
        """Test map_to_height with array input."""
        # Set the base input and output values
        in_args = [60, 15, 100, 100]
        test_mapped = np.full(shape=(2, 3),
                              fill_value=[60, 15.00000381, 0.0]).transpose()

        # Update inputs for one vectorized value
        in_args[ivec] = [in_args[ivec], in_args[ivec]]

        # Calculate and test function
        mapped = self.apex_out.map_to_height(*in_args)
        np.testing.assert_allclose(mapped, test_mapped, rtol=1e-5)
        return

    def test_map_to_height_ApexHeightError(self):
        """Test map_to_height raises ApexHeightError."""

        with pytest.raises(ApexHeightError) as aerr:
            self.apex_out.map_to_height(0, 15, 100, 10000)

        assert aerr.match("newheight is > apex height")
        return

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
    np.testing.assert_allclose(
        apex_out.map_E_to_height(60, 15, 100, 500, [1, 2, 3]),
        out_60_15_100_500, rtol=1e-5)
    np.testing.assert_allclose(
        apex_out.map_E_to_height(60, 15, 100, 500, [2, 3, 4]),
        out_60_15_100_500_234, rtol=1e-5)
    np.testing.assert_allclose(
        apex_out.map_E_to_height(60, 15, 100, 1000, [1, 2, 3]),
        out_60_15_100_1000, rtol=1e-5)
    np.testing.assert_allclose(
        apex_out.map_E_to_height(60, 15, 200, 500, [1, 2, 3]),
        out_60_15_200_500, rtol=1e-5)
    np.testing.assert_allclose(
        apex_out.map_E_to_height(60, 30, 100, 500, [1, 2, 3]),
        out_60_30_100_500, rtol=1e-5)
    np.testing.assert_allclose(
        apex_out.map_E_to_height(70, 15, 100, 500, [1, 2, 3]),
        out_70_15_100_500, rtol=1e-5)

    # vectorize lat
    np.testing.assert_allclose(
        apex_out.map_E_to_height([60, 70], 15, 100, 500,
                                 np.array([[1, 2, 3]] * 2).T),
        np.array([out_60_15_100_500, out_70_15_100_500]).T,
        rtol=1e-5)

    # vectorize lon
    np.testing.assert_allclose(
        apex_out.map_E_to_height(60, [15, 30], 100, 500,
                                 np.array([[1, 2, 3]] * 2).T),
        np.array([out_60_15_100_500, out_60_30_100_500]).T,
        rtol=1e-5)

    # vectorize height
    np.testing.assert_allclose(
        apex_out.map_E_to_height(60, 15, [100, 200], 500,
                                 np.array([[1, 2, 3]] * 2).T),
        np.array([out_60_15_100_500, out_60_15_200_500]).T,
        rtol=1e-5)

    # vectorize newheight
    np.testing.assert_allclose(
        apex_out.map_E_to_height(60, 15, 100, [500, 1000],
                                 np.array([[1, 2, 3]] * 2).T),
        np.array([out_60_15_100_500, out_60_15_100_1000]).T,
        rtol=1e-5)

    # vectorize E
    np.testing.assert_allclose(
        apex_out.map_E_to_height(60, 15, 100, 500,
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
    np.testing.assert_allclose(
        apex_out.map_V_to_height(60, 15, 100, 500, [1, 2, 3]),
        out_60_15_100_500, rtol=1e-5)
    np.testing.assert_allclose(
        apex_out.map_V_to_height(60, 15, 100, 500, [2, 3, 4]),
        out_60_15_100_500_234, rtol=1e-5)
    np.testing.assert_allclose(
        apex_out.map_V_to_height(60, 15, 100, 1000, [1, 2, 3]),
        out_60_15_100_1000, rtol=1e-5)
    np.testing.assert_allclose(
        apex_out.map_V_to_height(60, 15, 200, 500, [1, 2, 3]),
        out_60_15_200_500, rtol=1e-5)
    np.testing.assert_allclose(
        apex_out.map_V_to_height(60, 30, 100, 500, [1, 2, 3]),
        out_60_30_100_500, rtol=1e-5)
    np.testing.assert_allclose(
        apex_out.map_V_to_height(70, 15, 100, 500, [1, 2, 3]),
        out_70_15_100_500, rtol=1e-5)

    # vectorize lat
    np.testing.assert_allclose(
        apex_out.map_V_to_height([60, 70], 15, 100, 500,
                                 np.array([[1, 2, 3]] * 2).T),
        np.array([out_60_15_100_500, out_70_15_100_500]).T,
        rtol=1e-5)

    # vectorize lon
    np.testing.assert_allclose(
        apex_out.map_V_to_height(60, [15, 30], 100, 500,
                                 np.array([[1, 2, 3]] * 2).T),
        np.array([out_60_15_100_500, out_60_30_100_500]).T,
        rtol=1e-5)

    # vectorize height
    np.testing.assert_allclose(
        apex_out.map_V_to_height(60, 15, [100, 200], 500,
                                 np.array([[1, 2, 3]] * 2).T),
        np.array([out_60_15_100_500, out_60_15_200_500]).T,
        rtol=1e-5)

    # vectorize newheight
    np.testing.assert_allclose(
        apex_out.map_V_to_height(60, 15, 100, [500, 1000],
                                 np.array([[1, 2, 3]] * 2).T),
        np.array([out_60_15_100_500, out_60_15_100_1000]).T,
        rtol=1e-5)

    # vectorize E
    np.testing.assert_allclose(
        apex_out.map_V_to_height(60, 15, 100, 500,
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
    np.testing.assert_allclose(
        apex_out.basevectors_qd(60, 15, 100, coords='geo'),
        apex_out._basevec(60, 15, 100))


def test_basevectors_qd_scalar_apex():
    apex_out = Apex(date=2000, refh=300)
    glat, glon, _ = apex_out.apex2geo(60, 15, 100, precision=1e-2)
    np.testing.assert_allclose(
        apex_out.basevectors_qd(60, 15, 100, coords='apex',
                                precision=1e-2),
        apex_out._basevec(glat, glon, 100))


def test_basevectors_qd_scalar_qd():
    apex_out = Apex(date=2000, refh=300)
    glat, glon, _ = apex_out.qd2geo(60, 15, 100, precision=1e-2)
    np.testing.assert_allclose(
        apex_out.basevectors_qd(60, 15, 100, coords='qd',
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
    np.testing.assert_allclose(f1[:, 0], f1_lat0)
    np.testing.assert_allclose(f2[:, 0], f2_lat0)
    np.testing.assert_allclose(f1[:, 1], f1_lat30)
    np.testing.assert_allclose(f2[:, 1], f2_lat30)


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

    np.testing.assert_allclose(f1, f1_)
    np.testing.assert_allclose(f2, f2_)
    np.testing.assert_allclose(d1, d1_)
    np.testing.assert_allclose(d2, d2_)
    np.testing.assert_allclose(d3, d3_)
    np.testing.assert_allclose(e1, e1_)
    np.testing.assert_allclose(e2, e2_)
    np.testing.assert_allclose(e3, e3_)


def test_basevectors_apex_scalar_apex():
    apex_out = Apex(date=2000, refh=300)

    (f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2,
     e3) = apex_out.basevectors_apex(60, 15, 100, coords='apex', precision=1e-2)

    glat, glon, _ = apex_out.apex2geo(60, 15, 100, precision=1e-2)
    (_, _, _, _, f1_, f2_, _, d1_, d2_, d3_, _, e1_, e2_,
     e3_) = apex_out._geo2apexall(glat, glon, 100)

    np.testing.assert_allclose(f1, f1_)
    np.testing.assert_allclose(f2, f2_)
    np.testing.assert_allclose(d1, d1_)
    np.testing.assert_allclose(d2, d2_)
    np.testing.assert_allclose(d3, d3_)
    np.testing.assert_allclose(e1, e1_)
    np.testing.assert_allclose(e2, e2_)
    np.testing.assert_allclose(e3, e3_)


def test_basevectors_apex_scalar_qd():
    apex_out = Apex(date=2000, refh=300)

    (f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2,
     e3) = apex_out.basevectors_apex(60, 15, 100, coords='qd', precision=1e-2)

    glat, glon, _ = apex_out.qd2geo(60, 15, 100, precision=1e-2)
    (_, _, _, _, f1_, f2_, _, d1_, d2_, d3_, _, e1_, e2_,
     e3_) = apex_out._geo2apexall(glat, glon, 100)

    np.testing.assert_allclose(f1, f1_)
    np.testing.assert_allclose(f2, f2_)
    np.testing.assert_allclose(d1, d1_)
    np.testing.assert_allclose(d2, d2_)
    np.testing.assert_allclose(d3, d3_)
    np.testing.assert_allclose(e1, e1_)
    np.testing.assert_allclose(e2, e2_)
    np.testing.assert_allclose(e3, e3_)


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

    np.testing.assert_allclose(f1[:, 0], f1_1)
    np.testing.assert_allclose(f2[:, 0], f2_1)
    np.testing.assert_allclose(d1[:, 0], d1_1)
    np.testing.assert_allclose(d2[:, 0], d2_1)
    np.testing.assert_allclose(d3[:, 0], d3_1)
    np.testing.assert_allclose(e1[:, 0], e1_1)
    np.testing.assert_allclose(e2[:, 0], e2_1)
    np.testing.assert_allclose(e3[:, 0], e3_1)

    np.testing.assert_allclose(f3[:, 0],
                               np.array([-0.088671, -0.018272, 0.993576]),
                               rtol=1e-4)
    np.testing.assert_allclose(g1[:, 0],
                               np.array([0.903098, 0.245273, 0.085107]),
                               rtol=1e-4)
    np.testing.assert_allclose(g2[:, 0],
                               np.array([-0.103495, 1.072078, 0.01048]),
                               rtol=1e-4)
    np.testing.assert_allclose(g3[:, 0], np.array([0, 0, 1.006465]), rtol=1e-4)

    np.testing.assert_allclose(f1[:, 1], f1_2)
    np.testing.assert_allclose(f2[:, 1], f2_2)
    np.testing.assert_allclose(d1[:, 1], d1_2)
    np.testing.assert_allclose(d2[:, 1], d2_2)
    np.testing.assert_allclose(d3[:, 1], d3_2)
    np.testing.assert_allclose(e1[:, 1], e1_2)
    np.testing.assert_allclose(e2[:, 1], e2_2)
    np.testing.assert_allclose(e3[:, 1], e3_2)

    np.testing.assert_allclose(f3[:, 1],
                               np.array([-0.085415, -0.021176, 0.989645]),
                               rtol=1e-4)
    np.testing.assert_allclose(g1[:, 1],
                               np.array([0.902695, 0.246919, 0.083194]),
                               rtol=1e-4)
    np.testing.assert_allclose(g2[:, 1],
                               np.array([-0.11051, 1.066094, 0.013274]),
                               rtol=1e-4)
    np.testing.assert_allclose(g3[:, 1],
                               np.array([0, 0, 1.010463]), rtol=1e-4)


# test scalar return values

def test_basevectors_apex_scalar():
    apex_out = Apex(date=2000, refh=300)

    (f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2,
     e3) = apex_out.basevectors_apex(0, 15, 100, coords='geo')
    (_, _, _, _, f1_1, f2_1, _, d1_1, d2_1, d3_1, _, e1_1, e2_1,
     e3_1) = apex_out._geo2apexall(0, 15, 100)

    np.testing.assert_allclose(f1, f1_1)
    np.testing.assert_allclose(f2, f2_1)
    np.testing.assert_allclose(d1, d1_1)
    np.testing.assert_allclose(d2, d2_1)
    np.testing.assert_allclose(d3, d3_1)
    np.testing.assert_allclose(e1, e1_1)
    np.testing.assert_allclose(e2, e2_1)
    np.testing.assert_allclose(e3, e3_1)

    np.testing.assert_allclose(f3, np.array([0.092637, -0.245951, 0.938848]),
                               rtol=1e-4)
    np.testing.assert_allclose(g1, np.array([0.939012, 0.073416, -0.07342]),
                               rtol=1e-4)
    np.testing.assert_allclose(g2, np.array([0.055389, 1.004155, 0.257594]),
                               rtol=1e-4)
    np.testing.assert_allclose(g3, np.array([0, 0, 1.065135]), rtol=1e-4)


# test 1D array return values

def test_basevectors_apex_array():
    apex_out = Apex(date=2000, refh=300)
    (f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2,
     e3) = apex_out.basevectors_apex([0, 30], 15, 100, coords='geo')
    (_, _, _, _, f1_1, f2_1, _, d1_1, d2_1, d3_1, _, e1_1, e2_1,
     e3_1) = apex_out._geo2apexall(0, 15, 100)
    (_, _, _, _, f1_2, f2_2, _, d1_2, d2_2, d3_2, _, e1_2, e2_2,
     e3_2) = apex_out._geo2apexall(30, 15, 100)

    np.testing.assert_allclose(f1[:, 0], f1_1)
    np.testing.assert_allclose(f2[:, 0], f2_1)
    np.testing.assert_allclose(d1[:, 0], d1_1)
    np.testing.assert_allclose(d2[:, 0], d2_1)
    np.testing.assert_allclose(d3[:, 0], d3_1)
    np.testing.assert_allclose(e1[:, 0], e1_1)
    np.testing.assert_allclose(e2[:, 0], e2_1)
    np.testing.assert_allclose(e3[:, 0], e3_1)

    np.testing.assert_allclose(f3[:, 0],
                               np.array([0.092637, -0.245951, 0.938848]),
                               rtol=1e-4)
    np.testing.assert_allclose(g1[:, 0],
                               np.array([0.939012, 0.073416, -0.07342]),
                               rtol=1e-4)
    np.testing.assert_allclose(g2[:, 0],
                               np.array([0.055389, 1.004155, 0.257594]),
                               rtol=1e-4)
    np.testing.assert_allclose(g3[:, 0],
                               np.array([0, 0, 1.065135]), rtol=1e-4)

    np.testing.assert_allclose(f1[:, 1], f1_2)
    np.testing.assert_allclose(f2[:, 1], f2_2)
    np.testing.assert_allclose(d1[:, 1], d1_2)
    np.testing.assert_allclose(d2[:, 1], d2_2)
    np.testing.assert_allclose(d3[:, 1], d3_2)
    np.testing.assert_allclose(e1[:, 1], e1_2)
    np.testing.assert_allclose(e2[:, 1], e2_2)
    np.testing.assert_allclose(e3[:, 1], e3_2)

    np.testing.assert_allclose(f3[:, 1],
                               np.array([-0.036618, -0.071019, 0.861604]),
                               rtol=1e-4)
    np.testing.assert_allclose(g1[:, 1],
                               np.array([0.844391, 0.015353, 0.037152]),
                               rtol=1e-4)
    np.testing.assert_allclose(g2[:, 1],
                               np.array([0.050808, 1.02131, 0.086342]),
                               rtol=1e-4)
    np.testing.assert_allclose(g3[:, 1], np.array([0, 0, 1.160625]), rtol=1e-4)


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
                np.testing.assert_allclose(np.sum(f[i] * g[j]), delta,
                                           rtol=0, atol=1e-5)
                np.testing.assert_allclose(np.sum(d[i] * e[j]), delta,
                                           rtol=0, atol=1e-5)


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
            np.testing.assert_allclose(bvec, invalid)


# ============================================================================
#  Test the get_apex() method
# ============================================================================


def test_get_apex():
    apex_out = Apex(date=2000, refh=300)
    np.testing.assert_allclose(apex_out.get_apex(10), 507.409702543805)
    np.testing.assert_allclose(apex_out.get_apex(60), 20313.026999999987)


def test_get_apex_invalid_lat():
    apex_out = Apex(date=2000, refh=300)
    with pytest.raises(ValueError):
        apex_out.get_apex(91)
    with pytest.raises(ValueError):
        apex_out.get_apex(-91)
    apex_out.get_apex(90)
    apex_out.get_apex(-90)

    np.testing.assert_allclose(apex_out.get_apex(90 + 1e-5),
                               apex_out.get_apex(90),
                               rtol=0, atol=1e-8)


# ============================================================================
#  Test the set_epoch() method
# ============================================================================


def test_set_epoch():
    """Test successful setting of Apex epoch."""
    apex_out = Apex(date=2000.2, refh=300)
    np.testing.assert_allclose(apex_out.year, 2000.2)
    ret_2000_2_py = apex_out._geo2apex(60, 15, 100)
    apex_out.set_epoch(2000.8)
    np.testing.assert_allclose(apex_out.year, 2000.8)
    ret_2000_8_py = apex_out._geo2apex(60, 15, 100)

    assert ret_2000_2_py != ret_2000_8_py

    fa.loadapxsh(apex_out.datafile, 2000.2)
    ret_2000_2_apex = fa.apxg2all(60, 15, 100, 300, 0)[2:4]
    fa.loadapxsh(apex_out.datafile, 2000.8)
    ret_2000_8_apex = fa.apxg2all(60, 15, 100, 300, 0)[2:4]

    assert ret_2000_2_apex != ret_2000_8_apex

    np.testing.assert_allclose(ret_2000_2_py, ret_2000_2_apex)
    np.testing.assert_allclose(ret_2000_8_py, ret_2000_8_apex)


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

    np.testing.assert_allclose(ret_300, fa.apxg2all(60, 15, 100, 300, 0)[2:4])
    np.testing.assert_allclose(ret_500, fa.apxg2all(60, 15, 100, 500, 0)[2:4])


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
            np.testing.assert_allclose(output, expected[i][j], rtol=0,
                                       atol=1e-5)


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
            np.testing.assert_allclose(output.ravel()[j],
                                       expected[i].ravel()[j], rtol=0,
                                       atol=1e-5)
