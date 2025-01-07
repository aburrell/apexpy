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
import shutil
import warnings

import apexpy


class TestApexInit(object):
    """Test class for the Apex class object."""

    def setup_method(self):
        """Initialize all tests."""
        self.apex_out = None
        self.test_date = dt.datetime.now(tz=dt.timezone.utc)
        self.test_refh = 0
        self.bad_file = 'foo/path/to/datafile.blah'
        self.temp_file = None

    def teardown_method(self):
        """Clean up after each test."""
        del self.apex_out, self.test_date, self.test_refh
        del self.bad_file, self.temp_file

    def eval_date(self):
        """Evaluate the times in self.test_date and self.apex_out."""
        if isinstance(self.test_date, dt.datetime) \
           or isinstance(self.test_date, dt.date):
            self.test_date = apexpy.helpers.toYearFraction(self.test_date)

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
        self.apex_out = apexpy.Apex()
        self.eval_date()
        self.eval_refh()
        return

    def test_init_today(self):
        """Test Apex class initialization with today's date."""
        self.apex_out = apexpy.Apex(date=self.test_date)
        self.eval_date()
        self.eval_refh()
        return

    @pytest.mark.parametrize("in_date",
                             [2015, 2015.5, dt.date(2015, 1, 1),
                              dt.datetime(2015, 6, 1, 18, 23, 45)])
    def test_init_date(self, in_date):
        """Test Apex class with date initialization.

        Parameters
        ----------
        in_date : int, float, dt.date, or dt.datetime
            Input date in a variety of formats

        """
        self.test_date = in_date
        self.apex_out = apexpy.Apex(date=self.test_date)
        self.eval_date()
        self.eval_refh()
        return

    @pytest.mark.parametrize("new_date", [2015, 2015.5])
    def test_set_epoch(self, new_date):
        """Test successful setting of Apex epoch after initialization.

        Parameters
        ----------
        new_date : int or float
            New date for the Apex class

        """
        # Evaluate the default initialization
        self.apex_out = apexpy.Apex()
        self.eval_date()
        self.eval_refh()

        # Update the epoch
        ref_apex = eval(self.apex_out.__repr__())
        self.apex_out.set_epoch(new_date)
        assert ref_apex != self.apex_out
        self.test_date = new_date
        self.eval_date()
        return

    @pytest.mark.parametrize("in_refh", [0.0, 300.0, 30000.0, -1.0])
    def test_init_refh(self, in_refh):
        """Test Apex class with reference height initialization.

        Parameters
        ----------
        in_refh : float
            Input reference height in km

        """
        self.test_refh = in_refh
        self.apex_out = apexpy.Apex(refh=self.test_refh)
        self.eval_date()
        self.eval_refh()
        return

    @pytest.mark.parametrize("new_refh", [0.0, 300.0, 30000.0, -1.0])
    def test_set_refh(self, new_refh):
        """Test the method used to set the reference height after the init.

        Parameters
        ----------
        new_refh : float
            Reference height in km

        """
        # Verify the defaults are set
        self.apex_out = apexpy.Apex(date=self.test_date)
        self.eval_date()
        self.eval_refh()

        # Update to a new reference height and test
        ref_apex = eval(self.apex_out.__repr__())
        self.apex_out.set_refh(new_refh)

        if self.test_refh == new_refh:
            assert ref_apex == self.apex_out
        else:
            assert ref_apex != self.apex_out
            self.test_refh = new_refh
        self.eval_refh()
        return

    def copy_file(self, original, max_attempts=100):
        """Make a temporary copy of input file."""
        _, ext = os.path.splitext(original)
        self.temp_file = 'temp' + ext

        for _ in range(max_attempts):
            try:
                shutil.copy(original, self.temp_file)
                break
            except Exception:
                pass

        return self.temp_file

    def test_default_datafile(self):
        """Test that the class initializes with the default datafile."""
        self.apex_out = apexpy.Apex()
        assert os.path.isfile(self.apex_out.datafile)
        return

    def test_custom_datafile(self):
        """Test that the class initializes with a good datafile input."""

        # Get the original datafile name
        apex_out_orig = apexpy.Apex()
        original_file = apex_out_orig.datafile
        del apex_out_orig

        # Create copy of datafile
        custom_file = self.copy_file(original_file)

        apex_out = apexpy.Apex(datafile=custom_file)
        assert apex_out.datafile == custom_file

        os.remove(custom_file)
        return

    def test_init_with_bad_datafile(self):
        """Test raises IOError with non-existent datafile input."""
        with pytest.raises(IOError) as oerr:
            apexpy.Apex(datafile=self.bad_file)
        assert str(oerr.value).startswith('Data file does not exist')
        return

    def test_default_fortranlib(self):
        """Test that the class initializes with the default fortranlib."""
        apex_out = apexpy.Apex()
        assert os.path.isfile(apex_out.fortranlib)
        return

    def test_custom_fortranlib(self):
        """Test that the class initializes with a good fortranlib input."""

        # Get the original fortranlib name
        apex_out_orig = apexpy.Apex()
        original_lib = apex_out_orig.fortranlib
        del apex_out_orig

        # Create copy of datafile
        custom_lib = self.copy_file(original_lib)

        apex_out = apexpy.Apex(fortranlib=custom_lib)
        assert apex_out.fortranlib == custom_lib

        os.remove(custom_lib)
        return

    def test_init_with_bad_fortranlib(self):
        """Test raises IOError with non-existent fortranlib input."""
        with pytest.raises(IOError) as oerr:
            apexpy.Apex(fortranlib=self.bad_file)
        assert str(oerr.value).startswith('Fortran library does not exist')
        return

    def test_igrf_fn(self):
        """Test the default igrf_fn."""
        apex_out = apexpy.Apex()
        assert os.path.isfile(apex_out.igrf_fn)
        return

    def test_repr_eval(self):
        """Test the Apex.__repr__ results."""
        # Initialize the apex object
        self.apex_out = apexpy.Apex()
        self.eval_date()
        self.eval_refh()

        # Get and test the repr string
        out_str = self.apex_out.__repr__()
        assert out_str.find("apexpy.Apex(") == 0

        # Test the ability to re-create the apex object from the repr string
        new_apex = eval(out_str)
        assert new_apex == self.apex_out
        return

    def test_ne_other_class(self):
        """Test Apex class inequality to a different class."""
        self.apex_out = apexpy.Apex()
        self.eval_date()
        self.eval_refh()

        assert self.apex_out != self.test_date
        return

    def test_ne_missing_attr(self):
        """Test Apex class inequality when attributes are missing from one."""
        self.apex_out = apexpy.Apex()
        self.eval_date()
        self.eval_refh()
        ref_apex = eval(self.apex_out.__repr__())
        del ref_apex.RE

        assert ref_apex != self.apex_out
        assert self.apex_out != ref_apex
        return

    def test_eq_missing_attr(self):
        """Test Apex class equality when attributes are missing from both."""
        self.apex_out = apexpy.Apex()
        self.eval_date()
        self.eval_refh()
        ref_apex = eval(self.apex_out.__repr__())
        del ref_apex.RE, self.apex_out.RE

        assert ref_apex == self.apex_out
        return

    def test_str_eval(self):
        """Test the Apex.__str__ results."""
        # Initialize the apex object
        self.apex_out = apexpy.Apex()
        self.eval_date()
        self.eval_refh()

        # Get and test the printed string
        out_str = self.apex_out.__str__()
        assert out_str.find("Decimal year") > 0
        return


class TestApexMethod(object):
    """Test the Apex methods."""
    def setup_method(self):
        """Initialize all tests."""
        self.apex_out = apexpy.Apex(date=2000, refh=300)
        self.in_lat = 60
        self.in_lon = 15
        self.in_alt = 100

    def teardown_method(self):
        """Clean up after each test."""
        del self.apex_out, self.in_lat, self.in_lon, self.in_alt

    def get_input_args(self, method_name, precision=0.0):
        """Set the input arguments for the different Apex methods.

        Parameters
        ----------
        method_name : str
            Name of the Apex class method
        precision : float
            Value for the precision (default=0.0)

        Returns
        -------
        in_args : list
            List of the appropriate input arguments

        """
        in_args = [self.in_lat, self.in_lon, self.in_alt]

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

    def test_apex_conversion_today(self):
        """Test Apex class conversion with today's date."""
        self.apex_out = apexpy.Apex(date=dt.datetime.now(tz=dt.timezone.utc),
                                    refh=300)
        assert not np.isnan(self.apex_out.geo2apex(self.in_lat, self.in_lon,
                                                   self.in_alt)).any()
        return

    @pytest.mark.parametrize("apex_method,fortran_method,fslice",
                             [("_geo2qd", "apxg2q", slice(0, 2, 1)),
                              ("_geo2apex", "apxg2all", slice(2, 4, 1)),
                              ("_qd2geo", "apxq2g", slice(None)),
                              ("_basevec", "apxg2q", slice(2, 4, 1))])
    @pytest.mark.parametrize("lat", [0, 30, 60, 89])
    @pytest.mark.parametrize("lon", [-179, -90, 0, 90, 180])
    def test_fortran_scalar_input(self, apex_method, fortran_method, fslice,
                                  lat, lon):
        """Tests Apex/fortran interface consistency for scalars.

        Parameters
        ----------
        apex_method : str
            Name of the Apex class method to test
        fortran_method : str
            Name of the Fortran function to test
        fslice : slice
            Slice used select the appropriate Fortran outputs
        lat : int or float
            Latitude in degrees N
        lon : int or float
            Longitude in degrees E

        """
        # Set the input coordinates
        self.in_lat = lat
        self.in_lon = lon

        # Get the Apex class method and the fortran function call
        apex_func = getattr(self.apex_out, apex_method)
        fortran_func = getattr(apexpy.fortranapex, fortran_method)

        # Get the appropriate input arguments
        apex_args = self.get_input_args(apex_method)
        fortran_args = self.get_input_args(fortran_method)

        # Evaluate the equivalent function calls
        np.testing.assert_allclose(apex_func(*apex_args),
                                   fortran_func(*fortran_args)[fslice])
        return

    @pytest.mark.parametrize("apex_method,fortran_method,fslice",
                             [("_geo2qd", "apxg2q", slice(0, 2, 1)),
                              ("_geo2apex", "apxg2all", slice(2, 4, 1)),
                              ("_qd2geo", "apxq2g", slice(None)),
                              ("_basevec", "apxg2q", slice(2, 4, 1))])
    @pytest.mark.parametrize("lat", [0, 30, 60, 89])
    @pytest.mark.parametrize("lon1,lon2", [(180, 180), (-180, -180),
                                           (180, -180), (-180, 180),
                                           (-345, 15), (375, 15)])
    def test_fortran_longitude_rollover(self, apex_method, fortran_method,
                                        fslice, lat, lon1, lon2):
        """Tests Apex/fortran interface consistency for longitude rollover.

        Parameters
        ----------
        apex_method : str
            Name of the Apex class method to test
        fortran_method : str
            Name of the Fortran function to test
        fslice : slice
            Slice used select the appropriate Fortran outputs
        lat : int or float
            Latitude in degrees N
        lon1 : int or float
            Longitude in degrees E
        lon2 : int or float
            Equivalent longitude in degrees E

        """
        # Set the fixed input coordinate
        self.in_lat = lat

        # Get the Apex class method and the fortran function call
        apex_func = getattr(self.apex_out, apex_method)
        fortran_func = getattr(apexpy.fortranapex, fortran_method)

        # Get the appropriate input arguments
        self.in_lon = lon1
        apex_args = self.get_input_args(apex_method)

        self.in_lon = lon2
        fortran_args = self.get_input_args(fortran_method)

        # Evaluate the equivalent function calls
        np.testing.assert_allclose(apex_func(*apex_args),
                                   fortran_func(*fortran_args)[fslice])
        return

    @pytest.mark.parametrize("arr_shape", [(2, 2), (4,), (1, 4)])
    @pytest.mark.parametrize("apex_method,fortran_method,fslice",
                             [("_geo2qd", "apxg2q", slice(0, 2, 1)),
                              ("_geo2apex", "apxg2all", slice(2, 4, 1)),
                              ("_qd2geo", "apxq2g", slice(None)),
                              ("_basevec", "apxg2q", slice(2, 4, 1))])
    def test_fortran_array_input(self, arr_shape, apex_method, fortran_method,
                                 fslice):
        """Tests Apex/fortran interface consistency for array input.

        Parameters
        ----------
        arr_shape : tuple
            Expected output shape
        apex_method : str
            Name of the Apex class method to test
        fortran_method : str
            Name of the Fortran function to test
        fslice : slice
            Slice used select the appropriate Fortran outputs

        """
        # Get the Apex class method and the fortran function call
        apex_func = getattr(self.apex_out, apex_method)
        fortran_func = getattr(apexpy.fortranapex, fortran_method)

        # Set up the input arrays
        ref_lat = np.array([0, 30, 60, 90])
        ref_alt = np.array([100, 200, 300, 400])
        self.in_lat = ref_lat.reshape(arr_shape)
        self.in_alt = ref_alt.reshape(arr_shape)
        apex_args = self.get_input_args(apex_method)

        # Get the Apex class results
        aret = apex_func(*apex_args)

        # Get the fortran function results
        flats = list()
        flons = list()

        for i, lat in enumerate(ref_lat):
            self.in_lat = lat
            self.in_alt = ref_alt[i]
            fortran_args = self.get_input_args(fortran_method)
            fret = fortran_func(*fortran_args)[fslice]
            flats.append(fret[0])
            flons.append(fret[1])

        flats = np.array(flats)
        flons = np.array(flons)

        # Evaluate results
        try:
            # This returned value is array of floats
            np.testing.assert_allclose(aret[0].astype(float),
                                       flats.reshape(arr_shape).astype(float))
            np.testing.assert_allclose(aret[1].astype(float),
                                       flons.reshape(arr_shape).astype(float))
        except ValueError:
            # This returned value is array of arrays
            alats = aret[0].reshape((4,))
            alons = aret[1].reshape((4,))
            for i, flat in enumerate(flats):
                np.testing.assert_array_almost_equal(alats[i], flat, 2)
                np.testing.assert_array_almost_equal(alons[i], flons[i], 2)

        return

    @pytest.mark.parametrize("lat", [0, 30, 60, 89])
    @pytest.mark.parametrize("lon", [-179, -90, 0, 90, 180])
    def test_geo2apexall_scalar(self, lat, lon):
        """Test Apex/fortran geo2apexall interface consistency for scalars.

        Parameters
        ----------
        lat : int or float
            Latitude in degrees N
        long : int or float
            Longitude in degrees E

        """
        # Get the Apex and Fortran results
        aret = self.apex_out._geo2apexall(lat, lon, self.in_alt)
        fret = apexpy.fortranapex.apxg2all(lat, lon, self.in_alt, 300, 1)

        # Evaluate each element in the results
        for aval, fval in zip(aret, fret):
            np.testing.assert_allclose(aval, fval)

    @pytest.mark.parametrize("arr_shape", [(2, 2), (4,), (1, 4)])
    def test_geo2apexall_array(self, arr_shape):
        """Test Apex/fortran geo2apexall interface consistency for arrays.

        Parameters
        ----------
        arr_shape : tuple
            Expected output shape

        """
        # Set the input
        self.in_lat = np.array([0, 30, 60, 90])
        self.in_alt = np.array([100, 200, 300, 400])

        # Get the Apex class results
        aret = self.apex_out._geo2apexall(self.in_lat.reshape(arr_shape),
                                          self.in_lon,
                                          self.in_alt.reshape(arr_shape))

        # For each lat/alt pair, get the Fortran results
        fret = list()
        for i, lat in enumerate(self.in_lat):
            fret.append(apexpy.fortranapex.apxg2all(lat, self.in_lon,
                                                    self.in_alt[i], 300, 1))

        # Cycle through all returned values
        for i, ret in enumerate(aret):
            try:
                # This returned value is array of floats
                fret_test = np.array([fret[0][i], fret[1][i], fret[2][i],
                                      fret[3][i]]).reshape(arr_shape)
                np.testing.assert_allclose(ret.astype(float),
                                           fret_test.astype(float))
            except ValueError:
                # This returned value is array of arrays
                ret = ret.reshape((4,))
                for j, single_fret in enumerate(fret):
                    np.testing.assert_allclose(ret[j], single_fret[i])
        return

    @pytest.mark.parametrize("in_coord", ["geo", "apex", "qd"])
    @pytest.mark.parametrize("out_coord", ["geo", "apex", "qd"])
    def test_convert_consistency(self, in_coord, out_coord):
        """Test the self-consistency of the Apex convert method.

        Parameters
        ----------
        in_coord : str
            Input coordinate system
        out_coord : str
            Output coordinate system

        """
        if in_coord == out_coord:
            pytest.skip("Test not needed for same src and dest coordinates")

        # Define the method name
        method_name = "2".join([in_coord, out_coord])

        # Get the method and method inputs
        convert_kwargs = {'height': self.in_alt, 'precision': 0.0}
        apex_args = self.get_input_args(method_name)
        apex_method = getattr(self.apex_out, method_name)

        # Define the slice needed to get equivalent output from the named method
        mslice = slice(0, -1, 1) if out_coord == "geo" else slice(None)

        # Get output using convert and named method
        convert_out = self.apex_out.convert(self.in_lat, self.in_lon, in_coord,
                                            out_coord, **convert_kwargs)
        method_out = apex_method(*apex_args)[mslice]

        # Compare both outputs, should be identical
        np.testing.assert_allclose(convert_out, method_out)
        return

    @pytest.mark.parametrize("bound_lat", [90, -90])
    @pytest.mark.parametrize("in_coord", ["geo", "apex", "qd"])
    @pytest.mark.parametrize("out_coord", ["geo", "apex", "qd"])
    def test_convert_at_lat_boundary(self, bound_lat, in_coord, out_coord):
        """Test the conversion at the latitude boundary, with allowed excess.

        Parameters
        ----------
        bound_lat : int or float
            Boundary latitude in degrees N
        in_coord : str
            Input coordinate system
        out_coord : str
            Output coordinate system

        """
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
        return

    @pytest.mark.parametrize("src", ["geo", "apex", "qd"])
    @pytest.mark.parametrize("dest", ["geo", "apex", "qd"])
    def test_convert_withnan(self, src, dest):
        """Test Apex.convert success with NaN input.

        Parameters
        ----------
        src : str
            Input coordinate system
        dest : str
            Output coordinate system

        """
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

    @pytest.mark.parametrize("bad_lat", [91, -91])
    def test_convert_invalid_lat(self, bad_lat):
        """Test convert raises ValueError for invalid latitudes.

        Parameters
        ----------
        bad_lat : int or float
            Latitude ouside the supported range in degrees N

        """

        with pytest.raises(ValueError) as verr:
            self.apex_out.convert(bad_lat, 0, 'geo', 'geo')

        assert str(verr.value).find("must be in [-90, 90]") > 0
        return

    @pytest.mark.parametrize("coords", [("foobar", "geo"), ("geo", "foobar"),
                                        ("geo", "mlt")])
    def test_convert_invalid_transformation(self, coords):
        """Test raises NotImplementedError for bad coordinates.

        Parameters
        ----------
        coords : tuple
            Tuple specifying the input and output coordinate systems

        """
        if "mlt" in coords:
            estr = "datetime must be given for MLT calculations"
        else:
            estr = "Unknown coordinate transformation"

        with pytest.raises(ValueError) as verr:
            self.apex_out.convert(0, 0, *coords)

        assert str(verr).find(estr) >= 0
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
        """Test the user method against set values with scalars.

        Parameters
        ----------
        method_name : str
            Apex class method to be tested
        out_comp : tuple of floats
            Expected output values

        """
        # Get the desired methods
        user_method = getattr(self.apex_out, method_name)

        # Get the user output
        user_out = user_method(self.in_lat, self.in_lon, self.in_alt)

        # Evaluate the user output
        np.testing.assert_allclose(user_out, out_comp, rtol=1e-5, atol=1e-5)

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
        """Test the user method with inputs that require some broadcasting.

        Parameters
        ----------
        in_coord : str
            Input coordiante system
        out_coord : str
            Output coordiante system
        method_args : list
            List of input arguments
        out_shape : tuple
            Expected shape of output values

        """
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
    @pytest.mark.parametrize("bad_lat", [91, -91])
    def test_method_invalid_lat(self, in_coord, out_coord, bad_lat):
        """Test convert raises ValueError for invalid latitudes.

        Parameters
        ----------
        in_coord : str
            Input coordiante system
        out_coord : str
            Output coordiante system
        bad_lat : int
            Latitude in degrees N that is out of bounds

        """
        if in_coord == out_coord:
            pytest.skip("Test not needed for same src and dest coordinates")

        # Get the desired methods
        method_name = "2".join([in_coord, out_coord])
        user_method = getattr(self.apex_out, method_name)

        with pytest.raises(ValueError) as verr:
            user_method(bad_lat, 15, 100)

        assert str(verr.value).find("must be in [-90, 90]") > 0
        return

    @pytest.mark.parametrize("in_coord", ["geo", "apex", "qd"])
    @pytest.mark.parametrize("out_coord", ["geo", "apex", "qd"])
    @pytest.mark.parametrize("bound_lat", [90, -90])
    def test_method_at_lat_boundary(self, in_coord, out_coord, bound_lat):
        """Test user methods at the latitude boundary, with allowed excess.

        Parameters
        ----------
        in_coord : str
            Input coordiante system
        out_coord : str
            Output coordiante system
        bad_lat : int
            Latitude in degrees N that is at the limits of the boundary

        """
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
        self.apex_out = apexpy.Apex(date=2000, refh=10000)

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
    @pytest.mark.parametrize("delta_h", [1.0e-6, -1.0e-6])
    def test_quasidipole_apexheight_close(self, method_name, delta_h):
        """Test quasi-dipole success with a height close to the reference.

        Parameters
        ----------
        method_name : str
            Apex class method name to be tested
        delta_h : float
            tolerance for height in km

        """
        qd_method = getattr(self.apex_out, method_name)
        in_args = [0, 15, self.apex_out.refh + delta_h]
        out_coords = qd_method(*in_args)

        for i, out_val in enumerate(out_coords):
            np.testing.assert_almost_equal(out_val, in_args[i], decimal=3)
        return

    @pytest.mark.parametrize("method_name, hinc, msg",
                             [("apex2qd", 1.0, "is > apex height"),
                              ("qd2apex", -1.0, "is < reference height")])
    def test_quasidipole_raises_apexheight(self, method_name, hinc, msg):
        """Quasi-dipole raises ApexHeightError when height above reference.

        Parameters
        ----------
        method_name : str
            Apex class method name to be tested
        hinc : float
            Height increment in km
        msg : str
            Expected output message

        """
        qd_method = getattr(self.apex_out, method_name)

        with pytest.raises(apexpy.ApexHeightError) as aerr:
            qd_method(0, 15, self.apex_out.refh + hinc)

        assert str(aerr).find(msg) > 0
        return


class TestApexMLTMethods(object):
    """Test the Apex Magnetic Local Time (MLT) methods."""
    def setup_method(self):
        """Initialize all tests."""
        self.apex_out = apexpy.Apex(date=2000, refh=300)
        self.in_time = dt.datetime(2000, 2, 3, 4, 5, 6)

    def teardown_method(self):
        """Clean up after each test."""
        del self.apex_out, self.in_time

    @pytest.mark.parametrize("in_coord", ["geo", "apex", "qd"])
    def test_convert_to_mlt(self, in_coord):
        """Test the conversions to MLT using Apex convert.

        Parameters
        ----------
        in_coord : str
            Input coordinate system

        """

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
        """Test the conversions from MLT using Apex convert.

        Parameters
        ----------
        out_coord : str
            Output coordinate system

        """
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
        """Test mlon2mlt with scalar inputs.

        Parameters
        ----------
        mlon_kwargs : dict
            Input kwargs
        test_mlt : float
            Output MLT in hours

        """
        mlt = self.apex_out.mlon2mlt(0, self.in_time, **mlon_kwargs)

        np.testing.assert_allclose(mlt, test_mlt)
        assert np.asarray(mlt).shape == ()
        return

    @pytest.mark.parametrize("mlt_kwargs,test_mlon",
                             [({}, 14.705535888671875),
                              ({"ssheight": 100000}, 14.599319458007812)])
    def test_mlt2mlon_scalar_inputs(self, mlt_kwargs, test_mlon):
        """Test mlt2mlon with scalar inputs.

        Parameters
        ----------
        mlt_kwargs : dict
            Input kwargs
        test_mlon : float
            Output longitude in degrees E

        """
        mlon = self.apex_out.mlt2mlon(0, self.in_time, **mlt_kwargs)

        np.testing.assert_allclose(mlon, test_mlon)
        assert np.asarray(mlon).shape == ()
        return

    @pytest.mark.parametrize("mlon,test_mlt",
                             [([0, 180], [23.019261, 11.019261]),
                              (np.array([0, 180]), [23.019261, 11.019261]),
                              (np.array([[0], [180]]),
                               np.array([[23.019261], [11.019261]])),
                              ([[0, 180], [0, 180]], [[23.019261, 11.019261],
                                                      [23.019261, 11.019261]]),
                              (range(0, 361, 30),
                               [23.01963, 1.01963, 3.01963, 5.01963, 7.01963,
                                9.01963, 11.01963, 13.01963, 15.01963, 17.01963,
                                19.01963, 21.01963, 23.01963])])
    def test_mlon2mlt_array(self, mlon, test_mlt):
        """Test mlon2mlt with array inputs.

        Parameters
        ----------
        mlon : array-like
            Input longitudes in degrees E
        test_mlt : float
            Output MLT in hours

        """
        mlt = self.apex_out.mlon2mlt(mlon, self.in_time)

        assert mlt.shape == np.asarray(test_mlt).shape
        np.testing.assert_allclose(mlt, test_mlt, rtol=1e-4)
        return

    @pytest.mark.parametrize("mlt,test_mlon",
                             [([0, 12], [14.705551, 194.705551]),
                              (np.array([0, 12]), [14.705551, 194.705551]),
                              (np.array([[0], [12]]),
                               np.array([[14.705551], [194.705551]])),
                              ([[0, 12], [0, 12]], [[14.705551, 194.705551],
                                                    [14.705551, 194.705551]]),
                              (range(0, 25, 2),
                               [14.705551, 44.705551, 74.705551, 104.705551,
                                134.705551, 164.705551, 194.705551, 224.705551,
                                254.705551, 284.705551, 314.705551, 344.705551,
                                14.705551])])
    def test_mlt2mlon_array(self, mlt, test_mlon):
        """Test mlt2mlon with array inputs.

        Parameters
        ----------
        mlt : array-like
            Input MLT in hours
        test_mlon : float
            Output longitude in degrees E

        """
        mlon = self.apex_out.mlt2mlon(mlt, self.in_time)

        assert mlon.shape == np.asarray(test_mlon).shape
        np.testing.assert_allclose(mlon, test_mlon, rtol=1e-4)
        return

    @pytest.mark.parametrize("method_name", ["mlon2mlt", "mlt2mlon"])
    def test_mlon2mlt_diffdates(self, method_name):
        """Test that MLT varies with universal time.

        Parameters
        ----------
        method_name : str
            Name of Apex class method to be tested

        """
        apex_method = getattr(self.apex_out, method_name)
        mlt1 = apex_method(0, self.in_time)
        mlt2 = apex_method(0, self.in_time + dt.timedelta(hours=1))

        assert mlt1 != mlt2
        return

    @pytest.mark.parametrize("mlt_offset", [1.0, 10.0])
    def test_mlon2mlt_offset(self, mlt_offset):
        """Test the time wrapping logic for the MLT.

        Parameters
        ----------
        mlt_offset : float
            MLT offset in hours

        """
        mlt1 = self.apex_out.mlon2mlt(0.0, self.in_time)
        mlt2 = self.apex_out.mlon2mlt(-15.0 * mlt_offset,
                                      self.in_time) + mlt_offset

        np.testing.assert_allclose(mlt1, mlt2)
        return

    @pytest.mark.parametrize("mlon_offset", [15.0, 150.0])
    def test_mlt2mlon_offset(self, mlon_offset):
        """Test the time wrapping logic for the magnetic longitude.

        Parameters
        ----------
        mlt_offset : float
            MLT offset in hours

        """
        mlon1 = self.apex_out.mlt2mlon(0, self.in_time)
        mlon2 = self.apex_out.mlt2mlon(mlon_offset / 15.0,
                                       self.in_time) - mlon_offset

        np.testing.assert_allclose(mlon1, mlon2)
        return

    @pytest.mark.parametrize("order", [["mlt", "mlon"], ["mlon", "mlt"]])
    @pytest.mark.parametrize("start_val", [0, 6, 12, 18, 22])
    def test_convert_and_return(self, order, start_val):
        """Test the conversion to magnetic longitude or MLT and back again.

        Parameters
        ----------
        order : list
            List of strings specifying the order to run functions
        start_val : int or float
            Input value

        """
        first_method = getattr(self.apex_out, "2".join(order))
        second_method = getattr(self.apex_out, "2".join([order[1], order[0]]))

        middle_val = first_method(start_val, self.in_time)
        end_val = second_method(middle_val, self.in_time)

        np.testing.assert_allclose(start_val, end_val)
        return


class TestApexMapMethods(object):
    """Test the Apex height mapping methods."""
    def setup_method(self):
        """Initialize all tests."""
        self.apex_out = apexpy.Apex(date=2000, refh=300)

    def teardown_method(self):
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
        """Test the map_to_height function.

        Parameters
        ----------
        in_args : list
            List of input arguments
        test_mapped : list
            List of expected outputs

        """
        mapped = self.apex_out.map_to_height(*in_args)
        np.testing.assert_allclose(mapped, test_mapped, rtol=1e-5, atol=1e-5)
        return

    def test_map_to_height_same_height(self):
        """Test the map_to_height function when mapping to same height."""
        mapped = self.apex_out.map_to_height(60, 15, 100, 100, conjugate=False,
                                             precision=1e-10)
        np.testing.assert_allclose(mapped, (60.0, 15.000003814697266, 0.0),
                                   rtol=1e-5, atol=1e-5)
        return

    @pytest.mark.parametrize('arr_shape', [(2,), (2, 2), (1, 4)])
    @pytest.mark.parametrize('ivec', range(0, 4))
    def test_map_to_height_array_location(self, arr_shape, ivec):
        """Test map_to_height with array input.

        Parameters
        ----------
        arr_shape : tuple
            Expected array shape
        ivec : int
            Input argument index for vectorized input

        """
        # Set the base input and output values
        in_args = [60, 15, 100, 100]
        test_mapped = [60, 15.00000381, 0.0]

        # Update inputs for one vectorized value
        in_args[ivec] = np.full(shape=arr_shape, fill_value=in_args[ivec])

        # Calculate and test function
        mapped = self.apex_out.map_to_height(*in_args)
        for i, test_val in enumerate(test_mapped):
            assert mapped[i].shape == arr_shape
            np.testing.assert_allclose(mapped[i], test_val, rtol=1e-5,
                                       atol=1e-5)
        return

    @pytest.mark.parametrize("method_name,in_args",
                             [("map_to_height", [0, 15, 100, 10000]),
                              ("map_E_to_height",
                               [0, 15, 100, 10000, [1, 2, 3]]),
                              ("map_V_to_height",
                               [0, 15, 100, 10000, [1, 2, 3]])])
    def test_mapping_height_raises_ApexHeightError(self, method_name, in_args):
        """Test map_to_height raises ApexHeightError.

        Parameters
        ----------
        method_name : str
            Name of the Apex class method to test
        in_args : list
            List of input arguments

        """
        apex_method = getattr(self.apex_out, method_name)

        with pytest.raises(apexpy.ApexHeightError) as aerr:
            apex_method(*in_args)

        assert aerr.match("is > apex height")
        return

    @pytest.mark.parametrize("method_name",
                             ["map_E_to_height", "map_V_to_height"])
    @pytest.mark.parametrize("ev_input", [([1, 2, 3, 4, 5]),
                                          ([[1, 2], [3, 4], [5, 6], [7, 8]])])
    def test_mapping_EV_bad_shape(self, method_name, ev_input):
        """Test height mapping of E/V with baddly shaped input raises Error.

        Parameters
        ----------
        method_name : str
            Name of the Apex class method to test
        ev_input : list
            E/V input arguments

        """
        apex_method = getattr(self.apex_out, method_name)
        in_args = [60, 15, 100, 500, ev_input]
        with pytest.raises(ValueError) as verr:
            apex_method(*in_args)

        assert str(verr.value).find("must be (3, N) or (3,) ndarray") >= 0
        return

    def test_mapping_EV_bad_flag(self):
        """Test _map_EV_to_height raises error for bad data type flag."""
        with pytest.raises(ValueError) as verr:
            self.apex_out._map_EV_to_height(60, 15, 100, 500, [1, 2, 3], "P")

        assert str(verr.value).find("unknown electric field/drift flag") >= 0
        return

    @pytest.mark.parametrize("in_args,test_mapped",
                             [([60, 15, 100, 500, [1, 2, 3]],
                               [0.71152183, 2.35624876, 0.57260784]),
                              ([60, 15, 100, 500, [2, 3, 4]],
                               [1.56028502, 3.43916636, 0.78235384]),
                              ([60, 15, 100, 1000, [1, 2, 3]],
                               [0.67796492, 2.08982134, 0.55860785]),
                              ([60, 15, 200, 500, [1, 2, 3]],
                               [0.72377397, 2.42737471, 0.59083726]),
                              ([60, 30, 100, 500, [1, 2, 3]],
                               [0.68626344, 2.37530133, 0.60060124]),
                              ([70, 15, 100, 500, [1, 2, 3]],
                               [0.72760378, 2.18082305, 0.29141979])])
    def test_map_E_to_height_scalar_location(self, in_args, test_mapped):
        """Test mapping of E-field to a specified height.

        Parameters
        ----------
        in_args : list
            List of input arguments
        test_mapped : list
            List of expected outputs

        """
        mapped = self.apex_out.map_E_to_height(*in_args)
        np.testing.assert_allclose(mapped, test_mapped, rtol=1e-5)
        return

    @pytest.mark.parametrize('ev_flag, test_mapped',
                             [('E', [0.71152183, 2.35624876, 0.57260784]),
                              ('V', [0.81971957, 2.84512495, 0.69545001])])
    @pytest.mark.parametrize('arr_shape', [(2,), (5,)])
    @pytest.mark.parametrize('ivec', range(0, 5))
    def test_map_EV_to_height_array_location(self, ev_flag, test_mapped,
                                             arr_shape, ivec):
        """Test mapping of E-field/drift to a specified height with arrays.

        Parameters
        ----------
        ev_flag : str
            Character flag specifying whether to run 'E' or 'V' methods
        test_mapped : list
            List of expected outputs
        arr_shape : tuple
            Shape of the expected output
        ivec : int
            Index of the expected output

        """
        # Set the base input and output values
        eshape = list(arr_shape)
        eshape.insert(0, 3)
        edata = np.array([[1, 2, 3]] * np.prod(arr_shape)).transpose()
        in_args = [60, 15, 100, 500, edata.reshape(tuple(eshape))]

        # Update inputs for one vectorized value if this is a location input
        if ivec < 4:
            in_args[ivec] = np.full(shape=arr_shape, fill_value=in_args[ivec])

        # Get the mapped output
        apex_method = getattr(self.apex_out,
                              "map_{:s}_to_height".format(ev_flag))
        mapped = apex_method(*in_args)

        # Test the results
        for i, test_val in enumerate(test_mapped):
            assert mapped[i].shape == arr_shape
            np.testing.assert_allclose(mapped[i], test_val, rtol=1e-5)
        return

    @pytest.mark.parametrize("in_args,test_mapped",
                             [([60, 15, 100, 500, [1, 2, 3]],
                               [0.81971957, 2.84512495, 0.69545001]),
                              ([60, 15, 100, 500, [2, 3, 4]],
                               [1.83027746, 4.14346436, 0.94764179]),
                              ([60, 15, 100, 1000, [1, 2, 3]],
                               [0.92457698, 3.14997661, 0.85135187]),
                              ([60, 15, 200, 500, [1, 2, 3]],
                               [0.80388262, 2.79321504, 0.68285158]),
                              ([60, 30, 100, 500, [1, 2, 3]],
                               [0.76141245, 2.87884673, 0.73655941]),
                              ([70, 15, 100, 500, [1, 2, 3]],
                               [0.84681866, 2.5925821, 0.34792655])])
    def test_map_V_to_height_scalar_location(self, in_args, test_mapped):
        """Test mapping of velocity to a specified height.

        Parameters
        ----------
        in_args : list
            List of input arguments
        test_mapped : list
            List of expected outputs

        """
        mapped = self.apex_out.map_V_to_height(*in_args)
        np.testing.assert_allclose(mapped, test_mapped, rtol=1e-5)
        return


class TestApexBasevectorMethods(object):
    """Test the Apex height base vector methods."""
    def setup_method(self):
        """Initialize all tests."""
        self.apex_out = apexpy.Apex(date=2000, refh=300)
        self.lat = 60
        self.lon = 15
        self.height = 100
        self.test_basevec = None

    def teardown_method(self):
        """Clean up after each test."""
        del self.apex_out, self.test_basevec, self.lat, self.lon, self.height

    def get_comparison_results(self, bv_coord, coords, precision):
        """Get the base vector results using the hidden function for comparison.

        Parameters
        ----------
        bv_coord : str
            Basevector coordinate scheme, expects on of 'apex', 'qd',
            or 'bvectors_apex'
        coords : str
            Expects one of 'geo', 'apex', or 'qd'
        precision : float
            Float specifiying precision

        """
        if coords == "geo":
            glat = self.lat
            glon = self.lon
        else:
            apex_method = getattr(self.apex_out, "{:s}2geo".format(coords))
            glat, glon, _ = apex_method(self.lat, self.lon, self.height,
                                        precision=precision)

        if bv_coord == 'qd':
            self.test_basevec = self.apex_out._basevec(glat, glon, self.height)
        elif bv_coord == 'apex':
            (_, _, _, _, f1, f2, _, d1, d2, d3, _, e1, e2,
             e3) = self.apex_out._geo2apexall(glat, glon, 100)
            self.test_basevec = (f1, f2, d1, d2, d3, e1, e2, e3)
        else:
            # These are set results that need to be updated with IGRF
            if coords == "geo":
                self.test_basevec = (
                    np.array([4.42368795e-05, 4.42368795e-05]),
                    np.array([[0.01047826, 0.01047826],
                              [0.33089194, 0.33089194],
                              [-1.04941, -1.04941]]),
                    np.array([5.3564698e-05, 5.3564698e-05]),
                    np.array([[0.00865356, 0.00865356],
                              [0.27327004, 0.27327004],
                              [-0.8666646, -0.8666646]]))
            elif coords == "apex":
                self.test_basevec = (
                    np.array([4.48672735e-05, 4.48672735e-05]),
                    np.array([[-0.12510721, -0.12510721],
                              [0.28945938, 0.28945938],
                              [-1.1505738, -1.1505738]]),
                    np.array([6.38577444e-05, 6.38577444e-05]),
                    np.array([[-0.08790194, -0.08790194],
                              [0.2033779, 0.2033779],
                              [-0.808408, -0.808408]]))
            else:
                self.test_basevec = (
                    np.array([4.46348578e-05, 4.46348578e-05]),
                    np.array([[-0.12642345, -0.12642345],
                              [0.29695055, 0.29695055],
                              [-1.1517885, -1.1517885]]),
                    np.array([6.38626285e-05, 6.38626285e-05]),
                    np.array([[-0.08835986, -0.08835986],
                              [0.20754464, 0.20754464],
                              [-0.8050078, -0.8050078]]))

        return

    @pytest.mark.parametrize("bv_coord", ["qd", "apex"])
    @pytest.mark.parametrize("coords,precision",
                             [("geo", 1e-10), ("apex", 1.0e-2), ("qd", 1.0e-2)])
    def test_basevectors_scalar(self, bv_coord, coords, precision):
        """Test the base vector calculations with scalars.

        Parameters
        ----------
        bv_coord : str
            Name of the input coordinate system
        coords : str
            Name of the output coordinate system
        precision : float
            Level of run precision requested

        """
        # Get the base vectors
        base_method = getattr(self.apex_out,
                              "basevectors_{:s}".format(bv_coord))
        basevec = base_method(self.lat, self.lon, self.height, coords=coords,
                              precision=precision)
        self.get_comparison_results(bv_coord, coords, precision)
        if bv_coord == "apex":
            basevec = list(basevec)
            for i in range(4):
                # Not able to compare indices 2, 3, 4, and 5
                basevec.pop(2)

        # Test the results
        for i, vec in enumerate(basevec):
            np.testing.assert_allclose(vec, self.test_basevec[i])
        return

    @pytest.mark.parametrize("bv_coord", ["qd", "apex"])
    def test_basevectors_scalar_shape(self, bv_coord):
        """Test the shape of the scalar output.

        Parameters
        ----------
        bv_coord : str
            Name of the input coordinate system

        """
        base_method = getattr(self.apex_out,
                              "basevectors_{:s}".format(bv_coord))
        basevec = base_method(self.lat, self.lon, self.height)

        for i, vec in enumerate(basevec):
            if i < 2:
                assert vec.shape == (2,)
            else:
                assert vec.shape == (3,)
        return

    @pytest.mark.parametrize('arr_shape', [(2,), (5,)])
    @pytest.mark.parametrize("bv_coord", ["qd", "apex"])
    @pytest.mark.parametrize("ivec", range(3))
    def test_basevectors_array(self, arr_shape, bv_coord, ivec):
        """Test the output shape for array inputs.

        Parameters
        ----------
        arr_shape : tuple
            Expected output shape
        bv_coord : str
            Name of the input coordinate system
        ivec : int
            Index of the evaluated output value

        """
        # Define the input arguments
        in_args = [self.lat, self.lon, self.height]
        in_args[ivec] = np.full(shape=arr_shape, fill_value=in_args[ivec])

        # Get the basevectors
        base_method = getattr(self.apex_out,
                              "basevectors_{:s}".format(bv_coord))
        basevec = base_method(*in_args, coords='geo', precision=1e-10)
        self.get_comparison_results(bv_coord, "geo", 1e-10)
        if bv_coord == "apex":
            basevec = list(basevec)
            for i in range(4):
                # Not able to compare indices 2, 3, 4, and 5
                basevec.pop(2)

        # Evaluate the shape and the values
        for i, vec in enumerate(basevec):
            test_shape = list(arr_shape)
            test_shape.insert(0, 2 if i < 2 else 3)
            assert vec.shape == tuple(test_shape)
            assert np.all(self.test_basevec[i][0] == vec[0])
            assert np.all(self.test_basevec[i][1] == vec[1])
        return

    @pytest.mark.parametrize("coords", ["geo", "apex", "qd"])
    def test_bvectors_apex(self, coords):
        """Test the bvectors_apex method.

        Parameters
        ----------
        coords : str
            Name of the coordiante system

        """
        in_args = [[self.lat, self.lat], [self.lon, self.lon],
                   [self.height, self.height]]
        self.get_comparison_results("bvectors_apex", coords, 1e-10)

        basevec = self.apex_out.bvectors_apex(*in_args, coords=coords,
                                              precision=1e-10)
        for i, vec in enumerate(basevec):
            np.testing.assert_array_almost_equal(vec, self.test_basevec[i],
                                                 decimal=5)
        return

    def test_basevectors_apex_extra_values(self):
        """Test specific values in the apex base vector output."""
        # Set the testing arrays
        self.test_basevec = [np.array([0.092637, -0.245951, 0.938848]),
                             np.array([0.939012, 0.073416, -0.07342]),
                             np.array([0.055389, 1.004155, 0.257594]),
                             np.array([0, 0, 1.065135])]

        # Get the desired output
        basevec = self.apex_out.basevectors_apex(0, 15, 100, coords='geo')

        # Test the values not covered by `test_basevectors_scalar`
        for itest, ibase in enumerate(np.arange(2, 6, 1)):
            np.testing.assert_allclose(basevec[ibase],
                                       self.test_basevec[itest], rtol=1e-4)
        return

    @pytest.mark.parametrize("lat", range(0, 90, 10))
    @pytest.mark.parametrize("lon", range(0, 360, 15))
    def test_basevectors_apex_delta(self, lat, lon):
        """Test that vectors are calculated correctly.

        Parameters
        ----------
        lat : int or float
            Latitude in degrees N
        lon : int or float
            Longitude in degrees E

        """
        # Get the apex base vectors and sort them for easy testing
        (f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2,
         e3) = self.apex_out.basevectors_apex(lat, lon, 500)
        fvec = [np.append(f1, 0), np.append(f2, 0), f3]
        gvec = [g1, g2, g3]
        dvec = [d1, d2, d3]
        evec = [e1, e2, e3]

        for idelta, jdelta in [(i, j) for i in range(3) for j in range(3)]:
            delta = 1 if idelta == jdelta else 0
            np.testing.assert_allclose(np.sum(fvec[idelta] * gvec[jdelta]),
                                       delta, rtol=0, atol=1e-5)
            np.testing.assert_allclose(np.sum(dvec[idelta] * evec[jdelta]),
                                       delta, rtol=0, atol=1e-5)
        return

    def test_basevectors_apex_invalid_scalar(self):
        """Test warning and fill values for base vectors with bad inputs."""
        self.apex_out = apexpy.Apex(date=2000, refh=10000)
        invalid = np.full(shape=(3,), fill_value=np.nan)

        # Get the output and the warnings
        with warnings.catch_warnings(record=True) as warn_rec:
            basevec = self.apex_out.basevectors_apex(0, 0, 0)

        for i, bvec in enumerate(basevec):
            if i < 2:
                assert not np.allclose(bvec, invalid[:2])
            else:
                np.testing.assert_allclose(bvec, invalid)

        assert issubclass(warn_rec[-1].category, UserWarning)
        assert 'set to NaN where' in str(warn_rec[-1].message)
        return


class TestApexGetMethods(object):
    """Test the Apex `get` methods."""
    def setup_method(self):
        """Initialize all tests."""
        self.apex_out = apexpy.Apex(date=2000, refh=300)

    def teardown_method(self):
        """Clean up after each test."""
        del self.apex_out

    @pytest.mark.parametrize("alat, aheight",
                             [(10, 507.409702543805),
                              (60, 20313.026999999987),
                              ([10, 60],
                               [507.409702543805, 20313.026999999987]),
                              ([[10], [60]],
                               [[507.409702543805], [20313.026999999987]])])
    def test_get_apex(self, alat, aheight):
        """Test the apex height retrieval results.

        Parameters
        ----------
        alat : int or float
            Apex latitude in degrees N
        aheight : int or float
            Apex height in km

        """
        alt = self.apex_out.get_apex(alat)
        np.testing.assert_allclose(alt, aheight)
        return

    @pytest.mark.parametrize("glat,glon,height,test_bmag",
                             [([80], [100], [300], 5.100682377815247e-05),
                              ([80, 80], [100], [300],
                               [5.100682377815247e-05, 5.100682377815247e-05]),
                              ([[80], [80]], [100], [300],
                               [[5.100682377815247e-05],
                                [5.100682377815247e-05]]),
                              (range(50, 90, 8), range(0, 360, 80), [300] * 5,
                               np.array([4.18657154e-05, 5.11118114e-05,
                                         4.91969854e-05, 5.10519207e-05,
                                         4.90054816e-05])),
                              (90.0, 0, 1000, 3.7834718823432923e-05)])
    def test_get_babs(self, glat, glon, height, test_bmag):
        """Test the method to get the magnitude of the magnetic field.

        Parameters
        ----------
        glat : list
            List of latitudes in degrees N
        glon : list
            List of longitudes in degrees E
        height : list
            List of heights in km
        test_bmag : float
            Expected B field magnitude

        """
        bmag = self.apex_out.get_babs(glat, glon, height)
        np.testing.assert_allclose(bmag, test_bmag, rtol=0, atol=1e-5)
        return

    @pytest.mark.parametrize("bad_lat", [(91), (-91)])
    def test_get_apex_with_invalid_lat(self, bad_lat):
        """Test get methods raise ValueError for invalid latitudes.

        Parameters
        ----------
        bad_lat : int or float
            Bad input latitude in degrees N

        """

        with pytest.raises(ValueError) as verr:
            self.apex_out.get_apex(bad_lat)

        assert str(verr.value).find("must be in [-90, 90]") > 0
        return

    @pytest.mark.parametrize("bad_lat", [(91), (-91)])
    def test_get_babs_with_invalid_lat(self, bad_lat):
        """Test get methods raise ValueError for invalid latitudes.

        Parameters
        ----------
        bad_lat : int or float
            Bad input latitude in degrees N

        """

        with pytest.raises(ValueError) as verr:
            self.apex_out.get_babs(bad_lat, 15, 100)

        assert str(verr.value).find("must be in [-90, 90]") > 0
        return

    @pytest.mark.parametrize("bound_lat", [(90), (-90)])
    def test_get_at_lat_boundary(self, bound_lat):
        """Test get methods at the latitude boundary, with allowed excess.

        Parameters
        ----------
        bound_lat : int or float
            Boundary input latitude in degrees N

        """
        # Get a latitude just beyond the limit
        excess_lat = np.sign(bound_lat) * (abs(bound_lat) + 1.0e-5)

        # Get the two outputs, slight tolerance outside of boundary allowed
        bound_out = self.apex_out.get_apex(bound_lat)
        excess_out = self.apex_out.get_apex(excess_lat)

        # Test the outputs
        np.testing.assert_allclose(excess_out, bound_out, rtol=0, atol=1e-8)
        return

    @pytest.mark.parametrize("apex_height", [-100, 0, 300, 10000])
    def test_get_height_at_equator(self, apex_height):
        """Test that `get_height` returns apex height at equator.

        Parameters
        ----------
        apex_height : float
            Apex height

        """

        assert apex_height == self.apex_out.get_height(0.0, apex_height)
        return

    @pytest.mark.parametrize("lat, height", [
        (-90, -6371.009), (-80, -6088.438503309167), (-70, -5274.8091854339655),
        (-60, -4028.256749999999), (-50, -2499.1338178752017),
        (-40, -871.8751821247979), (-30, 657.2477500000014),
        (-20, 1903.8001854339655), (-10, 2717.4295033091657), (0, 3000.0),
        (10, 2717.4295033091657), (20, 1903.8001854339655),
        (30, 657.2477500000014), (40, -871.8751821247979),
        (50, -2499.1338178752017), (60, -4028.256749999999),
        (70, -5274.8091854339655), (80, -6088.438503309167)])
    def test_get_height_along_fieldline(self, lat, height):
        """Test that `get_height` returns expected height of field line.

        Parameters
        ----------
        lat : float
            Input latitude
        height : float
            Output field-line height for line with apex of 3000 km

        """

        fheight = self.apex_out.get_height(lat, 3000.0)
        assert abs(height - fheight) < 1.0e-7, \
            "bad height calculation: {:.7f} != {:.7f}".format(height, fheight)
        return


class TestApexMethodExtrapolateIGRF(object):
    """Test the Apex methods on a year when IGRF must be extrapolated.

    Notes
    -----
    Extrapolation should be done using a year within 5 years of the latest IGRF
    model epoch.

    """

    def setup_method(self):
        """Initialize all tests."""
        self.apex_out = apexpy.Apex(date=2030, refh=300)
        self.in_lat = 60
        self.in_lon = 15
        self.in_alt = 100
        self.in_time = dt.datetime(2029, 2, 3, 4, 5, 6)
        return

    def teardown_method(self):
        """Clean up after each test."""
        del self.apex_out, self.in_lat, self.in_lon, self.in_alt
        return

    @pytest.mark.parametrize("method_name, out_comp",
                             [("geo2apex",
                               (56.30344009399414, 91.99834442138672)),
                              ("apex2geo",
                               (54.22148513793945, -67.09271240234375,
                                4.268868451617891e-06)),
                              ("geo2qd",
                               (56.87860870361328, 91.99834442138672)),
                              ("apex2qd", (60.498401178276744, 15.0)),
                              ("qd2apex", (59.49138097045895, 15.0))])
    def test_method_scalar_input(self, method_name, out_comp):
        """Test the user method against set values with scalars.

        Parameters
        ----------
        method_name : str
            Apex class method to be tested
        out_comp : tuple of floats
            Expected output values

        """
        # Get the desired methods
        user_method = getattr(self.apex_out, method_name)

        # Get the user output
        user_out = user_method(self.in_lat, self.in_lon, self.in_alt)

        # Evaluate the user output
        np.testing.assert_allclose(user_out, out_comp, rtol=1e-5, atol=1e-5)

        for out_val in user_out:
            assert np.asarray(out_val).shape == (), "output is not a scalar"
        return

    def test_convert_to_mlt(self):
        """Test conversion from mlon to mlt with scalars."""

        # Get user output
        user_out = self.apex_out.mlon2mlt(self.in_lon, self.in_time)

        # Set comparison values
        out_comp = 23.933631388346353

        # Evaluate user output
        np.testing.assert_allclose(user_out, out_comp, rtol=1e-5, atol=1e-5)
        return
