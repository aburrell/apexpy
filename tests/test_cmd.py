# -*- coding: utf-8 -*-
"""Unit tests for command line execution."""

import numpy as np
import os
import pytest
import subprocess


class TestCommandLine(object):
    """Test class for the command-line apexpy interface."""

    def setup_method(self):
        """Runs before every test method to create a clean environment."""
        # Define the desired working paths
        self.startdir = os.path.abspath(os.path.curdir)
        split_dirs = os.path.split(os.path.dirname(os.path.abspath(__file__)))
        self.workdir = split_dirs[0]

        # Change directory, if needed
        if self.startdir != self.workdir:
            os.chdir(self.workdir)

        # Define the test filenames
        self.outfile = os.path.join(split_dirs[1], 'output.txt')
        self.infile = os.path.join(split_dirs[1], 'test_convert.txt')
        self.singlefile = os.path.join(split_dirs[1],
                                       'test_convert_single_line.txt')

    def teardown_method(self):
        """Runs after every method to clean up previous testing."""
        # Remove testing output
        if os.path.isfile(self.outfile):
            os.remove(self.outfile)

        # Return to starting directory
        if self.startdir != os.path.abspath(os.path.curdir):
            os.chdir(self.startdir)

        del self.outfile, self.infile, self.singlefile

    def execute_command_line(self, command, command_kwargs=None,
                             pipe_out=False):
        """Execute the command and load data from self.outfile

        Parameters
        ----------
        command : list or str
            List or string containing command to execute using subprocess
        command_kwargs : dict or NoneType
            Dict containing optional kwargs for subprocess.Popen command or
            None if using defaults. (default=None)
        pipe_out : bool
            Return pipe output instead of output from a data file if True
            (default=False)

        Returns
        -------
        data : np.array, NoneType, or subprocess.Popen attribute
            Numpy array of data from output file, None if no file was created,
            or the requested output from the pipe command.

        """
        data = None

        if pipe_out:
            if command_kwargs is None:
                command_kwargs = {}

            pipe = subprocess.Popen(command, **command_kwargs)
            data = pipe.communicate()
            pipe.wait()
        else:
            os.system(" ".join(command))
            if os.path.isfile(self.outfile):
                data = np.loadtxt(self.outfile)

        return data

    @pytest.mark.parametrize("date_str", [("2015"), ("201501"), ('20150101'),
                                          ('20150101000000')])
    def test_convert_w_datetime(self, date_str):
        """Test command line with different date and time specification.

        Parameters
        ----------
        date_str : str
           Input date string

        """
        # Build and execute the apexpy command line call
        cmd = ['python', '-m', 'apexpy', 'geo', 'apex', date_str, '--height',
               '300', '-i', self.infile, '-o', self.outfile]
        data = self.execute_command_line(cmd)

        # Test the outfile existance and values
        assert data is not None, 'error executing: {:s}'.format(' '.join(cmd))
        np.testing.assert_allclose(data, [[57.47145462, 93.62657928],
                                          [58.52458191, 94.03150177],
                                          [59.57331467, 94.46398163]],
                                   rtol=1e-4)
        return

    def test_convert_single_line(self):
        """Test command line with a single line of output."""
        # Build and execute the apexpy command line call
        cmd = ['python', '-m', 'apexpy', 'geo', 'apex', '20150101000000',
               '--height', '300', '-i', self.singlefile, '-o', self.outfile]
        data = self.execute_command_line(cmd)

        # Test the outfile existance and values
        assert data is not None, 'error executing: {:s}'.format(' '.join(cmd))
        np.testing.assert_allclose(data, [57.47145462, 93.62657928], rtol=1e-4)
        return

    @pytest.mark.parametrize("height, out_list",
                             [("300", [57.47145462, 93.62657928]),
                              ("100 --refh=300", [56.01779556, 93.35305023])])
    def test_convert_stdin_stdout_w_height_flags(self, height, out_list):
        """Test use of pipe input to command-line call with height flags.

        Parameters
        ----------
        height : str
            String specifying height with command line options
        out_list : list
            List of expected output values

        """
        # Build and execute the apexpy command line call
        cmd = ''.join(['echo 60 15 | python -m apexpy geo apex 2015 --height ',
                       '{:s}'.format(height)])
        cmd_kwargs = {'shell': True, 'stdout': subprocess.PIPE}
        stdout, _ = self.execute_command_line(cmd, cmd_kwargs, True)

        assert stdout is not None, 'error executing: {:s}'.format(' '.join(cmd))
        np.testing.assert_allclose(np.array(stdout.split(b' '), dtype=float),
                                   out_list, rtol=1e-4)
        return

    def test_convert_mlt(self):
        """Test magnetic local time conversion."""
        # Build and execute the apexpy command line call
        cmd = ['python', '-m', 'apexpy', 'geo', 'mlt', '20150101000000',
               '--height', '300', '-i', self.singlefile, '-o', self.outfile]
        data = self.execute_command_line(cmd)

        # Test the outfile existance and values
        assert data is not None, 'error executing: {:s}'.format(' '.join(cmd))
        np.testing.assert_allclose(data, [57.469547, 1.06324], rtol=1e-4)
        return

    @pytest.mark.parametrize("date_str", [("201501010"), ("2015010100000")])
    def test_invalid_date(self, date_str):
        """Test raises ValueError with an invalid input date.

        Parameters
        ----------
        date_str : str
           Input date string

        """
        # Build and execute the command
        cmd = 'echo 60 15 | python -m apexpy geo apex {:s}'.format(date_str)
        cmd_kwargs = {'shell': True, 'stderr': subprocess.PIPE}
        _, stderr = self.execute_command_line(cmd, cmd_kwargs, True)

        # Evaluate the error output
        assert stderr is not None, 'error executing: {:s}'.format(' '.join(cmd))
        assert b'ValueError' in stderr, 'invalid date error not raised'
        return

    def test_mlt_nodatetime(self):
        """Test raises ValueError when time not provided for MLT calc."""
        # Build and execute the command
        cmd = 'echo 60 15 | python -m apexpy geo mlt 20150101'
        cmd_kwargs = {'shell': True, 'stderr': subprocess.PIPE}
        _, stderr = self.execute_command_line(cmd, cmd_kwargs, True)

        # Evaluate the error output
        assert stderr is not None, 'error executing: {:s}'.format(' '.join(cmd))
        assert b'ValueError' in stderr, 'invalid time error not raised'
        return

    @pytest.mark.parametrize("coords", [("foobar apex"), ("geo foobar")])
    def test_invalid_coord(self, coords):
        """Test raises error when bad coordinate input provided.

        Parameters
        ----------
        coords : str
           Input/output coordinate pairs

        """
        # Build and execute the command
        cmd = 'echo 60 15 | python -m apexpy {:s} 2015'.format(coords)
        cmd_kwargs = {'shell': True, 'stderr': subprocess.PIPE}
        _, stderr = self.execute_command_line(cmd, cmd_kwargs, True)

        # Evaluate the error output
        assert stderr is not None, 'error executing: {:s}'.format(' '.join(cmd))
        assert b'invalid choice' in stderr, 'invalid coord error not raised'
        return
