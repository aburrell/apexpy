# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import
from __future__ import unicode_literals

import numpy as np
import os
import subprocess

outfile = os.path.abspath(os.path.join(os.path.dirname(__file__), 'output.txt'))
infile = os.path.abspath(os.path.join(os.path.dirname(__file__), 'test_convert.txt'))
singlefile = os.path.abspath(os.path.join(os.path.dirname(__file__), 'test_convert_single_line.txt'))


def test_convert_YYYY():
    subprocess.run(['python', '-m', 'apexpy', 'geo', 'apex', '2015',
                    '--height', '300', '-i', infile, '-o', outfile])

    assert os.path.isfile(outfile)
    data = np.loadtxt(outfile)
    np.testing.assert_allclose(data, [[57.469547, 93.639816],
                                      [58.522701, 94.044762],
                                      [59.571465, 94.477257]], rtol=1e-4)


def test_convert_YYYYMM():
    subprocess.run(['python', '-m', 'apexpy', 'geo', 'apex', '201501',
                    '--height', '300', '-i', infile, '-o', outfile])

    assert os.path.isfile(outfile)
    data = np.loadtxt(outfile)
    np.testing.assert_allclose(data, [[57.469547, 93.639816],
                                      [58.522701, 94.044762],
                                      [59.571465, 94.477257]], rtol=1e-4)


def test_convert_YYYYMMDD():
    subprocess.run(['python', '-m', 'apexpy', 'geo', 'apex', '20150101',
                    '--height', '300', '-i', infile, '-o', outfile])

    assert os.path.isfile(outfile)
    data = np.loadtxt(outfile)
    np.testing.assert_allclose(data, [[57.469547, 93.639816],
                                      [58.522701, 94.044762],
                                      [59.571465, 94.477257]], rtol=1e-4)


def test_convert_YYYYMMDDHHMMSS():
    subprocess.run(['python', '-m', 'apexpy', 'geo', 'apex', '20150101000000',
                    '--height', '300', '-i', infile, '-o', outfile])

    assert os.path.isfile(outfile)
    data = np.loadtxt(outfile)
    np.testing.assert_allclose(data, [[57.469547, 93.639816],
                                      [58.522701, 94.044762],
                                      [59.571465, 94.477257]], rtol=1e-4)


def test_convert_single_line():
    subprocess.run(['python', '-m', 'apexpy', 'geo', 'apex', '20150101000000',
                    '--height', '300', '-i', singlefile, '-o', outfile])

    assert os.path.isfile(outfile)
    data = np.loadtxt(outfile)
    np.testing.assert_allclose(data, [57.469547, 93.639816], rtol=1e-4)


def test_convert_stdin_stdout():
    """ Test use of pipe input to command-line call
    """
    pipe = subprocess.Popen('echo 60 15 | apexpy geo apex 2015 --height 300',
                            shell=True, stdout=subprocess.PIPE)
    stdout, _ = pipe.communicate()
    pipe.wait()
    np.testing.assert_allclose(np.array(stdout.split(b' '), dtype=float),
                               [57.469547, 93.639816], rtol=1e-4)


def test_convert_refh():
    pipe = subprocess.Popen(
        'echo 60 15 | apexpy geo apex 2000 --height 100 --refh=300',
        shell=True, stdout=subprocess.PIPE)
    stdout, _ = pipe.communicate()
    pipe.wait()
    np.testing.assert_allclose(np.array(stdout.split(b' '), dtype=float),
                               [55.94841766, 94.1068344], rtol=1e-4)


def test_convert_mlt():
    pipe = subprocess.Popen(['python', '-m', 'apexpy', 'geo', 'mlt',
                             '20150101000000', '--height', '300', '-i',
                             singlefile, '-o', outfile])
    pipe.communicate()
    pipe.wait()
    assert os.path.isfile(outfile)
    data = np.loadtxt(outfile)
    np.testing.assert_allclose(data, [57.469547, 1.06324], rtol=1e-4)


def test_invalid_date():
    pipe = subprocess.Popen('echo 60 15 | apexpy geo apex 201501010',
                            shell=True, stderr=subprocess.PIPE)
    _, stderr = pipe.communicate()
    pipe.wait()
    assert b'ValueError' in stderr

    pipe = subprocess.Popen('echo 60 15 | apexpy geo apex 2015010100000',
                            shell=True, stderr=subprocess.PIPE)
    _, stderr = pipe.communicate()
    pipe.wait()
    assert b'ValueError' in stderr


def test_mlt_nodatetime():
    pipe = subprocess.Popen('echo 60 15 | apexpy geo mlt 20150101', shell=True,
                            stderr=subprocess.PIPE)
    _, stderr = pipe.communicate()
    pipe.wait()
    assert b'ValueError' in stderr


def test_invalid_coord():
    pipe = subprocess.Popen('echo 60 15 | apexpy foobar apex 2015', shell=True,
                            stderr=subprocess.PIPE)
    _, stderr = pipe.communicate()
    pipe.wait()
    assert b'invalid choice' in stderr

    pipe = subprocess.Popen('echo 60 15 | apexpy geo foobar 2015', shell=True,
                            stderr=subprocess.PIPE)
    _, stderr = pipe.communicate()
    pipe.wait()
    assert b'invalid choice' in stderr
