# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals

import os
import subprocess

import numpy as np


def setup_function(function):
    try:
        os.remove('tests/output.txt')
    except:
        pass

teardown_function = setup_function


def test_module_invocation():
    p = subprocess.Popen(['python', '-m', 'apexpy', 'geo', 'apex', '20150224', '--height', '300',
                          '-i', 'tests/test_convert.txt', '-o', 'tests/output.txt'])
    p.communicate()
    p.wait()
    data = np.loadtxt('tests/output.txt')
    np.testing.assert_allclose(data, [[57.4761, 93.5572], [58.5332, 93.9607], [59.5852, 94.3897]], rtol=1e-4)


def test_convert_YYYY():
    p = subprocess.Popen(['apexpy', 'geo', 'apex', '2015', '--height', '300',
                          '-i', 'tests/test_convert.txt', '-o', 'tests/output.txt'])
    p.communicate()
    p.wait()
    data = np.loadtxt('tests/output.txt')
    np.testing.assert_allclose(data, [[57.4761, 93.5572], [58.5332, 93.9607], [59.5852, 94.3897]], rtol=1e-4)


def test_convert_YYYYMM():
    p = subprocess.Popen(['apexpy', 'geo', 'apex', '201501', '--height', '300',
                          '-i', 'tests/test_convert.txt', '-o', 'tests/output.txt'])
    p.communicate()
    p.wait()
    data = np.loadtxt('tests/output.txt')
    np.testing.assert_allclose(data, [[57.4761, 93.5572], [58.5332, 93.9607], [59.5852, 94.3897]], rtol=1e-4)


def test_convert_YYYYMMDD():
    p = subprocess.Popen(['apexpy', 'geo', 'apex', '20150101', '--height', '300',
                          '-i', 'tests/test_convert.txt', '-o', 'tests/output.txt'])
    p.communicate()
    p.wait()
    data = np.loadtxt('tests/output.txt')
    np.testing.assert_allclose(data, [[57.4761, 93.5572], [58.5332, 93.9607], [59.5852, 94.3897]], rtol=1e-4)


def test_convert_YYYYMMDDHHMMSS():
    p = subprocess.Popen(['apexpy', 'geo', 'apex', '20150101000000', '--height', '300',
                          '-i', 'tests/test_convert.txt', '-o', 'tests/output.txt'])
    p.communicate()
    p.wait()
    data = np.loadtxt('tests/output.txt')
    np.testing.assert_allclose(data, [[57.4761, 93.5572], [58.5332, 93.9607], [59.5852, 94.3897]], rtol=1e-4)


def test_convert_single_line():
    p = subprocess.Popen(['apexpy', 'geo', 'apex', '20150101000000', '--height', '300',
                          '-i', 'tests/test_convert_single_line.txt', '-o', 'tests/output.txt'])
    p.communicate()
    p.wait()
    data = np.loadtxt('tests/output.txt')
    np.testing.assert_allclose(data, [57.4761, 93.5572], rtol=1e-4)


def test_convert_stdin_stdout():
    p = subprocess.Popen('echo 60 15 | apexpy geo apex 2015', shell=True, stdout=subprocess.PIPE)
    stdout, _ = p.communicate()
    p.wait()
    assert b'57.47612194 93.55719875' in stdout


def test_convert_mlt_a2m():
    p = subprocess.Popen(['apexpy', 'geo', 'mlt', '20150101000000', '--height', '300',
                          '-i', 'tests/test_convert_single_line.txt', '-o', 'tests/output.txt'])
    p.communicate()
    p.wait()
    data = np.loadtxt('tests/output.txt')
    np.testing.assert_allclose(data, [9.057565, 9.790899, 10.524232], rtol=1e-4)
