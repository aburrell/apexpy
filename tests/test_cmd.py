#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest
from pytest import approx
import numpy as np
import subprocess
import sys
import os

OLDPY = sys.version_info < (3, 5)
PYEXE = sys.executable
TESTFN = 'output.txt'
INDIR = os.path.dirname(__file__)
INFN1 = os.path.join(INDIR, 'test_convert.txt')
INFN2 = os.path.join(INDIR, 'test_convert_single_line.txt')


@pytest.mark.skipif(OLDPY, reason='Python >= 3.5 for this test')
def test_module_invocation(tmp_path):
    outfile = str(tmp_path / TESTFN)
    subprocess.check_call([PYEXE, '-m', 'apexpy', 'geo', 'apex', '2015',
                             '--height', '300', '-i', INFN1,
                             '-o', outfile])

    data = np.loadtxt(outfile)
    assert data == approx(np.array([[57.469547, 93.639816],
                                  [58.522701, 94.044762],
                                  [59.571465, 94.477257]]), rel=1e-4)


@pytest.mark.skipif(OLDPY, reason='Python >= 3.5 for this test')
@pytest.mark.parametrize('dt', ['2015', '201501', '20150101', '20150101000000'])
def test_convert_dates(dt, tmp_path):
    outfile = str(tmp_path / TESTFN)
    subprocess.check_call([PYEXE, '-m', 'apexpy', 'geo', 'apex', dt,
                           '--height', '300',
                           '-i', INFN1, '-o', outfile])

    data = np.loadtxt(outfile)
    assert data == approx(np.array([[57.469547, 93.639816],
                                      [58.522701, 94.044762],
                                      [59.571465, 94.477257]]), rel=1e-4)


def test_convert_single_line(tmp_path):
    outfile = str(tmp_path / TESTFN)
    subprocess.check_call([PYEXE, '-m', 'apexpy', 'geo', 'apex',
                             '20150101000000', '--height', '300',
                             '-i', INFN2, '-o', outfile])

    data = np.loadtxt(outfile)
    assert data == approx([57.469547, 93.639816], rel=1e-4)


@pytest.mark.skipif(OLDPY, reason='Python >= 3.5 for this test')
def test_convert_stdin_stdout():
    ret = subprocess.check_output(['apexpy', 'geo', 'apex', '2015', '--height', '300'],
                          input='60 15',
                          universal_newlines=True)

    assert np.array(ret.split(' '), dtype=float) == approx([57.469547, 93.639816], rel=1e-4)


@pytest.mark.skipif(OLDPY, reason='Python >= 3.5 for this test')
def test_convert_refh():
    ret = subprocess.check_output(['apexpy', 'geo', 'apex', '2000', '--height', '100', '--refh=300'],
                          input='60 15',
                          universal_newlines=True)

    assert np.array(ret.split(' '), dtype=float) == approx([55.94841766, 94.1068344], rel=1e-4)


def test_convert_mlt(tmp_path):
    outfile = str(tmp_path / TESTFN)
    subprocess.check_call([PYEXE, '-m', 'apexpy', 'geo', 'mlt',
                             '20150101000000', '--height', '300',
                             '-i', INFN2, '-o', outfile])

    data = np.loadtxt(outfile)
    assert data == approx([57.469547, 1.06324], rel=1e-4)


@pytest.mark.skipif(OLDPY, reason='Python >= 3.5 for this test')
@pytest.mark.parametrize('bad', ['201501010', '2015010100000'])
def test_invalid_date(bad):
    ret = subprocess.run(['apexpy', 'geo', 'apex', bad],
                          input='60 15', stderr=subprocess.PIPE,
                          universal_newlines=True)

    assert 'ValueError' in ret.stderr


@pytest.mark.skipif(OLDPY, reason='Python >= 3.5 for this test')
def test_mlt_nodatetime():
    ret = subprocess.run(['apexpy', 'geo', 'mlt', '20150101'],
                          input='60 15', stderr=subprocess.PIPE,
                          universal_newlines=True)

    assert 'ValueError' in ret.stderr


@pytest.mark.skipif(OLDPY, reason='Python >= 3.5 for this test')
@pytest.mark.parametrize('bad1, bad2', [('foobar', 'apex'), ('geo', 'foobar')])
def test_invalid_coord(bad1, bad2):
    ret = subprocess.run(['apexpy', bad1, bad2, '2015'],
                          input='60 15', stderr=subprocess.PIPE,
                          universal_newlines=True)

    assert 'invalid choice' in ret.stderr



if __name__ == '__main__':
    pytest.main([__file__])
