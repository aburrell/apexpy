# -*- coding: utf-8 -*-

"""Entry point for the CLI"""

from __future__ import division, absolute_import

import sys
import argparse
import datetime as dt
import numpy as np

import apexpy

try:
    # Python 3
    STDIN = sys.stdin.buffer
    STDOUT = sys.stdout.buffer
except AttributeError:
    # Python 2
    STDIN = sys.stdin
    STDOUT = sys.stdout


def main():
    """Entry point for the script"""

    desc = 'Converts between geodetic, modified apex, quasi-dipole and MLT'
    parser = argparse.ArgumentParser(description=desc, prog='apexpy')

    parser.add_argument('source', metavar='SOURCE',
                        choices=['geo', 'apex', 'qd', 'mlt'],
                        help='Convert from {geo, apex, qd, mlt}')
    parser.add_argument('dest', metavar='DEST',
                        choices=['geo', 'apex', 'qd', 'mlt'],
                        help='Convert to {geo, apex, qd, mlt}')
    desc = 'YYYY[MM[DD[HHMMSS]]] date/time for IGRF coefficients, time part '
    desc += 'required for MLT calculations'
    parser.add_argument('date', metavar='DATE', help=desc)
    parser.add_argument('--height', dest='height', default=0, metavar='HEIGHT',
                        type=float, help='height for conversion')
    parser.add_argument('--refh', dest='refh', metavar='REFH', type=float,
                        default=0,
                        help='reference height for modified apex coordinates')

    parser.add_argument('-i', '--input', dest='file_in', metavar='FILE_IN',
                        type=argparse.FileType('r'), default=STDIN,
                        help='input file (stdin if none specified)')
    parser.add_argument('-o', '--output', dest='file_out', metavar='FILE_OUT',
                        type=argparse.FileType('wb'), default=STDOUT,
                        help='output file (stdout if none specified)')

    args = parser.parse_args()

    array = np.loadtxt(args.file_in, ndmin=2)

    if 'mlt' in [args.source, args.dest] and len(args.date) < 14:
        desc = 'full date/time YYYYMMDDHHMMSS required for MLT calculations'
        raise ValueError(desc)
    if 9 <= len(args.date) and len(args.date) <= 13:
        desc = 'full date/time must be given as YYYYMMDDHHMMSS, not ' \
               + 'YYYYMMDDHHMMSS'[:len(args.date)]
        raise ValueError(desc)
    datetime = dt.datetime.strptime(args.date,
                                    '%Y%m%d%H%M%S'[:len(args.date) - 2])
    A = apexpy.Apex(date=datetime, refh=args.refh)
    lats, lons = A.convert(array[:, 0], array[:, 1], args.source, args.dest,
                           args.height, datetime=datetime)
    np.savetxt(args.file_out, np.column_stack((lats, lons)), fmt='%.8f')


if __name__ == '__main__':
    sys.exit(main())
