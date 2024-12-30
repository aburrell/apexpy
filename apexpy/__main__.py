# -*- coding: utf-8 -*-

"""Entry point for the Command Line Interface"""

import argparse
import datetime as dt
import numpy as np
import sys

import apexpy

STDIN = sys.stdin.buffer
STDOUT = sys.stdout.buffer


def main():
    """Entry point for the script"""

    # Construct the description and parser for command-line arguments
    desc = 'Converts between geodetic, modified apex, quasi-dipole and MLT'
    parser = argparse.ArgumentParser(description=desc, prog='apexpy')

    parser.add_argument('source', metavar='SOURCE',
                        choices=['geo', 'apex', 'qd', 'mlt'],
                        help='Convert from {geo, apex, qd, mlt}')
    parser.add_argument('dest', metavar='DEST',
                        choices=['geo', 'apex', 'qd', 'mlt'],
                        help='Convert to {geo, apex, qd, mlt}')
    desc = ''.join(['YYYY[MM[DD[HHMMSS]]] date/time for IGRF coefficients, ',
                    'time part required for MLT calculations'])
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

    # Get the command line arguements
    args = parser.parse_args()
    arg_array = np.loadtxt(args.file_in, ndmin=2)

    # Test the input arguments
    if 'mlt' in [args.source, args.dest] and len(args.date) < 14:
        desc = 'full date/time YYYYMMDDHHMMSS required for MLT calculations'
        raise ValueError(desc)
    if 9 <= len(args.date) and len(args.date) <= 13:
        desc = 'full date/time must be given as YYYYMMDDHHMMSS, not ' \
               + 'YYYYMMDDHHMMSS'[:len(args.date)]
        raise ValueError(desc)

    # Format the time input
    in_time = dt.datetime.strptime(args.date,
                                   '%Y%m%d%H%M%S'[:len(args.date) - 2])

    # Run the desired apex conversion
    apex_obj = apexpy.Apex(date=in_time, refh=args.refh)
    lats, lons = apex_obj.convert(arg_array[:, 0], arg_array[:, 1], args.source,
                                  args.dest, args.height, datetime=in_time)

    # Save the output to a file.  Use the name for non-stdout inputs
    if args.file_out.name.lower().find('stdout') >= 0:
        fout_name = args.file_out
    else:
        fout_name = args.file_out.name
    np.savetxt(fout_name, np.column_stack((lats, lons)), fmt='%.8f')

    return


if __name__ == '__main__':
    sys.exit(main())
