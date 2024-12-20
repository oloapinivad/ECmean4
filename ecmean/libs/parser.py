#!/usr/bin/env python3
'''
Shared functions for Ecmean4: parsers arguments from command line
'''

import argparse
from ecmean import __version__


def parse_arguments(args, script):
    """Parse CLI arguments for global mean

    Args:
        args: command line argument
        script: script for which parsing is used, i.e. 'gm' or 'pi'
    Return:
        Parsed arguments
    """

    # common configuration to be parsed
    parser = argparse.ArgumentParser(
        description='ECmean global mean diagnostics for Global Climate models')
    parser.add_argument('exp', metavar='EXP', type=str, help='experiment ID')
    parser.add_argument('year1', metavar='Y1', type=int, help='starting year')
    parser.add_argument('year2', metavar='Y2', type=int, help='final year')
    parser.add_argument('-i', '--interface', type=str, default='',
                        help='interface (overrides config.yml)')
    parser.add_argument('-c', '--config', type=str, default='',
                        help='config file')
    parser.add_argument('-j', dest="numproc", type=int, default=1,
                        help='number of processors to use')
    parser.add_argument('-m', '--model', type=str, default='',
                        help='model name')
    parser.add_argument('-e', '--ensemble', type=str, default='r1i1p1f1',
                        help='variant label (ripf number for cmor)')
    parser.add_argument('-s', '--silent', action='store_true',
                        help='do not print anything to std output')
    parser.add_argument('--addnan', action='store_true',
                        help='provide figures also where observations are missing')
    parser.add_argument('-v', '--loglevel', type=str, default='WARNING',
                        help='define the level of logging.')
    parser.add_argument('-o', '--outputdir', type=str,
                        help='force the output directory')
    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + __version__)

    # specific to global mean
    if script == 'gm':
        parser.add_argument('-t', '--trend', action='store_true',
                            help='compute trends')
        parser.add_argument('-l', '--line', action='store_true',
                            help='appends also single line to a table')

    # specific to performance indices
    if script == 'pi':
        parser.add_argument('-k', '--climatology', type=str, default='EC23',
                            help='climatology to be compared. default: EC23. Options: [RK08, EC22, EC23]',
                            choices=['RK08', 'EC22', 'EC23'])
        parser.add_argument('-r', '--resolution', type=str, default='',
                            help='climatology resolution')

    return parser.parse_args(args)
