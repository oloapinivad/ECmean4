#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
 python3 version of ECmean global mean tool
 Using a reference file from yaml and cdo bindings

 @author Paolo Davini (p.davini@isac.cnr.it), March 2022
 @author Jost von Hardenberg (jost.hardenberg@polito.it), March 2022
'''

import os
import argparse
from statistics import mean
from pathlib import Path
import yaml
from tabulate import tabulate
import numpy as np
from cdo import Cdo

from functions import vars_are_there, is_number

def write_tuning_table(linefile, varstat, var_table, expname, year1, year2, ref):

    if not os.path.isfile(linefile):
        with open(linefile, 'w') as f:
            print('%exp from   to ', end='', file=f)
            for var in var_table:
                print('{:>12s}'.format(var), end=' ', file=f)
            print('\n%             ', end=' ', file=f)
            for var in var_table:
                print('{:>12s}'.format(ref[var]['units']), end=' ', file=f)
            print(file=f)

    with open(linefile, 'a') as f:
        print(expname,'{:4d} {:4d} '.format(year1, year2), end='', file=f)
        for var in var_table:
            print('{:12.5f}'.format(varstat[var] * ref[var].get('factor',1)), end=' ', file=f)
        print(file=f)


def make_filename(dr, var, expname, year, ref):
    if(ref[var].get('domain','atm')=='oce'):
        fname = dr / 'output/nemo' / \
                 f'{expname}_oce_1m_T_{year}-{year}.nc'
    else:
        fname = dr / 'output/oifs' / \
                 f'{expname}_atm_cmip6_1m_{year}-{year}.nc'
    return str(fname)

def main(args):

    expname = args.exp
    year1 = args.year1
    year2 = args.year2
    fverb = not args.silent
    ftable = args.line
    ftrend = args.trend

    cdo = Cdo()

    # config file (looks for it in the same dir as the .py program file
    INDIR = Path(os.path.dirname(os.path.abspath(__file__)))
    with open(INDIR / 'config.yml', 'r') as file:
        cfg = yaml.load(file, Loader=yaml.FullLoader)

    ECEDIR = Path(os.path.expandvars(cfg['dirs']['exp']), expname)
    TABDIR = Path(os.path.expandvars(cfg['dirs']['tab']))
    TMPDIR = Path(os.path.expandvars(cfg['dirs']['tmp']))
    os.makedirs(TABDIR, exist_ok=True)
    print(TMPDIR)

    # prepare grid description file
    GRIDFILE=str(TMPDIR / 'grid.txt')
    INIFILE=str(ECEDIR / f'ICMGG{expname}INIT')
    griddes = cdo.griddes(input=INIFILE)
    with open(GRIDFILE, 'w') as f:
        for line in griddes:
            print(line, file=f)

    # prepare ATM LSM
    LMFILE = cdo.selname('LSM', input=f'-setgridtype,regular {INIFILE}', options='-t ecmwf')
    SMFILE = cdo.mulc('-1', input=f'-subc,1 {LMFILE}')
    GAFILE = cdo.gridarea(input=f'-setgridtype,regular {LMFILE}')

    # prepare OCE areas
    OCEGAFILE = cdo.expr('area=e1t*e2t', input=cfg['areas']['oce']) 

    # reference data
    with open(INDIR / 'reference.yml', 'r') as file:
        ref = yaml.load(file, Loader=yaml.FullLoader)

    # list of vars on which to work
    var_field = cfg['global']['atm_vars']['field']
    var_radiation = cfg['global']['atm_vars']['radiation']
    var_ocean = cfg['global']['oce_vars']
    var_table = cfg['global']['tab_vars']

    # make sure all requested vars are available (use first year)
    # first find all needed variables (including those needed for derived ones)

    # create a list for all atm variables, avoid duplicates
    var_all = list(set(var_field + var_radiation + var_table))

    # Check if vars are available
    infile = make_filename(ECEDIR, var_all[0], expname, year1, ref)
    isavail = vars_are_there(infile, var_all, ref)
    infile = make_filename(ECEDIR, var_ocean[0], expname, year1, ref)
    isavail = {**isavail, **vars_are_there(infile, var_ocean, ref)} # python>=3.5 syntax for joining 2 dicts

    var_all = list(set(var_all + var_ocean))

    # loop
    varstat = {}
    trend = {}
    for var in var_all:
        if not isavail[var]:
            varstat[var] = float("NaN")
            trend[var] = float("NaN")
        else:
            a = []

            # check if var is derived
            if 'derived' in ref[var].keys():
                cmd = ref[var]['derived']
                der = f'-expr,{var}={cmd}'
            else:
                der = ''

            # ocean variables require specifying grid areas
            if(ref[var].get('domain','atm')=='oce'):
                pre = '-setgridarea,' + OCEGAFILE
            else:
                pre = ''

            # land/sea variables
            mask = ''
            op='-fldmean'
            if 'total' in ref[var].keys():
                mask_type = ref[var]['total']
                if mask_type == 'land':
                    mask = f'-mul {GAFILE} -mul {LMFILE}'
                    op='-fldsum'
                elif mask_type in ['sea', 'ocean']:
                    mask = f'-mul {GAFILE} -mul {SMFILE}'
                    op='-fldsum'

            # loop on years: call CDO to perform all the computations
            yrange = range(year1, year2+1)
            for year in yrange:
                infile =  make_filename(ECEDIR, var, expname, year, ref)
                cmd = f'-timmean {op} {mask} -setgridtype,regular ' \
                      f'-setgrid,{GRIDFILE} -selname,{var} {der} {pre} {infile}'
                x = float(cdo.output(input=cmd)[0])
                a.append(x)
            varstat[var] = mean(a)
            if ftrend:
                trend[var] = np.polyfit(yrange, a, 1)[0]
            if fverb:
                print('Average', var, mean(a))

    # define options for the output table
    if ftrend:
        head = ['Variable', 'Longname', 'Units', 'EC-Earth4', 'Trend', 'Obs.', 'Dataset', 'Years']
    else:
        head = ['Variable', 'Longname', 'Units', 'EC-Earth4', 'Obs.', 'Dataset', 'Years']
    global_table = []

    # loop on the variables to create the table
    for var in var_field + var_radiation + var_ocean:
        beta = ref[var]
        beta['value'] = varstat[var] * float(beta.get('factor', 1))
        if ftrend:
            beta['trend'] = trend[var] * float(beta.get('factor', 1))
            out_sequence = [var, beta['varname'], beta['units'], beta['value'],
                        beta['trend'],
                        float(beta['observations']['val']),
                        beta['observations'].get('data',''),
                        beta['observations'].get('years','')]
        else:
            out_sequence = [var, beta['varname'], beta['units'], beta['value'],
                        float(beta['observations']['val']),
                        beta['observations'].get('data',''),
                        beta['observations'].get('years','')]
        global_table.append(out_sequence)

    # write the file with tabulate: cool python feature
    tablefile = TABDIR / f'global_mean_{expname}_{year1}_{year2}.txt'
    print(tablefile)
    with open(tablefile, 'w') as f:
        f.write(tabulate(global_table, headers=head, tablefmt='orgtbl'))

    # Print appending one line to table (for tuning)
    if ftable:
        linefile = TABDIR / 'global_means.txt'
        if args.output:
            linefile = args.output
        write_tuning_table(linefile, varstat, var_table, expname, year1, year2, ref)

    # clean
    os.unlink(GRIDFILE)
    os.unlink(LMFILE)
    os.unlink(SMFILE)
    os.unlink(GAFILE)
    cdo.cleanTempDir()


if __name__ == "__main__":
    # arguments
    parser = argparse.ArgumentParser(
        description='ECmean global mean diagnostics for EC-Earth4')
    parser.add_argument('exp', metavar='EXP', type=str, help='experiment ID')
    parser.add_argument('year1', metavar='Y1', type=int, help='starting year')
    parser.add_argument('year2', metavar='Y2', type=int, help='final year')
    parser.add_argument('-s', '--silent', action='store_true',
                    help='do not print anything to std output')
    parser.add_argument('-t', '--trend', action='store_true',
                    help='compute trends')
    parser.add_argument('-l', '--line', action='store_true',
                    help='appends also single line to a table')
    parser.add_argument('-o', '--output', type=str, default='',
                    help='path of output one-line table')
    args = parser.parse_args()

    main(args)
