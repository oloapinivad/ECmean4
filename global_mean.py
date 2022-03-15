#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
 python3 version of ECmean global mean tool
 Using a reference file from yaml and cdo bindings

 @author Paolo Davini (p.davini@isac.cnr.it), March 2022
 @author Jost von Hardenberg (jost.hardenberg@polito.it), March 2022
'''

import os
import sys
import argparse
from statistics import mean
from pathlib import Path
import yaml
from tabulate import tabulate
import numpy as np
from cdo import Cdo
from functions import vars_are_there, is_number, load_config_file
from cdopipe import CdoPipe

cdo = Cdo()

def write_tuning_table(linefile, varmean, var_table, expname, year1, year2, ref):
    """Write results appending one line to a text file"""
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
            print('{:12.5f}'.format(varmean[var] * ref[var].get('factor',1)), end=' ', file=f)
        print(file=f)


def make_input_filename(dr, var, expname, year, ref):
    """Generate appropriate input filename for a variable"""
    if(ref[var].get('domain','atm')=='oce'):
        fname = dr / 'output/nemo' / \
                 f'{expname}_oce_1m_T_{year}-{year}.nc'
    else:
        fname = dr / 'output/oifs' / \
                 f'{expname}_atm_cmip6_1m_{year}-{year}.nc'
    return str(fname)


def main(args):

    assert sys.version_info >= (3, 5)

    expname = args.exp
    year1 = args.year1
    year2 = args.year2
    fverb = not args.silent
    ftable = args.line
    ftrend = args.trend
    modelname = args.model
    if year1 == year2: # Ignore if only one year requested
        ftrend = False

    # config file (looks for it in the same dir as the .py program file
    INDIR = Path(os.path.dirname(os.path.abspath(__file__)))
    cfg = load_config_file(INDIR)
    
    # define a few folders and create missing ones
    ECEDIR = Path(os.path.expandvars(cfg['dirs']['exp']), expname)
    TABDIR = Path(os.path.expandvars(cfg['dirs']['tab']))
    TMPDIR = Path(os.path.expandvars(cfg['dirs']['tmp']))
    os.makedirs(TABDIR, exist_ok=True)

    # prepare grid description file
    INIFILE=str(ECEDIR / f'ICMGG{expname}INIT')
    OCEINIFILE=cfg['areas']['oce']

    # Init CdoPipe object to use in the following, specifying the LM and SM files
    mycdo = CdoPipe(tempdir=TMPDIR)
    mycdo.make_grids(INIFILE, OCEINIFILE)

    # load reference data
    with open(INDIR / 'reference.yml', 'r') as file:
        ref = yaml.load(file, Loader=yaml.FullLoader)

    # list of vars on which to work
    var_atm = cfg['global']['atm_vars']
    var_oce = cfg['global']['oce_vars']
    var_table = cfg['global']['tab_vars']

    # make sure all requested vars are available (use first year)
    # first find all needed variables (including those needed for derived ones)

    # create a list for all atm variables, avoid duplicates
    var_all = list(set(var_atm + var_table))

    infile = make_input_filename(ECEDIR, var_atm[0], expname, year1, ref)

    # Check if vars are available
    infile = make_input_filename(ECEDIR, var_atm[0], expname, year1, ref)
    isavail = vars_are_there(infile, var_all, ref)
    infile = make_input_filename(ECEDIR, var_oce[0], expname, year1, ref)
    isavail = {**isavail, **vars_are_there(infile, var_oce, ref)} # python>=3.5 syntax for joining 2 dicts

    var_all = list(set(var_all + var_oce))
    
    # main loop
    varmean = {}
    vartrend = {}
    for var in var_all:
        if not isavail[var]:
            varmean[var] = float("NaN")
            vartrend[var] = float("NaN")
        else:
            # Refresh cdo pipe and specify type of file (atm or oce)
            mycdo.start(ref[var].get('domain','atm'))

            if 'derived' in ref[var].keys():
                mycdo.expr(var, ref[var]['derived'])

            mycdo.selname(var)

            # land/sea variables
            mycdo.masked_mean(ref[var].get('total','global'))
                
            mycdo.timmean()

            a = []
            # loop on years: call CDO to perform all the computations
            yrange = range(year1, year2+1)
            for year in yrange:
                infile =  make_input_filename(ECEDIR, var, expname, year, ref)
                x = mycdo.output(infile)
                a.append(x)
            varmean[var] = mean(a)
            if ftrend:
                vartrend[var] = np.polyfit(yrange, a, 1)[0]
            if fverb:
                print('Average', var, varmean[var])

    # loop on the variables to create the output table
    global_table = []
    for var in var_atm + var_oce:
        beta = ref[var]
        beta['value'] = varmean[var] * float(beta.get('factor', 1))
        if ftrend:
            beta['trend'] = vartrend[var] * float(beta.get('factor', 1))
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

    if ftrend:
        head = ['Variable', 'Longname', 'Units', modelname, 'Trend', 'Obs.', 'Dataset', 'Years']
    else:
        head = ['Variable', 'Longname', 'Units', modelname, 'Obs.', 'Dataset', 'Years']

    # write the file with tabulate: cool python feature
    tablefile = TABDIR / f'global_mean_{expname}_{year1}_{year2}.txt'
    if fverb: print(tablefile)
    with open(tablefile, 'w') as f:
        f.write(tabulate(global_table, headers=head, tablefmt='orgtbl'))

    # Print appending one line to table (for tuning)
    linefile = TABDIR / 'global_means.txt'
    if fverb: print(linefile)
    if args.output:
        linefile = args.output
        ftable = True
    if ftable:
        write_tuning_table(linefile, varmean, var_table, expname, year1, year2, ref)

    # clean
    os.unlink(mycdo.GRIDFILE)
    #cdo.cleanTempDir()


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
    parser.add_argument('-o', '--output', metavar='FILE', type=str, default='',
                    help='path of output one-line table')
    parser.add_argument('-m', '--model', type=str, default='EC-Earth4',
                    help='model name')
    args = parser.parse_args()

    main(args)
