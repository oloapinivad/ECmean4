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
from functions import *
from cdopipe import CdoPipe

cdo = Cdo()

def main(args):

    assert sys.version_info >= (3, 7)

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
    cdop = CdoPipe()
    cdop.make_grids(INIFILE, OCEINIFILE)

    # load reference data
    filename = 'reference.yml'
    with open(INDIR / filename, 'r') as file:
        ref = yaml.load(file, Loader=yaml.FullLoader)
    
    # loading the var-to-file interface
    filename = 'interface_ece4.yml'
    with open(INDIR / filename, 'r') as file:
        face = yaml.load(file, Loader=yaml.FullLoader)

    # list of vars on which to work
    var_atm = cfg['global']['atm_vars']
    var_oce = cfg['global']['oce_vars']
    var_table = cfg['global']['tab_vars']
    #var_all = list(set(var_atm + var_table + var_oce))
    var_all = list(dict.fromkeys(var_atm + var_table + var_oce)) # python 3.7+, preserver order

    # make sure all requested vars are available (use first year)
    # first find all needed variables (including those needed for derived ones)
    isavail = {}
    for var in var_all:
        infile = make_input_filename(
            ECEDIR, var, expname, year1, year1, face)
        isavail = {**isavail, **vars_are_there(infile, [var], face)}

    # main loop
    varmean = {}
    vartrend = {}
    for var in var_all:
        if not isavail[var]:
            varmean[var] = float("NaN")
            vartrend[var] = float("NaN")
        else:
            # Refresh cdo pipe
            cdop.start()

            if 'derived' in face[var].keys():
                cmd = face[var]['derived']
                dervars = (",".join(re.findall("[a-zA-Z]+", cmd)))
                cdop.selectname(dervars)
                cdop.expr(var, cmd)
            else:
                cdop.selectname(var)

            # Introduce grid fixes specifying type of file (atm or oce)
            cdop.fixgrid(domain=ref[var].get('domain','atm'))

            # land/sea variables
            cdop.masked_mean(ref[var].get('total','global'))
            cdop.timmean()

            a = []
            # loop on years: call CDO to perform all the computations
            yrange = range(year1, year2+1)
            for year in yrange:
                infile =  make_input_filename(ECEDIR, var, expname, year, year, face)
                x = cdop.output(infile, keep=True)
                a.append(x)

            varmean[var] = mean(a)
            if ftrend:
                vartrend[var] = np.polyfit(yrange, a, 1)[0]
            if fverb:
                print('Average', var, varmean[var])

    # loop on the variables to create the output table
    global_table = []
    for var in var_atm + var_oce:
        beta = face[var]
        gamma = ref[var]
        beta['value'] = varmean[var] * float(gamma.get('factor', 1))
        if ftrend:
            beta['trend'] = vartrend[var] * float(gamma.get('factor', 1))
            out_sequence = [var, beta['varname'], beta['units'], beta['value'],
                        beta['trend'],
                        float(gamma['observations']['val']),
                        gamma['observations'].get('data',''),
                        gamma['observations'].get('years','')]
        else:
            out_sequence = [var, beta['varname'], beta['units'], beta['value'],
                        float(gamma['observations']['val']),
                        gamma['observations'].get('data',''),
                        gamma['observations'].get('years','')]
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
    cdop.cdo.cleanTempDir()


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
