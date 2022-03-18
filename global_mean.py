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
import logging

cdo = Cdo()

def main(args):
    """The main EC-mean4 code"""

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
    filename = 'gm_reference.yml'
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

    # add missing unit definition
    units_extra_definition(units)

    # check if required variables are there: use interface file
    # check into first file, and load also model variable units
    isavail = {}
    varunit = {}
    for var in var_all:
        infile = make_input_filename(
            ECEDIR, var, expname, year1, year1, face)
        retavail, retunit = vars_are_there(infile, [var], face)
        isavail = {**isavail, **retavail}
        varunit = {**varunit, **retunit}


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

            # conversion debug
            logging.debug(var)
            logging.debug(varunit[var] + ' ---> ' + ref[var]['units'])

            # adjust integrated quantities
            new_units = units_are_integrals(varunit[var], ref[var])
            
            # unit conversion
            units_conversion = units_converter(new_units, ref[var]['units'])

            # sign adjustment (for heat fluxes)
            units_conversion['factor'] = units_conversion['factor'] * units_are_down(ref[var]) 

            # conversion debug
            logging.debug(units_conversion)

            if 'derived' in face[var].keys():
                cmd = face[var]['derived']
                dervars = (",".join(re.findall("[a-zA-Z]+", cmd)))
                cdop.selectname(dervars)
                cdop.expr(var, cmd)
            else:
                cdop.selectname(var)

            # Introduce grid fixes specifying type of file (atm or oce)
            cdop.setdomain(domain=face[var]['component'])
            cdop.fixgrid()

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

            varmean[var] = (mean(a) + units_conversion['offset']) * units_conversion['factor']
            if ftrend:
                vartrend[var] = np.polyfit(yrange, a, 1)[0]
            if fverb:
                print('Average', var, varmean[var])

    # loop on the variables to create the output table
    global_table = []
    for var in var_atm + var_oce:
        beta = face[var]
        gamma = ref[var]
        beta['value'] = varmean[var]
        if ftrend:
            beta['trend'] = vartrend[var]
            out_sequence = [var, beta['varname'], gamma['units'], beta['value'],
                        beta['trend'],
                        float(gamma['observations']['val']),
                        gamma['observations'].get('data',''),
                        gamma['observations'].get('years','')]
        else:
            out_sequence = [var, beta['varname'], gamma['units'], beta['value'],
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
    parser.add_argument('-v', '--loglevel', type=str, default='ERROR',
                    help='define the level of logging. default: error')
    args = parser.parse_args()

    # log level with logging
    # currently basic definition trought the text
    loglevel = args.loglevel.upper()
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level)

    main(args)
