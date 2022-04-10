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
import re
import argparse
from statistics import mean
from pathlib import Path
import logging
from multiprocessing import Process, Manager
import copy
from time import time
from tabulate import tabulate
import numpy as np
from ecmean import var_is_there, load_yaml, \
                      make_input_filename, write_tuning_table, \
                      units_extra_definition, units_are_integrals, \
                      units_converter, directions_match, chunks, \
                      Diagnostic, getinifiles, getdomain
from cdopipe import CdoPipe


def worker(cdopin, ref, face, diag, varmean, vartrend, varlist):
    """Main parallel diagnostic worker"""

    cdop = copy.copy(cdopin)  # Create a new local instance to avoid overlap of the cdo pipe
    for var in varlist:

        # check if required variables are there: use interface file
        # check into first file, and load also model variable units
        # with this implementation, files are accessed multiple times for each variables
        # this is simpler but might slow down the code
        infile = make_input_filename(var, diag.year1, diag.year1, face, diag)
        isavail, varunit = var_is_there(infile, var, face['variables'])

        if not isavail:
            varmean[var] = float("NaN")
            vartrend[var] = float("NaN")
        else:
            # Refresh cdo pipe
            cdop.start()

            # conversion debug
            logging.debug(var)
            logging.debug(varunit + ' ---> ' + ref[var]['units'])

            # adjust integrated quantities
            adjusted_units = units_are_integrals(varunit, ref[var])

            # unit conversion
            offset, factor = units_converter(adjusted_units, ref[var]['units'])

            # sign adjustment (for heat fluxes)
            factor = factor * directions_match(face['variables'][var], ref[var])

            # conversion debug
            logging.debug(offset, factor)

            if 'derived' in face['variables'][var].keys():
                cmd = face['variables'][var]['derived']
                dervars = (",".join(re.findall("[a-zA-Z]+", cmd)))
                cdop.selectname(dervars)
                cdop.expr(var, cmd)
            else:
                cdop.selectname(var)

            # Introduce grid fixes specifying type of file (atm or oce)
            domain = getdomain(var, face)
            cdop.fixgrid(domain=domain)

            # land/sea variables
            cdop.masked_meansum(ref[var].get('total', 'global'))
            cdop.timmean()

            a = []
            # loop on years: call CDO to perform all the computations
            yrange = range(diag.year1, diag.year2+1)
            for year in yrange:
                infile = make_input_filename(var, year, year, face, diag)
                x = cdop.output(infile, keep=True)
                a.append(x)

            varmean[var] = (mean(a) + offset) * factor
            if diag.ftrend:
                vartrend[var] = np.polyfit(yrange, a, 1)[0]
            if diag.fverb:
                print('Average', var, varmean[var])


def main(args):
    """The main EC-mean4 code"""

    assert sys.version_info >= (3, 7)

    # config file (looks for it in the same dir as the .py program file
    INDIR = Path(os.path.dirname(os.path.abspath(__file__)))
    cfg = load_yaml(INDIR / 'config.yml')

    # Setup all common variables, directories from arguments and config files
    diag = Diagnostic(args, cfg)

    # Create missing folders
    os.makedirs(diag.TABDIR, exist_ok=True)

    # Init CdoPipe object to use in the following
    cdop = CdoPipe()

    # load reference data
    ref = load_yaml(INDIR / 'gm_reference.yml')

    # loading the var-to-file interface
    face = load_yaml(INDIR / Path('interfaces', f'interface_{diag.interface}.yml'))

    # New bunch of functions to set grids, create correction command, masks and areas
    # Can probably be cleaned up further
    comp = face['model']['component']  # Get component for each domain
    atminifile, ocegridfile, oceareafile = getinifiles(face, diag)
    cdop.set_gridfixes(atminifile, ocegridfile, oceareafile, comp['atm'], comp['oce'])
    cdop.make_atm_masks(comp['atm'], atminifile)

    # list of vars on which to work
    var_atm = cfg['global']['atm_vars']
    var_oce = cfg['global']['oce_vars']
    var_table = cfg['global']['tab_vars']
    # var_all = list(set(var_atm + var_table + var_oce))
    var_all = list(dict.fromkeys(var_atm + var_table + var_oce))  # python 3.7+, preserve order

    # add missing unit definition
    units_extra_definition()

    # main loop: manager is required for shared variables
    mgr = Manager()

    # dictionaries are shared, so they have to be passed as functions
    varmean = mgr.dict()
    vartrend = mgr.dict()
    processes = []
    tic = time()

    # loop on the variables, create the parallel process
    for varlist in chunks(var_all, diag.numproc):
        p = Process(target=worker, args=(cdop, ref, face, diag,
                                         varmean, vartrend, varlist))
        p.start()
        processes.append(p)

    # wait for the processes to finish
    for proc in processes:
        proc.join()

    toc = time()

    # evaluate tic-toc time  of execution
    if diag.fverb:
        print('Done in {:.4f} seconds'.format(toc-tic))

    # loop on the variables to create the output table
    global_table = []
    for var in var_atm + var_oce:
        beta = face['variables'][var]
        gamma = ref[var]

        out_sequence = [var, beta['varname'], gamma['units'], varmean[var]]
        if diag.ftrend:
            out_sequence = out_sequence + [vartrend[var]]
        out_sequence = out_sequence + [float(gamma['val']),
                                       gamma.get('data', ''),
                                       gamma.get('years', '')]
        global_table.append(out_sequence)

    # prepare the header for the table
    head = ['Variable', 'Longname', 'Units', diag.modelname]
    if diag.ftrend:
        head = head + ['Trend']
    head = head + ['Obs.', 'Dataset', 'Years']

    # write the file with tabulate: cool python feature
    tablefile = diag.TABDIR / f'global_mean_{diag.expname}_{diag.year1}_{diag.year2}.txt'
    if diag.fverb:
        print(tablefile)
    with open(tablefile, 'w', encoding='utf-8') as f:
        f.write(tabulate(global_table, headers=head, tablefmt='orgtbl'))

    # Print appending one line to table (for tuning)
    if diag.fverb:
        print(diag.linefile)
    if diag.ftable:
        write_tuning_table(diag.linefile, varmean, var_table, diag.expname,
                           diag.year1, diag.year2, face, ref)

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
    parser.add_argument('-m', '--model', type=str, default='',
                        help='model name')
    parser.add_argument('-v', '--loglevel', type=str, default='ERROR',
                        help='define the level of logging.')
    parser.add_argument('-j', dest="numproc", type=int, default=1,
                        help='number of processors to use')

    args = parser.parse_args()

    # log level with logging
    # currently basic definition trought the text
    loglevel = args.loglevel.upper()
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level)

    main(args)
