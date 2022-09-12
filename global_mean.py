#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
   python3 version of ECmean global mean tool.
   Using a reference file from yaml and Xarray

   @author Paolo Davini (p.davini@isac.cnr.it), Sep 2022.
   @author Jost von Hardenberg (jost.hardenberg@polito.it), Sep 2022
"""

__author__ = "Paolo Davini (p.davini@isac.cnr.it), Sep 2022."

import os
import sys
import re
import argparse
from pathlib import Path
import logging
from multiprocessing import Process, Manager
from statistics import mean
from time import time
from tabulate import tabulate
import numpy as np
import xarray as xr

from ecmean import var_is_there, eval_formula, \
    masks_dictionary, areas_dictionary, get_inifiles, load_yaml, \
    units_extra_definition, units_are_integrals, \
    units_converter, directions_match, chunks, write_tuning_table, \
    Diagnostic, getdomain, make_input_filename, masked_meansum, \
    xr_preproc
import dask
dask.config.set(scheduler="synchronous")


def parse_arguments(args):
    """Parse CLI arguments for global mean"""

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
    parser.add_argument('-c', '--config', type=str, default='',
                        help='config file')
    parser.add_argument('-v', '--loglevel', type=str, default='ERROR',
                        help='define the level of logging.')
    parser.add_argument('-j', dest="numproc", type=int, default=1,
                        help='number of processors to use')
    parser.add_argument('-e', '--ensemble', type=str, default='r1i1p1f1',
                        help='variant label (ripf number for cmor)')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='activate cdo debugging')
    parser.add_argument('-i', '--interface', type=str, default='',
                        help='interface (overrides config.yml)')

    return parser.parse_args(args)

def gm_worker(util, ref, face, diag, varmean, vartrend, varlist):
    """Main parallel diagnostic worker for global mean
   
    Args: 
	util: the utility dictionary, including mask and weights
	ref: the reference dictionary for the global mean
	face: the interface to be used to access the data
	diag: the diagnostic class object
	varmean: the dictionary for the global mean (empty)
	vartrend: the dictionary for the trends (empty)
	varlist: the variable on which compute the global mean

    Returns:
	vartrend and varmean under the form of a dictionaries

    """

    for var in varlist:

        # compute weights
        vdom =  getdomain(var, face)
        weights = util[vdom + '_areas'] 

        if 'derived' in face['variables'][var].keys():
            cmd = face['variables'][var]['derived']
            dervars = re.findall("[a-zA-Z]+", cmd)
        else:
            dervars = [var]

        infile = make_input_filename(var, dervars, diag.year1, diag.year1, face, diag)
        isavail, varunit = var_is_there(infile, var, face['variables'])

        if not isavail:
                varmean[var] = float("NaN")
                vartrend[var] = float("NaN")
        else:
        
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
            logging.debug('Offset %f, Factor %f', offset, factor)

            a = []
            # loop on years: using xarray to perform the computations
            # could be replaced by open_mfdataset but need to handle the time-mean
            yrange = range(diag.year1, diag.year2+1)
            for year in yrange:

                infile = make_input_filename(var, dervars, year, year, face, diag)
                xfield = xr.open_mfdataset(infile, preprocess=xr_preproc)

                # time selection for longer dataset!
                xfield = xfield.sel(time=(xfield.time.dt.year == year))
                if 'derived' in face['variables'][var].keys():
                    cmd = face['variables'][var]['derived']
                    outfield = eval_formula(cmd, xfield)
                else:
                    outfield = xfield[var]
         
                tfield = outfield.mean(dim = 'time')
                x = masked_meansum(tfield, var, weights, ref[var].get('total', 'global'), util['atm_mask'])
                a.append(x)

            varmean[var] = (mean(a) + offset) * factor
            if diag.ftrend:
                vartrend[var] = np.polyfit(yrange, a, 1)[0]
            if diag.fverb:
                print('Average', var, varmean[var])

def main(argv):
    """The main ECmean4 global mean code"""

    #assert sys.version_info >= (3, 7)

    args = parse_arguments(argv)
    # log level with logging
    # currently basic definition trought the text
    numeric_level = getattr(logging, args.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % args.loglevel)
    logging.basicConfig(level=numeric_level)

    INDIR = Path(os.path.dirname(os.path.abspath(__file__)))
    # config file (looks for it in the same dir as the .py program file
    if args.config:
        cfg = load_yaml(args.config)
    else:
        cfg = load_yaml(INDIR / 'config.yml')

    # Setup all common variables, directories from arguments and config files
    diag = Diagnostic(args, cfg)

    # Create missing folders
    os.makedirs(diag.TABDIR, exist_ok=True)

    # load reference data
    ref = load_yaml(INDIR / 'gm_reference.yml')

    # loading the var-to-file interface
    face = load_yaml(INDIR / Path('interfaces', f'interface_{diag.interface}.yml'))

    # list of vars on which to work
    var_atm = cfg['global']['atm_vars']
    var_oce = cfg['global']['oce_vars']
    var_table = cfg['global']['tab_vars']
    # var_all = list(set(var_atm + var_table + var_oce))
    var_all = list(dict.fromkeys(var_atm + var_table + var_oce))  # python 3.7+, preserve order

    # We need 
    # Can probably be cleaned up further
    comp = face['model']['component']  # Get component for each domain

    # this required a change from the original file requirements of CDO version
    # now we have a mask file and two area files: need to be fixed and organized in the
    # config file in order to be more portable
    maskatmfile, atmareafile, oceareafile = get_inifiles(face, diag)

    # create util dictionary including mask and weights for both atmosphere and ocean grids
    areas = areas_dictionary(comp, atmareafile, oceareafile)
    masks = masks_dictionary(comp, maskatmfile)
    util_dictionary = {**areas, **masks}

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
        p = Process(target=gm_worker, args=(util_dictionary, ref, face, diag,
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
    tablefile = diag.TABDIR / f'global_mean_{diag.expname}_{diag.modelname}_{diag.ensemble}_{diag.year1}_{diag.year2}.txt'
    if diag.fverb:
        print(tablefile)
    with open(tablefile, 'w', encoding='utf-8') as f:
        f.write(tabulate(global_table, headers=head, tablefmt='orgtbl'))

    # Print appending one line to table (for tuning)
    if diag.ftable:
        if diag.fverb:
            print(diag.linefile)
        write_tuning_table(diag.linefile, varmean, var_table, diag, ref)

if __name__ == "__main__":

    sys.exit(main(sys.argv[1:]))
