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

import argparse
from pathlib import Path
import logging
from multiprocessing import Process, Manager
from statistics import mean
from time import time
from tabulate import tabulate
import numpy as np
import xarray as xr
from ecmean.libs.general import weight_split, write_tuning_table, Diagnostic, getdomain, numeric_loglevel, get_variables_to_load
from ecmean.libs.files import var_is_there, get_inifiles, load_yaml, make_input_filename, config_diagnostic
from ecmean.libs.formula import formula_wrapper
from ecmean.libs.masks import masks_dictionary, areas_dictionary, masked_meansum
from ecmean.libs.units import units_extra_definition, units_wrapper
from ecmean.libs.ncfixers import xr_preproc

import dask
dask.config.set(scheduler="synchronous")


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

        # get domain
        vdom = getdomain(var, face)

        # compute weights
        weights = util[vdom + '_areas']

        # get the list of the variables to be loaded
        dervars = get_variables_to_load(var, face)

        # create input filenames
        infile = make_input_filename(
            var, dervars, face, diag)

        # chck if variables are available
        isavail, varunit = var_is_there(infile, var, face['variables'])

        if not isavail:
            varmean[var] = float("NaN")
            vartrend[var] = float("NaN")
        else:

            # perform the unit conversion extracting offset and factor
            offset, factor = units_wrapper(var, varunit, ref, face)

            # load the object
            xfield = xr.open_mfdataset(infile, preprocess=xr_preproc, chunks={'time': 12})

            # get the data-array field for the required var
            cfield = formula_wrapper(var, face, xfield)

            a = []
            # loop on years: 
            for year in diag.years_joined:

                # time selection and mean
                tfield = cfield.sel(time=(cfield.time.dt.year == year)).mean(dim='time')

                # final operation on the field
                x = masked_meansum(
                    tfield, var, weights, ref[var].get(
                        'total', 'global'), util['atm_mask'])
                a.append(x)

            varmean[var] = (mean(a) + offset) * factor
            if diag.ftrend:
                vartrend[var] = np.polyfit(diag.years_joined, a, 1)[0]
            if diag.fverb:
                print('Average', var, varmean[var])


def global_mean(exp, year1, year2, 
            config = 'config.yml',
            loglevel = 'WARNING',
            numproc = 1, 
            interface = None, model = None, ensemble = 'r1i1p1f1', 
            silent = None, trend = None, line = None,
            output = None):
    
    """The main ECmean4 global mean function

    :param exp: Experiment name or ID
    :param year1: Initial year
    :param year2: Final year
    :param config: configuration file, optional (default 'config.yml')
    :param loglevel: level of logging, optional (default 'WARNING')
    :param numproc: number of multiprocessing cores, optional (default '1')
    :param interface: interface file to be used, optional (default as specifified in config file)
    :param model: model to be analyzed, optional (default as specifified in config file)
    :param ensemble: ensemble member to be analyzed, optional (default as 'r1i1p1f1')
    :param silent: do not print anything to std output, optional
    :param trend: compute yearly trends, optional
    :param line: appends also single line to a table, optional
    :param output: output directory for the single line output, optional

    :returns: the global mean txt table as defined in the output

    """

    # create a name space with all the arguments to feed the Diagnostic class
    # This is not the neatest option, but it is very compact
    argv = argparse.Namespace(**locals())

    # set loglevel
    logging.basicConfig(level=numeric_loglevel(argv.loglevel))

    # get local directory
    indir = Path(os.path.dirname(os.path.abspath(__file__)))
    logging.info(indir)

    # define config dictionary, interface dictionary and diagnostic class
    cfg, face, diag = config_diagnostic(indir, argv)

    # Create missing folders
    os.makedirs(diag.TABDIR, exist_ok=True)

    # load reference data
    ref = load_yaml(indir / 'reference/gm_reference.yml')

    # list of vars on which to work
    var_atm = cfg['global']['atm_vars']
    var_oce = cfg['global']['oce_vars']
    var_ice = cfg['global']['ice_vars']
    var_table = cfg['global']['tab_vars']
    # var_all = list(set(var_atm + var_table + var_oce))
    var_all = list(
        dict.fromkeys(
            var_atm +
            var_table +
            var_oce +
            var_ice))  # python 3.7+, preserve order

    # Can probably be cleaned up further
    comp = face['model']['component']  # Get component for each domain

    # this required a change from the original file requirements of CDO version
    # now we have a mask file and two area files:
    # needs to be fixed and organized in the
    # config file in order to be more portable
    maskatmfile, atmareafile, oceareafile = get_inifiles(face, diag)

    # create util dictionary including mask and weights for both atmosphere
    # and ocean grids
    areas = areas_dictionary(comp, atmareafile, oceareafile)
    masks = masks_dictionary(comp, maskatmfile)
    util_dictionary = {**areas, **masks}

    # add missing unit definition
    units_extra_definition()

    # We now use a list for all the years
    diag.years_joined = list(range(diag.year1, diag.year2 + 1))

    # main loop: manager is required for shared variables
    mgr = Manager()

    # dictionaries are shared, so they have to be passed as functions
    varmean = mgr.dict()
    vartrend = mgr.dict()
    processes = []
    tic = time()

    # loop on the variables, create the parallel process
    for varlist in weight_split(var_all, diag.numproc):
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
        print('Done in {:.4f} seconds'.format(toc - tic))

    # loop on the variables to create the output table
    global_table = []
    for var in var_atm + var_oce + var_ice:
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
    tablefile = diag.TABDIR / \
        f'global_mean_{diag.expname}_{diag.modelname}_{diag.ensemble}_{diag.year1}_{diag.year2}.txt'
    if diag.fverb:
        print(tablefile)
    with open(tablefile, 'w', encoding='utf-8') as f:
        f.write(tabulate(global_table, headers=head, tablefmt='orgtbl'))

    # Print appending one line to table (for tuning)
    if diag.ftable:
        if diag.fverb:
            print(diag.linefile)
        write_tuning_table(diag.linefile, varmean, var_table, diag, ref)


def gm_parse_arguments(args):
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
    parser.add_argument('-v', '--loglevel', type=str, default='WARNING',
                        help='define the level of logging.')
    parser.add_argument('-j', dest="numproc", type=int, default=1,
                        help='number of processors to use')
    parser.add_argument('-e', '--ensemble', type=str, default='r1i1p1f1',
                        help='variant label (ripf number for cmor)')
    parser.add_argument('-i', '--interface', type=str, default='',
                        help='interface (overrides config.yml)')

    return parser.parse_args(args)

def gm_entry_point():

    """
    Command line interface to run the global_mean function
    """

    # read arguments from command line
    args = gm_parse_arguments(sys.argv[1:])

    global_mean(exp = args.exp, year1 = args.year1, year2 = args.year2,
                numproc = args.numproc,
                silent = args.silent, trend = args.trend, line = args.line,
                output = args.output, loglevel = args.loglevel,
                interface = args.interface, config = args.config,
                model = args.model, ensemble = args.ensemble)


if __name__ == "__main__":

    sys.exit(gm_entry_point())
