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
from time import time
from tabulate import tabulate
import numpy as np
import xarray as xr
from ecmean.libs.general import weight_split, write_tuning_table, get_domain, numeric_loglevel, \
                                get_variables_to_load, check_time_axis
from ecmean.libs.files import var_is_there, get_inifiles, load_yaml, make_input_filename, init_diagnostic
from ecmean.libs.formula import formula_wrapper
from ecmean.libs.masks import masks_dictionary, masked_meansum, select_region
from ecmean.libs.areas import areas_dictionary
from ecmean.libs.units import units_extra_definition, units_wrapper
from ecmean.libs.ncfixers import xr_preproc
from ecmean.libs.parser import parse_arguments

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
        domain = get_domain(var, face)

        # compute weights
        weights = util[domain + '_areas']
        domain_mask = util[domain + '_mask']

        # get the list of the variables to be loaded
        dervars = get_variables_to_load(var, face)

        # create input filenames
        infile = make_input_filename(
            var, dervars, face, diag)

        # chck if variables are available
        isavail, varunit = var_is_there(infile, var, face['variables'])

        # store NaN in dict (can't use defaultdict due to multiprocessing)
        # result = defaultdict(lambda: defaultdict(lambda : float('NaN')))
        result = {}
        for season in diag.seasons:
            result[season] = {}
            for region in diag.regions:
                result[season][region] = float('NaN')

        if isavail:

            # perform the unit conversion extracting offset and factor
            offset, factor = units_wrapper(var, varunit, ref, face)

            # load the object
            xfield = xr.open_mfdataset(infile, preprocess=xr_preproc, chunks={'time': 12})

            # in case of big files with multi year, be sure of having opened the right records
            xfield = xfield.sel(time=xfield.time.dt.year.isin(diag.years_joined))

            # check time axis
            check_time_axis(xfield.time, diag.years_joined)

            # get the data-array field for the required var
            cfield = formula_wrapper(var, face, xfield).compute()

            for season in diag.seasons:

                if season != 'ALL':
                    tfield = cfield.sel(time=cfield.time.dt.season.isin(season)).mean(dim='time')
                else: 
                    tfield = cfield.mean(dim='time')

                for region in diag.regions:

                    slicefield = select_region(tfield, region)
                    sliceweights = select_region(weights, region)
                    if ref[var].get('mask', 'global') != 'Global':
                        slicemask = select_region(domain_mask, region)
                    else:
                        slicemask = 0.

                    #tfield = cfield.resample(time='1Y').mean('time')

                    # final operation on the field
                    a = masked_meansum(
                        xfield=slicefield, weights=sliceweights, mask=slicemask,
                        operation=ref[var].get('operation', 'mean'),
                        mask_type=ref[var].get('mask', 'global'),
                        domain=domain)
                    
                    try:
                        x = a.compute()
                    except:
                        x = a

                    result[season][region] = (np.nanmean(x) + offset) * factor

                    if diag.ftrend:
                        vartrend[var] = np.polyfit(diag.years_joined,x, 1)[0]
                    if diag.fverb and season == 'ALL' and region == 'Global':
                        print('Average', var, season, region, result[season][region])

        # nested dictionary, to be redifend as a dict to remove lambdas
        varmean[var] = result


def global_mean(exp, year1, year2,
                config='config.yml',
                loglevel='WARNING',
                numproc=1,
                interface=None, model=None, ensemble='r1i1p1f1',
                silent=None, trend=None, line=None,
                output=None):
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
    cfg, face, diag = init_diagnostic(indir, argv)

    # Create missing folders
    os.makedirs(diag.TABDIR, exist_ok=True)

    # load reference data
    ref = load_yaml(indir / 'reference/gm_reference_EC23.yml')

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

    # get file info
    inifiles = get_inifiles(face, diag)

    # add missing unit definition
    units_extra_definition()

    # create util dictionary including mask and weights for both atmosphere
    # and ocean grids
    areas = areas_dictionary(comp, inifiles['atm'], inifiles['oce'])
    masks = masks_dictionary(comp, inifiles['atm']['maskfile'], inifiles['oce']['maskfile'])
    util_dictionary = {**areas, **masks}

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
        # beta = face['variables'][var]
        gamma = ref[var]
        # get the predifined valuue or the ALL GLobal one
        if isinstance(gamma['obs'], dict):
            ff = gamma['obs']['ALL']['Global']
            outval = str(ff['mean']) + u'\u00B1' + str(ff['std'])
        else:
            outval = gamma['obs']

        if 'year1' in gamma.keys():
            years = str(gamma['year1']) + '-' + str(gamma['year2'])

        out_sequence = [var, gamma['longname'], gamma['units'], varmean[var]['ALL']['Global']]
        if diag.ftrend:
            out_sequence = out_sequence + [vartrend[var]]
        out_sequence = out_sequence + [outval, gamma.get('dataset', ''), years]
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
        f.write(tabulate(global_table, headers=head, stralign='center', tablefmt='orgtbl'))

    # Print appending one line to table (for tuning)
    if diag.ftable:
        if diag.fverb:
            print(diag.linefile)
        write_tuning_table(diag.linefile, varmean, var_table, diag, ref)


def gm_entry_point():
    """
    Command line interface to run the global_mean function
    """

    # read arguments from command line
    args = parse_arguments(sys.argv[1:], script='gm')

    global_mean(exp=args.exp, year1=args.year1, year2=args.year2,
                numproc=args.numproc,
                silent=args.silent, trend=args.trend, line=args.line,
                output=args.output, loglevel=args.loglevel,
                interface=args.interface, config=args.config,
                model=args.model, ensemble=args.ensemble)


if __name__ == "__main__":

    sys.exit(gm_entry_point())
