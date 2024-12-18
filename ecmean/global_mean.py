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
import logging
import argparse
from multiprocessing import Process, Manager
from time import time
from tabulate import tabulate
import numpy as np
import xarray as xr
import yaml
import dask

from ecmean import Diagnostic, Supporter, UnitsHandler
from ecmean.libs.general import weight_split, write_tuning_table, get_domain, \
    check_time_axis, init_mydict, \
        check_var_interface, check_var_climatology, set_multiprocessing_start_method
from ecmean.libs.files import var_is_there, get_inifiles, load_yaml, make_input_filename
from ecmean.libs.formula import formula_wrapper
from ecmean.libs.masks import masked_meansum, select_region
from ecmean.libs.units import units_extra_definition
from ecmean.libs.ncfixers import xr_preproc
from ecmean.libs.parser import parse_arguments
from ecmean.libs.plotting import heatmap_comparison_gm, prepare_clim_dictionaries_gm
from ecmean.libs.loggy import setup_logger

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

    loggy = logging.getLogger(__name__)

    for var in varlist:

        # create empty nested dictionaries
        result = init_mydict(diag.seasons, diag.regions)
        trend = init_mydict(diag.seasons, diag.regions)

        # check if the variable is in the interface file
        if check_var_interface(var, face):

            # get domain
            domain = get_domain(var, face)

            # compute weights
            weights = getattr(util, domain + 'area')
            domain_mask = getattr(util, domain + 'mask')

            # get input files/fielf
            infile = make_input_filename(var, face, diag)

            # check if variables are available
            isavail, varunit = var_is_there(infile, var, face)

            if isavail:

                # perform the unit conversion extracting offset and factor
                units_handler = UnitsHandler(var, org_units=varunit, clim=ref, face=face)
                offset, factor = units_handler.offset, units_handler.factor

                # load the object
                if not isinstance(infile, (xr.DataArray, xr.Dataset)):
                    xfield = xr.open_mfdataset(infile, preprocess=xr_preproc, chunks={'time': 12})
                else:
                    xfield = infile

                # in case of big files with multi year, be sure of having opened the right records
                xfield = xfield.sel(time=xfield.time.dt.year.isin(diag.years_joined))

                # check time axis
                check_time_axis(xfield.time, diag.years_joined)

                # get the data-array field for the required var
                cfield = formula_wrapper(var, face, xfield).compute()

                for season in diag.seasons:

                    # copy of the full field
                    tfield = cfield.copy(deep=True)

                    if season != 'ALL':
                        tfield = tfield.sel(time=cfield.time.dt.season.isin(season))

                    if diag.ftrend:
                        # this does not consider continuous seasons for DJF, but JF+D
                        tfield = tfield.groupby('time.year').mean('time')
                    else:
                        tfield = tfield.mean(dim='time')

                    for region in diag.regions:

                        slicefield = select_region(tfield, region)
                        sliceweights = select_region(weights, region)
                        if isinstance(domain_mask, xr.DataArray):
                            slicemask = select_region(domain_mask, region)
                        else:
                            slicemask = 0.

                        # final operation on the field
                        avg = masked_meansum(
                            xfield=slicefield, weights=sliceweights, mask=slicemask,
                            operation=ref[var].get('operation', 'mean'),
                            mask_type=ref[var].get('mask', 'global'),
                            domain=domain)

                        # if dask delayead object, compute
                        if isinstance(avg, dask.array.core.Array):
                            avg = avg.compute()

                        result[season][region] = float((np.nanmean(avg) + offset) * factor)

                        if diag.ftrend:
                            if len(avg) == len(diag.years_joined):
                                trend[season][region] = np.polyfit(diag.years_joined, avg, 1)[0]
                        if season == 'ALL' and region == 'Global':
                            loggy.info('Average: %s %s %s %s', var, season, region, result[season][region])

        # nested dictionary, to be redifend as a dict to remove lambdas
        varmean[var] = result
        vartrend[var] = trend


def global_mean(exp, year1, year2,
                config='config.yml',
                loglevel='WARNING',
                numproc=1,
                interface=None, model=None, ensemble='r1i1p1f1',
                addnan=False,
                silent=None, trend=None, line=None,
                outputdir=None, xdataset=None):
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
    :param add_nan: add to the final plots also fields which cannot be compared against observations
    :param silent: do not print anything to std output, optional
    :param trend: compute yearly trends, optional
    :param line: appends also single line to a table, optional
    :param outputdir: output directory for the single line output, optional
    :param xdataset: xarray dataset - already open - to be used without looking for files

    :returns: the global mean txt table as defined in the output

    """

    # create a name space with all the arguments to feed the Diagnostic class
    # This is not the neatest option, but it is very compact
    funcname = __name__
    argv = argparse.Namespace(**locals())

    # set loglevel
    loggy = setup_logger(level=argv.loglevel)

    # set dask and multiprocessing fork
    plat, mprocmethod = set_multiprocessing_start_method()
    loggy.info('Running on %s and multiprocessing method set as "%s"', plat, mprocmethod)

    # start time
    tic = time()

    # initialize the diag class, load the inteface and the reference file
    diag = Diagnostic(argv)
    face = load_yaml(diag.interface)
    ref = load_yaml(diag.reffile)

    # check that everyhing is there
    check_var_climatology(diag.var_all, ref.keys())

    # Create missing folders
    os.makedirs(diag.tabdir, exist_ok=True)
    os.makedirs(diag.figdir, exist_ok=True)

    # Can probably be cleaned up further
    comp = face['model']['component']  # Get component for each domain

    # get file info
    inifiles = get_inifiles(face, diag)


    # add missing unit definition
    units_extra_definition()

    # create util dictionary including mask and weights for both atmosphere
    # and ocean grids
    util_dictionary = Supporter(comp, inifiles['atm'], inifiles['oce'],
                                areas=True, remap=False)

    # main loop: manager is required for shared variables
    mgr = Manager()

    # dictionaries are shared, so they have to be passed as functions
    varmean = mgr.dict()
    vartrend = mgr.dict()
    processes = []

    # loop on the variables, create the parallel process
    for varlist in weight_split(diag.var_all, diag.numproc):
        core = Process(
            target=gm_worker, args=(util_dictionary, ref, face, diag,
                                    varmean, vartrend, varlist))
        core.start()
        processes.append(core)

    # wait for the processes to finish
    for proc in processes:
        proc.join()

    toc = time()

    # evaluate tic-toc time  of execution
    loggy.info('Analysis done in %.4f seconds', toc - tic)

    tic = time()


    # loop on the variables to create the output table
    global_table = []
    for var in diag.var_all:
        gamma = ref[var]
        # get the predefined value or the ALL Global one
        if isinstance(gamma['obs'], dict):
            tabval = gamma['obs']['ALL']['Global']
            outval = str(tabval['mean']) + '\u00B1' + str(tabval['std'])
        else:
            outval = gamma['obs']

        if 'year1' in gamma.keys():
            years = str(gamma['year1']) + '-' + str(gamma['year2'])
        else:
            raise ValueError('Year1 and Year2 are not defined in the reference file')

        out_sequence = [var, gamma['longname'], gamma['units'], varmean[var]['ALL']['Global']]
        if diag.ftrend:
            out_sequence = out_sequence + [vartrend[var]['ALL']['Global']]
        out_sequence = out_sequence + [outval, gamma.get('dataset', ''), years]
        global_table.append(out_sequence)

    # prepare the header for the table
    head = ['Variable', 'Longname', 'Units', diag.modelname]
    if diag.ftrend:
        head = head + ['Trend']
    head = head + ['Obs.', 'Dataset', 'Years']

    # write the file with tabulate
    tablefile = os.path.join(diag.tabdir,
                             f'global_mean_{diag.expname}_{diag.modelname}_{diag.ensemble}_{diag.year1}_{diag.year2}.txt')
    loggy.info('Table file is: %s', tablefile)
    with open(tablefile, 'w', encoding='utf-8') as out:
        out.write(tabulate(global_table, headers=head, stralign='center', tablefmt='orgtbl'))

    # required to avoid problem with multiprocessing
    ordered = {}
    for var in diag.var_all:
        ordered[var] = varmean[var]

    # dump the yaml file for global_mean, including all the seasons
    yamlfile = os.path.join(diag.tabdir,
                            f'global_mean_{diag.expname}_{diag.modelname}_{diag.ensemble}_{diag.year1}_{diag.year2}.yml')
    loggy.info('YAML file is: %s', tablefile)
    with open(yamlfile, 'w', encoding='utf-8') as file:
        yaml.safe_dump(ordered, file, default_flow_style=False, sort_keys=False)

    # prepare dictionaries for global mean
    obsmean, obsstd, data2plot, units_list = prepare_clim_dictionaries_gm(ordered, ref, diag.var_all, diag.seasons, diag.regions)

    # call the heatmap routine for a plot
    mapfile = os.path.join(diag.figdir,
                           f'global_mean_{diag.expname}_{diag.modelname}_r1i1p1f1_{diag.year1}_{diag.year2}.pdf')
    loggy.info('Figure file is: %s', mapfile)

    diag_dict = {'modelname': diag.modelname, 'expname': diag.expname,
                 'year1': diag.year1, 'year2': diag.year2,
                 'units_list': units_list}

    heatmap_comparison_gm(data_dict=data2plot, mean_dict=obsmean, std_dict=obsstd,
                          diag=diag_dict, filemap=mapfile, addnan=diag.addnan)

    # Print appending one line to table (for tuning)
    if diag.ftable:
        loggy.info('Line file is: %s', diag.linefile)
        write_tuning_table(diag.linefile, varmean, diag.var_table, diag, ref)

    toc = time()
    # evaluate tic-toc time of postprocessing
    loggy.info("Postproc done in %.4f seconds", toc - tic)
    print('ECmean4 Global Mean succesfully computed!')


def gm_entry_point():
    """
    Command line interface to run the global_mean function
    """

    # read arguments from command line
    args = parse_arguments(sys.argv[1:], script='gm')

    global_mean(exp=args.exp, year1=args.year1, year2=args.year2,
                numproc=args.numproc,
                trend=args.trend, line=args.line,
                loglevel=args.loglevel,
                interface=args.interface, config=args.config,
                model=args.model, ensemble=args.ensemble,
                addnan=args.addnan,
                outputdir=args.outputdir)
