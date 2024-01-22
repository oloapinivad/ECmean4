#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
 python3 version of ECmean performance indices tool
 Using a reference file from yaml and cdo bindings

 @author Paolo Davini (p.davini@isac.cnr.it), 2022
 @author Jost von Hardenberg (jost.hardenberg@polito.it), 2022
'''

import sys
import os
import logging
import argparse
from time import time
from multiprocessing import Process, Manager
import numpy as np
import xarray as xr
import yaml
import dask
from ecmean import Diagnostic, Supporter, UnitsHandler
from ecmean.libs.general import weight_split, get_domain, dict_to_dataframe, \
   check_time_axis, init_mydict, check_var_interface, check_var_climatology, \
   set_multiprocessing_start_method
from ecmean.libs.files import var_is_there, get_inifiles, load_yaml, \
    make_input_filename, get_clim_files
from ecmean.libs.formula import formula_wrapper
from ecmean.libs.masks import mask_field, select_region
from ecmean.libs.areas import guess_bounds
from ecmean.libs.units import units_extra_definition
from ecmean.libs.ncfixers import xr_preproc, adjust_clim_file
from ecmean.libs.plotting import heatmap_comparison_pi
from ecmean.libs.parser import parse_arguments
from ecmean.libs.loggy import setup_logger

dask.config.set(scheduler="synchronous")

def pi_worker(util, piclim, face, diag, field_3d, varstat, varlist):
    """Main parallel diagnostic worker for performance indices

    Args:
        util: the utility dictionary, including mask, area and remap weights
        piclim: the reference climatology for the global mean
        face: the interface to be used to access the data
        diag: the diagnostic class object
        field_3d: the list of 3d vars to be handled differently
        varstat: the dictionary for the variable PI (empty)
        varlist: the variable on which compute the global mean

    Returns:
        vartrend under the form of a dictionaries

    """

    loggy = logging.getLogger(__name__)

    for var in varlist:

        # store NaN in dict (can't use defaultdict due to multiprocessing)
        result = init_mydict(diag.seasons, diag.regions)

        if check_var_interface(var, face):

            # get domain
            domain = get_domain(var, face)

            # get masks
            domain_mask = getattr(util, domain + 'mask')

            # check if required variables are there: use interface file
            # check into first file, and load also model variable units
            infile = make_input_filename(var, face, diag)

            # check if var is available
            isavail, varunit = var_is_there(infile, var, face)

            # if the variable is available
            if isavail:

                # perform the unit conversion extracting offset and factor
                units_handler = UnitsHandler(var, org_units=varunit, clim=piclim, face=face)
                offset, factor = units_handler.offset, units_handler.factor

                # open file: chunking on time only, might be improved
                if not isinstance(infile, (xr.DataArray, xr.Dataset)):
                    xfield = xr.open_mfdataset(infile, preprocess=xr_preproc, chunks={'time': 12})
                else:
                    xfield = infile

                # in case of big files with multi year, be sure of having opened the right records
                xfield = xfield.sel(time=xfield.time.dt.year.isin(diag.years_joined))

                # check time axis
                check_time_axis(xfield.time, diag.years_joined)

                # get the data-array field for the required var
                outfield = formula_wrapper(var, face, xfield)

                # mean over time and fixing of the units
                for season in diag.seasons:

                    loggy.debug(season)

                    # copy of the full field
                    tmean = outfield.copy(deep=True)

                    # get filenames for climatology
                    clim, vvvv = get_clim_files(piclim, var, diag, season)

                    # open climatology files, fix their metadata
                    cfield = adjust_clim_file(
                        xr.open_mfdataset(clim, preprocess=xr_preproc))
                    vfield = adjust_clim_file(
                        xr.open_mfdataset(vvvv, preprocess=xr_preproc), remove_zero=True)

                    # season selection
                    if season != 'ALL':
                        tmean = tmean.sel(time=tmean.time.dt.season.isin(season))
                        cfield = cfield.sel(time=cfield.time.dt.season.isin(season))
                        vfield = vfield.sel(time=vfield.time.dt.season.isin(season))

                    # averaging
                    tmean = tmean.mean(dim='time')
                    cfield = cfield.load()  # this ensure no crash with multiprocessing
                    vfield = vfield.load()

                    # safe check for old RK08 which has a different format
                    if diag.climatology != 'RK08':
                        cfield = cfield.mean(dim='time')
                        vfield = vfield.mean(dim='time')

                    tmean = tmean * factor + offset

                    # this computation is required to ensure that file access is
                    # done in a correct way. resource errors found if disabled in some
                    # specific configuration
                    tmean = tmean.compute()

                    # apply interpolation, if fixer is availble and with different
                    # grids
                    fix = getattr(util, f'{domain}fix')
                    remap = getattr(util, f'{domain}remap')

                    if fix:
                        tmean = fix(tmean, keep_attrs=True)
                    try:
                        final = remap(tmean, keep_attrs=True)
                    except ValueError:
                        loggy.error('Cannot interpolate %s with the current weights...', var)
                        continue

                    # vertical interpolation
                    if var in field_3d:

                        # xarray interpolation on plev, forcing to be in Pascal
                        final = final.metpy.convert_coordinate_units('plev', 'Pa')
                        if set(final['plev'].data) != set(cfield['plev'].data):
                            loggy.warning('%s: Need to interpolate vertical levels...', var)
                            final = final.interp(plev=cfield['plev'].data, method='linear')

                            # safety check for missing values
                            sample = final.isel(lon=0, lat=0)
                            if np.sum(np.isnan(sample)) != 0:
                                loggy.warning('%s: You have NaN after the interpolation, this will affect your PIs...', var)
                                levnan = cfield['plev'].where(np.isnan(sample))
                                loggy.warning(levnan[~np.isnan(levnan)].data)

                        # zonal mean
                        final = final.mean(dim='lon')

                        # compute PI
                        complete = (final - cfield)**2 / vfield

                        # compute vertical bounds as weights
                        bounds_lev = guess_bounds(complete['plev'], name='plev')
                        bounds = abs(bounds_lev[:, 0] - bounds_lev[:, 1])
                        www = xr.DataArray(
                            bounds, coords=[
                                complete['plev']], dims=['plev'])

                        # vertical mean
                        outarray = complete.weighted(www).mean(dim='plev')

                    # horizontal averaging with land-sea mask
                    else:
                        
                        complete = (final - cfield)**2 / vfield
                        outarray = mask_field(xfield=complete,
                                            mask_type=piclim[var]['mask'],
                                            dom=domain, mask=domain_mask)

                    # loop on different regions
                    for region in diag.regions:

                        slicearray = select_region(outarray, region)
 
                        # latitude-based averaging
                        weights = np.cos(np.deg2rad(slicearray.lat))
                        out = slicearray.weighted(weights).mean().data
                         # store the PI
                        result[season][region] = round(float(out), 3)

                        # diagnostic
                        if region == 'Global':
                            logging.info('PI for', region, season, var, result[season][region])

        # nested dictionary, to be redifend as a dict to remove lambdas
        varstat[var] = result


def performance_indices(exp, year1, year2,
                        config='config.yml',
                        loglevel='WARNING',
                        numproc=1,
                        climatology='EC23',
                        interface=None, model=None, ensemble='r1i1p1f1',
                        silent=None, xdataset=None,
                        outputdir=None):
    """Main performance indices calculation

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
    :param climatology: climatology to be compared. default: EC23. Options: [RK08, EC22, EC23]
    :param xdataset: xarray dataset - already open - to be used without looking for files
    :param outputdir: if specified, override the target destination in the configuration files for both tables and figures

    :returns: the performance indices yaml file and heatmap

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
    piclim = load_yaml(diag.climfile)

    # check that everyhing is there
    check_var_climatology(diag.field_all, piclim.keys())

    # Create missing folders
    os.makedirs(diag.tabdir, exist_ok=True)
    os.makedirs(diag.figdir, exist_ok=True)

    # new bunch of functions to set grids, create correction command, masks
    # and areas
    comp = face['model']['component']  # Get component for each domain

    # all clim have the same grid, read from the first clim available and get
    # target grid
    clim, _ = get_clim_files(piclim, 'tas', diag, 'ALL')
    target_remap_grid = xr.open_dataset(clim)

    # get file info files
    inifiles = get_inifiles(face, diag)

    # add missing unit definitions
    units_extra_definition()

    # create remap dictionary with atm and oce interpolators
    util_dictionary = Supporter(comp, inifiles['atm'], inifiles['oce'],
                                areas=False, remap=True,
                                targetgrid=target_remap_grid)

    # main loop: manager is required for shared variables
    mgr = Manager()

    # dictionaries are shared, so they have to be passed as functions
    varstat = mgr.dict()
    processes = []

    toc = time()
    loggy.info('Preproc in {:.4f} seconds'.format(toc - tic))
    tic = time()

  
    # loop on the variables, create the parallel process
    for varlist in weight_split(diag.field_all, diag.numproc):
        #print(varlist)


        core = Process(
            target=pi_worker, args=(util_dictionary, piclim,
                face, diag, diag.field_3d, varstat, varlist))
        core.start()
        processes.append(core)

    # wait for the processes to finish
    for proc in processes:
        proc.join()


    toc = time()
    # evaluate tic-toc time  of execution
    loggy.info('Done in {:.4f} seconds'.format(toc - tic) + ' with ' + str(diag.numproc) + ' processors')

    tic = time()

    # order according to the original request the fields in the yaml file
    ordered = {}
    for item in diag.field_all:
        ordered[item] = varstat[item]

    # dump the yaml file for PI, including all the seasons (need to copy to avoid mess)
    yamlfile = diag.tabdir / \
        f'PI4_{diag.climatology}_{diag.expname}_{diag.modelname}_r1i1p1f1_{diag.year1}_{diag.year2}.yml'
    with open(yamlfile, 'w', encoding='utf-8') as file:
        yaml.safe_dump(ordered, file, default_flow_style=False, sort_keys=False)

    # to this date, only EC23 support comparison with CMIP6 data
    if diag.climatology == 'EC23':

        # uniform dictionaries
        filt_piclim = {}
        for k in piclim.keys():
            filt_piclim[k] = piclim[k]['cmip6']
            for f in ['models', 'year1', 'year2']:
                del filt_piclim[k][f]

        # set longname, reorganize the dictionaries
        plotted = {}
        cmip6 = {}
        longnames =[]
        for var in diag.field_all:
            plotted[piclim[var]['longname']] = ordered[var]
            cmip6[piclim[var]['longname']] = filt_piclim[var]
            longnames = longnames + [piclim[var]['longname']]

        # convert output dictionary to pandas dataframe
        data_table = dict_to_dataframe(plotted)
        loggy.debug(data_table)


        # relative pi with re-ordering of rows
        cmip6_table = data_table.div(dict_to_dataframe(cmip6).reindex(longnames))

        # compute the total PI mean
        cmip6_table.loc['Total PI'] = cmip6_table.mean()

        # reordering columns
        lll = [(x, y) for x in diag.seasons for y in diag.regions]
        cmip6_table = cmip6_table[lll]
        loggy.debug(cmip6_table)


        # call the heatmap routine for a plot
        mapfile = diag.figdir / \
            f'PI4_{diag.climatology}_{diag.expname}_{diag.modelname}_r1i1p1f1_{diag.year1}_{diag.year2}.pdf'
        # heatmap_comparison_old(data_table, diag, mapfile)
        heatmap_comparison_pi(cmip6_table, diag, mapfile)

    toc = time()
    # evaluate tic-toc time of postprocessing
    loggy.info('Postproc done in {:.4f} seconds'.format(toc - tic))
    print('ECmean4 Performance Indices succesfully computed!')


def pi_entry_point():
    """
    Command line interface to run the global_mean function
    """

    # read arguments from command line
    args = parse_arguments(sys.argv[1:], script='pi')

    performance_indices(exp=args.exp, year1=args.year1, year2=args.year2,
                        numproc=args.numproc,
                        loglevel=args.loglevel,
                        climatology=args.climatology,
                        interface=args.interface, config=args.config,
                        model=args.model, ensemble=args.ensemble,
                        outputdir=args.outputdir)

