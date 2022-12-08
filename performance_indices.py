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
import re
import argparse
from pathlib import Path
import logging
from time import time
from multiprocessing import Process, Manager
import numpy as np
from tabulate import tabulate
import xarray as xr
import pandas as pd
import yaml
from collections import defaultdict
from ecmean import var_is_there, eval_formula, \
    get_inifiles, adjust_clim_file, get_clim_files, \
    areas_dictionary, masks_dictionary, remap_dictionary, guess_bounds, \
    load_yaml, units_extra_definition, units_are_integrals, \
    units_converter, directions_match, chunks, \
    Diagnostic, getdomain, make_input_filename, xr_preproc, mask_field, \
    heatmap_comparison, select_region

# temporary disabling the scheduler
import dask
dask.config.set(scheduler="synchronous")


def parse_arguments(args):
    """Parse CLI arguments"""

    parser = argparse.ArgumentParser(
        description='ECmean Performance Indices for EC-Earth4')
    parser.add_argument('exp', metavar='EXP', type=str, help='experiment ID')
    parser.add_argument('year1', metavar='Y1', type=int, help='starting year')
    parser.add_argument('year2', metavar='Y2', type=int, help='final year')
    parser.add_argument('-s', '--silent', action='store_true',
                        help='do not print anything to std output')
    parser.add_argument('-v', '--loglevel', type=str, default='WARNING',
                        help='define the level of logging. default: warning')
    parser.add_argument('-j', dest="numproc", type=int, default=1,
                        help='number of processors to use')
    parser.add_argument('-c', '--config', type=str, default='',
                        help='config file')
    parser.add_argument('-m', '--model', type=str, default='',
                        help='model name')
    parser.add_argument(
        '-k',
        '--climatology',
        type=str,
        default='RK08',
        help='climatology to be compared. default: RK08. Options: [RK08, EC22]')
    parser.add_argument('-r', '--resolution', type=str, default='',
                        help='climatology resolution')
    parser.add_argument('-e', '--ensemble', type=str, default='r1i1p1f1',
                        help='variant label (ripf number for cmor)')
    parser.add_argument('-i', '--interface', type=str, default='',
                        help='interface (overrides config.yml)')
    return parser.parse_args(args)


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

    for var in varlist:

        vdom = getdomain(var, face)

        if 'derived' in face['variables'][var].keys():
            cmd = face['variables'][var]['derived']
            dervars = re.findall("[a-zA-Z]+", cmd)
        else:
            dervars = [var]

        # check if required variables are there: use interface file
        # check into first file, and load also model variable units
        infile = make_input_filename(
            var, dervars, diag.year1, diag.year1, face, diag)
        isavail, varunit = var_is_there(infile, var, face['variables'])

        # store NaN in a default dict
        result = defaultdict(lambda: defaultdict(lambda : float('NaN')))

        # if the variable is available
        if isavail : 

            logging.debug(var)
            logging.debug(varunit + ' ---> ' + piclim[var]['units'])

            # adjust integrated quantities
            new_units = units_are_integrals(varunit, piclim[var])

            # unit conversion based on metpy
            offset, factor = units_converter(new_units, piclim[var]['units'])

            # sign adjustment (for heat fluxes)
            factor = factor * \
                directions_match(face['variables'][var], piclim[var])
            logging.debug('Offset %f, Factor %f', offset, factor)

            # create a file list using bash wildcards
            infile = make_input_filename(
                var, dervars, diag.years_joined, '????', face, diag)


            # open file: chunking on time only, might be improved
            xfield = xr.open_mfdataset(infile, preprocess=xr_preproc, chunks={'time': 12})

            # in case of big files with multi year, be sure of having opened the right records
            xfield = xfield.sel(time=xfield.time.dt.year.isin(diag.years_joined))

            if 'derived' in face['variables'][var].keys():
                cmd = face['variables'][var]['derived']
                outfield = eval_formula(cmd, xfield)
            else:
                outfield = xfield[var]

            # mean over time and fixing of the units
            for season in diag.seasons : 

                # get filenames for climatology
                clim, vvvv = get_clim_files(piclim, var, diag, season)

                logging.info(season)

                # open climatology files, fix their metadata
                cfield = adjust_clim_file(xr.open_dataset(clim))
                vfield = adjust_clim_file(xr.open_dataset(vvvv), remove_zero=True)
                tmean = outfield.copy(deep=True)

                # season selection
                if season != 'ALL' :
                    tmean = tmean.sel(time=tmean.time.dt.season.isin(season))
                    cfield = cfield.sel(time=cfield.time.dt.season.isin(season))
                    vfield = vfield.sel(time=vfield.time.dt.season.isin(season))

                # averaging 
                tmean = tmean.mean(dim='time')

                # safe check for old RK08 which has a different format
                if diag.climatology != 'RK08' :
                    cfield = cfield.mean(dim='time')
                    vfield = vfield.mean(dim='time')

                tmean = tmean * factor + offset

                # apply interpolation, if fixer is availble and with different
                # grids
                if vdom in 'atm':
                    if util['atm_fix']:
                        tmean = util['atm_fix'](tmean, keep_attrs=True)
                    final = util['atm_remap'](tmean, keep_attrs=True)
                if vdom in 'oce':
                    if util['oce_fix']:
                        tmean = util['oce_fix'](tmean, keep_attrs=True)
                    final = util['oce_remap'](tmean, keep_attrs=True)

                # vertical interpolation
                if var in field_3d:

                    # xarray interpolation on plev, forcing to be in Pascal
                    final = final.metpy.convert_coordinate_units('plev', 'Pa')
                    interped = final.interp(plev=cfield['plev'].values)
                    final = interped.mean(dim='lon')

                # compute PI
                complete = (final - cfield)**2 / vfield

                # vertical averageing
                if var in field_3d:
                    # compute vertical bounds as weights
                    bounds_lev = guess_bounds(complete['plev'], name='plev')
                    bounds = abs(bounds_lev[:, 0] - bounds_lev[:, 1])
                    weights = xr.DataArray(
                        bounds, coords=[
                            complete['plev']], dims=['plev'])

                    # vertical mean
                    outarray = complete.weighted(weights).mean(dim='plev')

                # horizontal averaging with land-sea mask
                else:
                    outarray = mask_field(
                        complete, var, piclim[var]['mask'], util['atm_mask'])

                # run the computation before doing each region
                computed = outarray.load()

                # loop on different regions
                for region in diag.regions : 

                    slicearray = select_region(computed, region)

                    # latitude-based averaging
                    weights = np.cos(np.deg2rad(slicearray.lat))
                    out = slicearray.weighted(weights).mean().values

                    # store the PI
                    result[season][region] = float(out)

                    # diagnostic
                    if diag.fverb and region == 'Global':
                        print('PI for', region, season, var, result[season][region])
                

        # nested dictionary, to be redifend as a dict to remove lambdas
        new_result = {k: dict(v) for k, v in result.items()}
        varstat[var] = new_result

def pi_main(argv):
    """Main performance indices calculation"""

    # assert sys.version_info >= (3, 7)

    tic = time()
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
    os.makedirs(diag.FIGDIR, exist_ok=True)

    # loading the var-to-file interface
    face = load_yaml(
        INDIR /
        Path(
            'interfaces',
            f'interface_{diag.interface}.yml'))

    # load the climatology reference data
    piclim = load_yaml(diag.CLMDIR / f'pi_climatology_{diag.climatology}.yml')

    # new bunch of functions to set grids, create correction command, masks
    # and areas
    comp = face['model']['component']  # Get component for each domain

    # all clim have the same grid, read from the first clim available and get
    # target grid
    clim, _ = get_clim_files(piclim, 'tas', diag, 'ALL')
    target_remap_grid = xr.open_dataset(clim)

    # get file info files
    maskatmfile, atmareafile, oceareafile = get_inifiles(face, diag)

    # create remap dictionary with atm and oce interpolators
    remap = remap_dictionary(comp, atmareafile, oceareafile, target_remap_grid)

    # create util dictionary including mask and weights for both
    # atmosphere and ocean grids
    # use the atmospheric remap dictionary to remap the mask file
    areas = areas_dictionary(comp, atmareafile, oceareafile)
    masks = masks_dictionary(comp, maskatmfile, remap_dictionary=remap)

    # join the two dictionaries
    util_dictionary = {**masks, **areas, **remap}

    # add missing unit definitions
    units_extra_definition()

    # defines the two varlist
    field_2d = cfg['PI']['2d_vars']['field']
    field_3d = cfg['PI']['3d_vars']['field']
    field_oce = cfg['PI']['oce_vars']['field']
    field_ice = cfg['PI']['ice_vars']['field']
    field_all = field_2d + field_3d + field_oce + field_ice

    # We now use a list
    diag.years_joined = list(range(diag.year1, diag.year2 + 1))

    # main loop: manager is required for shared variables
    mgr = Manager()

    # dictionaries are shared, so they have to be passed as functions
    varstat = mgr.dict()
    processes = []

    toc = time()
    if diag.fverb:
        print('Preproc in {:.4f} seconds'.format(toc - tic))
    tic = time()

    # loop on the variables, create the parallel process
    for varlist in chunks(field_all, diag.numproc):
        p = Process(
            target=pi_worker,
            args=(
                util_dictionary,
                piclim,
                face,
                diag,
                field_3d,
                varstat,
                varlist))
        p.start()
        processes.append(p)

    # wait for the processes to finish
    for proc in processes:
        proc.join()

    toc = time()
    # evaluate tic-toc time  of execution
    if diag.fverb:
        print('Done in {:.4f} seconds'.format(toc - tic) + ' with ' + str(diag.numproc) + ' processors')

    # # define options for the output table
    head = ['Variable', 'Domain', 'Dataset'] + diag.regions + [s + ' CMIP6 Ratio' for s in diag.regions]
    global_table = []

    # loop on the variables
    for var in field_all:
        out_sequence = [
            var,
            piclim[var]['mask'],
            piclim[var]['dataset']]

        for region in diag.regions : 
            out_sequence = out_sequence + [varstat[var]['ALL'][region]]
        for region in diag.regions : 
            #out_sequence = out_sequence + [float(varstat[var]['ALL'][region]) / float(piclim[var]['cmip6'][region])]
            out_sequence = out_sequence + [1]
           
        global_table.append(out_sequence)

    # nice loop on dictionary to get the partial and total pi
    partial_pi = np.nanmean([varstat[k]['ALL']['Global'] for k in field_2d + field_3d])
    total_pi = np.nanmean([varstat[k]['ALL']['Global']
                          for k in field_2d + field_3d + field_oce + field_ice])

    # dump the yaml file for PI, including all the seasons (need to copy to avoid mess)
    yamlfile = diag.TABDIR / \
        f'PI4_{diag.climatology}_{diag.expname}_{diag.modelname}_r1i1p1f1_{diag.year1}_{diag.year2}.yml'
    with open(yamlfile, 'w') as file:
        yaml.safe_dump(varstat.copy(), file, default_flow_style=False, sort_keys=False)

    # very clumsy conversion of the dictionary to a pd.dataframe: NEED TO BE IMPROVED
    data_table = {}
    for i in varstat.keys(): 
        pippo = {}
        for outerKey, innerDict in varstat[i].items():
            for innerKey, values in innerDict.items():
                pippo[(outerKey, innerKey)] = values
        data_table[i] = pippo
    data_table = pd.DataFrame(data_table).T

    # write the file with tabulate only for yearly mean
    tablefile = diag.TABDIR / \
        f'PI4_{diag.climatology}_{diag.expname}_{diag.modelname}_r1i1p1f1_{diag.year1}_{diag.year2}.txt'
    with open(tablefile, 'w', encoding='utf-8') as f:
        f.write(
            tabulate(
                global_table,
                headers=head,
                tablefmt='orgtbl',
                floatfmt=".2f"))
        f.write('\n\nPartial PI (atm only) is   : ' +
                str(round(partial_pi, 3)))
        f.write('\nTotal Performance Index is : ' + str(round(total_pi, 3)))
 
    # call the heatmap routine for a plot
    mapfile = diag.FIGDIR / \
        f'PI4_{diag.climatology}_{diag.expname}_{diag.modelname}_r1i1p1f1_{diag.year1}_{diag.year2}.pdf'
    heatmap_comparison(data_table, diag, mapfile)

if __name__ == '__main__':

    sys.exit(pi_main(sys.argv[1:]))
