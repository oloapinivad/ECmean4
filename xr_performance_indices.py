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
from ecmean import  load_yaml, units_extra_definition, units_are_integrals, \
    units_converter, directions_match, chunks, \
    Diagnostic, getdomain, make_input_filename
import xarray as xr
from xrpipe import var_is_there, eval_formula, \
    xr_get_inifiles, \
    util_dictionary, remap_dictionary, guess_bounds


def parse_arguments(args):
    """Parse CLI arguments"""

    parser = argparse.ArgumentParser(
        description='ECmean Performance Indices for EC-Earth4')
    parser.add_argument('exp', metavar='EXP', type=str, help='experiment ID')
    parser.add_argument('year1', metavar='Y1', type=int, help='starting year')
    parser.add_argument('year2', metavar='Y2', type=int, help='final year')
    parser.add_argument('-s', '--silent', action='store_true',
                        help='do not print anything to std output')
    parser.add_argument('-v', '--loglevel', type=str, default='ERROR',
                        help='define the level of logging. default: error')
    parser.add_argument('-j', dest="numproc", type=int, default=1,
                        help='number of processors to use')
    parser.add_argument('-c', '--config', type=str, default='',
                        help='config file')
    parser.add_argument('-m', '--model', type=str, default='',
                        help='model name')
    parser.add_argument('-k', '--climatology', type=str, default='RK08',
                        help='climatology to be compared. default: RK08. Options: [RK08, EC22]')
    parser.add_argument('-r', '--resolution', type=str, default='',
                        help='climatology resolution')
    parser.add_argument('-e', '--ensemble', type=str, default='r1i1p1f1',
                        help='variant label (ripf number for cmor)')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='activate cdo debugging')
    parser.add_argument('-i', '--interface', type=str, default='',
                        help='interface (overrides config.yml)')
    return parser.parse_args(args)


def pi_worker(util, piclim, face, diag, field_3d, varstat, varlist):

    """Main parallel diagnostic worker"""


    for var in varlist:

        vdom =  getdomain(var, face)
        #weights = util[vdom + '_weights']

        if 'derived' in face['variables'][var].keys():
            cmd = face['variables'][var]['derived']
            dervars = re.findall("[a-zA-Z]+", cmd)
        else:
            dervars = [var]

        # check if required variables are there: use interface file
        # check into first file, and load also model variable units
        infile = make_input_filename(var, dervars, diag.year1, diag.year1, face, diag)
        isavail, varunit = var_is_there(infile, var, face['variables'])

        # if var is not available, store a NaN for the table
        if not isavail:
            varstat[var] = float('NaN')
        else:
            # unit conversion: from original data to data required by PI
            # using metpy avoid the definition of operations inside the dataset
            # use offset and factor separately (e.g. will not work with Fahrenait)
            # now in functions.py
            logging.debug(var)
            logging.debug(varunit + ' ---> ' + piclim[var]['units'])
            # adjust integrated quantities
            new_units = units_are_integrals(varunit, piclim[var])

            # unit conversion
            offset, factor = units_converter(new_units, piclim[var]['units'])

            # sign adjustment (for heat fluxes)
            factor = factor * directions_match(face['variables'][var], piclim[var])
            logging.debug('Offset %f, Factor %f', offset, factor)

            # extract info from pi_climatology.yml
            # reference dataset and reference varname
            # as well as years when available
            dataref = piclim[var]['dataset']
            dataname = piclim[var]['dataname']
            datayear1 = piclim[var].get('year1', 'nan')
            datayear2 = piclim[var].get('year2', 'nan')

            # get files for climatology
            if diag.climatology == 'RK08':
                clim = str(diag.RESCLMDIR / f'climate_{dataref}_{dataname}.nc')
                vvvv = str(diag.RESCLMDIR / f'variance_{dataref}_{dataname}.nc')
            elif diag.climatology == 'EC22':
                clim = str(diag.RESCLMDIR / f'climate_{dataname}_{dataref}_{diag.resolution}_{datayear1}-{datayear2}.nc')
                vvvv = str(diag.RESCLMDIR / f'variance_{dataname}_{dataref}_{diag.resolution}_{datayear1}-{datayear2}.nc')

            # create a file list using bash wildcards
            infile = make_input_filename(var, dervars, diag.years_joined, '????', face, diag)

            xfield = xr.open_mfdataset(infile)
            if vdom in 'oce' : 
                xfield = xfield.rename({"nav_lon": "lon", "nav_lat": "lat"})

            if 'derived' in face['variables'][var].keys():
                cmd = face['variables'][var]['derived']
                outfield = eval_formula(cmd, xfield)
            else:
                outfield = xfield[var]

            tmean = outfield.mean(dim='time_counter')
            tmean = tmean * factor + offset

            if vdom in 'atm':
                if util['atm_fix'] : 
                    fixed = util['atm_fix'](tmean, keep_attrs=True)
                    final = util['atm_remap'](fixed, keep_attrs=True)
            if vdom in 'oce':
                final = util['oce_remap'](tmean, keep_attrs=True)

            cfield = xr.open_dataset(clim)
            vfield = xr.open_dataset(vvvv)

            if var in field_3d:
                
                vname = list(cfield.data_vars)[-1]
                cfield = cfield[vname].metpy.convert_coordinate_units('lev', 'Pa')
                vfield = vfield[vname].metpy.convert_coordinate_units('lev', 'Pa')
                cfield = cfield.rename({"lev": "pressure_levels"})
                vfield = vfield.rename({"lev": "pressure_levels"})
                interped = final.interp(pressure_levels = cfield.pressure_levels.values)
                
                zonal = interped.mean(dim = 'lon')
                complete = (zonal - cfield)**2 / vfield
                bounds_lev = guess_bounds(complete['pressure_levels'], name = 'lev')
                bounds = abs(bounds_lev[:,0] - bounds_lev[:,1])
                weights = xr.DataArray(bounds, coords=[complete.pressure_levels], dims=['pressure_levels'])
                outarray = complete.weighted(weights).mean(dim = 'pressure_levels')

            else :
                cfield = xr.open_dataset(clim)
                vfield = xr.open_dataset(vvvv)
                if vdom in 'oce': 
                    cfield = cfield.rename({"LONGITUDE": "lon", "LATITUDE": "lat"})
                    vfield = vfield.rename({"LONGITUDE": "lon", "LATITUDE": "lat"})
                    
                cfield = cfield[list(cfield.data_vars)[-1]]
                vfield = vfield[list(vfield.data_vars)[-1]]

                outarray = (final - cfield)**2/vfield
               
            weights = np.cos(np.deg2rad(outarray.lat))
            out = outarray.weighted(weights).mean().values

                #outarray.to_netcdf(var + '_field.nc')             
                #outarray = mask(outarray, var, piclim[var]['mask'], util['atm_mask'] )

                
            
            # store the PI
            varstat[var] = float(out)
            if diag.fverb:
                print('PI for ', var, varstat[var])
            
def main(argv):

    """Main performance indices calculation"""

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

    # loading the var-to-file interface
    face = load_yaml(INDIR / Path('interfaces', f'interface_{diag.interface}.yml'))

    # load the climatology reference data
    piclim = load_yaml(diag.CLMDIR / f'pi_climatology_{diag.climatology}.yml')

    # new bunch of functions to set grids, create correction command, masks and areas
    comp = face['model']['component']  # Get component for each domain
    
    maskatmfile, atmareafile, oceareafile = xr_get_inifiles(face, diag)

    dataname = piclim['tas']['dataname']
    dataref = piclim['tas']['dataset']
    clim = str(diag.RESCLMDIR / f'climate_{dataref}_{dataname}.nc')

    #file = xr.open_dataset(clim)
    #print(area_cell(file).sum())

    #file = xr.open_dataset(atmareafile)
    #print(area_cell(file).sum())
  

    # create util dictionary including mask and weights for both atmosphere and ocean grids
    util = util_dictionary(comp, maskatmfile, clim, clim)
    remap = remap_dictionary(comp, atmareafile, oceareafile, 2)
    
    
    
    #print(util['atm_mask'])
    #print(remap['atm_fix'](util['atm_mask']))
    #util['atm_mask'] = remap['atm_remap'](remap['atm_fix'](util['atm_mask']))

    utildict = {**util, **remap}

    # add missing unit definitions
    units_extra_definition()

    # defines the two varlist
    field_2d = cfg['PI']['2d_vars']['field']
    field_3d = cfg['PI']['3d_vars']['field']
    field_oce = cfg['PI']['oce_vars']['field']
    field_ice = cfg['PI']['ice_vars']['field']
    field_all = field_2d + field_3d + field_oce + field_ice

    # trick to avoid the loop on years
    # define required years with a {year1,year2} and then use cdo select feature
    # years_list = [str(element) for element in range(diag.year1, diag.year2+1)]
    # diag.years_joined = ','.join(years_list)
    # special treatment to exploit bash wild cards on multiple years
    # if len(years_list) > 1:
    #    diag.years_joined = '{' + diag.years_joined + '}'

    # We now use a list
    diag.years_joined = list(range(diag.year1, diag.year2+1))

    # main loop: manager is required for shared variables
    mgr = Manager()

    # dictionaries are shared, so they have to be passed as functions
    varstat = mgr.dict()
    processes = []
    tic = time()

    # loop on the variables, create the parallel process
    for varlist in chunks(field_all, diag.numproc):
        p = Process(target=pi_worker,
                    args=(utildict, piclim, face, diag, field_3d, varstat, varlist))
        p.start()
        processes.append(p)

    # wait for the processes to finish
    for proc in processes:
        proc.join()

    toc = time()
    # evaluate tic-toc time  of execution
    if diag.fverb:
        print('Done in {:.4f} seconds'.format(toc-tic))

    # # define options for the output table
    head = ['Var', 'PI', 'Domain', 'Dataset', 'CMIP3', 'Ratio to CMIP3']
    global_table = []

    # loop on the variables
    for var in field_all:
        out_sequence = [var, varstat[var], piclim[var]['mask'], piclim[var]
                        ['dataset'], piclim[var]['cmip3'], varstat[var]/float(piclim[var]['cmip3'])]
        global_table.append(out_sequence)

    # nice loop on dictionary to get the partial and total pi
    partial_pi = np.nanmean([varstat[k] for k in field_2d + field_3d])
    total_pi = np.nanmean([varstat[k] for k in field_2d + field_3d + field_oce + field_ice])

    # write the file  with tabulate: cool python feature
    tablefile = diag.TABDIR / \
        f'xr_PI4_{diag.climatology}_{diag.expname}_{diag.modelname}_{diag.ensemble}_{diag.year1}_{diag.year2}.txt'
    if diag.fverb:
        print(tablefile)
    with open(tablefile, 'w', encoding='utf-8') as f:
        f.write(tabulate(global_table, headers=head, tablefmt='orgtbl', floatfmt=".2f"))
        f.write('\n\nPartial PI (atm only) is   : ' + str(round(partial_pi, 3)))
        f.write('\nTotal Performance Index is : ' + str(round(total_pi, 3)))


if __name__ == '__main__':

    sys.exit(main(sys.argv[1:]))