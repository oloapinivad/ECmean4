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
import numpy as np
from time import time
from tabulate import tabulate
from functions import vars_are_there, load_yaml, make_input_filename, \
                      get_levels, units_extra_definition, units_are_integrals, \
                      units_converter, directions_match, chunks
from cdopipe import CdoPipe
import copy
from multiprocessing import Process, Manager


def worker(cdopin, piclim, face, ECEDIR, CLMDIR, resolution, field_3d, year1, years_joined,
           expname, fverb, varstat, varlist):

    cdop = copy.copy(cdopin)  # Create a new local instance
    for var in varlist:

        # check if required variables are there: use interface file
        # check into first file, and load also model variable units
        infile = make_input_filename(ECEDIR, var, expname, year1, year1, face)
        isavail, varunit = vars_are_there(infile, [var], face)
        #varunit = {**varunit, **retunit}

        # if var is not available, store a NaN for the table
        if not isavail[var]:
            varstat[var] = float('NaN')
        else:
            # unit conversion: from original data to data required by PI
            # using metpy avoid the definition of operations inside the dataset
            # use offset and factor separately (e.g. will not work with Fahrenait)
            # now in functions.py
            logging.debug(var)
            logging.debug(varunit[var] + ' ---> ' + piclim[var]['units'])
            # adjust integrated quantities
            new_units = units_are_integrals(varunit[var], piclim[var])

            # unit conversion
            units_conversion = units_converter(new_units, piclim[var]['units'])

            # sign adjustment (for heat fluxes)
            units_conversion['factor'] = units_conversion['factor'] * \
                                         directions_match(face[var], piclim[var])
            logging.debug(units_conversion)

            # extract info from pi_climatology.yml
            # reference dataset and reference varname
            dataref = piclim[var]['dataset']
            dataname = piclim[var]['dataname']

            # get files for climatology
            clim = str(CLMDIR / f'climate_{dataref}_{dataname}.nc')
            vvvv = str(CLMDIR / f'variance_{dataref}_{dataname}.nc')

            # create a file list using bash wildcards
            infile = make_input_filename(
                ECEDIR, var, expname, years_joined, '????', face)

            # Start fresh pipe
            # This leaves the input file undefined for now. It can be set later with
            # cdop.set_infile(infile) or by specifying input=infile cdop.execute
            cdop.start()

            # set input file
            cdop.set_infile(infile)

            # check if var is derived
            # if this is the case, get the derived expression and select
            # the set of variables you need
            # otherwise, use only select (this avoid loop)
            # WARNING: it may scale badly with high-resolution centennial runs
            if 'derived' in face[var].keys():
                cmd = face[var]['derived']
                dervars = (",".join(re.findall("[a-zA-Z]+", cmd)))
                cdop.selectname(dervars)
                cdop.expr(var, cmd)
            else:
                cdop.selectname(var)

            # fix grids and set domain making use component key from interface file
            cdop.fixgrid(domain=face[var]['component'])
            cdop.timmean()

            # use convert() of cdopipe class to convert units
            cdop.convert(units_conversion['offset'], units_conversion['factor'])

            # temporarily using remapbil instead of remapcon due to NEMO grid missing corner
            outfile = cdop.execute('remapbil', resolution)

            # special treatment which includes vertical interpolation
            if var in field_3d:

                # extract the vertical levels from the file
                format_vlevels = get_levels(clim)

                cdop.chain(f'intlevelx,{format_vlevels}')
                cdop.zonmean()
                cdop.invertlat()
                cdop.sub(clim)
                cdop.sqr()
                cdop.div(vvvv)
                cdop.chain('vertmean -genlevelbounds,zbot=0,ztop=100000')

            else:
                cdop.invertlat()
                cdop.sub(clim)
                cdop.sqr()
                cdop.div(vvvv)
                # apply same mask as in climatology
                cdop.mask(piclim[var]['mask'])

            # execute command
            x = np.squeeze(cdop.execute('fldmean', input=outfile,
                                        returnCdf=True).variables[var])

            # store the PI
            varstat[var] = float(x)
            if fverb:
                print('PI for ', var, varstat[var])


def main(args):
    """Main performance indices calculation"""

    assert sys.version_info >= (3, 7)

    expname = args.exp
    year1 = args.year1
    year2 = args.year2
    fverb = not args.silent
    numproc = args.numproc

    # config file (looks for it in the same dir as the .py program file
    INDIR = Path(os.path.dirname(os.path.abspath(__file__)))

    # load - if exists - config file
    cfg = load_yaml(INDIR / 'config.yml')

    # hard-coded resolution (due to climatological dataset)
    resolution = cfg['PI']['resolution']

    # folder definition
    ECEDIR = Path(os.path.expandvars(cfg['dirs']['exp']), expname)
    TABDIR = Path(os.path.expandvars(cfg['dirs']['tab']))
    CLMDIR = Path(os.path.expandvars(cfg['dirs']['clm']), resolution)
    os.makedirs(TABDIR, exist_ok=True)

    # prepare grid description file
    ATMINIFILE = ECEDIR / f'ICMGG{expname}INIT'
    OCEINIFILE = cfg['areas']['oce']

    # Init CdoPipe object to use in the following
    # cdop = CdoPipe(debug=True)
    cdop = CdoPipe()

    # new bunch of functions to set grids, create correction command, masks and areas
    cdop.set_gridfixes(ATMINIFILE, OCEINIFILE, 'oifs', 'nemo')
    cdop.make_atm_masks(ATMINIFILE, extra=f'-invertlat -remapcon2,{resolution}')

    # add missing unit definitions
    units_extra_definition()

    # trick to avoid the loop on years
    # define required years with a {year1,year2} and then use cdo select feature
    years_list = [str(element) for element in range(year1, year2+1)]
    years_joined = ','.join(years_list)

    # special treatment to exploit bash wild cards on multiple years
    if len(years_list) > 1:
        years_joined = '{' + years_joined + '}'

    if fverb:
        print(years_joined)

    # loading the var-to-file interface
    face = load_yaml(INDIR / 'interface_ece4.yml')

    # reference data: it is badly written but it can be implemented in a much more intelligent
    # and modular way
    piclim = load_yaml('pi_climatology.yml')

    # defines the two varlist
    field_2d = cfg['PI']['2d_vars']['field']
    field_3d = cfg['PI']['3d_vars']['field']
    field_oce = cfg['PI']['oce_vars']['field']
    field_ice = cfg['PI']['ice_vars']['field']
    field_all = field_2d + field_3d + field_oce + field_ice

    # main loop
    mgr = Manager()
    varstat = mgr.dict()
    processes = []
    tic = time()

    for varlist in chunks(field_all, numproc):
        p = Process(target = worker, args=(cdop, piclim, face, ECEDIR, CLMDIR, resolution, field_3d,
                                           year1, years_joined, expname, fverb,
                                           varstat, varlist))
        p.start()
        processes.append(p)

    for proc in processes:
        proc.join()

    toc = time()
    if fverb:
        print('Done in {:.4f} seconds'.format(toc-tic))

    # define options for the output table
    head = ['Var', 'PI', 'Domain', 'Dataset', 'CMIP3', 'Ratio to CMIP3']
    global_table = []

    # loop on the variables
    for var in field_all:
        out_sequence = [var, varstat[var], piclim[var]['mask'], piclim[var]
                        ['dataset'], piclim[var]['cmip3'], varstat[var]/piclim[var]['cmip3']]
        global_table.append(out_sequence)

    # nice loop on dictionary to get the partial and total pi
    partial_pi = np.mean([varstat[k] for k in field_2d + field_3d])
    total_pi = np.mean([varstat[k] for k in field_2d + field_3d + field_oce + field_ice])

    # write the file  with tabulate: cool python feature
    tablefile = TABDIR / f'PI4_RK08_{expname}_{year1}_{year2}.txt'
    if fverb:
        print(tablefile)
    with open(tablefile, 'w', encoding='utf-8') as f:
        f.write(tabulate(global_table, headers=head, tablefmt='orgtbl'))
        f.write('\n\nPartial PI (atm only) is   : ' + str(partial_pi))
        f.write('\nTotal Performance Index is : ' + str(total_pi))

    # Make sure al temp files have been removed
    cdop.cdo.cleanTempDir()


if __name__ == '__main__':

    # arguments
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
    args = parser.parse_args()

    # log level with logging
    # currently basic definition trought the text
    loglevel = args.loglevel.upper()
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level)

    main(args)
