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
import copy
from time import time
from multiprocessing import Process, Manager
import numpy as np
from tabulate import tabulate
from ecmean import var_is_there, load_yaml, make_input_filename, \
    get_levels, units_extra_definition, units_are_integrals, \
    units_converter, directions_match, chunks, \
    Diagnostic, getinifiles, getdomain
from cdopipe import CdoPipe


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
    parser.add_argument('-e', '--ensemble', type=str, default='r1i1p1f1',
                        help='variant label (ripf number for cmor)')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='activate cdo debugging')
    parser.add_argument('-i', '--interface', type=str, default='',
                        help='interface (overrides config.yml)')
    return parser.parse_args(args)


def worker(cdopin, piclim, face, diag, field_3d, varstat, varlist):
    """Main parallel diagnostic worker"""

    cdop = copy.copy(cdopin)  # Create a new local instance

    for var in varlist:

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
            dataref = piclim[var]['dataset']
            dataname = piclim[var]['dataname']

            # get files for climatology
            clim = str(diag.CLMDIR / f'climate_{dataref}_{dataname}.nc')
            vvvv = str(diag.CLMDIR / f'variance_{dataref}_{dataname}.nc')

            # create a file list using bash wildcards
            infile = make_input_filename(var, dervars, diag.years_joined, '????', face, diag)

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
            if 'derived' in face['variables'][var].keys():
                cmd = face['variables'][var]['derived']
#                dervars = (",".join(re.findall("[a-zA-Z]+", cmd)))
                cdop.selectname(",".join(dervars))
                cdop.expr(var, cmd)
            else:
                cdop.selectname(var)

            # fix grids and set domain making use component key from interface file
            domain = getdomain(var, face)
            cdop.fixgrid(domain=domain)
            cdop.timmean()

            # use convert() of cdopipe class to convert units
            cdop.convert(offset, factor)

            # temporarily using remapbil instead of remapcon due to NEMO grid missing corner
            # outfile = cdop.execute('remapbil', diag.resolution)
            if getdomain(var, face) in 'atm':
                outfile = cdop.execute('remap', diag.resolution, cdop.ATMWEIGHTS)
            elif getdomain(var, face) in 'oce' + 'ice':
                outfile = cdop.execute('remap', diag.resolution, cdop.OCEWEIGHTS)

            # special treatment which includes vertical interpolation
            if var in field_3d:

                # extract the vertical levels from the file
                format_vlevels = get_levels(clim)

                cdop.chain(f'intlevelx,{format_vlevels}')
                cdop.zonmean()
                # cdop.invertlat()
                cdop.sub(clim)
                cdop.sqr()
                cdop.div(vvvv)
                cdop.chain('vertmean -genlevelbounds,zbot=0,ztop=100000')

            else:
                # cdop.invertlat()
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
            if diag.fverb:
                print('PI for ', var, varstat[var])


def main(argv):
    """Main performance indices calculation"""

    assert sys.version_info >= (3, 7)

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

    # Init CdoPipe object to use in the following
    cdop = CdoPipe(debug=diag.debug)

    # loading the var-to-file interface
    face = load_yaml(INDIR / Path('interfaces', f'interface_{diag.interface}.yml'))

    # load the climatology reference data
    piclim = load_yaml('pi_climatology.yml')

    # new bunch of functions to set grids, create correction command, masks and areas
    comp = face['model']['component']  # Get component for each domain
    atminifile, ocegridfile, oceareafile = getinifiles(face, diag)
    cdop.set_gridfixes(atminifile, ocegridfile, oceareafile, comp['atm'], comp['oce'])
    cdop.make_atm_masks(comp['atm'], atminifile, extra=f'-remapcon2,{diag.resolution}')

    # create interpolation weights
    cdop.make_atm_remap_weights(atminifile, 'remapcon', diag.resolution)
    cdop.make_oce_remap_weights(ocegridfile, 'remapbil', diag.resolution)

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
        p = Process(target=worker,
                    args=(cdop, piclim, face, diag, field_3d, varstat, varlist))
        p.start()
        processes.append(p)

    # wait for the processes to finish
    for proc in processes:
        proc.join()

    toc = time()
    # evaluate tic-toc time  of execution
    if diag.fverb:
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
    tablefile = diag.TABDIR / f'PI4_RK08_{diag.expname}_{diag.modelname}_{diag.ensemble}_{diag.year1}_{diag.year2}.txt'
    if diag.fverb:
        print(tablefile)
    with open(tablefile, 'w', encoding='utf-8') as f:
        f.write(tabulate(global_table, headers=head, tablefmt='orgtbl'))
        f.write('\n\nPartial PI (atm only) is   : ' + str(round(partial_pi, 3)))
        f.write('\nTotal Performance Index is : ' + str(round(total_pi, 3)))

    # Make sure al temp files have been removed
    cdop.cdo.cleanTempDir()


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
