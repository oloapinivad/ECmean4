#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
 python3 version of ECmean computes corraltion pattern (CP)
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
from tabulate import tabulate
from ecmean import var_is_there, load_yaml, make_input_filename, \
                      units_extra_definition, units_are_integrals, \
                      units_converter, directions_match, chunks, \
                      Diagnostic, getinifiles, getdomain
from cdopipe import CdoPipe


def worker(cdopin, piclim, face, diag, varstat, varlist):
    """Main parallel diagnostic worker"""

    cdop = copy.copy(cdopin)  # Create a new local instance

    for var in varlist:
        # check if required variables are there: use interface file
        # check into first file, and load also model variable units
        infile = make_input_filename(var, diag.year1, diag.year1, face, diag)
        isavail, varunit = var_is_there(infile, var, face['variables'])

        # if var is not available, store a NaN for the table
        if not isavail:
            varstat[var] = float('NaN')
        else:
            # unit conversion: from original data to data required by CP
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
            logging.debug(offset, factor)

            # extract info from pi_climatology.yml
            # reference dataset and reference varname
            dataref = piclim[var]['dataset']
            dataname = piclim[var]['dataname']

            # get files for climatology
            clim = str(diag.CLMDIR / f'climate_{dataref}_{dataname}.nc')

            # create a file list using bash wildcards
            infile = make_input_filename(var, diag.years_joined, '????', face, diag)

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
                dervars = (",".join(re.findall("[a-zA-Z]+", cmd)))
                cdop.selectname(dervars)
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
            outfile = cdop.execute('remapbil', diag.resolution)

            # avoid grid mismatch
            cdop.invertlat()

            # apply same mask as in climatology
            cdop.mask(piclim[var]['mask'])

            # compute spatial correlation
            cdop.fldcor(clim)

            # execute command
            x = cdop.execute('output', input=outfile)[0]

            # store the PI
            varstat[var] = float(x)
            if diag.fverb:
                print('Pattern Correlation for ', var, varstat[var])


def main(args):
    """Main pattern correlation calculation"""

    assert sys.version_info >= (3, 7)

    # config file (looks for it in the same dir as the .py program file
    INDIR = Path(os.path.dirname(os.path.abspath(__file__)))
    cfg = load_yaml(INDIR / 'config.yml')

    # Setup all common variables, directories from arguments and config files
    diag = Diagnostic(args, cfg)

    # Create missing folders
    os.makedirs(diag.TABDIR, exist_ok=True)

    # Init CdoPipe object to use in the following
    # cdop = CdoPipe(debug=True)
    cdop = CdoPipe()

    # loading the var-to-file interface
    face = load_yaml(INDIR / Path('interfaces', f'interface_{diag.modelname}.yml'))

    # load the climatology reference data
    piclim = load_yaml('pi_climatology.yml')

    # new bunch of functions to set grids, create correction command, masks and areas
    comp = face['model']['component']  # Get component for each domain
    atminifile, oceinifile = getinifiles(face, diag)
    cdop.set_gridfixes(atminifile, oceinifile, comp['atm'], comp['oce'])
    cdop.make_atm_masks(atminifile, extra=f'-invertlat -remapcon2,{diag.resolution}')

    # add missing unit definitions
    units_extra_definition()

    # defines the two varlist
    atm_vars = cfg['correlation']['atm_vars']
    all_vars = atm_vars

    # We now use a list
    diag.years_joined = list(range(diag.year1, diag.year2+1))

    # main loop: manager is required for shared variables
    mgr = Manager()

    # dictionaries are shared, so they have to be passed as functions
    varstat = mgr.dict()
    processes = []
    tic = time()

    # loop on the variables, create the parallel process
    for varlist in chunks(all_vars, diag.numproc):
        p = Process(target=worker,
                    args=(cdop, piclim, face, diag, varstat, varlist))
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
    head = ['Var', 'Correlation', 'Domain', 'Dataset']
    global_table = []

    # loop on the variables
    for var in all_vars:
        out_sequence = [var, varstat[var], piclim[var]['mask'], piclim[var]['dataset']]
        global_table.append(out_sequence)

    # nice loop on dictionary to get the partial and total pi
    #partial_pi = np.mean([varstat[k] for k in field_2d + field_3d])
    #total_pi = np.mean([varstat[k] for k in field_2d + field_3d + field_oce + field_ice])

    # write the file  with tabulate: cool python feature
    tablefile = diag.TABDIR / f'CP4_RK08_{diag.expname}_{diag.year1}_{diag.year2}.txt'
    if diag.fverb:
        print(tablefile)
    with open(tablefile, 'w', encoding='utf-8') as f:
        f.write(tabulate(global_table, headers=head, tablefmt='orgtbl'))
        #f.write('\n\nPartial PI (atm only) is   : ' + str(partial_pi))
        #f.write('\nTotal Performance Index is : ' + str(total_pi))

    # Make sure al temp files have been removed
    cdop.cdo.cleanTempDir()


if __name__ == '__main__':

    # arguments
    parser = argparse.ArgumentParser(
        description='ECmean Pattern Correlation for EC-Earth4')
    parser.add_argument('exp', metavar='EXP', type=str, help='experiment ID')
    parser.add_argument('year1', metavar='Y1', type=int, help='starting year')
    parser.add_argument('year2', metavar='Y2', type=int, help='final year')
    parser.add_argument('-s', '--silent', action='store_true',
                        help='do not print anything to std output')
    parser.add_argument('-v', '--loglevel', type=str, default='ERROR',
                        help='define the level of logging. default: error')
    parser.add_argument('-j', dest="numproc", type=int, default=1,
                        help='number of processors to use')
    parser.add_argument('-m', '--model', type=str, default='',
                        help='model name')
    args = parser.parse_args()

    # log level with logging
    # currently basic definition trought the text
    loglevel = args.loglevel.upper()
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level)

    main(args)
