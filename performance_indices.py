#!/usr/bin/env python3
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
import yaml
import numpy as np
from tabulate import tabulate
from cdo import *

from functions import *
from cdopipe import CdoPipe

cdo = Cdo()

def get_levels(infile):
    """Extract vertical levels from file,
       including Pa conversion"""

    # extract the vertical levels from the file
    vlevels = cdo.showlevel(input=infile)

    # perform multiple string manipulation to produce a Pa list of levels
    v0 = ''.join(vlevels).split()
    v1 = [int(x) for x in v0]

    # if the grid is hPa, move to Pa (there might be a better solution)
    if np.max(v1) < 10000:
        v1 = [x * 100 for x in v1]

    # format for CDO, converting to string
    return ' '.join(str(x) for x in v1).replace(' ', ',')


def main(args):
    """Main performance indices calculation"""

    assert sys.version_info >= (3, 7)

    expname = args.exp
    year1 = args.year1
    year2 = args.year2
    fverb = not args.silent

    # config file (looks for it in the same dir as the .py program file
    INDIR = Path(os.path.dirname(os.path.abspath(__file__)))

    # load - if exists - config file
    cfg = load_config_file(INDIR)

    # hard-coded resolution (due to climatological dataset)
    resolution = cfg['PI']['resolution']

    # folder definition
    ECEDIR = Path(os.path.expandvars(cfg['dirs']['exp']), expname)
    TABDIR = Path(os.path.expandvars(cfg['dirs']['tab']), 'ECmean4', 'table')
    CLMDIR = Path(os.path.expandvars(cfg['dirs']['clm']), resolution)
    TMPDIR = Path(os.path.expandvars(cfg['dirs']['tmp']))
    os.makedirs(TABDIR, exist_ok=True)

    #cdo.forceOutput = True

    # prepare grid description file
    INIFILE = ECEDIR / f'ICMGG{expname}INIT'
    OCEINIFILE=cfg['areas']['oce']

    # Init CdoPipe object to use in the following, specifying the LM and SM files
    #cdop = CdoPipe(debug=True)
    cdop = CdoPipe()
    cdop.make_grids(INIFILE, OCEINIFILE)

    # land-sea masks on regular 2x2 grid
    ocean_mask = cdo.setctomiss(0,
                                input=f'-ltc,0.5 -invertlat -remapcon2,{resolution} '
                                f'-setgridtype,regular -setgrid,{cdop.GRIDFILE} '
                                f'-selcode,172 {INIFILE}', options='-f nc')
    land_mask = cdo.addc(1, input=f'-setctomiss,1 -setmisstoc,0 {ocean_mask}')

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
    filename = 'interface_ece4.yml'
    with open(INDIR / filename, 'r') as file:
        face = yaml.load(file, Loader=yaml.FullLoader)

    # reference data: it is badly written but it can be implemented in a much more intelligent
    # and modular way
    piclim = load_yaml('pi_climatology.yml')

    # defines the two varlist
    field_2d = cfg['PI']['2d_vars']['field']
    field_3d = cfg['PI']['3d_vars']['field']
    field_oce = cfg['PI']['oce_vars']['field']
    field_ice = cfg['PI']['ice_vars']['field']
    field_all = field_2d + field_3d + field_oce + field_ice

    # check if required variables are there: use interface file
    # check into first file, and load also model variable units
    isavail = {}
    varunit = {}
    for field in field_all:
        infile = make_input_filename(
            ECEDIR, field, expname, year1, year1, face)
        retavail, retunit = vars_are_there(infile, [field], face)
        isavail = {**isavail, **retavail}
        varunit = {**varunit, **retunit}

    # main loop
    varstat = {}
    for var in field_all:

        # if var is not available, store a NaN for the table
        if not isavail[var]:
            varstat[var] = float('NaN')
        else:

            # unit conversion: from original data to data required by PI
            # using metpy avoid the definition of operations inside the dataset
            # use offset and factor separately (e.g. will not work with Fahrenait)
            # now in functions.py
            #print(var)
            #print(varunit[var] + ' ---> ' + piclim[var]['units'])
            # adjust integrated quantities
            new_units = units_are_integrals(varunit[var], piclim[var])

            # unit conversion
            units_conversion = units_converter(new_units, piclim[var]['units'])

            # sign adjustment (for heat fluxes)
            units_conversion['factor'] = units_conversion['factor'] * units_are_down(piclim[var])
            #print(units_conversion)
            #print('--------------')


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

             # set domain making use component key from interface file
            cdop.setdomain(face[var]['component'])

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

            cdop.fixgrid()
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

                # apply masks when needed
                # this is a hack for now
                if piclim[var]['mask'] == 'land':
                    cdop.mul(land_mask)
                elif piclim[var]['mask'] == 'ocean':
                    cdop.mul(ocean_mask)

                cdop.setname(var)

            # execute command
            x = np.squeeze(cdop.execute('fldmean', input=outfile,
                                        returnCdf=True).variables[var])

            # store the PI
            varstat[var] = float(x)
            if fverb:
                print('PI for ', var, varstat[var])

    # define options for the output table
    head = ['Var', 'PI', 'Domain', 'Dataset', 'CMIP3', 'Ratio to CMIP3']
    global_table = list()

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
    with open(tablefile, 'w') as f:
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
    args = parser.parse_args()

    main(args)
