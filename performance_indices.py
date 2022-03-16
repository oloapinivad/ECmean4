#!/usr/bin/env python3

# This is a tentative python script to convert ECmean global mean operation to python3
# It uses a reference file from yaml and cdo bindings (not very efficient)

import sys
import os
import yaml
import numpy as np
import argparse
from tabulate import tabulate
from cdo import *
from pathlib import Path
from pint import UnitRegistry
from metpy.units import units

# import functions from functions.py
import functions as fn


def main(args):

    assert sys.version_info >= (3, 5)

    cdo = Cdo()

    # special units definition, need to be moved in another placce
    units.define('fraction = [] = frac')
    units.define('percent = 1e-2 frac = pct')
    units.define('psu = 1e-3 frac')
    units.define('ppm = 1e-6 fraction')

    expname = args.exp
    year1 = args.year1
    year2 = args.year2
    fverb = not args.silent

    # config file (looks for it in the same dir as the .py program file
    INDIR = Path(os.path.dirname(os.path.abspath(__file__)))

    # load - if exists - config file
    cfg = fn.load_config_file(INDIR)

    # hard-coded resolution (due to climatological dataset)
    resolution = cfg['PI']['resolution']

    # folder definition
    ECEDIR = Path(os.path.expandvars(cfg['dirs']['exp']), expname)
    TABDIR = Path(os.path.expandvars(cfg['dirs']['tab']), 'ECmean4', 'table')
    CLMDIR = Path(os.path.expandvars(cfg['dirs']['clm']), resolution)
    TMPDIR = Path(os.path.expandvars(cfg['dirs']['tmp']))
    os.makedirs(TABDIR, exist_ok=True)

    #cdo.forceOutput = True
    #cdo.debug = True

    # prepare grid description file
    icmgg_file = ECEDIR / f'ICMGG{expname}INIT'
    gridfile = str(TMPDIR / f'grid_{expname}.txt')
    griddes = cdo.griddes(input=f'{icmgg_file}')
    with open(gridfile, 'w') as f:
        for line in griddes:
            print(line, file=f)

    # land-sea masks
    ocean_mask = cdo.setctomiss(0,
                                input=f'-ltc,0.5 -invertlat -remapcon2,{resolution} '
                                f'-setgridtype,regular -setgrid,{gridfile} '
                                f'-selcode,172 {icmgg_file}', options='-f nc')
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
    filename = 'pi_climatology.yml'
    with open(INDIR / filename, 'r') as file:
        piclim = yaml.load(file, Loader=yaml.FullLoader)

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
        infile = fn.make_input_filename(
            ECEDIR, field, expname, year1, year1, face)
        retavail, retunit = fn.vars_are_there(infile, [field], face)
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
            piunits = piclim[var]['units']

            if units(varunit[var]) != units(piunits) : 
                print(var+'... Unit converson required')
                offset_standard = 0 * units(varunit[var])
                factor_standard = 1 * units(varunit[var])
                try: 
                    offset = offset_standard.to(piunits).magnitude
                    if offset != 0:
                        print(offset)
                        oper = f'-addc,{offset}'

                    # for all the others
                    else :
                        factor = factor_standard.to(piunits).magnitude
                        print(factor)
                        oper = f'-mulc,{factor}'

                except : 
                    print("Assuming this as a precipitation field! Am I correct?")
                    density_water = units('kg / m^3') * 1000
                    factor = (factor_standard/density_water).to(piunits).magnitude
                    oper = f'-mulc,{factor}'
            else:
                oper = ''

            # extract info from pi_climatology.yml
            # reference dataset and reference varname
            dataref = piclim[var]['dataset']
            dataname = piclim[var]['dataname']

            # check if var is derived
            # if this is the case, get the derived expression and select
            # the set of variables you need
            # otherwise, use only select (this avoid loop)
            # WARNING: it may scale badly with high-resolution centennial runs
            if 'derived' in face[var].keys():
                cmd = face[var]['derived']
                dervars = (",".join(re.findall("[a-zA-Z]+", cmd)))
                cmd_select = f'-expr,{var}={cmd} -select,name={dervars}'
            else:
                cmd_select = f'-select,name={var}'

            # field atm cubic grid defition to be handled by cdo
            if var in field_2d + field_3d:
                cmd_grid = f'-setgridtype,regular -setgrid,{gridfile}'
            else:
                cmd_grid = ''

            # create a file list using bash wildcards
            infile = fn.make_input_filename(
                ECEDIR, var, expname, years_joined, '????', face)

            # get files for climatology
            clim = str(CLMDIR / f'climate_{dataref}_{dataname}.nc')
            vvvv = str(CLMDIR / f'variance_{dataref}_{dataname}.nc')

            # apply masks when needed (PI-dependent feature)
            if piclim[var]['domain'] == 'land':
                mask = f'-mul {land_mask}'
            elif piclim[var]['domain'] == 'ocean':
                mask = f'-mul {ocean_mask}'
            elif piclim[var]['domain'] == 'global':
                mask = ''

            # timmean and remap
            cmd1 = f'-timmean {cmd_grid} {cmd_select} {infile}'

            # temporarily using remapbil instead of remapcon due to NEMO grid missing corner
            outfile = cdo.remapbil(resolution, input=cmd1)

            # special treatment which includes vertical interpolation
            if var in field_3d:

                # extract the vertical levels from the file
                vlevels = cdo.showlevel(input=f'{clim}')

                # perform multiple string manipulation to produce a Pa list of levels
                v0 = ''.join(vlevels).split()
                v1 = [int(x) for x in v0]

                # if the grid is hPa, move to Pa (there might be a better solution)
                if np.max(v1) < 10000:
                    v1 = [x * 100 for x in v1]

                # format for CDO, converting to string
                format_vlevels = ' '.join(str(x) for x in v1).replace(' ', ',')

                # assign the vertical command for interpolation and zonal mean
                cmd_vertinterp = f'-zonmean -intlevelx,{format_vlevels}'
                cmd_vertmean = f'-vertmean -genlevelbounds,zbot=0,ztop=100000'
            else:
                cmd_vertmean = ''
                cmd_vertinterp = ''

            # defining the PI computation: RMS normalized by dataset interannual variance
            cmd2 = f'-setname,{var} {mask} {cmd_vertmean} -div -sqr -sub -invertlat {cmd_vertinterp} {oper} {outfile} {clim} {vvvv}'
        
            # calling PI computation
            x = np.squeeze(cdo.fldmean(input=cmd2, returnCdf=True).variables[var])

            # store the PI
            varstat[var] = float(x)
            if fverb:
                print('PI for ', var, varstat[var])

    # define options for the output table
    head = ['Var', 'PI', 'Domain', 'Dataset', 'CMIP3', 'Ratio to CMIP3']
    global_table = list()

    # loop on the variables
    for var in field_all:
        out_sequence = [var, varstat[var], piclim[var]['domain'], piclim[var]
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

    # cleaning
    os.unlink(gridfile)
    cdo.cleanTempDir()


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
