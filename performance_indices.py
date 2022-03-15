#!/usr/bin/env python3

# This is a tentative python script to convert ECmean global mean operation to python3
# It uses a reference file from yaml and cdo bindings (not very efficient)

import sys
import os
import yaml
import numpy as np
import argparse
from tabulate import tabulate
from statistics import mean
from cdo import *
from pathlib import Path

# import functions from functions.py
import functions as fn
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

    assert sys.version_info >= (3, 5)

    expname = args.exp
    year1 = args.year1
    year2 = args.year2
    fverb = not args.silent

    # config file (looks for it in the same dir as the .py program file
    INDIR = Path(os.path.dirname(os.path.abspath(__file__)))
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

    # prepare grid description file
    INIFILE = ECEDIR / f'ICMGG{expname}INIT'
    OCEINIFILE=cfg['areas']['oce']
 
    # Init CdoPipe object to use in the following, specifying the LM and SM files
    cdop = CdoPipe()
    #cdop.debug=True
    cdop.make_grids(INIFILE, OCEINIFILE)

    # land-sea masks
    ocean_mask = cdo.setctomiss( 0,
        input=f'-ltc,0.5 -invertlat -remapcon2,{resolution} ' \
              f'-setgridtype,regular -setgrid,{cdop.GRIDFILE} ' \
              f'-selcode,172 {INIFILE}', options='-f nc')
    land_mask = cdo.addc( 1, input=f'-setctomiss,1 -setmisstoc,0 {ocean_mask}') 

    # trick to avoid the loop on years
    # define required years with a {year1,year2} and then use cdo select feature
    years_list = [str(element) for element in range(year1, year2+1)]
    years_joined = ','.join(years_list)

    # special treatment to exploit bash wild cards on multiple years
    if len(years_list) > 1:
        years_joined = '{' + years_joined + '}'

    if fverb: print(years_joined)

    # Load reference data
    ref = fn.load_yaml('pi_climatology.yml')

    # defines the two varlist
    field_2d = cfg['PI']['2d_vars']['field']
    field_3d = cfg['PI']['3d_vars']['field']
    field_oce = cfg['PI']['oce_vars']['field']
    field_ice = cfg['PI']['ice_vars']['field']
    field_all = field_2d + field_3d + field_oce + field_ice

    # check if required vars are available in the output
    # create a filename from the first year
    INFILE_2D = str(ECEDIR / 'output/oifs' / f'{expname}_atm_cmip6_1m_{year1}-{year1}.nc')
    INFILE_3D = str(ECEDIR / 'output/oifs' / f'{expname}_atm_cmip6_pl_1m_{year1}-{year1}.nc')
    INFILE_OCE = str(ECEDIR / 'output/nemo' / f'{expname}_oce_1m_T_{year1}-{year1}.nc')
    INFILE_ICE = str(ECEDIR / 'output/nemo' / f'{expname}_ice_1m_{year1}-{year1}.nc')

    # alternative method with loop
    isavail={}
    for a,b in zip([INFILE_2D, INFILE_3D, INFILE_OCE, INFILE_ICE],
                   [field_2d, field_3d, field_oce, field_ice]) : 
        isavail={**isavail, **fn.vars_are_there(a,b,ref)}

    # main loop
    varstat = {}
    for var in field_all :

        if not isavail[var] : 
            varstat[var] = float('NaN')
        else:

            # Refresh cdo pipe and specify type of file (atm or oce)
            # Temporary hack ... wil need to be homeginized with global_mean
            # The type of model should be read from ref
            if var in field_oce + field_ice : 
                cdop.domain = 'oce'
                model = 'nemo'
            else : 
                cdop.domain = 'atm'
                model = 'oifs' 

            cdop.start() # Start fresh pipe

            # This leaves the input file undefined for now. It can be set later with 
            # cdop.set_infile(infile) or by specifying input=infile cdop.execute

            # Select variable
            cdop.selectname(var)

            # extract info from reference.yml
            dataref = ref[var]['dataset']
            dataname = ref[var]['dataname']
            filetype = ref[var]['filetype']

            infile = str(ECEDIR / 'output' / model  / f'{expname}_{filetype}_{years_joined}-????.nc')
            clim = str(CLMDIR / f'climate_{dataref}_{dataname}.nc')
            vvvv = str(CLMDIR / f'variance_{dataref}_{dataname}.nc')

            cdop.timmean()

            oper = ref[var]['oper']
            if oper:
                cdop.chain(oper[1:]) # Hack to remove leading '-' - to be fixed
    
            # temporarily using remapbil instead of remapcon due to NEMO grid missing corner
            outfile = cdop.execute('remapbil', resolution, input=infile)

            if var in field_3d:

                # special treatment which includes vertical interpolation

                # extract the vertical levels from the file
                format_vlevels = get_levels(clim)

                cdop.chain(f'intlevelx,{format_vlevels}')
                cdop.zonmean()
                cdop.invertlat()
                cdop.sub(clim)
                cdop.sqr()
                cdop.div(vvvv)
                cdop.chain(f'vertmean -genlevelbounds,zbot=0,ztop=100000')

            else:
                cdop.invertlat()
                cdop.sub(clim)
                cdop.sqr()
                cdop.div(vvvv)

                # apply masks when needed
                # this is a hack for now
                if ref[var]['domain'] == 'land':
                    cdop.mul(land_mask)
                elif ref[var]['domain'] == 'ocean':
                    cdop.mul(ocean_mask)

                cdop.setname(var)

            #cmd2 = f'-setname,{var} {mask} {cmd_vertmean} -div -sqr -sub -invertlat {cmd_vertinterp} {oper} {outfile} {clim} {vvvv}'
            x = np.squeeze(cdop.execute('fldmean', input=outfile,
                                        returnCdf=True).variables[var])

            # store the PI
            varstat[var] = float(x)
            if fverb: print('PI for ', var, varstat[var])

    # define options for the output table
    head = ['Var', 'PI', 'Domain', 'Dataset', 'CMIP3', 'Ratio to CMIP3']
    global_table = list()

    # loop on the variables
    for var in field_all:
        out_sequence = [var, varstat[var], ref[var]['domain'], ref[var]
                        ['dataset'], ref[var]['cmip3'], varstat[var]/ref[var]['cmip3']]
        global_table.append(out_sequence)

    # nice loop on dictionary to get the partial and total pi
    partial_pi = np.mean([varstat[k] for k in field_2d + field_3d])
    total_pi = np.mean([varstat[k] for k in field_2d + field_3d + field_oce])

    # write the file  with tabulate: cool python feature
    tablefile = TABDIR / f'PI4_RK08_{expname}_{year1}_{year2}.txt'
    if fverb: print(tablefile)
    with open(tablefile, 'w') as f:
        f.write(tabulate(global_table, headers=head, tablefmt='orgtbl'))
        f.write('\n\nPartial PI (atm only) is   : ' + str(partial_pi))
        f.write('\nTotal Performance Index is : ' + str(total_pi))

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
