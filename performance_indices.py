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

cdo = Cdo()

# arguments
parser = argparse.ArgumentParser(
    description='ECmean Performance Indices for EC-Earth4')
parser.add_argument('exp', metavar='EXP', type=str, help='experiment ID')
parser.add_argument('year1', metavar='Y1', type=int, help='starting year')
parser.add_argument('year2', metavar='Y2', type=int, help='final year')
args = parser.parse_args()
expname = args.exp
year1 = args.year1
year2 = args.year2

# config file (looks for it in the same dir as the .py program file
indir = Path(os.path.dirname(os.path.abspath(__file__)))
with open(indir / "config.yml", "r") as file:
    cfg = yaml.load(file, Loader=yaml.FullLoader)

# hard-coded resolution (due to climatological dataset)
resolution = "r180x91"

# folder definition
ECEDIR = Path(os.path.expandvars(cfg["dirs"]["exp"]), expname)
TABDIR = Path(os.path.expandvars(cfg["dirs"]["tab"]), "ECmean4", "table")
TMPDIR = Path(os.path.expandvars(cfg["dirs"]["tmp"]), "tmp", expname)
CLMDIR = Path(os.path.expandvars(cfg["dirs"]["clm"]))
os.makedirs(TABDIR, exist_ok=True)
os.makedirs(TMPDIR, exist_ok=True)

#cdo.forceOutput = True
# cdo.debug = True

# prepare grid description file
icmgg_file = ECEDIR / f"ICMGG{expname}INIT"
gridfile = str(TABDIR / "grid.txt")
griddes = cdo.griddes(input=f"{icmgg_file}")
with open(gridfile, "w") as f:
    for line in griddes:
        print(line, file=f)

# land-sea masks
land_mask = TMPDIR / "land_mask_2x2.nc"
ocean_mask = TMPDIR / "ocean_mask_2x2.nc"
cdo.setctomiss(
    0, input=f"-ltc,0.5 -invertlat -remapcon2,{resolution} -setgridtype,regular -setgrid,{gridfile} -selcode,172 {icmgg_file}", output=f"{ocean_mask}", options="-f nc")
cdo.addc(
    1, input=f"-setctomiss,1 -setmisstoc,0 {ocean_mask}", output=f"{land_mask}")

# trick to avoid the loop on years
# define required years with a {year1,year2} and then use cdo select feature
years_list = [str(element) for element in range(year1, year2+1)]
years_joined = ",".join(years_list)

# special treatment to exploit bash wild cards on multiple years
if len(years_list) > 1:
    years_joined = '{' + years_joined + '}'

print(years_joined)

# reference data: it is badly written but it can be implemented in a much more intelligent
# and modular way
filename = "pi_climatology.yml"
with open(filename, 'r') as file:
    ref = yaml.load(file, Loader=yaml.FullLoader)

# defines the two varlist
field_2d = cfg["PI"]["2d_vars"]["field"]
field_3d = cfg["PI"]["3d_vars"]["field"]
field_oce = cfg["PI"]["oce_vars"]["field"]
field_all = field_2d + field_3d + field_oce

# check if required vars are available in the output
# create a filename from the first year
INFILE_2D = str(ECEDIR / 'output/oifs' / f'{expname}_atm_cmip6_1m_{year1}-{year1}.nc')
INFILE_3D = str(ECEDIR / 'output/oifs' / f'{expname}_atm_cmip6_pl_1m_{year1}-{year1}.nc')
INFILE_OCE = str(ECEDIR / 'output/nemo' / f'{expname}_oce_cmip6_1m_{year1}-{year1}.nc')
#isavail_2d=fn.vars_are_there(INFILE_2D, field_2d, ref)
#isavail_3d=fn.vars_are_there(INFILE_3D, field_3d, ref)
#isavail_oce=fn.vars_are_there(INFILE_OCE, field_oce, ref)
#isavail={**isavail_2d, **isavail_3d, **isavail_oce}

# alternative method with loop
isavail={}
for a,b in zip([INFILE_2D, INFILE_3D, INFILE_OCE], [field_2d, field_3d, field_oce]) : 
    isavail={**isavail, **fn.vars_are_there(a,b,ref)}


# loop
varstat = {}
for var in field_all :

    if not isavail[var] : 
        varstat[var] = float("NaN")
    else:

        # extract info from reference.yml
        oper = ref[var]["oper"]
        dataref = ref[var]["dataset"]
        dataname = ref[var]["dataname"]
        filetype = ref[var]["filetype"]

        # file names
        infile = str(ECEDIR / "output/oifs" / \
            f"{expname}_atm_cmip6_{filetype}_{years_joined}-????.nc")
        outfile = str(TMPDIR / f"tmp_{expname}_{var}_2x2grid.nc")
        clim = str(CLMDIR / f"climate_{dataref}_{dataname}.nc")
        vvvv = str(CLMDIR / f"variance_{dataref}_{dataname}.nc")

        # apply masks when needed
        if ref[var]["domain"] == "land":
            mask = f"-mul {land_mask}"
        elif ref[var]["domain"] == "ocean":
            mask = f"-mul {ocean_mask}"
        elif ref[var]["domain"] == "global":
            mask = ""

        # timmean and remap
        cmd1 = f"-timmean -setgridtype,regular -setgrid,{gridfile} -select,name={var} {infile}"
        cdo.remapcon2(resolution, input=cmd1, output=outfile)

        if var in field_3d:

            # special treatment which includes vertical interpolation
            # and extraction of pressure weights

            # extract the vertical levels from the file
            vlevels = cdo.showlevel(input=f"{clim}")

            # perform multiple string manipulation to produce a Pa list of levels
            v0 = ''.join(vlevels).split()
            v1 = [int(x) for x in v0]

            # if the grid is hPa, move to Pa (there might be a better solution)
            if np.max(v1) < 10000:
                v1 = [x * 100 for x in v1]

            # compute level weights (numpy is not that smart, or it's me?)
            half_levels = np.convolve(v1, np.ones(2)/2, mode='valid')
            args = (np.array([0]), half_levels, np.array([100000]))
            level_weights = np.diff(np.concatenate(args))

            # format for CDO, converting to string
            format_vlevels = ' '.join(str(x) for x in v1).replace(" ", ",")

            # assign the vertical command for interpolation and zonal mean
            cmd_vertical = f"-zonmean -intlevelx,{format_vlevels}"
        else:
            cmd_vertical = ""

        # compute the PI
        cmd2 = f"-setname,{var} {mask} -div -sqr -sub -invertlat {cmd_vertical} {oper} {outfile} {clim} {vvvv}"
        x = np.squeeze(cdo.fldmean(input=cmd2, returnCdf=True).variables[var])

        # deprecated: pre-estimated weights from previous version of ECmean (climatology dependent!)
        #level_weights = np.array([30, 45, 75, 100, 100, 100, 150, 175, 112.5, 75, 37.5])

        # weighted average
        if var in field_3d:
            x = np.average(x, weights=level_weights)

        # store the PI
        varstat[var] = float(x)
        print("PI for ", var, varstat[var])

        # clean
        os.unlink(outfile)

# define options for the output table
head = ["Var", "PI", "Domain", "Dataset", "CMIP3", "Ratio to CMIP3"]
global_table = list()

# loop on the variables
for var in field_all:
    out_sequence = [var, varstat[var], ref[var]["domain"], ref[var]
                    ["dataset"], ref[var]["cmip3"], varstat[var]/ref[var]["cmip3"]]
    global_table.append(out_sequence)

# write the file  with tabulate: cool python feature
tablefile = TABDIR / f"PI4_RK08_{expname}_{year1}_{year2}.txt"
print(tablefile)
with open(tablefile, 'w') as f:
    f.write(tabulate(global_table, headers=head, tablefmt="orgtbl"))
