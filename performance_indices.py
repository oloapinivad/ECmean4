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
ECEDIR = Path(cfg["dirs"]["exp"], expname)
TABDIR = Path(cfg["dirs"]["tab"], "ECmean4", "table")
TMPDIR = Path(cfg["dirs"]["tmp"], "tmp", expname)
CLMDIR = Path(cfg["dirs"]["clm"])
os.makedirs(TABDIR, exist_ok=True)
os.makedirs(TMPDIR, exist_ok=True)

#cdo.forceOutput = True
#cdo.debug = True

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

# reference data: it is badly written but it can be implemented in a much more intelligent
# and modular way
filename = "pi_climatology.yml"
with open(filename, 'r') as file:
    ref = yaml.load(file, Loader=yaml.FullLoader)

# loop
varstat = {}
field_2d = cfg["PI"]["2d_vars"]["field"]
field_3d = cfg["PI"]["3d_vars"]["field"]
for var in field_2d :

    # extract info from reference.yml
    oper = ref[var]["oper"]
    dataref = ref[var]["dataset"]
    dataname = ref[var]["dataname"]

    # file names
    infile = ECEDIR / "output/oifs" / \
        f"{expname}_atm_cmip6_1m_{{{years_joined}}}-????.nc"
    outfile = TMPDIR / f"tmp_{expname}_{var}_2x2grid.nc"
    clim = CLMDIR / f"climate_{dataref}_{dataname}.nc"
    vvvv = CLMDIR / f"variance_{dataref}_{dataname}.nc"

    # apply masks when needed
    if ref[var]["domain"] == "land":
        mask = f"-mul {land_mask}"
    elif ref[var]["domain"] == "ocean":
        mask = f"-mul {ocean_mask}"
    elif ref[var]["domain"] == "global":
        mask = ""

    # timmean and remap
    cmd1 = f"-timmean -setgridtype,regular -setgrid,{gridfile} -select,name={var} {infile}"
    cdo.remapcon2(resolution, input=cmd1, output=f"{outfile}")

    # cannot test it, missing 3d vars in ECE output still
    # ERA40 data ha no weights, so pressure level are not weighted
    if var in field_3d :
        cmd_vertical = f"-zonmean -intllevel,{clim}"
        mask="-vertmean"
    else :
        cmd_vertical = ""

    # compute the PI
    cmd2 = f"-setname,{var} {mask} -div -sqr -sub -invertlat {cmd_vertical} {oper} {outfile} {clim} {vvvv}"
    varstat[var] = float(cdo.fldmean(input=cmd2, returnCdf=True).variables[var][:])
    print("PI for ", var, varstat[var])

    # clean
    os.remove(f"{outfile}")

# define options for the output table
head = ["Var", "PI", "Domain", "Dataset", "CMIP3", "Ratio to CMIP3"]
global_table = list()

# loop on the variables
for var in field_2d:
    out_sequence = [var, varstat[var], ref[var]["domain"], ref[var]
                    ["dataset"], ref[var]["cmip3"], varstat[var]/ref[var]["cmip3"]]
    global_table.append(out_sequence)

# write the file  with tabulate: cool python feature
tablefile = TABDIR / f"PI4_RK08_{expname}_{year1}_{year2}.txt"
print(tablefile)
with open(tablefile, 'w') as f:
    f.write(tabulate(global_table, headers=head, tablefmt="orgtbl"))
