#!/usr/bin/env python

# this is a tentative python script to convert ECmean global mean operation to python3
# uses a reference file from yaml and cdo bindings (not very efficient)

import sys
import os
import glob
import yaml
import netCDF4 as nc
import argparse
from tabulate import tabulate
from cdo import *
from pathlib import Path
cdo = Cdo()

# arguments
parser = argparse.ArgumentParser(description='ECmean global mean diagnostics for EC-Earth4')
parser.add_argument('exp', metavar='EXP', type=str, help='experiment ID')
parser.add_argument('year1', metavar='Y1', type=int, help='starting year')
parser.add_argument('year2', metavar='Y2', type=int, help='final year')
args = parser.parse_args()
expname = args.exp
year1 = args.year1
year2 = args.year2

# config file
with open("config.yml", "r") as file:
    cfg = yaml.load(file, Loader=yaml.FullLoader)

ECEDIR = Path(cfg["dirs"]["exp"], expname)
TABDIR = Path(cfg["dirs"]["tab"], "ECmean4", "table")
TMPDIR = Path(cfg["dirs"]["tmp"], "ECmean4", "tmp")
#ECEDIR = os.path.join(cfg["dirs"]["exp"], expname)
#TABDIR = os.path.join(cfg["dirs"]["tab"], "ECmean4", "table")
#TMPDIR = os.path.join(cfg["dirs"]["tmp"], "ECmean4", "tmp")
os.makedirs(TABDIR, exist_ok=True)
os.makedirs(TMPDIR, exist_ok=True)
eceinitfile = ECEDIR / f"ICMGG{expname}INIT"

#cdo.forceOutput = True
#cdo.debug = True
vardict = {}

# loop
var_field = cfg["atm_vars"]["field"]
var_radiation = cfg["atm_vars"]["radiation"]
for var in var_field + var_radiation:
    for year in range(year1, year2+1):
        print(var, year)
        infile = ECEDIR / "output/oifs" / f"{expname}_atm_cmip6_1m_{year}-{year}.nc"
        outfile = TMPDIR / f"tmp_{year}.nc"
        # regrid and select variable
        cdo.fldmean(input=f"-setgridtype,regular -setgrid,{eceinitfile} -selname,{var} {infile}", output=str(outfile))
    # cat and timmean
    cdo.timmean(input=f"-cat {TMPDIR}/tmp*.nc", output=f"{TMPDIR}/{var}_mean.nc")
    # store mean value in a local dictionary
    vardict[var] = float(nc.Dataset( TMPDIR / f"{var}_mean.nc").variables[var][:])

# extra radiative variables
extra_radiation = ["net_toa", "net_sfc"]
vardict["net_toa"] = vardict["rsnt"] + vardict["rlnt"]
vardict["net_sfc"] = vardict["rsns"] + \
    vardict["rlns"] - vardict["hfls"] - vardict["hfss"]

# reference data: it is badly written but it can be implemented in a much more intelligent
# and modulable way
filename = "reference.yml"
with open(filename, 'r') as file:
    yummy = yaml.load(file, Loader=yaml.FullLoader)

# define options for the output table
head = ["Var", "Longname", "Units", "ECE4", "OBS", "Obs Dataset"]
global_table = list()

# loop on the variables
for var in var_field + var_radiation + extra_radiation:
    print(var)
    beta = yummy[var]
    beta['value'] = vardict[var] * float(beta['factor'])
    out_sequence = [var, beta['varname'], beta['units'], beta['value'],
                    float(beta['observations']['val']), beta['observations']['data']]
    global_table.append(out_sequence)

# write the file  with tabulate: cool python feature
tablefile = TABDIR / f"Global_Mean_{expname}_{year1}_{year2}.txt"
with open(tablefile, 'w') as f:
    f.write(tabulate(global_table, headers=head, tablefmt="orgtbl"))

# clean
for f in TMPDIR.glob("*.nc"):
    os.remove(f)
