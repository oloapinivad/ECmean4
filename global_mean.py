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
    description='ECmean global mean diagnostics for EC-Earth4')
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

ECEDIR = Path(cfg["dirs"]["exp"], expname)
TABDIR = Path(cfg["dirs"]["tab"], "ECmean4", "table")
os.makedirs(TABDIR, exist_ok=True)

#cdo.forceOutput = True
#cdo.debug = True

# prepare grid description file
gridfile = str(TABDIR / "grid.txt")
griddes = cdo.griddes(input=str(ECEDIR / f"ICMGG{expname}INIT"))
with open(gridfile, "w") as f:
    for line in griddes:
        print(line, file=f)

# reference data: it is badly written but it can be implemented in a much more intelligent
# and modular way
filename = "reference.yml"
with open(filename, 'r') as file:
    ref = yaml.load(file, Loader=yaml.FullLoader)

# loop
varstat = {}
var_field = cfg["global"]["atm_vars"]["field"]
var_radiation = cfg["global"]["atm_vars"]["radiation"]
for var in var_field + var_radiation:
    a = []
    if 'derived' in ref[var].keys():
        cmd = ref[var]['derived']
        der = f" -expr,{var}={cmd} "
    else:
        der = ""
    for year in range(year1, year2+1):
        infile = ECEDIR / "output/oifs" / \
            f"{expname}_atm_cmip6_1m_{year}-{year}.nc"
        #cmd = f"-timmean -zonmean -setgrid,{gridfile} -selname,{var} {der} {infile}"
        #x=cdo.fldmean(input=cmd, returnCdf = True).variables[var][:]
        #a.append(x.item())
        # simpler with cdo output
        cmd = f"-timmean -fldmean -setgridtype,regular -setgrid,{gridfile} -selname,{var} {der} {infile}"
        x = float(cdo.output(input=cmd)[0])
        a.append(x)
    varstat[var] = mean(a)
    print("Average", var, mean(a))

# define options for the output table
head = ["Var", "Longname", "Units", "ECE4", "OBS", "Obs Dataset"]
global_table = list()

# loop on the variables
for var in var_field + var_radiation:
    print(var)
    beta = ref[var]
    beta['value'] = varstat[var] * float(beta['factor'])
    out_sequence = [var, beta['varname'], beta['units'], beta['value'],
                    float(beta['observations']['val']), beta['observations']['data']]
    global_table.append(out_sequence)

# write the file  with tabulate: cool python feature
tablefile = TABDIR / f"Global_Mean_{expname}_{year1}_{year2}.txt"
print(tablefile)
with open(tablefile, 'w') as f:
    f.write(tabulate(global_table, headers=head, tablefmt="orgtbl"))
