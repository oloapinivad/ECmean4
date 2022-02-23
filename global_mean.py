#!/usr/bin/env python3

# python3 version of ECmean global mean tool
# It uses a reference file from yaml and cdo bindings

import sys
import os
import yaml
import numpy as np
import argparse
from tabulate import tabulate
from statistics import mean
from cdo import *
from pathlib import Path
import re

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# arguments
parser = argparse.ArgumentParser(
    description='ECmean global mean diagnostics for EC-Earth4')
parser.add_argument('exp', metavar='EXP', type=str, help='experiment ID')
parser.add_argument('year1', metavar='Y1', type=int, help='starting year')
parser.add_argument('year2', metavar='Y2', type=int, help='final year')
parser.add_argument('-s', '--silent', action='store_true', help='do not print anything to std output')
parser.add_argument('-t', '--table', action='store_true', help='appends also single line to a table')
parser.add_argument('-o', '--output', type=str, default='',  help='path of output one-line table')
args = parser.parse_args()
expname = args.exp
year1 = args.year1
year2 = args.year2
fverb = not args.silent
ftable = args.table

cdo = Cdo()

# config file (looks for it in the same dir as the .py program file
indir = Path(os.path.dirname(os.path.abspath(__file__)))
with open(indir / 'config.yml', 'r') as file:
    cfg = yaml.load(file, Loader=yaml.FullLoader)

#cdo.forceOutput = True
#cdo.debug = True

ECEDIR = Path(cfg['dirs']['exp'], expname)
TABDIR = Path(cfg['dirs']['tab'], 'ECmean4', 'table')
os.makedirs(TABDIR, exist_ok=True)

# prepare grid description file
gridfile=str(TABDIR / 'grid.txt')
inifile=str(ECEDIR / f'ICMGG{expname}INIT')
griddes = cdo.griddes(input=inifile)
with open(gridfile, 'w') as f:
    for line in griddes:
        print(line, file=f)

# prepare LSM
lmfile=str(TABDIR / 'lmask.nc')
smfile=str(TABDIR / 'smask.nc')
gafile=str(TABDIR / 'ga.nc')
cdo.selname('LSM', input=f'-setgridtype,regular {inifile}', output=lmfile, options='-t ecmwf')
cdo.mulc('-1', input=f'-subc,1 {lmfile}', output=smfile)  
cdo.gridarea(input=f'-setgridtype,regular {lmfile}', output=gafile)

# reference data
filename = 'reference.yml'
with open(filename, 'r') as file:
    ref = yaml.load(file, Loader=yaml.FullLoader)

# list of vars on which to work
var_field = cfg['global']['atm_vars']['field']
var_radiation = cfg['global']['atm_vars']['radiation']
var_table = cfg['global']['tab_vars']

# make sure all requested vars are available (use first year)
# first find all needed variables (including those needed for derived ones)

# create a list for all variable, avoid duplicates
var_all = list(set(var_field + var_radiation + var_table))

infile = str(ECEDIR / 'output/oifs' / f'{expname}_atm_cmip6_1m_{year1}-{year1}.nc')
var_avail = [v.split()[1] for v in cdo.pardes(input=infile)]

#var_der = []
isavail = {}
for v in var_all:
     isavail[v] = True
     d = ref[v].get('derived')
     if d:
#         var_der.append(v)
         var_req = re.split('[\*+-]', d)
         for x in var_req:
             if is_number(x):
                 var_req.remove(x)
     else:
         var_req = [v] 

     for x in var_req:
         if x not in var_avail:
             isavail[v] = False
             print(f"Variable {x} needed by {v} is not available in the model output!")

# loop
varstat = {}
for var in var_all:
    if not isavail[var]:
        varstat[var] = float("NaN")
    else:
        a = []

        # check if var is derived
        if 'derived' in ref[var].keys():
            cmd = ref[var]['derived']
            der = f' -expr,{var}={cmd} '
        else:
            der = ''
    
        # land/sea variables
        mask = ''
        op='-fldmean'
        if 'total' in ref[var].keys():
            mask_type = ref[var]['total']
            if mask_type == 'land':
                mask = f'-mul {gafile} -mul {lmfile}'
                op='-fldsum'
            elif mask_type in ['sea', 'ocean']:
                mask = f'-mul {gafile} -mul {smfile}'
                op='-fldsum'

        # loop on years: call CDO to perform all the computations
        for year in range(year1, year2+1):
            infile = ECEDIR / 'output/oifs' / \
               f'{expname}_atm_cmip6_1m_{year}-{year}.nc'
            cmd = f'-timmean {op} {mask} -setgridtype,regular -setgrid,{gridfile} -selname,{var} {der} {infile}'
            x = float(cdo.output(input=cmd)[0])
            a.append(x)
        varstat[var] = mean(a)
        if fverb: print('Average', var, mean(a))

# define options for the output table
head = ['Variable', 'Longname', 'Units', 'EC-Earth4', 'Obs.', 'Dataset', 'Years']
global_table = list()

# loop on the variables to create the table
for var in var_field + var_radiation:
    beta = ref[var]
    beta['value'] = varstat[var] * float(beta['factor'])
    out_sequence = [var, beta['varname'], beta['units'], beta['value'],
                    float(beta['observations']['val']),
                    beta['observations'].get('data',''),  beta['observations'].get('years','')]
    global_table.append(out_sequence)

# write the file with tabulate: cool python feature
tablefile = TABDIR / f'Global_Mean_{expname}_{year1}_{year2}.txt'
print(tablefile)
with open(tablefile, 'w') as f:
    f.write(tabulate(global_table, headers=head, tablefmt='orgtbl'))

# Print appending one line to table (for tuning)
if ftable:
    linefile = TABDIR / 'global_means.txt'
    if args.output: 
        linefile = args.output
    if not os.path.isfile(linefile):
        with open(linefile, 'w') as f:
            print('%exp from   to ', end='', file=f)
            for var in var_table:
                print('{:>12s}'.format(var), end=' ', file=f)
            print('\n%             ', end=' ', file=f)
            for var in var_table:
                print('{:>12s}'.format(ref[var]['units']), end=' ', file=f)
            print(file=f)

    with open(linefile, 'a') as f:
        print(expname,'{:4d} {:4d} '.format(year1, year2), end='', file=f)
        for var in var_table: 
            print('{:12.5f}'.format(varstat[var] * ref[var]['factor']), end=' ', file=f)
        print(file=f)

# clean
os.unlink(gridfile)
os.unlink(lmfile)
os.unlink(smfile)
os.unlink(gafile)
