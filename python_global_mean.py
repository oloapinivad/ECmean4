#!/usr/bin/env python3

# python3 version of ECmean global mean tool
# It uses a reference file from yaml and cdo bindings

import sys
import os
import yaml
import numpy as np
#import netCDF4 as nc
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
parser = argparse.ArgumentParser(description='ECmean global mean diagnostics for EC-Earth4')
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

ECEDIR = Path(cfg['dirs']['exp'], expname)
TABDIR = Path(cfg['dirs']['tab'], 'ECmean4', 'table')
#TMPDIR = Path(cfg['dirs']['tmp'], 'ECmean4', 'tmp')
os.makedirs(TABDIR, exist_ok=True)

# prepare grid description file
gridfile=str(TABDIR / 'grid.txt')
griddes = cdo.griddes(input=str(ECEDIR / f'ICMGG{expname}INIT'))
with open(gridfile, 'w') as f:
    for line in griddes:
        print(line, file=f)

# reference data
filename = 'reference.yml'
with open(filename, 'r') as file:
    ref = yaml.load(file, Loader=yaml.FullLoader)

# list of vars on which to work
var_field = cfg['atm_vars']['field']
var_radiation = cfg['atm_vars']['radiation']
var_table = cfg['tab_vars']
var_all = list(set(var_field + var_radiation + var_table)) # Extract all variables, avoid duplicates

# make sure all requested vars are available (use first year)
# first find all needed variables (including those needed for derived ones)
var_req = []
var_der = []
for v in var_all:
     d = ref[v].get('derived')
     if d:
         var_req.extend(re.split('[\*+-]', d))
         var_der.append(v)

for v in var_req:
    if is_number(v):
        var_req.remove(v)
# print("Vars requested", var_req)

var_req = list(set(var_req + var_all)) # this avoids duplicates
for v in var_der:
    var_req.remove(v)
# find available vars and check
infile = str(ECEDIR / 'output/oifs' / f'{expname}_atm_cmip6_1m_{year1}-{year1}.nc')
var_avail = [v.split()[1] for v in cdo.pardes(input=infile)]
for v in var_req:
    if v not in var_avail:
        sys.exit(f'Variable {v} needed but not available in model output!')

# loop
varstat = {}
for var in var_all:
    a = []
    if 'derived' in ref[var].keys():
       cmd = ref[var]['derived']
       der = f' -expr,{var}={cmd} '
    else: 
       der = ''
    for year in range(year1, year2+1):
        infile = ECEDIR / 'output/oifs' / f'{expname}_atm_cmip6_1m_{year}-{year}.nc'
#        cmd = f'-timmean -setgridtype,regular -setgrid,{gridfile} -selname,{var} {infile}'
        cmd = f'-timmean -zonmean -setgrid,{gridfile} -selname,{var} {der} {infile}' # Equivalent, faster
        x=cdo.fldmean(input=cmd, returnCdf = True).variables[var][:]
        a.append(x.item())
    varstat[var] = mean(a)
    if fverb: print('Average', var, mean(a))

# define options for the output table
head = ['Var', 'Longname', 'Units', 'ECE4', 'OBS', 'Obs Dataset']
global_table = list()

# loop on the variables
for var in var_field + var_radiation:
    beta = ref[var]
    beta['value'] = varstat[var] * float(beta['factor'])
    out_sequence = [var, beta['varname'], beta['units'], beta['value'],
                    float(beta['observations']['val']), beta['observations']['data']]
    global_table.append(out_sequence)

# write the file  with tabulate: cool python feature
tablefile = TABDIR / f'Global_Mean_{expname}_{year1}_{year2}.txt'
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
