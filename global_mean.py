#!/usr/bin/env python3

# python3 version of ECmean global mean tool
# It uses a reference file from yaml and cdo bindings

import os
import argparse
from statistics import mean
from pathlib import Path
import yaml
from tabulate import tabulate
from cdo import Cdo

# import functions from functions.py
import functions as fn

# arguments
parser = argparse.ArgumentParser(
    description='ECmean global mean diagnostics for EC-Earth4')
parser.add_argument('exp', metavar='EXP', type=str, help='experiment ID')
parser.add_argument('year1', metavar='Y1', type=int, help='starting year')
parser.add_argument('year2', metavar='Y2', type=int, help='final year')
parser.add_argument('-s', '--silent', action='store_true',
                    help='do not print anything to std output')
parser.add_argument('-t', '--table', action='store_true',
                    help='appends also single line to a table')
parser.add_argument('-o', '--output', type=str, default='',
                    help='path of output one-line table')
args = parser.parse_args()
expname = args.exp
year1 = args.year1
year2 = args.year2
fverb = not args.silent
ftable = args.table

cdo = Cdo()

# config file (looks for it in the same dir as the .py program file
INDIR = Path(os.path.dirname(os.path.abspath(__file__)))
with open(INDIR / 'config.yml', 'r') as file:
    cfg = yaml.load(file, Loader=yaml.FullLoader)

ECEDIR = Path(os.path.expandvars(cfg['dirs']['exp']), expname)
TABDIR = Path(os.path.expandvars(cfg['dirs']['tab']))
TMPDIR = Path(os.path.expandvars(cfg['dirs']['tmp']))
os.makedirs(TABDIR, exist_ok=True)

print(TMPDIR)

# prepare grid description file
GRIDFILE=str(TMPDIR / 'grid.txt')
INIFILE=str(ECEDIR / f'ICMGG{expname}INIT')
griddes = cdo.griddes(input=INIFILE)
with open(GRIDFILE, 'w') as f:
    for line in griddes:
        print(line, file=f)

# prepare LSM
LMFILE=str(TMPDIR / 'lmask.nc')
SMFILE=str(TMPDIR / 'smask.nc')
GAFILE=str(TMPDIR / 'ga.nc')
cdo.selname('LSM', input=f'-setgridtype,regular {INIFILE}', output=LMFILE, options='-t ecmwf')
cdo.mulc('-1', input=f'-subc,1 {LMFILE}', output=SMFILE)
cdo.gridarea(input=f'-setgridtype,regular {LMFILE}', output=GAFILE)

# reference data
with open(INDIR / 'reference.yml', 'r') as file:
    ref = yaml.load(file, Loader=yaml.FullLoader)

# list of vars on which to work
var_field = cfg['global']['atm_vars']['field']
var_radiation = cfg['global']['atm_vars']['radiation']
var_table = cfg['global']['tab_vars']

# make sure all requested vars are available (use first year)
# first find all needed variables (including those needed for derived ones)

# create a list for all variable, avoid duplicates
var_all = list(set(var_field + var_radiation + var_table))

# create a filename
INFILE = str(ECEDIR / 'output/oifs' / f'{expname}_atm_cmip6_1m_{year1}-{year1}.nc')

# which vars are available
isavail=fn.vars_are_there(INFILE, var_all, ref)

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
            der = f'-expr,{var}={cmd}'
        else:
            der = ''

        # land/sea variables
        mask = ''
        op='-fldmean'
        if 'total' in ref[var].keys():
            mask_type = ref[var]['total']
            if mask_type == 'land':
                mask = f'-mul {GAFILE} -mul {LMFILE}'
                op='-fldsum'
            elif mask_type in ['sea', 'ocean']:
                mask = f'-mul {GAFILE} -mul {SMFILE}'
                op='-fldsum'

        # loop on years: call CDO to perform all the computations
        for year in range(year1, year2+1):
            INFILE = ECEDIR / 'output/oifs' / \
               f'{expname}_atm_cmip6_1m_{year}-{year}.nc'
            cmd = f'-timmean {op} {mask} -setgridtype,regular ' \
                  f'-setgrid,{GRIDFILE} -selname,{var} {der} {INFILE}'
            x = float(cdo.output(input=cmd)[0])
            a.append(x)
        varstat[var] = mean(a)
        if fverb:
            print('Average', var, mean(a))

# define options for the output table
head = ['Variable', 'Longname', 'Units', 'EC-Earth4', 'Obs.', 'Dataset', 'Years']
global_table = []

# loop on the variables to create the table
for var in var_field + var_radiation:
    beta = ref[var]
    beta['value'] = varstat[var] * float(beta['factor'])
    out_sequence = [var, beta['varname'], beta['units'], beta['value'],
                    float(beta['observations']['val']),
                    beta['observations'].get('data',''),
                    beta['observations'].get('years','')]
    global_table.append(out_sequence)

# write the file with tabulate: cool python feature
tablefile = TABDIR / f'global_mean_{expname}_{year1}_{year2}.txt'
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
os.unlink(GRIDFILE)
os.unlink(LMFILE)
os.unlink(SMFILE)
os.unlink(GAFILE)
