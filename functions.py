#!/usr/bin/env python3
'''
Shared functions for ECmean4
'''
import re
import os.path
import sys
import yaml
import numpy as np
from cdo import Cdo

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

def is_number(s):
    """Check if input is a float type"""
    try:
        float(s)
        return True
    except ValueError:
        return False

def vars_are_there(infile, var_needed, reference):
    """Check if a list of variables is available in the input file.
       Make sure all requested vars are available (use first year)
       first find all needed variables (including those needed for derived ones)
       added extra if to check if file exists"""

    isavail = {}
    # if file exists, check which variables are inside
    if os.path.isfile(infile):
        var_avail = [v.split()[1] for v in cdo.pardes(input=infile)]

        # loop on vars
        for v in var_needed:
            isavail[v] = True
            d = reference[v].get('derived')
            # if variable is derived, extract required vars
            if d:
                var_req = re.split('[*+-]', d)
                for x in var_req:
                    if is_number(x):
                        var_req.remove(x)
            else:
                var_req = [v]

            # check if required varialbes are in model output
            for x in var_req:
                if x not in var_avail:
                    isavail[v] = False
                    print(f"Variable {x} needed by {v} is not available in the model output!")
    else:
        for v in var_needed :
            isavail[v] = False
    return isavail

def load_yaml(infile):
    """Load generic yaml file"""
    try:
        with open(infile, 'r', encoding='utf-8') as file:
            cfg = yaml.load(file, Loader=yaml.FullLoader)
    except IOError:
        sys.exit(f'{infile} not found: you need to have this configuration file!')
    return cfg

def make_input_filename(dr, var, expname, year1, year2, face):
    """Create input filenames for the required variable and a given year"""

    filetype = face[var]['filetype']
    fname = dr / 'output' / face[var]['component'] / \
                f'{expname}_{filetype}_{year1}-{year2}.nc'
    return str(fname)

def write_tuning_table(linefile, varmean, var_table, expname, year1, year2, face, ref):
    """Write results appending one line to a text file.
       Write a tuning table: need to fix reference to face/ref"""

    if not os.path.isfile(linefile):
        with open(linefile, 'w', encoding='utf-8') as f:
            print('%exp from   to ', end='', file=f)
            for var in var_table:
                print('{:>12s}'.format(var), end=' ', file=f)
            print('\n%             ', end=' ', file=f)
            for var in var_table:
                print('{:>12s}'.format(face[var]['units']), end=' ', file=f)
            print(file=f)

    with open(linefile, 'a', encoding='utf-8') as f:
        print(expname,'{:4d} {:4d} '.format(year1, year2), end='', file=f)
        for var in var_table:
            print('{:12.5f}'.format(varmean[var] * ref[var].get('factor',1)), end=' ', file=f)
        print(file=f)
