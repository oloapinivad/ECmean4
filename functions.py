#!/usr/bin/env python3
import re
import os.path
import sys
import yaml
from cdo import Cdo
cdo = Cdo()

# small function to define if a number is a float 
def is_number(s):
    """Check if input is a float type"""
    try:
        float(s)
        return True
    except ValueError:
        return False


# make sure all requested vars are available (use first year)
# first find all needed variables (including those needed for derived ones)
# added extra if to check if file exists
def vars_are_there(infile, var_needed, reference):
    """Check if a list of variables is available in the input file"""

    isavail = {}
    # if file exists, check which variables are inside
    if os.path.isfile(infile) : 
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

# given a folder, verify that the config.yml exists and open it
def load_config_file(indir): 
    """Load configuration file, once you have it!"""
    CONFIGFILE = str(indir / 'config.yml')
    if os.path.exists(CONFIGFILE):
        with open(CONFIGFILE, 'r') as file:
            cfg = yaml.load(file, Loader=yaml.FullLoader)
    else:
        sys.exit('config.yml not found: you need to have a configuration file!')

    return cfg

# create input filenames for the required variable and a given year
def make_input_filename(dr, var, expname, year1, year2, face):
    """Generate appropriate input filename for a variable"""
    filetype = face[var]['filetype']
    fname = dr / 'output' / face[var]['component'] / \
                f'{expname}_{filetype}_{year1}-{year2}.nc'
    return str(fname)

# write a tuning table: need to fix reference to face/ref
def write_tuning_table(linefile, varmean, var_table, expname, year1, year2, face, ref):
    """Write results appending one line to a text file"""
    if not os.path.isfile(linefile):
        with open(linefile, 'w') as f:
            print('%exp from   to ', end='', file=f)
            for var in var_table:
                print('{:>12s}'.format(var), end=' ', file=f)
            print('\n%             ', end=' ', file=f)
            for var in var_table:
                print('{:>12s}'.format(face[var]['units']), end=' ', file=f)
            print(file=f)

    with open(linefile, 'a') as f:
        print(expname,'{:4d} {:4d} '.format(year1, year2), end='', file=f)
        for var in var_table:
            print('{:12.5f}'.format(varmean[var] * ref[var].get('factor',1)), end=' ', file=f)
        print(file=f)


