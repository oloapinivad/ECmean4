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
# extract also units from files
def vars_are_there(infile, var_needed, reference):
    """Check if a list of variables is available in the input file"""

    # if file exists, check which variables are inside
    isavail = {}
    isunit = {}
    if os.path.isfile(infile) : 

        # extract units from attributes: messy string manipulation to get it
        #units_avail_list = re.findall('"([^"]*)"', ",".join(cdo.showattribute('*@units', input=infile)[1::2]))
        units_avail_list = cdo.showattribute('*@units', input=infile)
        
        # extract variables
        var_avail_list = [v.split()[1] for v in cdo.pardes(input=infile)]

        # create a dictionary, taking care when units are missing
        # not really pythonic, need to be improved
        var_avail = {} 
        for u in units_avail_list : 
            if u.replace(':','') in var_avail_list : 
                ind = units_avail_list.index(u)
                if 'units' in units_avail_list[ind+1] : 
                    found_unit = re.findall('"([^"]*)"', units_avail_list[ind+1])[0]
                else : 
                    found_unit = ''
                var_avail[u.replace(':','')] = found_unit

        # loop on vars
        for v in var_needed:
            isavail[v] = True
            
            d = reference[v].get('derived')
            # if variable is derived, extract required vars
            if d:
                var_req = re.split('[*+-]', d)
                isunit[v] = var_avail[var_req[0]]
                for x in var_req:
                    if is_number(x):
                        var_req.remove(x)
            else:
                var_req = [v]
                isunit[v] = var_avail[v]

            # check if required varialbes are in model output
            for x in var_req:
                if x not in var_avail:
                    isavail[v] = False
                    isunit[v] = None
                    print(f"Variable {x} needed by {v} is not available in the model output!")
    else: 
        for v in var_needed : 
            isavail[v] = False
            isunit[v] = None
    return isavail, isunit

# given a folder, verify that the config.yml exists and open it
def load_config_file(indir): 

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


