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

    if os.path.isfile(infile) : 

        var_avail = [v.split()[1] for v in cdo.pardes(input=infile)]

        for v in var_needed:
            isavail[v] = True
            d = reference[v].get('derived')
            if d:
                var_req = re.split('[*+-]', d)
                for x in var_req:
                    if is_number(x):
                        var_req.remove(x)
            else:
                var_req = [v]

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

    CONFIGFILE = str(indir / 'config.yml')
    if os.path.exists(CONFIGFILE):
        with open(CONFIGFILE, 'r') as file:
            cfg = yaml.load(file, Loader=yaml.FullLoader)
    else:
        sys.exit('config.yml not found: you need to have a configuration file!')

    return cfg


