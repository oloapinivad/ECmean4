#!/usr/bin/env python3
'''
Shared functions for ECmean4
'''
import re
import os.path
import sys
import logging
import yaml
import numpy as np
from cdo import Cdo
from metpy.units import units

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

    # if file exists, check which variables are inside
    isavail = {}
    isunit = {}
    if os.path.isfile(infile):

        # extract units from attributes: messy string manipulation to get it
        units_avail_list = cdo.showattribute('*@units', input=infile)

        # extract variables
        var_avail_list = [v.split()[1] for v in cdo.pardes(input=infile)]

        # create a dictionary, taking care when units are missing
        # not really pythonic, need to be improved
        var_avail = {}
        for u in units_avail_list:
            if u.replace(':', '') in var_avail_list:
                ind = units_avail_list.index(u)
                if 'units' in units_avail_list[ind+1]:
                    found_unit = re.findall('"([^"]*)"', units_avail_list[ind+1])[0]
                else:
                    found_unit = ''
                var_avail[u.replace(':', '')] = found_unit

        # loop on vars
        for v in var_needed:
            isavail[v] = True

            d = reference[v].get('derived')
            # if variable is derived, extract required vars
            if d:
                var_req = re.split('[*+-]', d)

                # check of unit is specified in the interface file
                u = reference[v].get('units')
                if u:
                    isunit[v] = u
                else:
                    logging.warning('%s is a derived var, assuming unit as the first of its term', v)
                    isunit[v] = var_avail[var_req[0]]

                # remove numbers
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
                    logging.warning("Variable %s needed by %s is not available in the model output!", x, v)
    else:
        for v in var_needed:
            isavail[v] = False
            isunit[v] = None
    return isavail, isunit


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


def units_extra_definition():
    """Add units to the pint registry required by ECMean4"""

    # special units definition, need to be moved in another placce
    units.define('fraction = [] = frac')
    units.define('psu = 1e-3 frac')
    units.define('Sv = 1e+6 m^3/s')  # Replace Sievert with Sverdrup


# use metpy/pint to provide factors for correction of units
def units_converter(org_units, tgt_units):
    """Units conversion using metpy and pint.
    From a org_units convert to tgt_units providing offset and factor.
    Some assumptions are done for precipitation field: must be extended to other vars.
    It will not work if BOTH factor and offset are required"""

    units_relation = (units(org_units)/units(tgt_units)).to_base_units()
    logging.debug(units_relation)
    if units_relation.magnitude != 1:
        logging.info('Unit conversion required...')
        offset_standard = 0 * units(org_units)
        factor_standard = 1 * units(org_units)
        if units_relation.units == units('dimensionless'):
            offset = offset_standard.to(tgt_units).magnitude
            if offset == 0:
                factor = factor_standard.to(tgt_units).magnitude
            else:
                factor = 1.

        elif units_relation.units == units('kg / m^3'):
            logging.debug("Assuming this as a water flux! Am I correct?")
            logging.debug("Dividing by water density...")
            density_water = units('kg / m^3') * 1000
            offset = 0.
            factor = (factor_standard/density_water).to(tgt_units).magnitude

        else:
            logging.error(units_relation)
            sys.exit("Units mismatch, this cannot be handled!")
    else:
        offset = 0.
        factor = 1.

    return {'offset': offset, 'factor': factor}


def units_are_integrals(org_units, ref_var):
    """Check functions for spatially integrated variables"""
    if 'total' in ref_var.keys():
        new_units = str((units(org_units) * units('m^2')).units)
    else:
        new_units = org_units
    return new_units


def directions_match(org, dst):
    """Check function for fluxes direction: they should match. Default is down"""
    direction_org = org.get('direction', 'down')
    direction_dst = dst.get('direction', 'down')
    if direction_org != direction_dst:
        factor = -1.
    else:
        factor = 1.
    return factor


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
        print(expname, '{:4d} {:4d} '.format(year1, year2), end='', file=f)
        for var in var_table:
            print('{:12.5f}'.format(varmean[var] * ref[var].get('factor', 1)), end=' ', file=f)
        print(file=f)
