#!/usr/bin/env python3
'''
Shared functions and classes for ECmean4
'''
import re
import os.path
import sys
import logging
import itertools
from pathlib import Path
import yaml
import numpy as np
from cdo import Cdo
from metpy.units import units

cdo = Cdo()


class Diagnostic():
    """General container class for common variables"""
    def __init__(self, args, cfg):
        self.expname = args.exp
        self.year1 = args.year1
        self.year2 = args.year2
        self.fverb = not args.silent
        self.ftable = getattr(args, 'line', False)
        self.ftrend = getattr(args, 'trend', False)
        self.numproc = args.numproc
        self.modelname = getattr(args, 'model', '')
        if not self.modelname:
            self.modelname = cfg['model']['name']
        if self.year1 == self.year2:  # Ignore if only one year requested
            self.ftrend = False

        # hard-coded resolution (due to climatological dataset)
        self.resolution = cfg['PI']['resolution']

        # Various input and output directories
        self.ECEDIR = Path(os.path.expandvars(cfg['dirs']['exp']))
        self.TABDIR = Path(os.path.expandvars(cfg['dirs']['tab']))
        self.CLMDIR = Path(os.path.expandvars(cfg['dirs']['clm']), self.resolution)
        self.years_joined = ''

        self.linefile = self.TABDIR / 'global_means.txt'

        # check if output attribute exists
        if hasattr(self, 'output') : 
           self.linefile = args.output
           self.ftable = True
          

def chunks(iterable, num):
    """Generate num adjacent chunks of data from a list iterable"""
    """Split lists in a convienet way for a parallel process"""
    size = int(np.ceil(len(iterable) / num))
    it = iter(iterable)
    return iter(lambda: tuple(itertools.islice(it, size)), ())


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
                    logging.warning('%s is a derived var, assuming unit '
                                    'as the first of its term', v)
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
                    logging.warning("Variable %s needed by %s is not "
                                    "available in the model output!", x, v)
    else:
        for v in var_needed:
            isavail[v] = False
            isunit[v] = None
        print(f'Not available: {var_needed} File: {infile}')
    return isavail, isunit


def load_yaml(infile):
    """Load generic yaml file"""
    try:
        with open(infile, 'r', encoding='utf-8') as file:
            cfg = yaml.load(file, Loader=yaml.FullLoader)
    except IOError:
        sys.exit(f'{infile} not found: you need to have this configuration file!')
    return cfg


def make_input_filename(var, year1, year2, face, diag):
    """Create input filenames for the required variable and a given year"""

    filetype = face[var]['filetype']
    filemask = face['model']['filetype'][filetype]['filename']
    filename = Path(os.path.expandvars(filemask.format(expname=diag.expname, 
                                                       year1=year1,
                                                       year2=year2,
                                                       var=var)))
    filedir = Path(os.path.expandvars(face['model']['basedir'].format(expname=diag.expname)),
                   os.path.expandvars(face['model']['filetype'][filetype]['dir'].format(expname=diag.expname)))
    return str(diag.ECEDIR / filedir / filename)


def units_extra_definition():
    """Add units to the pint registry required by ECMean4"""

    # special units definition, need to be moved in another placce
    units.define('fraction = [] = frac')
    units.define('psu = 1e-3 frac')
    units.define('Sv = 1e+6 m^3/s')  # Replace Sievert with Sverdrup


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

    return offset, factor


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


def getdomain(var, face):
    """Given a variable var extract its domain (ace or atm) from the interface"""
    component = face['model']['filetype'][face[var]['filetype']]['component']
    domain = face['model']['component'][component]['domain']
    return domain


def getcomponent(face):
    """Return a dictionary providing the component associated with each domain
       (the interface file specifies the domain for each component instead)"""
    d = face['model']['component']
    p = dict(zip([list(d.values())[x]['domain'] for x in range(len(d.values()))], d.keys()))
    return p


def getinifiles(face, comp, diag):
    """Return the inifiles from the interface, needs the component dictionary"""
    atminifile = os.path.expandvars(face['model']['component']
                                        [comp['atm']]['inifile'].format(expname=diag.expname))
    oceinifile = os.path.expandvars(face['model']['component']
                                        [comp['oce']]['inifile'].format(expname=diag.expname))
    if not atminifile[0] == '/':
        atminifile = str(diag.ECEDIR /
                         Path(os.path.expandvars(face['model']['basedir'].format(expname=diag.expname))) /
                         Path(atminifile))
    if not oceinifile[0] == '/':
        oceinifile = str(diag.ECEDIR /
                         Path(os.path.expandvars(face['model']['basedir'].format(expname=diag.expname))) /
                         Path(oceinifile))
    if not os.path.exists(oceinifile):
        oceinifile = ''
    return atminifile, oceinifile
