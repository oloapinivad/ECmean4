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
from glob import glob
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
        self.debug = getattr(args, 'debug', False)
        self.numproc = args.numproc
        self.modelname = getattr(args, 'model', '')
        if not self.modelname:
            self.modelname = cfg['model']['name']
        if self.year1 == self.year2:  # Ignore if only one year requested
            self.ftrend = False
        #  These are here in prevision of future expansion to CMOR
        self.interface = cfg['interface']
        self.frequency = '*mon'
        self.ensemble = getattr(args, 'ensemble', 'r1i1p1f1')
        self.grid = '*'
        self.version = '*'

        # hard-coded resolution (due to climatological dataset)
        self.resolution = cfg['PI']['resolution']

        # Various input and output directories
        self.ECEDIR = Path(os.path.expandvars(cfg['dirs']['exp']))
        self.TABDIR = Path(os.path.expandvars(cfg['dirs']['tab']))
        self.CLMDIR = Path(os.path.expandvars(cfg['dirs']['clm']), self.resolution)
        self.years_joined = ''

        self.linefile = self.TABDIR / 'global_means.txt'

        # check if output attribute exists
        if hasattr(self, 'output'):
            self.linefile = args.output
            self.ftable = True


def chunks(iterable, num):
    """Generate num adjacent chunks of data from a list iterable
       Split lists in a convenient way for a parallel process"""
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


def var_is_there(infile, var, reference):
    """Check if a variable is available in the input file and provide its units."""

    # Use only first file if list is passed
    if infile:
        if isinstance(infile, list):
            ffile = infile[0]
        else:
            ffile = infile
    else:
        ffile = ''

    # File could still be a merge: find all files
    flist = [x for x in ffile.split() if x not in ['[', ']', '-merge']]

    isavail = True
    for f in flist:
        isavail = isavail and os.path.isfile(f)
    isavail = isavail and (len(flist) > 0)

    # if file exists, check which variables are inside
    if isavail:
        params = cdo.pardes(input=ffile)
        # Extract list of vars and of units in infile
        vars_avail = [v.split()[1] for v in params]
        # The following is a trick to obtain the units list (placing '' for missing ones)
        units_avail = [(v.replace('[', ']').split(']')[1:2] or [''])[0] for v in params]
        units_avail = dict(zip(vars_avail, units_avail))  # Transform into dict

        # if variable is derived, extract required vars
        d = reference[var].get('derived')
        if d:
            var_req = re.split('[*+-]', d)
            # remove numbers
            for x in var_req:
                if is_number(x):
                    var_req.remove(x)

            # check of unit is specified in the interface file
            varunit = reference[var].get('units')
            if not varunit:
                logging.warning('%s is a derived var, assuming unit '
                                'as the first of its term', var)
                varunit = units_avail.get(var_req[0])
        else:
            var_req = [var]
            varunit = units_avail.get(var)

        # check if all required variables are in model output
        isavail = True
        for x in var_req:
            if x not in vars_avail:
                isavail = False
                logging.warning("Variable %s needed by %s is not "
                                "available in the model output!", x, var)
    else:
        varunit = None
        print(f'Not available: {var} File: {infile}')
        logging.warning("Requested file %s is not available.", infile)

    return isavail, varunit


def load_yaml(infile):
    """Load generic yaml file"""
    try:
        with open(infile, 'r', encoding='utf-8') as file:
            cfg = yaml.load(file, Loader=yaml.FullLoader)
    except IOError:
        sys.exit(f'{infile} not found: you need to have this configuration file!')
    return cfg


def _expand_filename(fn, var, year1, year2, diag):
    """Expands a path (filename or dir) for var, expname, frequency, ensemble etc. and
       environment variables."""
    return Path(str(os.path.expandvars(fn)).format(
        expname=diag.expname,
        year1=year1,
        year2=year2,
        var=var,
        frequency=diag.frequency,
        ensemble=diag.ensemble,
        grid=diag.grid,
        model=diag.modelname,
        version=diag.version
    ))


def _filter_filename_by_year(fname, year):
    """Find filename containing a given year in a list of filenames"""
    filenames = glob(str(fname))
    # Assumes that the file name ends with 199001-199012.nc or 1990-1991.nc
    year1 = [int(x.split('_')[-1].split('-')[0][0:4]) for x in filenames]
    year2 = [int(x.split('_')[-1].split('-')[1][0:4]) for x in filenames]
    return [filenames[i] for i in range(len(year1)) if year >= year1[i] and year <= year2[i]]


def make_input_filename(var0, varlist, year1, year2, face, diag):
    """Create full input filepaths for the required variable and a given year"""

    filetype = face['variables'][var0]['filetype']
    filepath = Path(diag.ECEDIR) / \
        Path(face['model']['basedir']) / \
        Path(face['filetype'][filetype]['dir']) / \
        Path(face['filetype'][filetype]['filename'])
    # if year1 is a list, loop over it (we cannot use curly brackets anymore, now we pass a list)
    filename = []
    # Make an iterable even if year1 is not a list
    yy = year1
    if not isinstance(year1, list):
        yy = [year1]
    for year in yy:
        filename1 = []
        for var in varlist:
            fname = _expand_filename(filepath, var, '*', '*', diag)
            fname = _filter_filename_by_year(fname, year)
            filename1 = filename1 + fname
        filename1 = list(dict.fromkeys(filename1))  # Filter unique ones
        if len(filename1) <= 1:
            filename = filename + filename1
        else:
            filename = filename + ['-merge [ ' + ' '.join(filename1) + ' ]']
    if len(filename) == 1:  # glob always returns a list, return str if only one
        filename = filename[0]
    logging.debug("Filenames: %s", filename)
    return filename


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


def write_tuning_table(linefile, varmean, var_table, diag, ref):
    """Write results appending one line to a text file.
       Write a tuning table: need to fix reference to face/ref"""
    if not os.path.isfile(linefile):
        with open(linefile, 'w', encoding='utf-8') as f:
            print('%model  ens  exp from   to ', end='', file=f)
            for var in var_table:
                print('{:>12s}'.format(var), end=' ', file=f)
            print('\n%                         ', end=' ', file=f)
            for var in var_table:
                print('{:>12s}'.format(ref[var]['units']), end=' ', file=f)
            print(file=f)

    with open(linefile, 'a', encoding='utf-8') as f:
        print(f'{diag.modelname} {diag.ensemble} {diag.expname}',
              '{:4d} {:4d} '.format(diag.year1, diag.year2), end='', file=f)
        for var in var_table:
            print('{:12.5f}'.format(varmean[var] * ref[var].get('factor', 1)), end=' ', file=f)
        print(file=f)


def getdomain(var, face):
    """Given a variable var extract its domain (ace or atm) from the interface.
       To do so it creates a dictionary providing the domain associated with a component.
       (the interface file specifies the component for each domain instead)"""
    comp = face['filetype'][face['variables'][var]['filetype']]['component']
    d = face['model']['component']
    domain = dict(zip([list(d.values())[x] for x in range(len(d.values()))], d.keys()))
    return domain[comp]


def getcomponent(face):  # unused function
    """Return a dictionary providing the domain associated with a variable
       (the interface file specifies the domain for each component instead)"""
    d = face['component']
    p = dict(zip([list(d.values())[x]['domain'] for x in range(len(d.values()))], d.keys()))
    return p


def getinifiles(face, diag):
    """
    Return the inifiles from the interface, needs the component dictionary
    Check if inifiles exist.
    """
    dictcomp = face['model']['component']

    # use a dictionary to create the list of initial files
    inifiles = {}
    for comp, filename, filein in zip(['atm', 'oce', 'oce'],
                                      ['atminifile', 'ocegridfile', 'oceareafile'],
                                      ['inifile', 'gridfile', 'areafile']):

        inifile = face['component'][dictcomp[comp]].get(filein, '')

        # add the full path if missing
        inifiles[filename] = ''
        if inifile:
            if inifile[0] == '/':
                inifiles[filename] = str(_expand_filename(inifile,
                                                          '', diag.year1, diag.year1, diag))
            else:
                inifiles[filename] = Path(diag.ECEDIR) / \
                    Path(face['model']['basedir']) / \
                    Path(inifile)
                inifiles[filename] = str(_expand_filename(inifiles[filename],
                                                          '', diag.year1, diag.year1, diag))

            # safe check if inifile exist in the experiment folder
            if not glob(inifiles[filename]):
                inifiles[filename] = ''
        else:
            inifiles[filename] = ''

    # return dictionary values only
    return inifiles.values()
