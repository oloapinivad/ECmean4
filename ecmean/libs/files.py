#!/usr/bin/env python3
'''
Shared functions for XArray ECmean4
'''

import os
import re
import logging
from pathlib import Path
from glob import glob
import yaml
import xarray as xr

loggy = logging.getLogger(__name__)

##################
# FILE FUNCTIONS #
##################


def inifiles_priority(inidict):
    """
    For areas dictionary and remap dictionary, provides a priority of which
    files to be used for interpolation and area computation
    Areas files comes first, then gridfile and finally land-sea mask.
    Provides flexibility for multiple models with different data access
    """
    file = None
    for key in ['areafile', 'gridfile', 'maskfile']:
        if inidict.get(key):
            file = inidict[key]
            break
    return file

def var_is_there(flist, var, face):
    """
    Check if a variable is available in the input file and provide its units.
    Returns:
        isavail (bool): if the variable is found or not
        varunit (string):  if the variable is there, its unit (None otherwise)
    """

    # we expect a list obtained by glob
    if not isinstance(flist, (xr.DataArray, xr.Dataset)):
        isavail = all(os.path.isfile(f) for f in flist) and len(flist) > 0
    else:
        # isavail = True if var in flist else False
        isavail = True

    if isavail:
        # no need of preproc here
        if not isinstance(flist, (xr.DataArray, xr.Dataset)):
            xfield = xr.open_mfdataset(flist)
        else:
            xfield = flist

        # if variable is derived, extract required vars
        var_req = _get_variables_to_load(var, face)

        # check if all required variables are in model output
        if set(var_req).issubset(set(xfield.data_vars)):
            units_avail = {}
            # if I don't know the unit, assuming is a fraction
            for i in xfield.data_vars:
                units_avail[i] = getattr(xfield[i], 'units', 'frac')
                # this is because I can't get rid of this unit
                if units_avail[i] == '(0 - 1)':
                    units_avail[i] = 'frac'
        else:
            isavail = False
            varunit = None
            x = [e for e in var_req if e not in xfield.data_vars]
            loggy.warning("Variable %s requires %s which is not "
                            "available in the model output. Ignoring it.", var, ' '.join(x))
            return isavail, varunit

        # get units
        varunit = face['variables'][var].get('units', None)
        if not varunit:
            varunit = units_avail.get(var_req[0])
            if len(var_req) > 1:
                loggy.debug('%s is a derived var, assuming unit '
                             'as the first of its term', var)

    else:
        varunit = None
        # print(f'Not available: {var} File: {flist}')
        loggy.error("No file found for variable %s. Ignoring it.", var)

    return isavail, varunit


def get_clim_files(piclim, var, diag, season):
    """Function to extra names for the climatology files"""

    # extract info from pi_climatology.yml
    # reference dataset and reference varname
    # as well as years when available
    dataref = piclim[var]['dataset']
    datayear1 = piclim[var].get('year1', None)
    datayear2 = piclim[var].get('year2', None)

    # get files for climatology
    if diag.climatology == 'RK08':
        dataname = piclim[var]['dataname']
        clim = str(diag.resclmdir / f'climate_{dataref}_{dataname}.nc')
        vvvv = str(diag.resclmdir / f'variance_{dataref}_{dataname}.nc')
    elif diag.climatology in 'EC23':
        if season == 'ALL':
            stringname = 'climate'
        else:
            stringname = 'seasons'
        clim = str(
            diag.resclmdir /
            f'{stringname}_average_{var}_{dataref}_{diag.resolution}_{datayear1}-{datayear2}.nc')
        vvvv = str(
            diag.resclmdir /
            f'{stringname}_variance_{var}_{dataref}_{diag.resolution}_{datayear1}-{datayear2}.nc')

    return clim, vvvv


def get_inifiles(face, diag):
    """Return the inifiles from the interface, needs the component dictionary.
    Check if inifiles exist.

    Args:
        face: the interface dictionary
        diag: the diagnostic object

    Returns:
        a dictionary with the different initial files

    """

    dictcomp = face['model']['component']

    ifiles = {}
    for comp in dictcomp.keys():
        ifiles[comp] = {}
        dictnames = face['component'][dictcomp[comp]]
        for name in dictnames.keys():
            inifile = dictnames[name]
            if inifile:
                if inifile[0] != '/':
                    inifile = Path(diag.ecedir) / \
                        Path(face['model']['basedir']) / \
                        Path(inifile)

                expandfile = _expand_filename(inifile, '', diag)
                filenames = glob(str(expandfile))
                if not filenames:
                    loggy.warning('Inifile %s cannot be found!', str(expandfile))
                    ifiles[comp][name] = ''
                else:
                    ifiles[comp][name] = str(_filter_filename_by_year(str(inifile), filenames, diag.year1)[0])

                loggy.info('%s for component %s is: %s', name, comp, ifiles[comp][name])

            else:
                ifiles[comp][name] = ''

    return ifiles


def _expand_filename(filenames, var, diag):
    """Expands a path (filename or dir) for var, expname, frequency, ensemble etc.
    and environment variables. Years are set as a wildcard and filtered by _filter_by_year"""

    return Path(str(os.path.expandvars(filenames)).format(
        expname=diag.expname,
        year1='*',
        year2='*',
        var=var,
        frequency=diag.frequency,
        ensemble=diag.ensemble,
        grid=diag.grid,
        model=diag.modelname,
        version=diag.version
    ))


def _filter_filename_by_year(template, filenames, year):
    """Find filename containing a given year in a list of filenames"""

    # if year1 is used in the file template
    if 'year1' in template:
        # Assumes that the file name ends with 199001-199012.nc or 1990-1991.nc
        year1 = [int(x.split('_')[-1].split('-')[0][0:4]) for x in filenames]
        # if year2 is used in the file template
        if 'year2' in template:
            year2 = [int(x.split('_')[-1].split('-')[1][0:4]) for x in filenames]
        else:
            year2 = year1
        # filter names
        filternames = [filenames[i] for i in range(len(year1))
                       if year >= year1[i] and year <= year2[i]]
    else:
        # this is introduced for file that does not have year in their filename
        filternames = filenames

    # safety warning if something is missing
    if not filternames and len(filenames) > 0:
        loggy.warning('Data for year %s has not been found!', str(year))

    loggy.debug('Filtered filenames: %s', filternames)
    return filternames


def load_yaml(infile):
    """Load generic yaml file"""

    if not os.path.isfile(infile):
        raise FileNotFoundError(f'ERROR: {infile} not found: you need to have this configuration file!')

    with open(infile, 'r', encoding='utf-8') as file:
        cfg = yaml.safe_load(file)

    if cfg is None:
        raise yaml.YAMLError(f'ERROR: An error occurred while parsing the file {infile}')

    return cfg


def _create_filepath(cmorname, face, diag):
    """Create filepath with wildcards"""

    filetype = face['variables'][cmorname]['filetype']
    filepath = Path(diag.ecedir) / \
        Path(face['model']['basedir']) / \
        Path(face['filetype'][filetype]['dir']) / \
        Path(face['filetype'][filetype]['filename'])
    loggy.debug('Filepath: %s', filepath)

    return filepath


def make_input_filename(cmorname, face, diag):
    """Create full input filepaths for the required variable and a given year

    Returns:
        a list of files to be loaded by xarray
    """

    # if a dataarray is provided
    if diag.xdataset is not None:
        return diag.xdataset

    # detect if it is a derived vars and get the list of var to be loaded
    varname = _get_variables_to_load(cmorname, face)

    # create filepath
    filepath = _create_filepath(cmorname, face, diag)

    # create the file structure according to the interface file
    filename = []
    for year in diag.years_joined:
        filename1 = []
        for var in varname:
            expandfile = _expand_filename(filepath, var, diag)
            filenames = glob(str(expandfile))
            fname = _filter_filename_by_year(str(filepath), filenames, year)
            filename1 = filename1 + fname
        filename = filename + filename1
    filename = list(dict.fromkeys(filename))  # Filter unique ones
    loggy.debug("Filenames: %s", filename)
    return filename


def _get_variables_to_load(var, face):
    """Function to extract from the interface file the list of derived variable,
    i.e. the real variables to be loaded, for each of the cmorname introduced in the
    interface file

    Args:
        var: the cmorname variable of the data to be loaded
        face: the interface file
    """

    if 'derived' in face['variables'][var].keys():
        cmd = face['variables'][var]['derived']
        dervars = [x for x in re.split('[*+-/]', cmd) if not x.isdigit()]
    else:
        dervars = [var]

    return dervars
