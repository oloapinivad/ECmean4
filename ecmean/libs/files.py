#!/usr/bin/env python3
'''
Shared functions for XArray ECmean4
'''

import os
import re
import xarray as xr
from pathlib import Path
import logging
from glob import glob
import yaml
import sys
from ecmean.libs.general import is_number
from ecmean.libs.ncfixers import xr_preproc


##################
# FILE FUNCTIONS #
##################


def var_is_there(flist, var, reference):
    """Check if a variable is available in the input file and provide its units"""

    # we expect a list obtained by glob
    isavail = True
    for f in flist:
        isavail = isavail and os.path.isfile(f)
    isavail = isavail and (len(flist) > 0)

    if isavail:
        xfield = xr.open_mfdataset(flist, preprocess=xr_preproc)
        # vars_avail = [i for i in xfield.data_vars]
        vars_avail = xfield.data_vars
        units_avail = {}
        # if I don't know the unit, assuming is a fraction
        for i in vars_avail:
            try:
                k = xfield[i].units
            except BaseException:
                k = 'frac'
            units_avail[i] = k

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
                logging.info('%s is a derived var, assuming unit '
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
                logging.warning("Variable %s requires %s which is not "
                                "available in the model output. Ignoring it.", var, x)
    else:
        varunit = None
        # print(f'Not available: {var} File: {flist}')
        logging.warning("No data found for variable %s. Ignoring it.", var)

    return isavail, varunit


def get_clim_files(piclim, var, diag, season):
    """Function to extra names for the climatology files"""

    # extract info from pi_climatology.yml
    # reference dataset and reference varname
    # as well as years when available
    dataref = piclim[var]['dataset']
    datayear1 = piclim[var].get('year1', 'nan')
    datayear2 = piclim[var].get('year2', 'nan')

    # get files for climatology
    if diag.climatology == 'RK08':
        dataname = piclim[var]['dataname']
        clim = str(diag.RESCLMDIR / f'climate_{dataref}_{dataname}.nc')
        vvvv = str(diag.RESCLMDIR / f'variance_{dataref}_{dataname}.nc')
    elif diag.climatology in 'EC22':
        dataname = piclim[var]['dataname']
        clim = str(
            diag.RESCLMDIR /
            f'climate_{dataname}_{dataref}_{diag.resolution}_{datayear1}-{datayear2}.nc')
        vvvv = str(
            diag.RESCLMDIR /
            f'variance_{dataname}_{dataref}_{diag.resolution}_{datayear1}-{datayear2}.nc')
    elif diag.climatology in 'EC23':
        if season == 'ALL':
            clim = str(
                diag.RESCLMDIR /
                f'climate_average_{var}_{dataref}_{diag.resolution}_{datayear1}-{datayear2}.nc')
            vvvv = str(
                diag.RESCLMDIR /
                f'climate_variance_{var}_{dataref}_{diag.resolution}_{datayear1}-{datayear2}.nc')
        else:
            clim = str(
                diag.RESCLMDIR /
                f'seasons_average_{var}_{dataref}_{diag.resolution}_{datayear1}-{datayear2}.nc')
            vvvv = str(
                diag.RESCLMDIR /
                f'seasons_variance_{var}_{dataref}_{diag.resolution}_{datayear1}-{datayear2}.nc')

    return clim, vvvv


def get_inifiles(face, diag):
    """Return the inifiles from the interface, needs the component dictionary.
    Check if inifiles exist."""

    dictcomp = face['model']['component']

    # use a dictionary to create the list of initial files
    inifiles = {}
    for comp, filename, filein in zip(['atm', 'atm', 'oce'],
                                      ['maskatmfile', 'atmareafile', 'oceareafile'],
                                      ['inifile', 'atmfile', 'areafile']):

        inifile = face['component'][dictcomp[comp]].get(filein, '')

        # add the full path if missing
        inifiles[filename] = ''
        if inifile:
            if inifile[0] == '/':
                inifiles[filename] = str(
                    _expand_filename(
                        inifile,
                        '',
                        diag.year1,
                        diag.year1,
                        diag))
            else:
                inifiles[filename] = Path(diag.ECEDIR) / \
                    Path(face['model']['basedir']) / \
                    Path(inifile)
                inifiles[filename] = str(
                    _expand_filename(
                        inifiles[filename],
                        '',
                        diag.year1,
                        diag.year1,
                        diag))

            # safe check if inifile exist in the experiment folder
            if not glob(inifiles[filename]):
                inifiles[filename] = ''
        else:
            inifiles[filename] = ''

    # return dictionary values only
    return inifiles.values()


def _expand_filename(fn, var, year1, year2, diag):
    """Expands a path (filename or dir) for var, expname, frequency, ensemble etc.
    and environment variables."""

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
    try:
        year2 = [int(x.split('_')[-1].split('-')[1][0:4]) for x in filenames]
    except IndexError:
        # this is introduced to handle files which have only one year in their filename
        year2 = year1
    return [filenames[i]
            for i in range(len(year1)) if year >= year1[i] and year <= year2[i]]


def load_yaml(infile):
    """Load generic yaml file"""

    try:
        with open(infile, 'r', encoding='utf-8') as file:
            cfg = yaml.load(file, Loader=yaml.FullLoader)
    except IOError:
        sys.exit(f'ERROR: {infile} not found: you need to have this configuration file!')
    return cfg


def make_input_filename(var0, varlist, year1, year2, face, diag):
    """Create full input filepaths for the required variable and a given year"""

    filetype = face['variables'][var0]['filetype']
    filepath = Path(diag.ECEDIR) / \
        Path(face['model']['basedir']) / \
        Path(face['filetype'][filetype]['dir']) / \
        Path(face['filetype'][filetype]['filename'])
    # if year1 is a list, loop over it (we cannot use curly brackets anymore,
    # now we pass a list)
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
        # filename1 = list(dict.fromkeys(filename1))
        filename = filename + filename1
    filename = list(dict.fromkeys(filename))  # Filter unique ones
    # print(filename)
    logging.debug("Filenames: %s", filename)
    return filename