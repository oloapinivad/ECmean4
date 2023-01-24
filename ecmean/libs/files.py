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
from ecmean.libs.general import is_number, Diagnostic
from ecmean.libs.ncfixers import xr_preproc


##################
# FILE FUNCTIONS #
##################

def init_diagnostic(indir, argv) : 
    """
    configuration function to load config and interface files
    and to initilized the diagnotic object
    """
    
    # config file (looks for it in the same dir as the .py program file
    if argv.config:
        cfg = load_yaml(argv.config)
    else:
        cfg = load_yaml(indir / 'config.yml')

    # Setup all common variables, directories from arguments and config files
    logging.info(argv)
    diag = Diagnostic(argv, cfg)

    # loading the var-to-file interface
    # allow for both interface name or interface file
    fff, ext = os.path.splitext(diag.interface)
    if ext: 
        faceload = diag.interface
    else :  
        faceload = indir / Path(
            'interfaces',
            f'interface_{diag.interface}.yml')   
    face = load_yaml(faceload)

    return cfg, face, diag

def inifiles_priority(inidict) :

    """
    For areas dictionary and remap dictionary, provides a priority of which 
    files to be used for interpolation and area computation
    Areas files comes first, then gridfile and finally land-sea mask. 
    Provides flexibility for multiple models with different data access
    """ 

    if inidict['areafile']: 
        file = inidict['areafile']
    elif inidict['gridfile']: 
        file = inidict['gridfile']
    else:
        file = inidict['maskfile'] 

    return file

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
        ifiles[comp] ={}
        dictnames = face['component'][dictcomp[comp]]
        for name in dictnames.keys():
            inifile = dictnames[name]
            if inifile : 
                if inifile[0] != '/':
                    inifile = Path(diag.ECEDIR) / \
                        Path(face['model']['basedir']) / \
                        Path(inifile)
                    
                ifiles[comp][name] = str(_expand_filename(inifile, '', diag))
                logging.info(f'{name} for component {comp} is: {ifiles[comp][name]}')
               
                # safe check if inifile exist 
                if not glob(ifiles[comp][name]):
                    logging.error('Inifile %s cannot be found!', ifiles[comp][name])
                    ifiles[comp][name] = ''

            else: 
                ifiles[comp][name] = ''

    return ifiles

def get_inifiles_old(face, diag):
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
                        diag))
            else:
                inifiles[filename] = Path(diag.ECEDIR) / \
                    Path(face['model']['basedir']) / \
                    Path(inifile)
                inifiles[filename] = str(
                    _expand_filename(
                        inifiles[filename],
                        '',
                        diag))

            # safe check if inifile exist in the experiment folder
            if not glob(inifiles[filename]):
                inifiles[filename] = ''
        else:
            inifiles[filename] = ''

    # return dictionary values only
    return inifiles.values()


def _expand_filename(fn, var, diag):
    """Expands a path (filename or dir) for var, expname, frequency, ensemble etc.
    and environment variables. Years are set as a wildcard and filtered by _filter_by_year"""

    return Path(str(os.path.expandvars(fn)).format(
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
        filternames = [filenames[i] for i in range(len(year1)) if year >= year1[i] and year <= year2[i]]
    else :
        # this is introduced for file that does not have year in their filename
        filternames = filenames
    
    # safety warning if something is missing
    if not filternames and len(filenames)>0: 
        logging.warning('Data for year ' + str(year) + ' has not been found!')

    logging.info('Filtered filenames: %s', filternames)
    return filternames
        
  
def load_yaml(infile):
    """Load generic yaml file"""

    try:
        with open(infile, 'r', encoding='utf-8') as file:
            cfg = yaml.load(file, Loader=yaml.FullLoader)
    except IOError:
        sys.exit(f'ERROR: {infile} not found: you need to have this configuration file!')
    return cfg

def _create_filepath(cmorname, face, diag):
    """Create filepath with wildcards"""

    filetype = face['variables'][cmorname]['filetype']
    filepath = Path(diag.ECEDIR) / \
        Path(face['model']['basedir']) / \
        Path(face['filetype'][filetype]['dir']) / \
        Path(face['filetype'][filetype]['filename'])
    logging.info('Filepath: %s', filepath)

    return filepath

def make_input_filename(cmorname, varname, face, diag):
    """Create full input filepaths for the required variable and a given year
   
    Returns:
        a list of files to be loaded by xarray
    """

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
    logging.info("Filenames: %s", filename)
    return filename