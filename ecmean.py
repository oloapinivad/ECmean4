#!/usr/bin/env python3
'''
Shared functions for XArray ECmean4
'''

import os
import re
import logging
import operator
import sys
from pathlib import Path
from glob import glob
import itertools
import numpy as np
import xarray as xr
import xesmf as xe
from metpy.units import units
import yaml
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.colors import TwoSlopeNorm


####################
# DIAGNOSTIC CLASS #
####################


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
        self.climatology = getattr(args, 'climatology', 'RK08')
        self.interface = getattr(args, 'interface', '')
        self.resolution = getattr(args, 'resolution', '')
        self.regions = cfg['PI']['regions']
        self.seasons = cfg['PI']['seasons']
        if not self.modelname:
            self.modelname = cfg['model']['name']
        if self.year1 == self.year2:  # Ignore if only one year requested
            self.ftrend = False
        #  These are here in prevision of future expansion to CMOR
        if not self.interface:
            self.interface = cfg['interface']
        self.frequency = '*mon'
        self.ensemble = getattr(args, 'ensemble', 'r1i1p1f1')
        self.grid = '*'
        self.version = '*'

        # hard-coded resolution (due to climatological dataset)
        if self.climatology == 'RK08':
            logging.error('RK08 can work only with r180x91 grid')
            self.resolution = 'r180x91'
        else:
            if not self.resolution:
                self.resolution = cfg['PI']['resolution']

        # hard-coded seasons (due to climatological dataset)
        if self.climatology in ['EC22', 'RK08']:
            logging.error('only EC23 climatology support multiple seasons! Keeping only yearly seasons!')
            self.seasons = ['ALL']

        # Various input and output directories
        self.ECEDIR = Path(os.path.expandvars(cfg['dirs']['exp']))
        self.TABDIR = Path(os.path.expandvars(cfg['dirs']['tab']))
        self.FIGDIR = Path(os.path.expandvars(cfg['dirs']['fig']))
        self.CLMDIR = Path(
            os.path.expandvars(
                cfg['dirs']['clm']),
            self.climatology)
        self.RESCLMDIR = Path(self.CLMDIR, self.resolution)
        self.years_joined = ''

        if hasattr(args, 'output') and args.output:
            self.linefile = args.output
            self.ftable = True
        else:
            self.linefile = self.TABDIR / 'global_means.txt'


##################
# HELP FUNCTIONS #
##################


def is_number(s):
    """Check if input is a float type"""

    try:
        float(s)
        return True
    except ValueError:
        return False


# def chunks(iterable, num):
#     """Generate num adjacent chunks of data from a list iterable
#     Split lists in a convenient way for a parallel process"""

#     size = int(np.ceil(len(iterable) / num))
#     it = iter(iterable)
#     return iter(lambda: tuple(itertools.islice(it, size)), ())

# def split(iterable, num):
#     k, m = divmod(len(iterable), num)
#     return (iterable[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(num))

def runtime_weights(varlist) : 
    """Define the weights to estimate the best repartition of the cores
    This is done a-priori, considering that 1) compound variables are more difficult to compute 
    2) 3d variables requires more evaluation"""

    w = {}
    for k in varlist : 
        if k in ['ua', 'ta', 'va', 'hus'] :
            t = 8
        elif k in ['pme', 'net_sfc_nosn', 'net_sfc', 'toamsfc_nosn', 'toamsfc', 
            'pr_oce', 'pme_oce', 'pr_land', 'pme_land', 'net_toa'] :
            t = 3
        else :
            t = 1
        w[k] = t 

    return w


def weight_split(a, n) : 
    """use the weights by runtime_weights to provide a number of n chunks based on the
    minimum possible value of each chunk computing time. Then, provide the list of variables
    to be used in the multiprocessing routine"""

    ws = sorted(runtime_weights(a).items(), key=lambda x:x[1], reverse=True)
    ordered = dict(ws)

    elists = [ [0] for _ in range(n) ]
    olists = [ [] for _ in range(n) ]
    
    count=0
    for f in ordered.keys() : 
        elists[count].append(ordered[f])
        olists[count].append(f)
        count=elists.index(min(elists)) 
        
    return olists


def getdomain(var, face):
    """Given a variable var extract its domain (oce or atm) from the interface.
    To do so it creates a dictionary providing the domain associated with a component.
    (the interface file specifies the component for each domain instead)"""

    comp = face['filetype'][face['variables'][var]['filetype']]['component']
    d = face['model']['component']
    domain = dict(zip([list(d.values())[x]
                  for x in range(len(d.values()))], d.keys()))
    return domain[comp]


def getcomponent(face):  # unused function
    """Return a dictionary providing the domain associated with a variable
    (the interface file specifies the domain for each component instead)"""

    d = face['component']
    p = dict(zip([list(d.values())[x]['domain']
             for x in range(len(d.values()))], d.keys()))
    return p


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
        #print(f'Not available: {var} File: {flist}')
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


##########################
# MASK-RELATED FUNCTIONS #
##########################


def masks_dictionary(component, maskatmfile, remap_dictionary=None):
    """Create a dictionary with atmospheric land-sea mask"""

    mask = {
        'atm_mask': _make_atm_masks(
            component['atm'],
            maskatmfile,
            remap_dictionary=remap_dictionary),
    }

    return mask


def _make_atm_masks(component, maskatmfile, remap_dictionary=None):
    """Create land-sea masks for atmosphere model"""

    # prepare ATM LSM: this needs to be improved, since it is clearly model
    # dependent
    logging.debug('maskatmfile is' + maskatmfile)
    if not maskatmfile: 
        sys.exit("ERROR: maskatmfile cannot be found")

    if component == 'oifs':
        # create mask: opening a grib and loading only lsm to avoid
        # inconsistencies # in the grib structure ->
        # see here https://github.com/ecmwf/cfgrib/issues/13
        mask = xr.open_mfdataset(
            maskatmfile,
            engine="cfgrib",
            filter_by_keys={
                'shortName': 'lsm'},
            preprocess=xr_preproc)
        mask = mask['lsm']
    elif component == 'cmoratm':
        mask = xr.open_mfdataset(maskatmfile, preprocess=xr_preproc)
        mask = mask['sftlf']
    elif component == 'globo':
        mask = xr.open_mfdataset(maskatmfile, preprocess=xr_preproc)
        mask = mask['lsm'].mean(dim='time')
        mask = abs(1-mask)
    else:
        sys.exit("ERROR: Mask undefined yet mismatch, this cannot be handled!")

    if remap_dictionary is not None:
        if remap_dictionary['atm_fix']:
            mask = remap_dictionary['atm_fix'](mask, keep_attrs=True)
        mask = remap_dictionary['atm_remap'](mask, keep_attrs=True)

    return mask


def masked_meansum(xfield, var, weights, mask_type, mask):
    """For global variables rvaluate the weighted averaged
    or weighted integral when required by the variable properties"""

    # call the mask_field to mask where necessary
    masked = mask_field(xfield, var, mask_type, mask)

    if mask_type in ['global']:
        out = masked.weighted(weights.fillna(0)).mean().values
    elif mask_type in ['land', 'ocean', 'sea']:
        out = masked.weighted(weights.fillna(0)).sum().values
    else:
        sys.exit("ERROR: Mask undefined, this cannot be handled!")

    return float(out)


def mask_field(xfield, var, mask_type, mask):
    """Apply a land/sea mask on a xarray variable var"""

    # nothing to be done
    if mask_type == 'global':
        out = xfield

    else:

        # check that we are receiving a dataset and not a datarray
        if isinstance(xfield, xr.DataArray):
            xfield = xfield.to_dataset(name=var)

        # convert from datarray to dataset and merge
        mask = mask.to_dataset(name='mask')

        # the compat='override' option forces the merging. some CMIP6 data might
        # have different float type, this simplies the handling
        bfield = xr.merge([xfield, mask], compat='override')

        # conditions
        if mask_type == 'land':
            out = bfield[var].where(bfield['mask'] >= 0.5)
        elif mask_type in ['sea', 'ocean']:
            out = bfield[var].where(bfield['mask'] < 0.5)
        else:
            sys.exit("ERROR: Mask undefined, this cannot be handled!")

    return out


def select_region(xfield, region):
    """Trivial function to convert region definition to xarray
    sliced array to compute the PIs or global means on selected regions"""

    if region == 'Global':
        slicearray = xfield
    elif region == 'North Midlat':
        slicearray = xfield.sel(lat=slice(30, 90))
    elif region == 'South Midlat':
        slicearray = xfield.sel(lat=slice(-90, -30))
    elif region == 'Tropical':
        slicearray = xfield.sel(lat=slice(-30, 30))
    else:
        sys.exit(region + "region not supported!!!")

    return slicearray


##################################
# AREA-WEIGHT AND MASK FUNCTIONS #
##################################


def areas_dictionary(component, atmareafile, oceareafile):
    """Create a dictionary with atmospheric and oceanic area weights"""

    areas = {
        'atm_areas': _make_atm_areas(component['atm'], atmareafile),
        'oce_areas': _make_oce_areas(component['oce'], oceareafile)
    }

    return areas


def _make_atm_areas(component, atmareafile):
    "Create atmospheric weights for area operations"

    logging.debug('Atmareafile is ' + atmareafile)
    if not atmareafile : 
        sys.exit("ERROR: Atmareafile cannot be found")

    if component == 'oifs':
        xfield = xr.open_mfdataset(atmareafile, preprocess=xr_preproc)
        area = _area_cell(xfield)
    elif component == 'cmoratm':
        xfield = xr.open_mfdataset(atmareafile, preprocess=xr_preproc)
        area = _area_cell(xfield)
    elif component == 'globo':
        xfield = xr.open_mfdataset(atmareafile, preprocess=xr_preproc)
        area = _area_cell(xfield)
    else:
        sys.exit("ERROR: Area for this configuration cannot be handled!")
    return area


def _make_oce_areas(component, oceareafile):
    "Create atmospheric weights for area operations"

    logging.debug('Oceareafile is ' + oceareafile)
    if not oceareafile : 
        logging.warning("Ocereafile cannot be found, assuming this is an AMIP run")

    if oceareafile:
        if component == 'nemo':
            xfield = xr.open_mfdataset(oceareafile, preprocess=xr_preproc)
            if 'e1t' in xfield.data_vars:
                area = xfield['e1t'] * xfield['e2t']
            else:
                area = _area_cell(xfield)
        elif component == 'cmoroce':
            xfield = xr.open_mfdataset(oceareafile, preprocess=xr_preproc)
            if 'areacello' in xfield.data_vars:
                area = xfield['areacello']
            else:
                area = _area_cell(xfield)
        else:
            sys.exit("ERROR: Area for this configuration cannot be handled!")
    else:
        area = None
    return area


def guess_bounds(axis, name='lon'):
    """Basic function that estimates the boundaries for lon and lat if they are not
    available. Works only with regular grids.
    It also avoid having values larger than 90N/90S as well as 0 and 10^5 Pa"""
    # inspired by
    # https://gist.github.com/dennissergeev/60bf7b03443f1b2c8eb96ce0b1880150

    # this define the proportion of the bounds, assumed to half of the levels
    bound_position = 0.5
    diffs = np.diff(axis)
    diffs = np.insert(diffs, 0, diffs[0])
    diffs = np.append(diffs, diffs[-1])

    # pair of dounds
    min_bounds = axis - diffs[:-1] * bound_position
    max_bounds = axis + diffs[1:] * (1 - bound_position)

    # safety check, to be generalized
    if name in 'lat':
        max_bounds[-1] = 90
        min_bounds[0] = (-90)
    if name in 'plev':
        min_bounds[0] = 0
        max_bounds[-1] = 100000

    # should we use a xarray object instead of a numpy?
    bounds = np.array([min_bounds, max_bounds]).transpose()
    return bounds


def _lonlat_to_sphere(lon, lat):
    """Convert from lon lat coordinates to a 3d sphere of unity radius"""

    vec = np.array([
        np.cos(np.deg2rad(lon)) * np.cos(np.deg2rad(lat)),
        np.sin(np.deg2rad(lon)) * np.cos(np.deg2rad(lat)),
        np.sin(np.deg2rad(lat))
    ])
    return vec


def _huilier(a, b, c):
    """Apply the L'Huilier theorem from the three side of the spherical triangle
    obtaining the spherical excess, i.e. the solid angle of the triangle,
    i.e. the area of the spherical surface"""
    # More info at https://mathworld.wolfram.com/LHuiliersTheorem.html

    s = (a + b + c) * 0.5
    t = np.tan(s * 0.5) * np.tan((s - a) * 0.5) * \
        np.tan((s - b) * 0.5) * np.tan((s - c) * 0.5)
    area = abs(4. * np.arctan(np.sqrt(abs(t))))
    return area


def _vector_spherical_triangle(p1, p2, p3):
    """Given the coordinates of three points on a sphere, estimate the length
    of their vectors a,b,c connecting to the centre of the spere.
    Then using L'Huilier formula derive the solid angle among the three
    vectors, which is multiplied by squared Earth Radius is
    exactly the surface of the corresponding spherical triangle. """
    # This is inspired by CDO code found at
    # https://code.mpimet.mpg.de/projects/cdo/repository/cdo/revisions/331ab3f7fd18295cf6a433fb799034c7589a4a61/entry/src/grid_area.cc

    a = np.arcsin(np.linalg.norm(np.cross(p1, p2), axis=1))
    b = np.arcsin(np.linalg.norm(np.cross(p1, p3), axis=1))
    c = np.arcsin(np.linalg.norm(np.cross(p3, p2), axis=1))

    area = _huilier(a, b, c)
    return area


def _area_cell(xfield, formula='triangles'):
    """
    Function which estimate the area cell from bounds. This is done assuming
    making use of spherical triangels.
    Working also on regular grids which does not have lon/lat bounds
    via the guess_bounds function. Curvilinear/unstructured grids are not supported,
    especially if with more with more than 4 vertices are not supported.

    Args:
    xfield: a generic xarray dataset
    formula: 'squares' or 'trapezoids' or 'triangles' equation for the area cell
        'triangles' is the default, uses the spherical triangles - same as
        used by CDO - and it is very accurate

    Returns:
    An xarray dataarray with the area for each grid point
    """

    earth_radius = 6371000.

    # some check to starts
    if all(x in xfield.dims for x in ['lon', 'lat']):
        logging.debug('Regulard grid recognized..')
        regular_grid = True
    else:
        regular_grid = False

    if all(x in xfield.data_vars for x in ['lon_bnds', 'lat_bnds']):
        logging.debug('cmor lon/lat_bounds found...')
        cmor_bounds = True
    else:
        cmor_bounds = False

    # this is a nightmare, so far working only for ECE4 gaussian reduced
    if not regular_grid:

        logging.debug('Curvilinear/Unstructured grid, tryin to get grid info...')

        blondim = None
        blatdim = None
        # trying to find bounderies
        for g in (
            t for t in list(
                xfield.data_vars) if t in [
                'lon_bnds',
                'bounds_lon']):
            blondim = g
        for g in (
            t for t in list(
                xfield.data_vars) if t in [
                'lat_bnds',
                'bounds_lat']):
            blatdim = g

        # checking
        if blondim is None and blatdim is None:
            sys.exit(
                "ERROR: Can't find any lon/lat boundaries and grid is unstructured, need some help!")

        logging.debug('Unstructured grid, special ECE4 treatment...')
        # ATTENTION: this is a very specific ECE4 definition, it will not work
        # with other unstructured grids. The assumption of the vertex position
        # is absolutely random. Needs to be generalized.
        bounds_lon = np.column_stack((xfield[blondim].isel(nvertex=1),
                                      xfield[blondim].isel(nvertex=2)))
        bounds_lat = np.column_stack((xfield[blatdim].isel(nvertex=2),
                                      xfield[blatdim].isel(nvertex=3)))
        area_dims = 'cell'
        # set full lat
        full_lat = xfield['lat'].values

    # if we are dealing with a regular grid
    if regular_grid:

        # if we have bounds, just check they have the right dimensions names
        if cmor_bounds:

            # if dimension is not called bnds, rename it
            if 'bnds' not in list(xfield.dims):
                logging.debug('bnds not found, trying to rename it...')
                for g in (
                    t for t in list(
                        xfield.dims) if t not in [
                        'lon',
                        'lat',
                        'time']):
                    bdim = g
                xfield = xfield.rename_dims({bdim: "bnds"})

        # else use guess_bounds() and expand the xarray dataset including them
        if not cmor_bounds:

            logging.debug('Bounds estimation from lon/lat...')
            # create and xarray dataset which the boundaries
            xbounds = xr.Dataset(
                data_vars=dict(
                    lat_bnds=(
                        ('lat', 'bnds'), guess_bounds(
                            xfield['lat'], name='lat')), lon_bnds=(
                        ('lon', 'bnds'), guess_bounds(
                            xfield['lon'], name='lon'))), coords=dict(
                    lat=(
                        'lat', xfield['lat'].values), lon=(
                        'lon', xfield['lon'].values)))
            xfield = xfield.merge(xbounds)

        # create numpy array
        blon = np.column_stack((xfield['lon_bnds'].isel(bnds=0),
                                xfield['lon_bnds'].isel(bnds=1)))
        blat = np.column_stack((xfield['lat_bnds'].isel(bnds=0),
                                xfield['lat_bnds'].isel(bnds=1)))
        full_lat = np.repeat(xfield['lat'].values, len(xfield['lon']), axis=0)

        # 2d matrix of bounds
        expansion = np.array([(y, x) for y in blat for x in blon])
        bounds_lon = expansion[:, 1, :]
        bounds_lat = expansion[:, 0, :]
        area_dims = ('lat', 'lon')

    # cell dimension
    if formula == "triangles":
        p1 = _lonlat_to_sphere(bounds_lon[:, 0], bounds_lat[:, 0]).transpose()
        p2 = _lonlat_to_sphere(bounds_lon[:, 0], bounds_lat[:, 1]).transpose()
        p3 = _lonlat_to_sphere(bounds_lon[:, 1], bounds_lat[:, 1]).transpose()
        p4 = _lonlat_to_sphere(bounds_lon[:, 1], bounds_lat[:, 0]).transpose()
        area_cell = _vector_spherical_triangle(
            p1, p2, p3) + _vector_spherical_triangle(p1, p4, p3)
        area_cell = area_cell * earth_radius**2
    else:
        # cell dimension
        dlon = abs(bounds_lon[:, 0] - bounds_lon[:, 1])
        dlat = abs(bounds_lat[:, 0] - bounds_lat[:, 1])

        # safe check on cosine of 90 included:
        # assume a trapezoid or a squared cell
        if formula == 'trapezoids':
            arclon1 = earth_radius * \
                abs(np.cos(abs(np.deg2rad(bounds_lat[:, 0])))) * np.deg2rad(dlon)
            arclon2 = earth_radius * \
                abs(np.cos(abs(np.deg2rad(bounds_lat[:, 1])))) * np.deg2rad(dlon)
        if formula == 'squares':
            full_lat = np.repeat(
                xfield['lat'].values, len(
                    xfield['lon']), axis=0)
            arclon1 = arclon2 = earth_radius * \
                abs(np.cos(abs(np.deg2rad(full_lat)))) * np.deg2rad(dlon)

        arclat = earth_radius * np.deg2rad(dlat)

        # trapezoid area
        area_cell = (arclon1 + arclon2) * arclat / 2

    if regular_grid:
        area_cell = area_cell.reshape([len(xfield['lat']), len(xfield['lon'])])

    # since we are using numpy need to bring them back into xarray dataset
    xfield['area'] = (area_dims, area_cell)

    # check the total area
    logging.debug('Total Earth Surface: %s Km2',
                  str(xfield['area'].sum().values / 10**6))

    return xfield['area']


#####################
# FORMULA FUNCTIONS #
#####################


# this is a tool to parse a CDO-based formula into mathematical operatos
# there might exists something more intelligent such as the pyparsing package

def eval_formula(mystring, xdataset):
    """Evaluate the cmd string provided by the yaml file
    producing a parsing for the derived variables"""

    # Tokenize the original string
    token = [i for i in re.split('(\\W+)', mystring) if i]
    if len(token) > 1:
        # Use order of operations
        out = _operation(token, xdataset)
    else:
        out = xdataset[token[0]]
    return out


def _operation(token, xdataset):
    """Parsing of the CDO-based commands using operator package
    and an ad-hoc dictionary. Could be improved, working with four basic
    operations only."""

    # define math operators: order is important, since defines
    # which operation is done at first!
    ops = {
        '/': operator.truediv,
        "*": operator.mul,
        "-": operator.sub,
        "+": operator.add
    }

    # use a dictionary to store xarray field and call them easily
    dct = {}
    for k in token:
        if k not in ops:
            if not is_number(k):
                dct[k] = xdataset[k]
            else:
                dct[k] = float(k)

    # apply operators to all occurrences, from top priority
    # so far this is not parsing parenthesis
    code = 0
    for p in ops:
        while p in token:
            code += 1
            # print(token)
            x = token.index(p)
            name = 'op' + str(code)
            replacer = ops.get(p)(dct[token[x - 1]], dct[token[x + 1]])
            dct[name] = replacer
            token[x - 1] = name
            del token[x:x + 2]
    return replacer


###########################
# INTERPOLATION FUNCTIONS #
###########################


def remap_dictionary(component, atmareafile, oceareafile, target_grid):
    """Create a dicitionary with atmospheric and oceanic weights for
    interpolation. There is an option of fix grid before the real
    interpolation: this is used for Gaussian reduced grids"""

    atmfix, atmremap = _make_atm_interp_weights(
        component['atm'], atmareafile, target_grid)
    if oceareafile:
        ocefix, oceremap = _make_oce_interp_weights(
            component['oce'], oceareafile, target_grid)
    else:
        ocefix = None
        oceremap = None

    remap = {
        'atm_fix': atmfix,
        'atm_remap': atmremap,
        'oce_fix': ocefix,
        'oce_remap': oceremap,
    }

    return remap


def _make_atm_interp_weights(component, atmareafile, target_grid):
    """Create atmospheric interpolator"""

    logging.debug('Atmareafile is ' + atmareafile)
    if not atmareafile : 
        sys.exit("ERROR: Atmareafile cannot be found")

    if component == 'oifs':

        # this is to get lon and lat from the Equator
        xfield = xr.open_mfdataset(atmareafile, preprocess=xr_preproc).load()
        xname = list(xfield.data_vars)[-1]
        m = xfield[xname].isel(time=0).load()
        g = sorted(list(set(m.lat.values)))
        f = sorted(list(m.sel(cell=m.lat == g[int(len(g) / 2)]).lon.values))

        # this creates a a gaussian non reduced grid
        ds_out = xr.Dataset({"lon": (["lon"], f), "lat": (["lat"], g)})

        # use nearest neighbour to remap to gaussian regular
        fix = xe.Regridder(
            xfield[xname],
            ds_out,
            method="nearest_s2d",
            locstream_in=True,
            periodic=True)

        # create bilinear interpolator
        interp = xe.Regridder(
            fix(xfield[xname]), target_grid, periodic=True, method="bilinear")

    elif component == 'cmoratm':

        fix = None
        xfield = xr.open_mfdataset(atmareafile, preprocess=xr_preproc).load()
        interp = xe.Regridder(
            xfield,
            target_grid,
            periodic=True,
            method="bilinear")

    elif component == 'globo':

        fix = None
        xfield = xr.open_mfdataset(atmareafile, preprocess=xr_preproc).load()
        interp = xe.Regridder(
            xfield,
            target_grid,
            periodic=True,
            ignore_degenerate=True,
            method="bilinear")

    else:
        sys.exit(
            "ERROR: Atm weights not defined for this component, this cannot be handled!")

    return fix, interp


def _make_oce_interp_weights(component, oceareafile, target_grid):
    """Create oceanic interpolator weights"""

    logging.debug('Oceareafile is ' + oceareafile)
    if not oceareafile : 
        sys.exit("ERROR: Oceareafile cannot be found")

    if component == 'nemo':
        fix = None
        xfield = xr.open_mfdataset(oceareafile, preprocess=xr_preproc).load()

        # set coordinates which are missing
        for cl in ['nav_lon', 'nav_lat', 'nav_lev', 'time_counter', 'x', 'y']:
            if cl in xfield.data_vars:
                xfield = xfield.set_coords([cl])

        # rename dimensions and coordinates
        xfield = xfield.rename(
            {"nav_lon": "lon", "nav_lat": "lat", "nav_lev": "deptht"})

        # use grid distance as generic variable
        interp = xe.Regridder(
            xfield['e1t'],
            target_grid,
            method="bilinear",
            periodic=True,
            ignore_degenerate=True)

    elif component == 'cmoroce':

        fix = None
        xfield = xr.open_mfdataset(oceareafile, preprocess=xr_preproc)
        xname = list(xfield.data_vars)[-1]
        # print(len(xfield.coords['lon'].shape))

        # check if oceanic grid is regular: lon/lat dims should be 1d
        # if not all(x in xfield.dims for x in ['lon', 'lat']) and (len(xfield.dims) < 3) :
        if len(xfield.coords['lon'].shape) == 1 and len(xfield.coords['lat'].shape) == 1:

            print("Detecting a unstructured grid, using nearest neighbour!")
            interp = xe.Regridder(
                xfield[xname].load(),
                target_grid,
                method="nearest_s2d",
                locstream_in=True,
                periodic=True)
        else:
            print("Detecting regular or curvilinear grid, using bilinear!")
            interp = xe.Regridder(
                xfield[xname].load(),
                target_grid,
                method="bilinear",
                ignore_degenerate=True,
                periodic=True)

    else:
        sys.exit(
            "ERROR: Oce weights not defined for this component, this cannot be handled!")

    return fix, interp


#########################
# FILE FORMAT FUNCTIONS #
#########################

def xr_preproc(ds):
    """Preprocessing functuon to adjust coordinate and dimensions
    names to a common format. To be called by xr.open_mf_dataset()"""

    # print(ds)
    if 'time_counter' in list(ds.dims):
        ds = ds.rename({"time_counter": "time"})

    if 'time_counter' in list(ds.coords):
        ds = ds.rename({"time_counter": "time"})

    if 'pressure_levels' in list(ds.coords):
        ds = ds.rename({"pressure_levels": "plev"})

    if 'plevel' in list(ds.dims):
        ds = ds.rename({"plevel": "plev"})

    # fix for NEMO eORCA grid (nav_lon, nav_lat)
    for h in ['lon', 'lat'] : 
        for f in ['', 'grid_T'] :
            g = 'nav_'+ h + '_' + f 
            if g in list(ds.coords):
                ds = ds.rename({g: h})

    # fix for NEMO eORCA grid (x_grid_T, etc.)
    for h in ['x', 'y'] :
        for f in ['grid_T'] :
            g = h + '_'+ f 
            if g in list(ds.dims):
                ds = ds.rename({g: h})

    if 'longitude' in list(ds.dims):
        ds = ds.rename({"longitude": "lon"})

    if 'latitude' in list(ds.dims):
        ds = ds.rename({"latitude": "lat"})

    if 'longitude' in list(ds.coords):
        ds = ds.rename({"longitude": "lon"})

    if 'nav_lon' in list(ds.coords):
        ds = ds.rename({"nav_lon": "lon"})

    if 'latitude' in list(ds.coords):
        ds = ds.rename({"latitude": "lat"})

    if 'nav_lat' in list(ds.coords):
        ds = ds.rename({"nav_lat": "lat"})

    if 'values' in list(ds.dims):
        ds = ds.rename({"values": "cell"})

    return ds


def adjust_clim_file(cfield, remove_zero=False):
    """Routine to fix file format of climatology"""

    # fix coordinates
    org = ['LONGITUDE', 'LATITUDE', 'lev']
    new = ['lon', 'lat', 'plev']
    for o, n in zip(org, new):
        if o in cfield.coords:
            cfield = cfield.rename({o: n})

    # extract data_array
    cname = list(cfield.data_vars)[-1]
    field = cfield[cname]

    if remove_zero:
        field = field.where(field != 0)

    # convert vertical levels
    if 'plev' in cfield.coords:
        field = field.metpy.convert_coordinate_units('plev', 'Pa')

    return field


#############################
# UNIT ADJUSTMENT FUNCTIONS #
#############################


def units_extra_definition():
    """Add units to the pint registry required by ECMean4"""

    # special units definition
    units.define('fraction = [] = frac')
    units.define('psu = 1e-3 frac')
    units.define('PSU = 1e-3 frac')
    units.define('Sv = 1e+6 m^3/s')  # Replace Sievert with Sverdrup


def units_converter(org_units, tgt_units):
    """Units conversion using metpy and pint.
    From a org_units convert to tgt_units providing offset and factor.
    Some assumptions are done for precipitation field: must be extended
    to other vars. It will not work if BOTH factor and offset are required"""

    units_relation = (units(org_units) / units(tgt_units)).to_base_units()
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
            logging.info("Assuming this as a water flux! Am I correct?")
            logging.info("Dividing by water density...")
            density_water = units('kg / m^3') * 1000
            offset = 0.
            factor = (factor_standard / density_water).to(tgt_units).magnitude

        else:
            logging.error(units_relation)
            sys.exit("ERROR: Units mismatch, this cannot be handled!")
    else:
        offset = 0.
        factor = 1.

    logging.info('Offset is ' + str(offset))
    logging.info('Factor is ' + str(factor))
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


####################
# OUTPUT FUNCTIONS #
####################

def dict_to_dataframe(varstat):
    """very clumsy conversion of the nested 3-level dictionary
    to a pd.dataframe: NEED TO BE IMPROVED"""
    data_table = {}
    for i in varstat.keys():
        pippo = {}
        for outerKey, innerDict in varstat[i].items():
            for innerKey, values in innerDict.items():
                pippo[(outerKey, innerKey)] = values
        data_table[i] = pippo
    data_table = pd.DataFrame(data_table).T
    return data_table


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
            print(
                '{:12.5f}'.format(
                    varmean[var] *
                    ref[var].get(
                        'factor',
                        1)),
                end=' ',
                file=f)
        print(file=f)

##################
# PLOT FUNCTIONS #
##################


def heatmap_comparison(absolute_table, relative_table, diag, filemap):
    """Function to produce a heatmap - seaborn based - for Performance Indices
    based on CMIP6 ratio"""

    nplots = 1
    # real plot
    fig, axs = plt.subplots(nplots, 1, sharey=True, tight_layout=True, figsize=(15, nplots*8))

    for k in range(nplots):
        if k == 0:
            thr = [0, 1, 5]
            tictoc = [0, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5]
            title = 'CMIP6 RELATIVE PI'
            myfield = relative_table
        elif k == 1:
            thr = [0, 10, 20]
            tictoc = list(range(0, 20, 2))
            title = 'ABSOLUTE PI'
            myfield = absolute_table

        # axs.subplots_adjust(bottom=0.2)
        # pal = sns.diverging_palette(h_neg=130, h_pos=10, s=99, l=55, sep=3, as_cmap=True)
        tot = (len(myfield.columns))
        sss = (len(set([tup[1] for tup in myfield.columns])))
        divnorm = TwoSlopeNorm(vmin=thr[0], vcenter=thr[1], vmax=thr[2])
        pal = sns.color_palette("Spectral_r", as_cmap=True)
        # pal = sns.diverging_palette(220, 20, as_cmap=True)
        chart = sns.heatmap(myfield, norm=divnorm, cmap=pal,
                            cbar_kws={"ticks": tictoc, 'label': title},
                            ax=axs, annot=True, linewidth=0.5,
                            annot_kws={'fontsize': 11, 'fontweight': 'bold'})
        chart = chart.set_facecolor('whitesmoke')
        axs.set_title(f'{title} {diag.modelname} {diag.year1} {diag.year2}', fontsize=25)
        axs.vlines(list(range(sss, tot+sss, sss)), ymin=-1, ymax=len(myfield.index), colors='k')
        axs.hlines(len(myfield.index)-1, xmin=-1, xmax=len(myfield.columns), colors='purple', lw=0.8)
        names = [' '.join(x) for x in myfield.columns]
        if (k == (nplots-1)):
            axs.set_xticks([x+.5 for x in range(len(names))], names, rotation=45, ha='right', fontsize=15)
        else:
            axs.set_xticks([])
        axs.set_yticks([x+.5 for x in range(len(myfield.index))], myfield.index, rotation=0, fontsize=18)
        axs.set(xlabel=None)

    # save and close
    plt.savefig(filemap)
    plt.cla()
