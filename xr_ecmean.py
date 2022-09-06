#!/usr/bin/env python3
'''
Shared functions for XArray ECmean4
'''

import numpy as np
import xarray as xr
import os
import re
import logging
import operator
import sys
from pathlib import Path
from glob import glob
import xesmf as xe
from metpy.units import units
import yaml
import itertools


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
            self.resolution = 'r180x91'
        else:
            if not self.resolution:
                self.resolution = cfg['PI']['resolution']

        # Various input and output directories
        self.ECEDIR = Path(os.path.expandvars(cfg['dirs']['exp']))
        self.TABDIR = Path(os.path.expandvars(cfg['dirs']['tab']))
        self.CLMDIR = Path(os.path.expandvars(cfg['dirs']['clm']), self.climatology)
        self.RESCLMDIR = Path(self.CLMDIR, self.resolution)
        self.years_joined = ''

        self.linefile = self.TABDIR / 'global_means.txt'

        # check if output attribute exists
        if hasattr(self, 'output'):
            self.linefile = args.output
            self.ftable = True


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

def chunks(iterable, num):
    """Generate num adjacent chunks of data from a list iterable
       Split lists in a convenient way for a parallel process"""
    size = int(np.ceil(len(iterable) / num))
    it = iter(iterable)
    return iter(lambda: tuple(itertools.islice(it, size)), ())

def getdomain(var, face):
    """Given a variable var extract its domain (oce or atm) from the interface.
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

##################
# FILE FUNCTIONS #
##################

def var_is_there(flist, var, reference):
    """Check if a variable is available in the input file and provide its units."""

    # we expect a list obtained by glob
    isavail = True
    for f in flist:
        isavail = isavail and os.path.isfile(f)
    isavail = isavail and (len(flist) > 0)

    if isavail:
        xfield = xr.open_mfdataset(flist)
        vars_avail = [i for i in xfield.data_vars]
        units_avail ={}
        for i in vars_avail :
            try: 
                k = xfield[i].units
            except : 
                k = "None"
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
        print(f'Not available: {var} File: {flist}')
        logging.warning("Requested file %s is not available.", flist)

    return isavail, varunit

def get_clim_files(piclim, var, diag) : 

    """Function to extra names for the climatology files"""

    # extract info from pi_climatology.yml
    # reference dataset and reference varname
    # as well as years when available
    dataref = piclim[var]['dataset']
    dataname = piclim[var]['dataname']
    datayear1 = piclim[var].get('year1', 'nan')
    datayear2 = piclim[var].get('year2', 'nan')

    # get files for climatology
    if diag.climatology == 'RK08':
        clim = str(diag.RESCLMDIR / f'climate_{dataref}_{dataname}.nc')
        vvvv = str(diag.RESCLMDIR / f'variance_{dataref}_{dataname}.nc')
    elif diag.climatology == 'EC22':
        clim = str(diag.RESCLMDIR / f'climate_{dataname}_{dataref}_{diag.resolution}_{datayear1}-{datayear2}.nc')
        vvvv = str(diag.RESCLMDIR / f'variance_{dataname}_{dataref}_{diag.resolution}_{datayear1}-{datayear2}.nc')

    return clim, vvvv

def get_inifiles(face, diag):
    """
    Return the inifiles from the interface, needs the component dictionary
    Check if inifiles exist.
    """
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

def load_yaml(infile):
    """Load generic yaml file"""
    try:
        with open(infile, 'r', encoding='utf-8') as file:
            cfg = yaml.load(file, Loader=yaml.FullLoader)
    except IOError:
        sys.exit(f'{infile} not found: you need to have this configuration file!')
    return cfg

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
            filename = filename + filename1 
    #if len(filename) == 1:  # glob always returns a list, return str if only one
    #    filename = filename[0]
    logging.debug("Filenames: %s", filename)
    return filename


#####################
# AVERAGE FUNCTIONS #
#####################

def masked_meansum(xfield, var, weights, mask_type, mask):
    """Evaluate the weighted averaged for global varialbes and for 
    required vars estimate the land-only or ocean only surface integral"""

    #use generator expression
    for g in (t for t in ['time', 'time_counter'] if t in list(xfield.dims)) : 
        tfield = xfield.mean(dim=g).to_dataset(name = var)

    if mask_type == 'land':   
        tfield['mask'] = (tuple(tfield.coords), mask.values)
        out = tfield[var].where(tfield['mask'] >= 0.5).weighted(weights).sum().values
    elif mask_type in ['sea', 'ocean']:
        tfield['mask'] = (tuple(tfield.coords), mask.values)
        out = tfield[var].where(tfield['mask'] < 0.5).weighted(weights).sum().values
    else:
        out = tfield[var].weighted(weights.fillna(0)).mean().values
        
    return float(out)

##################################
# AREA-WEIGHT AND MASK FUNCTIONS #
##################################

def util_dictionary(component, maskatmfile, atmareafile, oceareafile) : 
    """Create a dictionary with atmospheric mask and 
    atmospheric and oceanic area weights"""

    util = {
        'atm_mask' : _make_atm_masks(component['atm'], maskatmfile),
        'atm_weights': _make_atm_areas(component['atm'], atmareafile),
        'oce_weights': _make_oce_areas(component['oce'], oceareafile)
    }
    
    return util

def _make_atm_masks(component, maskatmfile):
    """Create land-sea masks for atmosphere model"""
    # prepare ATM LSM: this need to be improved, since it is clearly model dependent
    if component == 'oifs':
        # create mask: opening a grib and loading only lsm to avoid inconsistencies in the grib 
        # structure -> see here https://github.com/ecmwf/cfgrib/issues/13
        mask = xr.open_dataset(maskatmfile, engine="cfgrib", filter_by_keys={'shortName': 'lsm'})
        mask = mask['lsm']
    elif component == 'cmoratm':
        mask = xr.open_mfdataset(maskatmfile)
        mask = mask['sftlf']
    else : 
        sys.exit("Mask undefined yet mismatch, this cannot be handled!")

    return mask

def _make_atm_areas(component, atmareafile) : 
    "Create atmospheric weights for area operations"
    if component == 'oifs' : 
        xfield = xr.open_dataset(atmareafile)
        area = area_cell(xfield)
    elif component == 'cmoratm' :
        xfield = xr.open_mfdataset(atmareafile)
        area = area_cell(xfield)
    else :
        sys.exit("Area for this configuration cannot be handled!")
    return area

def _make_oce_areas(component, oceareafile) : 
    "Create atmospheric weights for area operations"
    if oceareafile: 
        if component == 'nemo' : 
            xfield = xr.open_dataset(oceareafile)
            if 'e1t' in xfield.data_vars : 
                area = xfield['e1t']*xfield['e2t']
            else : 
                area = area_cell(xfield)
        elif component == 'cmoroce' :
            area = xr.open_mfdataset(oceareafile)['areacello']
        else :
            sys.exit("Area for this configuration cannot be handled!")
    else :
        area = None
    return area

def guess_bounds(axis, name = 'lon') : 
    #inspired by https://gist.github.com/dennissergeev/60bf7b03443f1b2c8eb96ce0b1880150
    """Basic function that estimates the boundaries for lon and lat if they are not
    available. Works only with regular grids. 
    It also avoid having values larger than 90N/90S as well as 0 and 10^5 Pa"""

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
    if name in 'lev' :
        min_bounds[0] = 0
        max_bounds[-1] = 100000

    # should we a xarray object instead of a numpy?
    bounds = np.array([min_bounds, max_bounds]).transpose()
    return(bounds)

def area_cell(xfield): 
    """Function which estimate the area cell from bounds. This is done assuming 
    trapezoidal shape of the grids - useful for reduced grids. 
    Working also on regular grids which does not have lon/lat bounds
    via the guess_bounds function"""    

    Earth_Radius = 6371000.
    
    # bounds available: EC-Earth4 definition
    if 'bounds_lon' in xfield.data_vars and 'bounds_lat' in xfield.data_vars :

        # logging 
        logging.debug('bounds_lon/lat found...')

        # the assumption of the vertex position is absolutely random, need to generalized
        bounds_lon = np.column_stack((xfield['bounds_lon'].isel(nvertex=1), 
            xfield['bounds_lon'].isel(nvertex=2)))
        bounds_lat = np.column_stack((xfield['bounds_lat'].isel(nvertex=2), 
            xfield['bounds_lat'].isel(nvertex=3)))
        #dlon = xfield['bounds_lon'].isel(nvertex=1) - xfield['bounds_lon'].isel(nvertex=2)
        #dlat = xfield['bounds_lat'].isel(nvertex=2) - xfield['bounds_lat'].isel(nvertex=3)

    # bounds available: cmor definition
    # direction of tiling is very weird, this has to be rearranged
    if 'lon_bnds' in xfield.data_vars and 'lat_bnds' in xfield.data_vars :

         # logging 
        logging.debug('lon/lat_bounds found...')

        # get bounds dimension
        dimslist = list(xfield.dims) 

        # if it is not called bnds, rename it
        if 'bnds' not in  dimslist : 
            logging.debug('bnds not found, renaming it...')         
            for g in (t for t in list(xfield.dims) if t not in ['lon', 'lat', 'time']) : 
                bdim = g
            xfield = xfield.rename_dims({bdim: "bnds"})
    
        blon = np.column_stack((xfield['lon_bnds'].isel(bnds=0), 
            xfield['lon_bnds'].isel(bnds=1)))
        blat = np.column_stack((xfield['lat_bnds'].isel(bnds=0), 
            xfield['lat_bnds'].isel(bnds=1)))

        bounds_lon = np.repeat(blon, len(xfield['lat']), axis = 0)
        bounds_lat = np.tile(blat.transpose(), len(xfield['lon'])).transpose()

    # no bounds defined, using lon/lat estimation
    # direction of tiling is very weird, this has to be rearranged
    else :

        # logging 
        logging.debug('no bounds found, guessing area cell...')
        
        blon = guess_bounds(xfield['lon'], name = 'lon')
        blat =  guess_bounds(xfield['lat'], name = 'lat')

        # all this is made with numpy, perhaps better to use xarray? 
        # should be improved, it is very clumsy
        if list(xfield.coords)[0] == 'lat' :
            bounds_lon = np.repeat(blon, len(xfield['lat']), axis = 0)
            bounds_lat = np.tile(blat.transpose(), len(xfield['lon'])).transpose()
        elif list(xfield.coords)[0] == 'lon' :
            
            bounds_lon = np.tile(blon.transpose(), len(xfield['lat'])).transpose()
            bounds_lat = np.repeat(blat, len(xfield['lon']), axis = 0)

    # cell dimension
    dlon = abs(bounds_lon[:,0] - bounds_lon[:,1])
    dlat = abs(bounds_lat[:,0] - bounds_lat[:,1])

    # safe check on cosine of 90: estimate the trapezoid
    arclon1 =  Earth_Radius * abs(np.cos(abs(np.deg2rad(bounds_lat[:,0])))) * np.deg2rad(dlon)
    arclon2 =  Earth_Radius * abs(np.cos(abs(np.deg2rad(bounds_lat[:,1])))) * np.deg2rad(dlon)
    arclat = Earth_Radius * np.deg2rad(dlat)

    # trapezoid area
    area_cell = (arclon1 + arclon2) * arclat / 2
    
    #use generator expression to remove time dependency
    for g in (t for t in ['time', 'time_counter'] if t in list(xfield.dims)) : 
        xfield['area'] = xfield['area'].mean(dim=g)

    # if we are using a lon/lat regular grid reshape area_cell
    if 'lon' in list(xfield.dims) : 
        area_cell = area_cell.reshape([len(xfield['lon']), len(xfield['lat'])]).transpose() 

    # since we are using numpy need to bring them back into xarray dataset
    #xfield['area'].values = np.squeeze(area_cell)
    xfield['area'] = (('lat', 'lon'), area_cell)

    # check the total area
    logging.debug('Total Earth Surface: ' + str(xfield['area'].sum().values/10**6) + ' Km2')

    return xfield['area']


#####################
# FORMULA FUNCTIONS #
#####################

# this is a tool to parse CDO-based formula into mathematical operatos
# there might exists something more intelligent as pyparsing package
def eval_formula(mystring, xdataset):
    """Evaluate the cmd string provided by the yaml file
    producing a parsing for the derived variables""" 
    # Tokenize the original string
    token = [i for i in re.split('(\W+)', mystring) if i ]
    if (len(token)>1) :
        # Use order of operations 
        out = operation(token, xdataset)
    else :
        out = xdataset[token[0]]
    return out
    
# core of the parsing operation, using dictionaries and operator package
def operation(token, xdataset) : 
    """Parsing of the CDO-based function using operator package
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
    dict = {}
    for k in token :
        if k not in ops : 
            if not is_number(k): 
                dict[k] = xdataset[k]
            else : 
                dict[k] = float(k)
    
    # apply operators to all occurrences, from top priority
    # so far this is not parsing parenthesis
    code = 0
    for p in ops:
        #print('Operation:' + p)   
        while p in token:
            code += 1
            #print(token) 
            x = token.index(p)
            name = 'op' + str(code)
            replacer = ops.get(p)(dict[token[x-1]], dict[token[x+1]])
            #print(replacer)
            dict[name] = replacer
            token[x-1] = name
            del token[x:x+2]
            #print(token)
    return replacer

###########################
# INTERPOLATION FUNCTIONS #
###########################

def remap_dictionary(component, atmareafile, oceareafile, target_grid) : 
    """Create a dicitionary with atmospheric and oceanic weights for 
    interpolation. There is an option of fix grid before the real interpolation
    this is used for gaussian reduced grids"""

    atmfix, atmremap = _make_atm_interp_weights(component['atm'], atmareafile, target_grid)
    if oceareafile : 
        ocefix, oceremap = _make_oce_interp_weights(component['oce'], oceareafile, target_grid)
    else :
        ocefix = None
        oceremap = None

    remap = {
        'atm_fix' : atmfix, 
        'atm_remap': atmremap,
        'oce_fix' : ocefix, 
        'oce_remap': oceremap,
    }

    return remap

def _make_atm_interp_weights(component, atmareafile, target_grid) :
    """"Create atmospheric interpolator"""
    if component == 'oifs':

        # this is to get lon and lat from the Equator
        xfield = xr.open_dataset(atmareafile)
        m = xfield['tas'].isel(time_counter=0).load()
        g = sorted(list(set(m.lat.values)))
        f = sorted(list(m.sel(cell=m.lat==g[int(len(g)/2)]).lon.values))

        # this creates a a gaussian non reduced grid
        ds_out = xr.Dataset({"lon": (["lon"], f), "lat": (["lat"], g)})

        # use nearest neighbour to remap to gaussian regular
        fix = xe.Regridder(xfield['tas'], ds_out, 
            method = "nearest_s2d", locstream_in=True, periodic = True)

        # create bilinear interpolator
        interp = xe.Regridder(fix(xfield['tas']), target_grid, periodic = True, method = "bilinear")

    elif component == 'cmoratm':
        sys.exit("Mask from cmor non defined yet mismatch, this cannot be handled!")
    
    return fix, interp

def _make_oce_interp_weights(component, oceareafile, target_grid) :
    """"Create atmospheric interpolator"""
    if component == 'nemo':
        fix = None
        xfield = xr.open_dataset(oceareafile)
        # set coordinates which are missing
        xfield = xfield.set_coords(['nav_lon', 'nav_lat', 'nav_lev', 'time_counter'])
        # rename lon and lat for interpolation
        xfield = xfield.rename_dims({"z": "deptht"})
        xfield = xfield.rename({"nav_lon": "lon", "nav_lat": "lat",  "nav_lev": "deptht" })

        #final = xe.util.grid_global(target, target)
        # use grid distance as generic variable
        interp = xe.Regridder(xfield['e1t'], 
            target_grid, method = "bilinear", ignore_degenerate=True)
    elif component == 'cmoratm':
        sys.exit("Mask from cmor non defined yet mismatch, this cannot be handled!")

    return fix, interp


#########################
# FILE FORMAT FUNCTIONS #
#########################

def adjust_clim_file(cfield, remove_zero = False) : 
    """Routine to fix file format of climatology"""

    # fix coordinates
    org = ['LONGITUDE', 'LATITUDE', 'lev', 'plev']
    new = ['lon', 'lat', 'pressure_levels', 'pressure_levels']
    for o, n in zip(org, new) : 
        if o in cfield.coords : 
            cfield = cfield.rename({o: n})

    # extract data_array
    cname = list(cfield.data_vars)[-1]
    field = cfield[cname]

    #print(field)
    if remove_zero : 
        field = field.where(field!=0)
    #print(field)
    
    # convert vertical levels 
    if 'pressure_levels' in cfield.coords :
        field = field.metpy.convert_coordinate_units('pressure_levels', 'Pa')

    return field

#############################
# UNIT ADJUSTMENT FUNCTIONS #
#############################

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

####################
# OUTPUT FUNCTIONS #
####################

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