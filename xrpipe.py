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

    isavail = True
    isavail = isavail and os.path.isfile(ffile)
    #isavail = isavail and (len(flist) > 0)

    if isavail:
        xfield = xr.open_dataset(ffile)
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
        print(f'Not available: {var} File: {infile}')
        logging.warning("Requested file %s is not available.", infile)

    return isavail, varunit

def masked_meansum(xfield, var, weights, mask_type, mask):
    """Evaluate the weighted averaged for global varialbes and for 
    required vars estimate the land-only or ocean only surface integral"""

    tfield = xfield.mean(dim='time_counter').to_dataset(name = var)
    if mask_type == 'land':   
        tfield['mask'] = (('cell'), mask.values)
        out = tfield[var].where(tfield['mask'] >= 0.5).weighted(weights).sum().item()
    elif mask_type in ['sea', 'ocean']:
        tfield['mask'] = (('cell'), mask.values)
        out = tfield[var].where(tfield['mask'] < 0.5).weighted(weights).sum().item()
    else:
        out = tfield[var].weighted(weights).mean().item()
        
    return out

# def mask(xfield, var, mask_type, mask) : 
#     tfield = xfield.to_dataset(name = var)
#     if mask_type == 'land':
#         tfield['mask'] = (('cell'), mask.values)
#         out = tfield[var].where(tfield['mask'] >= 0.5)
#     elif mask_type in ['sea', 'ocean']:
#         tfield['mask'] = (('cell'), mask.values)
#         out = tfield[var].where(tfield['mask'] >= 0.5)
#     return out

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
    xfield['area'] = xfield[list(xfield.data_vars)[-1]]

    # bounds available
    if 'bounds_lon' in xfield.data_vars and 'bounds_lat' in xfield.data_vars :

        # the assumption of the vertex position is absolutely random, need to generalized
        bounds_lon = np.column_stack((xfield['bounds_lon'].isel(nvertex=1), 
            xfield['bounds_lon'].isel(nvertex=2)))
        bounds_lat = np.column_stack((xfield['bounds_lat'].isel(nvertex=2), 
            xfield['bounds_lat'].isel(nvertex=3)))
        #dlon = xfield['bounds_lon'].isel(nvertex=1) - xfield['bounds_lon'].isel(nvertex=2)
        #dlat = xfield['bounds_lat'].isel(nvertex=2) - xfield['bounds_lat'].isel(nvertex=3)

    # no bounds defined, using lon/lat estimation
    else :

        # all this is made with numpy, perhaps better to use xarray?
        blon = guess_bounds(xfield['lon'], name = 'lon')
        bounds_lon = np.tile(blon.transpose(), len(xfield['lat'])).transpose()
        blat =  guess_bounds(xfield['lat'], name = 'lat')
        bounds_lat = np.repeat(blat, len(xfield['lon']), axis = 0)

    dlon = abs(bounds_lon[:,0] - bounds_lon[:,1])
    dlat = abs(bounds_lat[:,0] - bounds_lat[:,1])

    # safe check on cosine of 90: estimate the trapezoid
    arclon1 =  Earth_Radius * abs(np.cos(abs(np.deg2rad(bounds_lat[:,0])))) * np.deg2rad(dlon)
    arclon2 =  Earth_Radius * abs(np.cos(abs(np.deg2rad(bounds_lat[:,0])))) * np.deg2rad(dlon)
    arclat = Earth_Radius * np.deg2rad(dlat)

    # trapezoid area
    area_cell = (arclon1 + arclon2) * arclat / 2

    # we want a mask which is not dependent on time
    for t in ['time', 'time_counter'] :
        if t in list(xfield.dims) : 
            xfield['area'] = xfield['area'].mean(dim=t)

    # if we are using a lon/lat regular grid reshape area_cell
    if 'lon' in list(xfield.dims) : 
        area_cell = area_cell.reshape([len(xfield['lon']), len(xfield['lat'])]).transpose() 
        
    
    # since we are using numpy need to bring them back into xarray dataset
    xfield['area'].values = np.squeeze(area_cell)

    return xfield['area']


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

def util_dictionary(component, maskatmfile, atmareafile, oceareafile) : 
    """Create a dictionary with atmospheric mask and 
    atmospheric and oceanic area weights"""

    util = {
        'atm_mask' : _make_atm_masks(component['atm'], maskatmfile),
        'atm_weights': _make_atm_areas(component['atm'], atmareafile),
        'oce_weights': _make_oce_areas(component['oce'], oceareafile)
    }
    
    return util

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
    #if remove_zero : 
    #    field = field.where(field==0, 1, cfield)
    #print(field)
    
    # convert vertical levels 
    if 'pressure_levels' in cfield.coords :
        field = field.metpy.convert_coordinate_units('pressure_levels', 'Pa')

    return field

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
    

def _make_atm_masks(component, maskatmfile):
    """Create land-sea masks for atmosphere model"""
    # prepare ATM LSM: this need to be improved, since it is clearly model dependent
    if component == 'oifs':
        # create mask: opening a grib and loading only lsm to avoid inconsistencies in the grib 
        # structure -> see here https://github.com/ecmwf/cfgrib/issues/13
        mask = xr.open_dataset(maskatmfile, engine="cfgrib", filter_by_keys={'shortName': 'lsm'})
        mask = mask['lsm']
    elif component == 'cmoratm':
        sys.exit("Mask from cmor non defined yet mismatch, this cannot be handled!")
        #self.LANDMASK = self.cdo.selname('sftlf',
        #                                     input=f'-gec,50 {extra} {self.atmfix} {atminifile}')
        #    self.SEAMASK = self.cdo.mulc('-1', input=f'-subc,1 {self.LANDMASK}')
    return mask

def xr_get_clim_files(piclim, var, diag) : 

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



def _make_atm_areas(component, atmareafile) : 
    "Create atmospheric weights for area operations"
    if component == 'oifs' : 
        xfield = xr.open_dataset(atmareafile)
        area = area_cell(xfield)
    else :
        sys.exit("Area this cannot be handled!")

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
        else :
            sys.exit("Area this cannot be handled!")
    else :
        area = None
    return area


def xr_get_inifiles(face, diag):
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




