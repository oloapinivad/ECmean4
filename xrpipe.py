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
import cf2cdm

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

def area_cell(xfield): 
    """Function which estimate the area cell from bounds. This is done assuming 
    trapezoidal shape of the grids - useful for reduced grids"""    

    dlon = xfield['bounds_lon'].isel(nvertex=1) - xfield['bounds_lon'].isel(nvertex=2)
    dlat = xfield['bounds_lat'].isel(nvertex=2) - xfield['bounds_lat'].isel(nvertex=3)
    Earth_Radius = 6371000.
    #arclon = Earth_Radius * np.cos(np.deg2rad(xfield.lat)) * np.deg2rad(dlon)
    arclon1 =  Earth_Radius * np.cos(np.deg2rad(abs((xfield['bounds_lat'].isel(nvertex=2))))) * np.deg2rad(dlon)
    arclon2 = Earth_Radius * np.cos(np.deg2rad(abs(xfield['bounds_lat'].isel(nvertex=3)))) * np.deg2rad(dlon)
    arclat = Earth_Radius * np.deg2rad(dlat)

    #area = arclon * arclat
    area_cell = (arclon1 + arclon2) * arclat / 2
    return area_cell


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
    and an ad-hoc dictionary. Could be improved"""
    
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
    
    return (util)


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
            area = xfield['e1t']*xfield['e2t']
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




