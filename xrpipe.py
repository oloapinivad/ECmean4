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



def make_atm_masks(component, atminifile):
    """Create land-sea masks for atmosphere model"""
    # prepare ATM LSM: this need to be improved, since it is clearly model dependent
    if component == 'oifs':
        # create mask
        mask = xr.open_dataset(atminifile, engine="cfgrib")
        print(mask)
        mask = mask['lsm']
    elif component == 'cmoratm':
        sys.exit("Mask from cmor non defined yet mismatch, this cannot be handled!")
        #self.LANDMASK = self.cdo.selname('sftlf',
        #                                     input=f'-gec,50 {extra} {self.atmfix} {atminifile}')
        #    self.SEAMASK = self.cdo.mulc('-1', input=f'-subc,1 {self.LANDMASK}')
    return mask




