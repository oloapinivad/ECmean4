#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""
    Tool to create a new ECmean4 climatology.
    It requires to have cdo and cdo-bindings installed
"""

__author__ = "Paolo Davini (p.davini@isac.cnr.it), Gen 2023."


from ecmean.libs.files import load_yaml
from ecmean.libs.units import units_extra_definition, units_converter
from ecmean.libs.ncfixers import xr_preproc
from ecmean.libs.masks import mask_field, _area_cell, masked_meansum
import logging
from time import time
import numpy as np
import xarray as xr
import os


variables = ['pr']

# to set: time period (default, can be shorter if data are missing) #WARNING MISSING
year1 = 1990
year2 = 1990

# yml file to get information on dataset on some machine
clim_info = 'create-reference-wilma.yml'

# add other units
units_extra_definition()

# open the clim info file
info = load_yaml(clim_info)
years = list(range(year1, year2 + 1))

# skip NaN: if False, yearly/season average require that all the points are defined in the correspondent
# time window.
nanskipper = False

 # loop on variables to be processed
for var in variables:

    print(var)
    tic = time()
    # get the directory
    filedata = str(os.path.expandvars(info[var]['dir'])).format(
        datadir=info['dirs']['datadir'],
        dataset=info[var]['dataset'],
        varname=info[var]['varname'])
    logging.info(filedata)

     # load data and time select
    print("Loading multiple files...")
    # unable to operate with Parallel=True
    xfield = xr.open_mfdataset(filedata, chunks='auto',
                                parallel=False, preprocess=xr_preproc, engine='netcdf4')
    xfield = xfield.rename({info[var]['varname']: var})
    cfield = xfield[var].sel(time=xfield.time.dt.year.isin(years))

    real_year1 = np.min(cfield.time.dt.year.values)
    real_year2 = np.max(cfield.time.dt.year.values)
    if (real_year1 != year1):
        logging.warning("Initial year different from what expected: " + str(real_year1))
    if (real_year2 != year2):
        logging.warning("Final year different from what expected: " + str(real_year2))

    # get a single record and exploit of ecmean function to estimate areas
    first = cfield.isel(time=0).to_dataset(name=var)
    weights = _area_cell(first)
    #print(weights)
    
    # global mean
    mask_type = 'global'

    #alfa = masked_meansum(cfield, var, weights, mask_type, 0)
    #print(alfa)

    masked = mask_field(cfield, var, mask_type, 0)
    print(masked)
    # average
    if mask_type in ['global']:
        out = masked[var].weighted(weights).mean(('lon', 'lat'))
    # global integrals
    elif mask_type in ['land', 'ocean', 'sea', 'north', 'south']:
        out = masked.weighted(weights).sum(dim=['lon', 'lat']).values

    print(out)

    offset, factor = units_converter(info[var]['org_units'], info[var]['tgt_units'])

    final = out.mean() * factor + offset
    print(final)

    #zfield = cfield.resample(time='1MS', skipna=nanskipper).mean('time', skipna=nanskipper)

