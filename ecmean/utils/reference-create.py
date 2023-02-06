#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""
    Tool to create a new ECmean4 reference climatology for global mean.
    It requires to have cdo and cdo-bindings installed
"""

__author__ = "Paolo Davini (p.davini@isac.cnr.it), Feb 2023."


from ecmean.libs.files import load_yaml
from ecmean.libs.units import units_extra_definition, units_converter, _units_are_integrals
from ecmean.libs.ncfixers import xr_preproc
from ecmean.libs.masks import masked_meansum, select_region
from ecmean.libs.areas import _area_cell
import logging
import yaml
import numpy as np
import xarray as xr
import os
from cdo import *
from glob import glob
cdo = Cdo()

# set default logging
logging.basicConfig(level=logging.WARNING)

# variable list
variables = ['pr', 'pr_land']

# to set: time period (default, can be shorter if data are missing) #WARNING MISSING
year1 = 1990
year2 = 1994

# yml file to get information on dataset on some machine
clim_info = '/home/paolo/ECmean4/ecmean/utils/create-reference-wilma.yml'

# add other units
units_extra_definition()

# open the clim info file
info = load_yaml(clim_info)
years = list(range(year1, year2 + 1))

# skip NaN: if False, yearly/season average require that all the points are defined in the correspondent
# time window.
nanskipper = False

# set the land-sea mask (ERA5)
maskfile = info['mask']


# climatology yaml output
clim_name = 'EC23'
clim_file = os.path.join('gm_reference_' + clim_name + '.yml')

 # loop on variables to be processed
for var in variables:

    print(var)
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
    print("Compute cell area for weights...")
    first = cfield.to_dataset(name=var)
    weights = _area_cell(first).load()

    mask_type = info[var].get('total', 'global')
    domain = info[var].get('domain', 'atm')

    # compute land sea mask
    if mask_type != 'global':
        print("Mask...")
        xmask = cdo.remapbil(glob(filedata)[0], input = maskfile, returnXArray =  'lsm')
        xmask.load()
    else: 
        xmask = 0.

    # yearly and season averages
    print("Time averages...")
    gfield1 = cfield.resample(time='AS', skipna=nanskipper).mean('time', skipna=nanskipper).load()
    gfield2 = cfield.resample(time='Q-NOV', skipna=nanskipper).mean('time', skipna=nanskipper).load()

    print("Season loop...")
    mf = {}
    sf = {}
    for season in ['ALL', 'DJF', 'MAM', 'JJA', 'SON']:
    # select the season
        if season == 'ALL':
            gfield = gfield1
        else:
            gfield = gfield2.sel(time=gfield2.time.dt.season.isin(season))
        # for winter, we drop first and last to have only complete season.
        # this reduces the sample by one but it is safer for variance
        if season == 'DJF':
            gfield = gfield.drop_isel(time=[0, gfield.sizes['time'] - 1])

        mr = {}
        ms = {}
        print("Region loop...")
        for region in ['Global', 'North Midlat', 'Tropical', 'South Midlat']:
            
            slicefield = select_region(gfield, region)
            sliceweights = select_region(weights, region)
            if mask_type != 'global':
                slicemask = select_region(xmask, region)
            else:
                slicemask = 0.

            out = masked_meansum(xfield=slicefield, weights=sliceweights,
                    mask_type=info[var].get('total', 'global'),
                    dom=domain, mask=slicemask)
            
            new_units = _units_are_integrals(info[var]['org_units'], info[var])
            offset, factor = units_converter(new_units, info[var]['tgt_units'])
            final = out * factor + offset

            omean = np.mean(final)
            ostd = np.std(final)
            mr[region] = str(round(omean,3))
            ms[region] = str(round(ostd, 3))
            print('{} {} {} mean is: {:.2f} +- {:.2f}'.format(var, season, region, omean, ostd))

        # store everything
        mf[season] = mr
        sf[season] = ms

    # log output
    logging.info(mf)
    logging.info(sf)

    # preparing clim file
    if os.path.isfile(clim_file):
        dclim = load_yaml(clim_file)
    else:
        dclim = {}

    # initialize variable if not exists
    if var not in dclim:
        dclim[var] = {}

    dclim[var]['dataset'] = info[var]['dataset']
    # dclim[var]['dataname'] = info[var]['varname']
    if mask_type != 'global':
        dclim[var]['total'] = mask_type
    dclim[var]['units'] = info[var]['tgt_units']
    dclim[var]['year1'] = int(real_year1)
    dclim[var]['year2'] = int(real_year2)
    dclim[var]['mean'] =  mf
    dclim[var]['std'] = sf

    # dump the yaml file
    with open(clim_file, 'w') as file:
        yaml.safe_dump(dclim, file, sort_keys=False)
    

    #zfield = cfield.resample(time='1MS', skipna=nanskipper).mean('time', skipna=nanskipper)

