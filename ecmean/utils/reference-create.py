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
from ecmean.libs.formula import _eval_formula
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
variables = ['tas', 'psl', 'pr', 'evspsbl', 'pme', 'clt', 'cll', 'clm', 'clh',
             'pr_oce', 'pme_oce', 'pr_land', 'pme_land']
variables = ['tos', 'siconc']
# clh does not work, double check


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

def expand_filedata(directory, var, info) : 

    return os.path.expandvars(directory).format(datadir=info['dirs']['datadir'],
            eradir=info['dirs']['eradir'], esadir=info['dirs']['esadir'],
            dataset=info[var]['dataset'],filevar=info[var]['filevar'])

# climatology yaml output
clim_name = 'EC23'
clim_file = os.path.join('gm_reference_' + clim_name + '.yml')

 # loop on variables to be processed
for var in variables:

    print(var)
    # get the directory (specific treatment if a list, use glob to expand)
    if isinstance(info[var]['dir'], list):
        temp_list = [glob(expand_filedata(ss, var, info)) for ss in info[var]['dir']]
        filedata = [item for sublist in temp_list for item in sublist]
    # use wildcards from xarray
    else:
        filedata = str(expand_filedata(info[var]['dir'], var, info))
    logging.warning(filedata)

     # load data and time select
    print("Loading multiple files...")
    xfield = xr.open_mfdataset(filedata, chunks = 'auto', preprocess=xr_preproc, engine='netcdf4')

    # if derived, use the formula skill (or just rename)
    if 'derived' in info[var].keys():
        cmd = info[var]['derived']
        xfield = _eval_formula(cmd, xfield).to_dataset(name=var)
    else:
        xfield = xfield.rename({info[var]['filevar']: var})

    # select time
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

    # get infos on domain, operation and masks
    mask_type = info[var].get('mask', 'global')
    domain = info[var].get('domain', 'atm')
    operation = info[var].get('operation', 'mean')

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

        mf[season] = {}
        print("Region loop...")
        for region in ['Global', 'North Midlat', 'Tropical', 'South Midlat']:
            
            # slice everything
            slicefield = select_region(gfield, region)
            sliceweights = select_region(weights, region)
            if mask_type != 'global':
                slicemask = select_region(xmask, region)
            else:
                slicemask = 0.

            # get the masked-mean-sum
            out = masked_meansum(xfield=slicefield, weights=sliceweights,
                    mask=slicemask, operation=operation,
                    mask_type=mask_type, domain=domain)
            
            # set the units
            new_units = _units_are_integrals(info[var]['org_units'], info[var])
            offset, factor = units_converter(new_units, info[var]['tgt_units'])
            final = out * factor + offset

            omean = np.mean(final)
            ostd = np.std(final)
            mf[season][region] = {}
            mf[season][region]['mean'] = float(str(round(omean,3)))
            mf[season][region]['std'] = float(str(round(ostd, 3)))
            print('{} {} {} mean is: {:.2f} +- {:.2f}'.format(var, season, region, omean, ostd))


    # log output
    logging.info(mf)

    # preparing clim file
    if os.path.isfile(clim_file):
        dclim = load_yaml(clim_file)
    else:
        dclim = {}

    # initialize variable if not exists
    if var not in dclim:
        dclim[var] = {}

    dclim[var]['longname'] = info[var]['longname']
    dclim[var]['dataset'] = info[var]['dataset']
    if mask_type != 'global':
        dclim[var]['mask'] = mask_type
    if operation != 'mean':
        dclim[var]['operation'] = operation
    if operation != 'atm':
        dclim[var]['domain'] = domain
    dclim[var]['units'] = info[var]['tgt_units']
    dclim[var]['year1'] = int(real_year1)
    dclim[var]['year2'] = int(real_year2)
    dclim[var]['obs'] = mf


    # dump the yaml file
    with open(clim_file, 'w') as file:
        yaml.safe_dump(dclim, file, sort_keys=False)
    

    #zfield = cfield.resample(time='1MS', skipna=nanskipper).mean('time', skipna=nanskipper)

