#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""
    Tool to create a new ECmean4 climatology
"""

__author__ = "Paolo Davini (p.davini@isac.cnr.it), Sep 2022."

import yaml
import os
import xarray as xr
from cdo import * 
import numpy as np
import pandas as pd
import logging
from time import time
from dask.distributed import Client, LocalCluster, progress
from ecmean import xr_preproc, load_yaml
cdo = Cdo()

# set default logging
logging.basicConfig(level=logging.INFO)

# variable list
variables = ['tas', 'pr', 'net_sfc', 'tauu', 'tauv',
    'ua', 'va', 'ta', 'hus', 'tos', 'sos', 'sic']

# to set: time period (default, can be shorter if data are missing) #WARNING MISSING
year1 = 1990
year2 = 2019

# yml file to get information on dataset on some machine
clim_info = '/home/paolo/ECmean4/climatology/create-clim-wilma-EC23.yml'

# targets resolution and years
grids = ['r360x180']

# number of dask workes
workers = 4
threads = 2

# some dataset show very low variance in some grid point: this might create
# irrealistic high values of PI due to the  division by variance performend
# a hack is to use 5 sigma from the mean of the log10 distribution of variance
# define a couple of threshold to remove variance outliers
def variance_threshold(xvariance) : 
    f = np.log10(xvariance.where(xvariance>0))
    m = float(np.mean(f).values)
    s = float(np.std(f).values)
    low = 10**(m-5*s)
    high = 10**(m+5*s)
    return low, high

# get domain of the variable from the fraction of NaN: UNDER TESTING
def mask_from_field(xfield) :
    ratio = float(xfield.count() / np.prod(np.array(xfield.shape)))
    logging.info(ratio)
    if ratio < 0.2 : # this is a special case for ice, need to be double checked
        mask = 'ocean'
    elif 0.2 < ratio < 0.3 :
        mask = 'land'
    elif 0.6 < ratio < 0.7 : 
        mask = 'ocean'
    elif ratio > 0.95 : 
        mask = 'global'
    else : 
        mask = 'undefined'

    logging.info(mask)
    return mask

# to exploit of dask we need a main function 
def main() : 

    # always keep the attributes along the xarray
    xr.set_options(keep_attrs=True)

    # open the clim info file
    info = load_yaml(clim_info)

    # set few parameters
    clim_name = info['clim']
    years = list(range(year1, year2 + 1))

    # directory definitions and creations
    tmpdir = info['dirs']['tmpdir']
    tgtdir = info['dirs']['tgtdir'].format(clim = clim_name)
    datadir = info['dirs']['datadir']
    for dir in [tmpdir, tgtdir, datadir] :
        os.makedirs(dir, exist_ok=True)

    # climatology yaml output
    clim_file = os.path.join(tgtdir, 'pi_climatology_' + clim_name + '.yml')

    # loop on variables to be processed
    for var in variables : 

        print(var)
        tic = time()
        # get the directory
        filedata = str(os.path.expandvars(info[var]['dir'])).format(
            datadir = info['dirs']['datadir'],
            dataset = info[var]['dataset'],
            varname = info[var]['varname'])
        logging.info(filedata)


        # load data and time select
        print("Loading multiple files...")
        xfield = xr.open_mfdataset(filedata, chunks='auto', preprocess=xr_preproc)
        xfield = xfield.rename({info[var]['varname']: var})
        cfield = xfield[var].sel(time=xfield.time.dt.year.isin(years))

        real_year1 = np.min(cfield.time.dt.year.values)
        real_year2 = np.max(cfield.time.dt.year.values)
        if (real_year1 != year1) :
            logging.warning("Initial year different from what expected: " + str(real_year1))
        if (real_year2 != year2) :
            logging.warning("Final year different from what expected: " + str(real_year2))

        # check existence of unit, then apply from file
        if not hasattr(cfield, 'units') :
            if info[var]['units'] : 
                cfield.attrs['units'] = info[var]['units']
            else : 
                sys.exit('ERROR: no unit found or defined!')

        print(cfield)

        # season average using resample/pandas
        print("resampling...")
        zfield = cfield.resample(time="Q-NOV").mean('time')

        print("computation...")
        zfield = zfield.persist()
        progress(zfield)
        zfield.compute()
        
        # dump the netcdf file to disk
        print("new file...") 
        tmpout = os.path.join(tmpdir, 'tmp.nc')
        zfield.to_netcdf(tmpout)

        # loop on grids
        for grid in grids : 

            # create target directory
            os.makedirs(os.path.join(tgtdir, grid), exist_ok=True)

            # use cdo to interpolate: call to attribute to exploit different interpolation
            print("interpolation..") 
            interpolator = getattr(cdo,  info[var]['remap'])
            yfield = interpolator(grid, input = tmpout, returnXArray = var)
            #yfield = cdo.remapbil(grid, input = tmpout, returnXArray = info[var]['varname'])
            os.remove(tmpout)

            # create empty lists
            d1 = []
            d2 = []

            # loop on seasons
            for season in ['DJF', 'MAM', 'JJA', 'SON'] :
                print(season)

                # select the season 
                gfield = yfield.sel(time=yfield.time.dt.season.isin(season))
      
                # zonal averaging for 3D fields
                if 'plev' in gfield.coords :
                    gfield = gfield.mean(dim = 'lon') 

                # create a reference time (average year, average month of the season)
                reftime = pd.to_datetime(str(int((year1+year2)/2)) + '-' + 
                    str(gfield.time.dt.month.values[0]) + '-15')

                # compute mean and variance
                ymean = gfield.mean('time', keepdims=True)
                yvar = gfield.var('time', skipna=True, keepdims=True)

                # define the variance threshold
                low, high = variance_threshold(yvar)
                logging.info(low)
                logging.info(high)

                # clean according to thresholds
                yvar = yvar.where((yvar >= low) & (yvar <= high))
                ymean = ymean.where((yvar >= low) & (yvar <= high))

                # add a reference time
                ymean = ymean.assign_coords({"time": ("time",  [reftime])})
                yvar = yvar.assign_coords({"time": ("time",  [reftime])})

                # append the dataset in the list
                d1.append(ymean)
                d2.append(yvar)

            # merge into a single dataarray
            mean = xr.concat(d1, dim = 'time')
            variance = xr.concat(d2, dim = 'time')
            
            # define compression and dtype for time
            compression = {var: {"zlib": True, '_FillValue': -999.0}, 'time': {'dtype': 'f8'}}

            # save
            suffix = var + '_' + info[var]['dataset'] + '_' + grid + '_' +  str(real_year1) + '_' + str(real_year2) + '.nc'
            variance.to_netcdf(os.path.join(tgtdir, grid, 'variance_' + suffix), encoding = compression)
            mean.to_netcdf(os.path.join(tgtdir, grid, 'climate_' + suffix), encoding = compression)

            toc = time()
            print('Processing in {:.4f} seconds'.format(toc - tic))

            # preparing clim file
            if os.path.isfile(clim_file) : 
                dclim = load_yaml(clim_file)
            else :
                dclim = {}
            
            # initialize variable if not exists
            if not var in dclim : 
                dclim[var] = {}

            # assign to the dictionary the required info
            dclim[var]['dataset'] = info[var]['dataset']
            dclim[var]['dataname'] = info[var]['varname']
            dclim[var]['remap'] = info[var]['remap']
            dclim[var]['mask'] = mask_from_field(mean)
            dclim[var]['units'] = mean.attrs['units']
            dclim[var]['year1'] = int(real_year1)
            dclim[var]['year2'] = int(real_year2)

            # dump the yaml file
            with open(clim_file, 'w') as file:
                yaml.safe_dump(dclim, file, sort_keys=False)
            
            print(dclim)

# setting up dask
if __name__ == "__main__":

    # set up clusters
    cluster = LocalCluster(threads_per_worker=threads, n_workers = workers)
    client = Client(cluster)
    logging.warning(client)
    main()
        
