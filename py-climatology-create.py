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
import matplotlib.pyplot as plt
from dask.distributed import Client, LocalCluster, progress
from ecmean import xr_preproc, load_yaml, units_extra_definition
import tempfile
cdo = Cdo()

# set default logging
logging.basicConfig(level=logging.INFO)

# variable list
variables = ['tas', 'pr', 'net_sfc', 'tauu', 'tauv', 'psl',
    'ua', 'va', 'ta', 'hus', 'tos', 'sos', 'siconc']

# to set: time period (default, can be shorter if data are missing) #WARNING MISSING
year1 = 1990
year2 = 2019

# yml file to get information on dataset on some machine
clim_info = '/home/paolo/ECmean4/climatology/create-clim-wilma-EC23.yml'

# figures : some diagnostic figures can be saved to check the consistency of mean and variance fields
do_figures = True
figdir = '/work/users/paolo/figures/ecmean-py-variances/'

# targets resolutio
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

# function to set absurd value from a specific dataset
def fix_specific_dataset(var, dataset, field) : 
    if var == 'tos' and dataset == 'ESA-CCI-L4' : 
        field = field.where(field > 5*10**-3)
    
    return field

def full_histogram(field, figname, n_bins = 100) : 

    fig, axs = plt.subplots(1,1, sharey=True, tight_layout=True, figsize=(15, 5))
    field.plot.hist(ax=axs, bins=n_bins, yscale = 'log')
    axs.title.set_text('Complete original values ' + field.name)
    plt.savefig(figname)

def check_histogram(ymean, yvar, yvar_filtered, figname, n_bins = 50) :

    fig, axs = plt.subplots(4,1, sharey=True, tight_layout=True, figsize=(20, 10))

    # log 10 fields
    f = np.log10(yvar.where(yvar>0))
    g = np.log10(yvar_filtered.where(yvar>0))

    # stats
    avg = f.mean()
    sss = 5*f.std()
    left = avg - sss
    right = avg + sss
    mmm = np.min(f)
    xxx = np.max(f)
    
    # mean and variance field 
    ymean.plot.hist(ax=axs[0], bins=n_bins, yscale = 'log')
    axs[0].title.set_text('Original Mean ' + yvar.name)
    yvar.plot.hist(ax=axs[1], bins=n_bins, yscale = 'log')
    axs[1].title.set_text('Original variance ' + yvar.name)

    # log10 plots
    f.plot.hist(ax=axs[2], bins=n_bins, yscale = 'log', xlim =[np.min([left, mmm]), np.max([right, xxx])])
    axs[2].title.set_text('Original variance log10 ' + yvar.name)
    g.plot.hist(ax=axs[3], bins=n_bins, yscale = 'log', color = 'red', xlim = [np.min([left, mmm]), np.max([right, xxx])])
    axs[3].title.set_text('Filtered variance log10 ' + yvar.name)
    for k in [2,3] : 
        axs[k].axvline(avg, color='k', linewidth=1)
        axs[k].axvline(left, color='k', linestyle='dashed', linewidth=1)
        axs[k].axvline(right, color='k', linestyle='dashed', linewidth=1)

    plt.savefig(figname)
 

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
        sys.exit('ERROR: cant recognize mask')

    logging.info(mask)
    return mask

# add other units
units_extra_definition()

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
        mfield = xfield[var].sel(time=xfield.time.dt.year.isin(years))

        # specific fixer for each varialbe
        cfield = fix_specific_dataset(var, info[var]['dataset'], mfield)

        real_year1 = np.min(cfield.time.dt.year.values)
        real_year2 = np.max(cfield.time.dt.year.values)
        if (real_year1 != year1) :
            logging.warning("Initial year different from what expected: " + str(real_year1))
        if (real_year2 != year2) :
            logging.warning("Final year different from what expected: " + str(real_year2))

        if do_figures: 
            figname = 'values_' + var + '_' + info[var]['dataset'] + '_' +  str(real_year1) + '_' + str(real_year2) + '.pdf'
            os.makedirs(os.path.join(figdir, var), exist_ok=True)
            file = os.path.join(figdir, var, figname)
            full_histogram(mfield, file)

        # check existence of unit, then apply from file
        if 'units' in info[var] : 
            cfield.attrs['units'] = info[var]['units']
        elif not hasattr(cfield, 'units') :
            sys.exit('ERROR: no unit found or defined!')

        logging.info(cfield)

        # monthly average using resample/pandas
        print("resampling...")
        zfield = cfield.resample(time='1MS').mean('time')

        print("computation...")
        zfield = zfield.persist()
        progress(zfield)
        zfield.compute()
        
        # dump the netcdf file to disk
        print("new file...") 
        temp_name = next(tempfile._get_candidate_names()) + '.nc'
        tmpout = os.path.join(tmpdir, temp_name)
        zfield.to_netcdf(tmpout)

        # loop on grids
        for grid in grids : 

            # create target directory
            os.makedirs(os.path.join(tgtdir, grid), exist_ok=True)

            # use cdo to interpolate: call to attribute to exploit different interpolation
            print("interpolation..") 
            interpolator = getattr(cdo,  info[var]['remap'])
            yfield = interpolator(grid, input = tmpout, returnXArray = var)
            os.remove(tmpout)

            # create empty lists
            d1 = []
            d2 = []

            # loop on seasons
            for season in ['ALL', 'DJF', 'MAM', 'JJA', 'SON'] :
                print(season)

                # select the season
                if season == 'ALL' :
                    # yearly average
                    gfield = yfield.resample(time='AS').mean('time')
                else :
                    # season average
                    gfield = yfield.resample(time="Q-NOV").mean('time')
                    gfield = gfield.sel(time=gfield.time.dt.season.isin(season))
                    # for winter, we drop first and last to have only complete season. 
                    # this reduces the sample by one but it is safer especially for variance
                    if season == 'DJF' :
                        gfield = gfield.drop_isel(time=[0, gfield.sizes['time']-1])
                
                logging.info(gfield.shape)

                # zonal averaging for 3D fields
                if 'plev' in gfield.coords :
                    gfield = gfield.mean(dim = 'lon') 

                # create a reference time (average year, average month of the season)
                reftime = pd.to_datetime(str(int((year1+year2)/2)) + '-' + 
                    str(gfield.time.dt.month.values[0]) + '-15')

                # compute mean and variance
                ymean0 = gfield.mean('time', keepdims=True)
                yvar0 = gfield.var('time', skipna=True, keepdims=True)
                
                # define the variance threshold
                low, high = variance_threshold(yvar0)
                logging.info(low)
                logging.info(high)

                # clean according to thresholds
                yvar = yvar0.where((yvar0 >= low) & (yvar0 <= high))
                ymean = ymean0.where((yvar0 >= low) & (yvar0 <= high))

                if do_figures: 
                    figname = var + '_' + info[var]['dataset'] + '_' + grid + '_' +  str(real_year1) + '_' + str(real_year2) + '_' + season + '.pdf'
                    os.makedirs(os.path.join(figdir, var), exist_ok=True)
                    file = os.path.join(figdir, var, figname)
                    check_histogram(ymean0, yvar0, yvar, file)


                # add a reference time
                ymean = ymean.assign_coords({"time": ("time",  [reftime])})
                yvar = yvar.assign_coords({"time": ("time",  [reftime])})

                # append the dataset in the list
                d1.append(ymean)
                d2.append(yvar)

            # merge into a single dataarray
            season_mean = xr.concat(d1[1:], dim = 'time')
            season_variance = xr.concat(d2[1:], dim = 'time')
            full_mean = d1[0]
            full_variance = d2[0]
            
            # define compression and dtype for time
            compression = {var: {"zlib": True, '_FillValue': -999.0}, 'time': {'dtype': 'f8'}}

            # define file suffix
            suffix = var + '_' + info[var]['dataset'] + '_' + grid + '_' +  str(real_year1) + '-' + str(real_year2) + '.nc'

            # save full - standard format
            full_variance.to_netcdf(os.path.join(tgtdir, grid, 'climate_variance_' + suffix), encoding = compression)
            full_mean.to_netcdf(os.path.join(tgtdir, grid, 'climate_average_' + suffix), encoding = compression)

             # save season - 4 season format
            season_variance.to_netcdf(os.path.join(tgtdir, grid, 'seasons_variance_' + suffix), encoding = compression)
            season_mean.to_netcdf(os.path.join(tgtdir, grid, 'seasons_average_' + suffix), encoding = compression)

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
            #dclim[var]['dataname'] = info[var]['varname']
            dclim[var]['remap'] = info[var]['remap']
            dclim[var]['mask'] = mask_from_field(full_mean)
            dclim[var]['units'] = full_mean.attrs['units']
            dclim[var]['year1'] = int(real_year1)
            dclim[var]['year2'] = int(real_year2)

            # dump the yaml file
            with open(clim_file, 'w') as file:
                yaml.safe_dump(dclim, file, sort_keys=False)
            
            logging.debug(dclim)

# setting up dask
if __name__ == "__main__":

    # set up clusters
    cluster = LocalCluster(threads_per_worker=threads, n_workers = workers)
    client = Client(cluster)
    logging.warning(client)
    main()
        
