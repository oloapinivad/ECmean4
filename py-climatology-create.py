#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""
    Tool to create a new ECmean4 climatology. 
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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dask.distributed import Client, LocalCluster, progress
import dask.array as da
from ecmean import xr_preproc, load_yaml, units_extra_definition
import tempfile
cdo = Cdo(logging=True)

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

# number of dask workers and threads
workers = 8
threads = 1

# skip NaN: if False, yearly/season average require that all the points are defined in the correspondent 
# time window. 
nanskipper = False

# some dataset show very low variance in some grid point: this might create
# irrealistic high values of PI due to the  division by variance performend
# a hack is to use 5 sigma from the mean of the log10 distribution of variance
# define a couple of threshold to remove variance outliers
def variance_threshold(xvariance) : 
    """this defines the two thresholds (high and low) for filtering the dataset
    So far it is done on the 5-sigma of the log10 distribution"""
    f = np.log10(xvariance.where(xvariance>0))
    m = float(np.mean(f).values)
    s = float(np.std(f).values)
    low = 10**(m-5*s)
    high = 10**(m+5*s)
    return low, high

# function to set absurd value from a specific dataset
def fix_specific_dataset(var, dataset, xfield) : 
    """in the case some dataset show unexpected values this can be filtered here"""
    #variance of SST under sea ice is almost zero. We need to get rid of those points
    if var == 'tos' and dataset == 'ESA-CCI-L4' : 
        #xfield = xfield.where(xfield > 271.15)
        xfield = xfield.where(xfield > 1*10**-2)   
    return xfield

def full_histogram(field, figname, n_bins = 100) : 
    """Compute the histogram of the full field before it is processed.
    this is done to check the presence of irrealistic values within the dataset
    dask.array.histogram is used to speed up the computation."""

    fig, axs = plt.subplots(1,1, sharey=True, tight_layout=True, figsize=(15, 5))
    
    # test using underlying dask array
    mmm = da.nanmin(field.data).compute()
    xxx = da.nanmax(field.data).compute()
    extra = (xxx - mmm) / 20
    hist, bins = da.histogram(field.data, bins= n_bins, range=[mmm - extra, xxx + extra])
    x = 0.5 * (bins[1:] + bins[:-1])
    width = np.diff(bins)
    axs.bar(x, hist.compute(), width,  log = True)
    axs.title.set_text('Complete original values ' + field.name)
    fig.savefig(figname)

def check_histogram(ymean, yvar, yvar_filtered, figname, n_bins = 100) :
    """Four histograms made for inspection of mean and variance of the field
    Mean field, variance and variance after filtering are passed and then plotted
    using histograms. log10 scales is used to highlight outliers."""

    fig, axs = plt.subplots(4,1, sharey=True, tight_layout=True, figsize=(20, 15))

    # log 10 fields
    f = np.log10(yvar.where(yvar>0))
    g = np.log10(yvar_filtered.where(yvar>0))

    # stats
    avg = f.mean()
    sss = 5*f.std()
    qqq = f.quantile([0.25, 0.75])
    iqr = qqq[1] - qqq[0]
    iqleft =  qqq[0] - 1.5 * iqr
    iqright = qqq[1] + 1.5 * iqr
    left = np.min([avg - sss, f.min(skipna=True)])
    right = np.max([avg + sss, f.max(skipna=True)])
    extra = abs(left - right) / 20
    #print([avg, sss, left, right, left - extra, right - extra])
    
    # mean and variance field 
    ymean.plot.hist(ax=axs[0], bins=n_bins, yscale = 'log', color = "goldenrod")
    axs[0].title.set_text('Original Mean ' + yvar.name)
    yvar.plot.hist(ax=axs[1], bins=n_bins, yscale = 'log')
    axs[1].title.set_text('Original variance ' + yvar.name)

    # log10 plots
    f.plot.hist(ax=axs[2], bins=n_bins, yscale = 'log', xlim =[left - extra, right + extra])
    axs[2].title.set_text('Original variance log10 ' + yvar.name)
    g.plot.hist(ax=axs[3], bins=n_bins, yscale = 'log', color = 'red', xlim = [left - extra, right + extra])
    axs[3].title.set_text('Filtered variance log10 ' + yvar.name)
    for k in [2,3] : 
        axs[k].axvline(avg, color='k', linewidth=1)
        axs[k].axvline(avg - sss, color='k', linestyle='dashed', linewidth=1)
        axs[k].axvline(avg + sss, color='k', linestyle='dashed', linewidth=1)
        axs[k].axvline(iqleft, color='r', linestyle='dashed', linewidth=1)
        axs[k].axvline(iqright, color='r', linestyle='dashed', linewidth=1)

    fig.savefig(figname)
 

# get domain of the variable from the fraction of NaN: UNDER TESTING
def mask_from_field(xfield) :
    """get the domain to be passed to the climatology .yml file from the number of 
    missing point. Special treatment for sea ice. Use with caution."""
    ratio = float(xfield.count() / np.prod(np.array(xfield.shape)))
    logging.info(ratio)
    if ratio < 0.2 : # this is a special case for ice, need to be double checked
        mask = 'ocean'
    elif 0.2 < ratio < 0.3 :
        mask = 'land'
    elif 0.55 < ratio < 0.7 : 
        mask = 'ocean'
    elif ratio > 0.95 : 
        mask = 'global'
    else : 
        mask = 'undefined'
        sys.exit('ERROR: cant recognize mask')

    logging.debug(mask)
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
        # unable to operate with Parallel=True
        xfield = xr.open_mfdataset(filedata, chunks='auto', 
            parallel=False, preprocess=xr_preproc, engine='netcdf4')
        xfield = xfield.rename({info[var]['varname']: var})
        cfield = xfield[var].sel(time=xfield.time.dt.year.isin(years))

        # keep original dtype
        #cfield = cfield.astype(cfield.encoding['dtype'])

        real_year1 = np.min(cfield.time.dt.year.values)
        real_year2 = np.max(cfield.time.dt.year.values)
        if (real_year1 != year1) :
            logging.warning("Initial year different from what expected: " + str(real_year1))
        if (real_year2 != year2) :
            logging.warning("Final year different from what expected: " + str(real_year2))

        # check existence of unit, then apply from file
        if 'units' in info[var] : 
            cfield.attrs['units'] = info[var]['units']
        elif not hasattr(cfield, 'units') :
            sys.exit('ERROR: no unit found or defined!')

        # cleaning
        #cfield = fix_specific_dataset(var, info[var]['dataset'], cfield)
        logging.info(cfield)

        # monthly average using resample/pandas
        print("resampling...")
        zfield = cfield.resample(time='1MS', skipna=nanskipper).mean('time', skipna=nanskipper)
        zfield = zfield.persist()
        progress(zfield)
        #zfield.compute()

        if do_figures: 
            print("Full histogram...")
            figname = 'values_' + var + '_' + info[var]['dataset'] + '_' +  str(real_year1) + '_' + str(real_year2) + '.pdf'
            os.makedirs(os.path.join(figdir, var), exist_ok=True)
            file = os.path.join(figdir, var, figname)
            full_histogram(zfield, file)
        
        # dump the netcdf file to disk
        print("new file...") 
        temp_name = next(tempfile._get_candidate_names()) + '.nc'
        tmpout = os.path.join(tmpdir, temp_name)

        # preserve dtype for numerical reasons
        codes = ['dtype', '_FillValue', 'scale_factor', 'add_offset', 'missing_value']
        ftype = { k:v for k,v in cfield.encoding.items() if k in codes }
        logging.info(ftype)
        zfield.to_netcdf(tmpout, encoding = {var : ftype})

        # loop on grids
        for grid in grids : 

            # create target directory
            os.makedirs(os.path.join(tgtdir, grid), exist_ok=True)

            # use cdo to interpolate: call to attribute to exploit different interpolation
            print("interpolation..") 
            interpolator = getattr(cdo,  info[var]['remap'])
            yfield = interpolator(grid, input = tmpout, returnXArray = var)

            # keep original dtype
            #yfield = yfield.astype(yfield.encoding['dtype'])
            #print(yfield)
            os.remove(tmpout)

            # create empty lists
            d1 = []
            d2 = []

            # compute the yearly mean and the season mean
            print("Averaging...")
            gfield1 = yfield.resample(time='AS', skipna=nanskipper).mean('time', skipna=nanskipper) 
            gfield2 = yfield.resample(time='Q-NOV', skipna=nanskipper).mean('time', skipna=nanskipper) 

            # loop on seasons
            for season in ['ALL', 'DJF', 'MAM', 'JJA', 'SON'] :
                print(season)

                # select the season
                if season == 'ALL' :
                    gfield = gfield1
                else :
                    gfield = gfield2.sel(time=gfield2.time.dt.season.isin(season))
                    # for winter, we drop first and last to have only complete season. 
                    # this reduces the sample by one but it is safer for variance
                    if season == 'DJF' :
                        gfield = gfield.drop_isel(time=[0, gfield.sizes['time']-1])
                
                logging.info(gfield.shape)

                # zonal averaging for 3D fields
                if 'plev' in gfield.coords :
                    gfield = gfield.mean(dim = 'lon') 

                # create a reference time (average year, average month of the season)
                reftime = pd.to_datetime(str(int((year1+year2)/2)) + '-' + 
                    str(gfield.time.dt.month.values[0]) + '-15')

                # compute mean and variance: remove NaN in this case only
                omean = gfield.mean('time', skipna=True, keepdims=True)
                ovar = gfield.var('time', skipna=True, keepdims=True)

                # special fix
                #ovar = fix_specific_dataset(var, info[var]['dataset'], ovar)
                
                # define the variance threshold
                low, high = variance_threshold(ovar)
                logging.info(low)
                logging.info(high)

                # clean according to thresholds
                fvar= ovar.where((ovar >= low) & (ovar <= high))
                fmean = omean.where((ovar >= low) & (ovar <= high))

                print(fvar)

                if do_figures: 
                    print("Mean and variance histograms...")
                    figname = var + '_' + info[var]['dataset'] + '_' + grid + '_' +  str(real_year1) + '_' + str(real_year2) + '_' + season + '.pdf'
                    os.makedirs(os.path.join(figdir, var), exist_ok=True)
                    file = os.path.join(figdir, var, figname)
                    check_histogram(omean, ovar, fvar, file)

                # add a reference time
                ymean = fmean.assign_coords({"time": ("time",  [reftime])})
                yvar = fvar.assign_coords({"time": ("time",  [reftime])})

                # append the dataset in the list
                d1.append(ymean)
                d2.append(yvar)

            # merge into a single dataarray
            season_mean = xr.concat(d1[1:], dim = 'time')
            season_variance = xr.concat(d2[1:], dim = 'time')
            full_mean = d1[0]
            full_variance = d2[0]
            
            # define compression and dtype for time, keep original dtype
            ftype["zlib"] = True
            compression = {var: ftype, 'time': {'dtype': 'f8'}}

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
        
