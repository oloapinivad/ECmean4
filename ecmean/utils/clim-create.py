#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""
    Tool to create a new ECmean4 climatology.
    It requires to have cdo and cdo-bindings installed
"""

__author__ = "Paolo Davini (paolo.davini@cnr.it), Sep 2022."

import logging
import os
import tempfile
from time import time
from dask.distributed import Client, LocalCluster

import matplotlib
import pandas as pd
import xarray as xr
import yaml
from cdo import *
#from dask.distributed import Client, LocalCluster, progress

#from ecmean.libs.climatology import full_histogram
from ecmean.libs.climatology import check_histogram, \
    mask_from_field, variance_threshold, variance_fraction, variance_iqr, \
    variance_iqr_adjusted, variance_combined,\
    select_time_period, timeframe_years, \
    parse_create_args, select_time_data, get_climatology_files, CLIMATOLOGY_PREFIXES
from ecmean.libs.files import load_yaml
from ecmean.libs.ncfixers import xr_preproc
from ecmean.libs.units import units_extra_definition

# activate CDO class
cdo = Cdo(logging=True)

# output for matplot lib
matplotlib.use('Agg')

# setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(name)s | %(levelname)8s -> %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# variable list
variables = ['tas', 'pr', 'net_sfc', 'tauu', 'tauv', 'psl',
             'ua', 'va', 'ta', 'hus', 'tos', 'sos', 'siconc']
#variables = ['tas']
#variables.reverse()

# target resolution
GRID = 'r360x180'

# method for variance filtering: "sigma", "fraction" or "iqr"
METHOD = "combined"
cutting = "clipping"
# method = "sigma"
SIGMA = 3
# fraction method: fraction of the median
FRACTION = 3

# skip NaN: if False, yearly/season average require that all
# the points are defined in the correspondent time window.
NANSKIP = False

# some dataset show very low variance in some grid point: this might create
# irrealistic high values of PI due to the  division by variance performend
# a hack is to use 5 sigma from the mean of the log10 distribution of variance
# define a couple of threshold to remove variance outliers


# add other units
units_extra_definition()


def main(climdata='EC26', timeframe='HIST', machine='wilma', do_figures=False, overwrite=False, method="sigma"):
    """Main function to create the climatology.
    
    Parameters
    ----------
    climdata : str
        Climate dataset name (default: 'EC26')
    timeframe : str
        Time period - HIST, PDAY, or CMIP (default: 'HIST')
    machine : str
        Machine name for data paths (default: 'wilma')
    do_figures : bool
        Generate diagnostic figures (default: False)
    overwrite : bool
        Overwrite existing climatology files (default: False)
    method : str
        Variance filtering method (default: 'sigma')   
    """
    FIGDIR = f'/scratch/users/paolo/ecmean-py-variances/{method}-{cutting}/'
    TMPDIR = '/scratch/users/paolo'

    # define the years from timeframe
    year1, year2 = timeframe_years(timeframe)
    climname = f'{climdata}-{timeframe}'

    # climatology yaml output
    tgtdir = os.path.join('..', 'climatology', climname)
    os.makedirs(tgtdir, exist_ok=True)
    clim_file = os.path.join(tgtdir, 'pi_climatology_' + climname + '.yml')
    # yml file to get information on dataset on some machine
    clim_info = f'create-clim-{machine}-{climdata}.yml'

    # figures directory
    if do_figures:
        figdir = os.path.join(FIGDIR, climname)
        os.makedirs(figdir, exist_ok=True)

    # always keep the attributes along the xarray
    xr.set_options(keep_attrs=True)

    # open the clim info file
    logging.debug("Loading climatology info from %s", clim_info)
    info = load_yaml(clim_info)

    # directory definitions and creations
    datadir = info['dirs']['datadir']
    archivedir = info['dirs']['archivedir']

    # loop on variables to be processed
    for var in variables:

        logging.warning('Processing variable: %s', var)

        #if var not in 'siconc':
        #    continue

        tic = time()
        # get the directory
        filedata = str(os.path.expandvars(info[var]['dir'])).format(
            datadir=datadir, archivedir=archivedir,
            dataset=info[var]['dataset'],
            varname=info[var]['varname'])
        logging.debug(filedata)

        # load data and time select
        logging.info("Loading multiple files...")
        # unable to operate with Parallel=True
        xfield = xr.open_mfdataset(filedata, chunks='auto',
                                   parallel=True, preprocess=xr_preproc, engine='netcdf4',
                                   data_vars='all', join='outer', compat='no_conflicts')
        xfield = xfield.rename({info[var]['varname']: var})

        # select time based on data availability
        cfield, real_year1, real_year2 = select_time_period(xfield, var, year1, year2)

        # check if files already exist and skip if overwrite is False
        if not overwrite:
            files_to_check = get_climatology_files(
                tgtdir, var, info[var]['dataset'], GRID, real_year1, real_year2
            ).values()
            if all(os.path.isfile(f) for f in files_to_check):
                logging.warning('All files for %s already exist. Skipping computation.', var)
                continue

        # check existence of unit, then apply from file
        if 'org_units' in info[var]:
            cfield.attrs['units'] = info[var]['org_units']
        elif 'units' in info[var]:
            # Backward compatibility
            cfield.attrs['units'] = info[var]['units']
        elif not hasattr(cfield, 'units'):
            raise ValueError('no unit found or defined!')

        logging.debug(cfield)

        # monthly average using resample/pandas
        logging.info("resampling...")
        zfield = cfield.resample(time='1MS', skipna=NANSKIP).mean('time', skipna=NANSKIP)

        #if do_figures:
        #    logging.info("Full histogram...")
        #    figname = f'values_{var}_{info[var]["dataset"]}_{real_year1}_{real_year2}_full.pdf'
        #    os.makedirs(os.path.join(figdir, var), exist_ok=True)
        #    file = os.path.join(figdir, var, figname)
        #    full_histogram(zfield, file)

        # dump the netcdf file to disk
        logging.info("new file...")
        tmpfile = tempfile.NamedTemporaryFile(suffix='.nc', dir=TMPDIR, delete=False)

        # preserve dtype for numerical reasons
        codes = ['dtype', '_FillValue', 'scale_factor', 'add_offset', 'missing_value']
        ftype = {k: v for k, v in cfield.encoding.items() if k in codes}
        logging.debug(ftype)
        zfield.to_netcdf(tmpfile.name, encoding={var: ftype})

        # create target directory
        os.makedirs(os.path.join(tgtdir, GRID), exist_ok=True)

        # use cdo to interpolate: call to attribute to exploit different interpolation
        logging.info("interpolation..")
        interpolator = getattr(cdo, info[var]['remap'])
        yfield = interpolator(GRID, input=tmpfile.name, returnXArray=var)

        # create empty lists
        d1 = []
        d2 = []

        # compute the yearly mean and the season mean
        logging.info("Averaging...")
        gfield1 = yfield.resample(time='YS', skipna=NANSKIP).mean('time', skipna=NANSKIP).persist()
        gfield2 = yfield.resample(time='QE-NOV', skipna=NANSKIP).mean('time', skipna=NANSKIP).persist()

        # loop on seasons
        for season in ['ALL', 'DJF', 'MAM', 'JJA', 'SON']:
            logging.info(season)

            gfield = select_time_data(gfield1, gfield2, season)

            logging.debug(gfield.shape)

            # zonal averaging for 3D fields
            if 'plev' in gfield.coords:
                gfield = gfield.mean(dim='lon')
                # select only up to 10hpa
                gfield = gfield.sel(plev=slice(100000, 1000))

            # create a reference time (average year, average month of the season)
            timestring = f'{int((year1 + year2) / 2)}-{str(gfield.time.dt.month.values[0])}-15'
            reftime = pd.to_datetime(timestring)

            # compute mean and variance: remove NaN in this case only
            omean = gfield.mean('time', skipna=True, keepdims=True)
            ovar = gfield.var('time', skipna=True, keepdims=True)

            #if do_figures:
            #    os.makedirs(os.path.join(figdir, var), exist_ok=True)
            #    omean.to_netcdf(os.path.join(figdir, var, f'mean_{season}.nc'))
            #    ovar.to_netcdf(os.path.join(figdir, var, f'var_{season}.nc'))

            # define the variance threshold
            if method == "sigma":
                low, center, high = variance_threshold(ovar, sigma=SIGMA)
                logging.info('Variance threshold: low = %s, center = %s, high = %s', low, center, high)
            elif method == "fraction":
                low, center, high = variance_fraction(ovar, fraction=FRACTION)
                logging.info('Variance fraction: low = %s, center = %s, high = %s', low, center, high)
            elif method == "iqr":
                low, center, high = variance_iqr(ovar)
                logging.info('Variance IQR: low = %s, center = %s, high = %s', low, center, high)
            elif method == "iqr_adjusted":
                low, center, high = variance_iqr_adjusted(ovar)
                logging.info('Variance IQR adjusted: low = %s, center = %s, high = %s', low, center, high)
            elif method == "combined":
                low, center, high = variance_combined(ovar, sigma=SIGMA, fraction=FRACTION)
                logging.info('Variance combined: low = %s, center = %s, high = %s', low, center, high)
            else:
                raise ValueError(f"Unknown method for variance filtering: {method}")

            # clipping
            if cutting == "clipping":
                logging.info('Applying variance clipping to minimum...')
                fvar = ovar.clip(min=low)
                fmean = omean
            elif cutting == "nan":
                logging.info('Applying variance sigma filtering...')
                fvar = ovar.where((ovar >= low) & (ovar <= high))
                fmean = omean.where((ovar >= low) & (ovar <= high))
            else:
                raise ValueError(f"Unknown method for variance filtering: {cutting}")

            if do_figures:
                logging.info("Mean and variance histograms...")
                figname = f'{var}_{info[var]["dataset"]}_{GRID}_{real_year1}_{real_year2}_{season}.pdf'
                os.makedirs(os.path.join(figdir, var), exist_ok=True)
                file = os.path.join(figdir, var, figname)
                check_histogram(
                    ovar, fvar, file, method=method, 
                    center=center, low=low, high=high)

            # add a reference time
            ymean = fmean.assign_coords({"time": ("time", [reftime])})
            yvar = fvar.assign_coords({"time": ("time", [reftime])})

            # append the dataset in the list
            d1.append(ymean)
            d2.append(yvar)
        
        # cleanup temporary file
        os.remove(tmpfile.name)

        # merge into a single dataarray
        season_mean = xr.concat(d1[1:], dim='time')
        season_variance = xr.concat(d2[1:], dim='time')
        full_mean = d1[0]
        full_variance = d2[0]

        # define compression and dtype for time, keep original dtype
        ftype["zlib"] = True
        compression = {var: ftype, 'time': {'dtype': 'f8'}}

        # save all climatology files (map climatology prefixes)
        climatology_data = [
            full_variance,
            full_mean,
            season_variance,
            season_mean
    ]
        
        climatology_files = get_climatology_files(
            tgtdir, var, info[var]['dataset'], GRID, real_year1, real_year2
        )

        for prefix, dataset in zip(CLIMATOLOGY_PREFIXES, climatology_data):
            logging.info('Saving %s...', prefix)
            dataset.to_netcdf(climatology_files[prefix], encoding=compression)

        toc = time()
        logging.warning('Processing in {:.4f} seconds'.format(toc - tic))

        # preparing clim file
        if os.path.isfile(clim_file):
            dclim = load_yaml(clim_file)
        else:
            dclim = {}

        # initialize variable if not exists
        if var not in dclim:
            dclim[var] = {}

        # assign to the dictionary the required info
        dclim[var]['dataset'] = info[var]['dataset']
        dclim[var]['description'] = info[var]['description']
        dclim[var]['longname'] = info[var]['longname']
        if 'version' in info[var].keys():
            dclim[var]['version'] = info[var]['version']
        # dclim[var]['dataname'] = info[var]['varname']
        dclim[var]['remap'] = info[var]['remap']
        dclim[var]['mask'] = mask_from_field(full_mean)
        dclim[var]['units'] = full_mean.attrs['units']
        dclim[var]['year1'] = int(real_year1)
        dclim[var]['year2'] = int(real_year2)

        # dump the yaml file
        with open(clim_file, 'w', encoding='utf8') as file:
            yaml.safe_dump(dclim, file, sort_keys=False)

        logging.debug(dclim)

    # Reorder climatology dictionary according to variable list order
    logging.info("Reordering climatology file according to variable list...")
    dclim = load_yaml(clim_file)
    dclim = {var: dclim[var] for var in variables if var in dclim}
    
    # dump the reordered yaml file
    with open(clim_file, 'w', encoding='utf8') as file:
        yaml.safe_dump(dclim, file, sort_keys=False)
    logging.warning("Climatology file reordered and saved: %s", clim_file)


# setting up dask
if __name__ == "__main__":
    parser = parse_create_args()
    parser.add_argument(
        '--overwrite',
        action='store_true',
        default=False,
        help='Not working yet: overwrite existing climatology files (default: False)'
    )
    parser.add_argument(
        '--figures',
        action='store_true',
        default=False,
        help='Generate diagnostic figures (default: False)'
    )
    parser.add_argument(
        '--method',
        type=str,
        default=METHOD,
        help='Variance filtering method (default: %(default)s)',
        choices=['sigma', 'fraction', 'iqr', 'iqr_adjusted', 'combined']
    )

    args = parser.parse_args()

    logging.getLogger().setLevel(args.loglevel.upper())

    if args.cores > 1:
        workers = args.cores
        cluster = LocalCluster(threads_per_worker=1, n_workers=workers)
        client = Client(cluster)
        logging.warning(client)

    main(args.climdata, args.timeframe, args.machine, args.figures, args.overwrite, args.method)

    if args.cores > 1:
        client.close()
        cluster.close()