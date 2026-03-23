#!/usr/bin/env python3
'''
Common functions for reference and climatology creation.
'''

import argparse
import os
import logging
import numpy as np
import dask.array as da
import matplotlib.pyplot as plt

loggy = logging.getLogger(__name__)


# Climatology file prefixes
CLIMATOLOGY_PREFIXES = [
    'climate_variance_',
    'climate_average_',
    'seasons_variance_',
    'seasons_average_'
]

# supported climatology
SUPPORTED_REFERENCE = ['EC23', 'EC26-PDAY', 'EC26-CMIP', 'EC26-HIST']
SUPPORTED_CLIMATOLOGY = ['EC23', 'EC24', 'EC26-HIST', 'EC26-CMIP']

def timeframe_years(timeframe):
    """Timeframe to years mapping."""
    if timeframe == 'HIST':
        year1 = 1981
        year2 = 2010
    elif timeframe == 'PDAY':
        year1 = 2000
        year2 = 2024
    elif timeframe == 'CMIP':
        year1 = 1985
        year2 = 2014
    else:
        raise ValueError(f"TIMEFRAME {timeframe} not recognized."
                         "Please choose from HIST, PDAY, CMIP.")
    return year1, year2

def select_time_data(yearly, seasonly, season):
    """Select time data based on season."""

    if season == 'ALL':
        return yearly

    gfield = seasonly.sel(time=seasonly.time.dt.season.isin(season))
    # for winter, we drop first and last to have only complete season.
    # this reduces the sample by one but it is safer for variance
    if season == 'DJF':
        gfield = gfield.drop_isel(time=[0, gfield.sizes['time'] - 1])

    return gfield


def get_climatology_files(tgtdir, var, dataset, grid, year1, year2):
    """Generate climatology file paths.
    """

    suffix = f'{var}_{dataset}_{grid}_{year1}-{year2}.nc'

    climatology_files = {
        prefix: os.path.join(tgtdir, grid, prefix + suffix)
        for prefix in CLIMATOLOGY_PREFIXES
    }

    return climatology_files


def expand_filedata(directory, var, info):
    """Expand filedata directory with environment variables and info dictionary."""

    return os.path.expandvars(directory).format(
        datadir=info['dirs']['datadir'], mswepdir=info['dirs']['mswepdir'],
        eradir=info['dirs']['eradir'], esadir=info['dirs']['esadir'],
        dataset=info[var]['dataset']
    )

def select_time_period(xfield, var, year1, year2):
    """
    Select time period based on data availability.
    """
    # Check available years
    valid_time = xfield[var].dropna("time", how="all").time
    avail_years = np.unique(valid_time.dt.year.values)
    
    # Calculate intersection between requested and available
    real_year1 = max(year1, int(avail_years.min()))
    real_year2 = min(year2, int(avail_years.max()))

    # Raise error if incompatible
    if real_year1 > real_year2:
        raise ValueError(
            f"{var}: requested period {year1}-{year2} "
            f"not compatible with data availability "
            f"{avail_years.min()}-{avail_years.max()}"
        )

    # Warn if years differ from requested
    if real_year1 != year1 or real_year2 != year2:
        logging.warning(
            "%s: using %s-%s instead of requested %s-%s",
            var, real_year1, real_year2, year1, year2
        )

    # Select time
    years_eff = list(range(real_year1, real_year2 + 1))
    cfield = xfield[var].sel(time=xfield[var].time.dt.year.isin(years_eff))

    return cfield, real_year1, real_year2

def variance_threshold(xvariance, sigma=5):
    """
    This defines the two thresholds (high and low) for filtering the dataset
    So far it is done on the 5-std of the log10 distribution
    """
    f = np.log10(xvariance.where(xvariance > 0))
    m = float(f.mean())
    loggy.debug('Median log10 variance: %s', m)
    s = float(f.std())
    low = 10**(m - sigma * s)
    high = 10**(m + sigma * s)
    return low, high

def variance_iqr(xvariance):
    """This defines the two thresholds (high and low) for filtering the dataset
    based on the interquartile range of the log10 distribution."""

    f = np.log10(xvariance.where(xvariance > 0))
    qqq = f.quantile([0.25, 0.75])
    iqr = qqq[1] - qqq[0]
    iqleft = 10**(qqq[0] - 1.5 * iqr)
    iqright = 10**(qqq[1] + 1.5 * iqr)
    return iqleft.values, iqright.values

def variance_fraction(xvariance, fraction=1e-3):
    """
    Alternative method for variance clipping, based on the median of the distribution. 
    The threshold is defined as a fraction of the median.
    """
    f = xvariance.where(xvariance > 0)
    m = float(f.median())
    loggy.debug('Median variance: %s', m)
    low = fraction * m
    high = m / fraction
    return low, high

def full_histogram(field, figname, n_bins=100):
    """Compute the histogram of the full field before it is processed.
    this is done to check the presence of irrealistic values within the dataset
    dask.array.histogram is used to speed up the computation."""

    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True, figsize=(15, 5))

    # test using underlying dask array
    mmm = da.nanmin(field.data).compute()
    xxx = da.nanmax(field.data).compute()
    extra = (xxx - mmm) / 20
    hist, bins = da.histogram(field.data, bins=n_bins, range=[mmm - extra, xxx + extra])
    x = 0.5 * (bins[1:] + bins[:-1])
    width = np.diff(bins)
    axs.bar(x, hist.compute(), width, log=True)
    axs.title.set_text('Complete original values ' + field.name)
    fig.savefig(figname)


def check_histogram(ymean, yvar, yvar_filtered, figname, n_bins=100, sigma=5, fraction=1e-2):
    """Four histograms made for inspection of mean and variance of the field
    Mean field, variance and variance after filtering are passed and then plotted
    using histograms. log10 scales is used to highlight outliers."""

    fig, axs = plt.subplots(4, 1, sharey=True, tight_layout=True, figsize=(20, 15))

    # log 10 fields
    f = np.log10(yvar.where(yvar > 0))
    g = np.log10(yvar_filtered.where(yvar > 0))

    # stats
    avg = f.mean()
    median = yvar.where(yvar > 0).median()

    sss = sigma * f.std()
    qqq = f.quantile([0.25, 0.75])
    iqr = qqq[1] - qqq[0]
    iqleft = qqq[0] - 1.5 * iqr
    iqright = qqq[1] + 1.5 * iqr
    left = np.min([avg - sss, f.min(skipna=True)])
    right = np.max([avg + sss, f.max(skipna=True)])
    extra = abs(left - right) / 20
    # print([avg, sss, left, right, left - extra, right - extra])

    # mean and variance field
    ymean.plot.hist(ax=axs[0], bins=n_bins, yscale='log', color="goldenrod")
    axs[0].title.set_text('Original Mean ' + yvar.name)
    yvar.plot.hist(ax=axs[1], bins=n_bins, yscale='log')
    axs[1].title.set_text('Original variance ' + yvar.name)

    # log10 plots
    f.plot.hist(ax=axs[2], bins=n_bins, yscale='log', xlim=[left - extra, right + extra])
    axs[2].title.set_text('Original variance log10 ' + yvar.name)
    g.plot.hist(ax=axs[3], bins=n_bins, yscale='log', color='red', xlim=[left - extra, right + extra])
    axs[3].title.set_text('Filtered variance log10 ' + yvar.name)
    for k in [2, 3]:
        axs[k].axvline(avg, color='k', linewidth=1)
        axs[k].axvline(avg - sss, color='k', linestyle='dashed', linewidth=1)
        axs[k].axvline(avg + sss, color='k', linestyle='dashed', linewidth=1)
        axs[k].axvline(np.log10(median), color='g', linewidth=1)
        axs[k].axvline(np.log10(median/fraction), color='g', linestyle='dashed', linewidth=1)
        axs[k].axvline(np.log10(median*fraction), color='g', linestyle='dashed', linewidth=1)

        axs[k].axvline(iqleft, color='r', linestyle='dashed', linewidth=1)
        axs[k].axvline(iqright, color='r', linestyle='dashed', linewidth=1)

    fig.savefig(figname)


# get domain of the variable from the fraction of NaN: UNDER TESTING
def mask_from_field(xfield):
    """get the domain to be passed to the climatology .yml file from the number of
    missing point. Special treatment for sea ice. Use with caution."""
    ratio = float(xfield.count() / np.prod(np.array(xfield.shape)))
    loggy.info(ratio)
    if ratio < 0.2:  # this is a special case for ice, need to be double checked
        mask = 'ocean'
    elif 0.2 < ratio < 0.3:
        mask = 'land'
    elif 0.55 < ratio < 0.7:
        mask = 'ocean'
    elif ratio > 0.95:
        mask = 'global'
    else:
        mask = 'undefined'
        raise ValueError('ERROR: cant recognize mask')

    loggy.debug(mask)
    return mask


def parse_create_args():
    """Parse command line arguments for reference/climatology creation."""
    parser = argparse.ArgumentParser(
        description='Create ecmean reference/climatology for global mean/performance indices.'
    )
    parser.add_argument(
        '-c', '--climdata',
        type=str,
        default='EC26',
        choices=['EC26', 'EC24'],
        help='Climatology/reference dataset name (default: EC26)'
    )
    parser.add_argument(
        '-t', '--timeframe',
        type=str,
        default='HIST',
        choices=['HIST', 'PDAY', 'CMIP'],
        help='Time period for the climatology (default: HIST)'
    )
    parser.add_argument(
        '--machine',
        type=str,
        default='wilma',
        help='Machine name for input data path (default: wilma)'
    )
    parser.add_argument(
        '-j', '--cores',
        type=int,
        default=1,
        help='Number of cores for parallel processing (default: 1). ' \
        'Two or more cores will activate Dask for parallel computation.'
    )
    parser.add_argument(
        '-l', '--loglevel',
        type=str.upper,
        default='WARNING',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        help='Logging level (default: WARNING)'
    )
    return parser
