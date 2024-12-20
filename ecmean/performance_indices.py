#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
 python3 version of ECmean performance indices tool
 Using a reference file from yaml and cdo bindings

 @author Paolo Davini (p.davini@isac.cnr.it), 2022
 @author Jost von Hardenberg (jost.hardenberg@polito.it), 2022
'''

import sys
import os
import logging
from time import time
from multiprocessing import Process, Manager
import numpy as np
import xarray as xr
import yaml
import dask
from ecmean import Diagnostic, Supporter, UnitsHandler
from ecmean.libs.general import weight_split, get_domain, \
    check_time_axis, init_mydict, check_var_interface, check_var_climatology, \
    set_multiprocessing_start_method
from ecmean.libs.files import var_is_there, get_inifiles, load_yaml, \
    make_input_filename, get_clim_files
from ecmean.libs.formula import formula_wrapper
from ecmean.libs.masks import mask_field, select_region
from ecmean.libs.areas import guess_bounds
from ecmean.libs.units import units_extra_definition
from ecmean.libs.ncfixers import xr_preproc, adjust_clim_file
from ecmean.libs.plotting import heatmap_comparison_pi, prepare_clim_dictionaries_pi
from ecmean.libs.parser import parse_arguments
from ecmean.libs.loggy import setup_logger

dask.config.set(scheduler="synchronous")


class PerformanceIndices:
    """
    Class to compute the performance indices for a given experiment and years.

    Attributes:
        exp (str): Experiment name.
        year1 (int): Start year of the experiment.
        year2 (int): End year of the experiment.
        config (str): Path to the configuration file. Default is 'config.yml'.
        loglevel (str): Logging level. Default is 'WARNING'.
        numproc (int): Number of processes to use. Default is 1.
        climatology (str): Climatology to use. Default is 'EC23'.
        interface (str): Path to the interface file.
        model (str): Model name.
        ensemble (str): Ensemble identifier. Default is 'r1i1p1f1'.
        silent (bool): If True, suppress output. Default is None.
        xdataset (xarray.Dataset): Dataset to use.
        outputdir (str): Directory to store output files.
        loggy (logging.Logger): Logger instance.
        diag (Diagnostic): Diagnostic instance.
        face (dict): Interface dictionary.
        piclim (dict): Climatology dictionary.
        util_dictionary (Supporter): Utility dictionary for remapping and masks.
        varstat (dict): Dictionary to store variable statistics.
        funcname (str): Name of the class.
        start_time (float): Start time for performance measurement.
    Methods:
        toc(message):
            Update the timer and log the elapsed time.
        prepare():
            Prepare the necessary components for performance indices calculation.
        run():
            Run the performance indices calculation.
        store(yamlfile=None):
            Store the performance indices in a yaml file.
        plot(mapfile=None, figformat='pdf'):
            Generate the heatmap for performance indices.
        pi_worker(util, piclim, face, diag, field_3d, varstat, varlist):
            Main parallel diagnostic worker for performance indices.
    """

    def __init__(self, exp, year1, year2, config='config.yml',
                 loglevel='WARNING', numproc=1, climatology='EC23',
                 interface=None, model=None, ensemble='r1i1p1f1',
                 silent=None, xdataset=None, outputdir=None):
        """Initialize the PerformanceIndices class with the given parameters."""

        self.loglevel = loglevel
        self.loggy = setup_logger(level=self.loglevel)
        self.diag = Diagnostic(exp=exp, year1=year1, year2=year2, config=config,
                               funcname=self.__class__.__name__,
                               numproc=numproc, climatology=climatology,
                               interface=interface,
                               modelname=model, ensemble=ensemble,
                               outputdir=outputdir, xdataset=xdataset)
        self.silent = silent
        self.face = None
        self.piclim = None
        self.util_dictionary = None
        self.varstat = None
        self.start_time = time()

    def toc(self, message):
        """Update the timer and log the elapsed time."""
        elapsed_time = time() - self.start_time
        self.start_time = time()
        self.loggy.info('%s time: %.2f seconds', message, elapsed_time)

    def prepare(self):
        """Prepare the necessary components for performance indices calculation."""
        # set dask and multiprocessing fork
        plat, mprocmethod = set_multiprocessing_start_method()
        self.loggy.info('Running on %s and multiprocessing method set as "%s"', plat, mprocmethod)

        # initialize the diag class, load the interface and the reference file
        self.face = load_yaml(self.diag.interface)
        self.piclim = load_yaml(self.diag.climfile)

        # check that everything is there
        check_var_climatology(self.diag.field_all, self.piclim.keys())

        # Create missing folders
        os.makedirs(self.diag.tabdir, exist_ok=True)
        os.makedirs(self.diag.figdir, exist_ok=True)

        # new bunch of functions to set grids, create correction command, masks and areas
        comp = self.face['model']['component']  # Get component for each domain

        # all clim have the same grid, read from the first clim available and get target grid
        clim, _ = get_clim_files(self.piclim, 'tas', self.diag, 'ALL')
        target_remap_grid = xr.open_dataset(clim)

        # get file info files
        inifiles = get_inifiles(self.face, self.diag)

        # add missing unit definitions
        units_extra_definition()

        # create remap dictionary with atm and oce interpolators
        self.util_dictionary = Supporter(comp, inifiles['atm'], inifiles['oce'],
                                         areas=False, remap=True, targetgrid=target_remap_grid)
        self.toc('Preparation')

    def run(self):
        """Run the performance indices calculation."""
        # main loop: manager is required for shared variables
        mgr = Manager()

        # dictionaries are shared, so they have to be passed as functions
        self.varstat = mgr.dict()
        processes = []

        # loop on the variables, create the parallel process
        for varlist in weight_split(self.diag.field_all, self.diag.numproc):
            core = Process(target=self.pi_worker, args=(self.util_dictionary, self.piclim,
                                                        self.face, self.diag, self.diag.field_3d,
                                                        self.varstat, varlist))
            core.start()
            processes.append(core)

        # wait for the processes to finish
        for proc in processes:
            proc.join()
        self.toc('Computation')

    def store(self, yamlfile=None):
        """Store the performance indices in a yaml file."""

        # order according to the original request the fields in the yaml file
        self.varstat = {var: self.varstat[var] for var in self.diag.field_all}

        # dump the yaml file for PI, including all the seasons (need to copy to avoid mess)
        if yamlfile is None:
            yamlfile = self.diag.filenames('yml')
        with open(yamlfile, 'w', encoding='utf-8') as file:
            yaml.safe_dump(self.varstat, file, default_flow_style=False, sort_keys=False)
        self.toc('Storing')

    def plot(self, mapfile=None, figformat='pdf'):
        """
        Generate the heatmap for performance indices.

        Args:
            mapfile (str): Path to the output file. If None, it will be defined automatically following ECmean syntax.
            figformat (str): Format of the output file. Default is 'pdf'.
        """
        # to this date, only EC23 support comparison with CMIP6 data
        if self.diag.climatology == 'EC23':

            # load yaml file if is missing
            if not self.varstat:
                yamlfile = self.diag.filenames('yml')
                self.loggy.info('Loading the stored data from the yaml file %s', yamlfile)
                if os.path.isfile(yamlfile):
                    with open(yamlfile, 'r', encoding='utf-8') as file:
                        self.varstat = yaml.safe_load(file)
                else:
                    raise FileNotFoundError(f'File {yamlfile} not found.')

            # prepare the data for the heatmap from the original yaml dictionaries
            data2plot, cmip6, longnames = prepare_clim_dictionaries_pi(data=self.varstat,
                                                                       clim=self.piclim,
                                                                       shortnames=self.diag.field_all)

            # call the heatmap routine for a plot
            if mapfile is None:
                mapfile = self.diag.filenames(figformat)

            heatmap_comparison_pi(data_dict=data2plot, cmip6_dict=cmip6, diag=self.diag, longnames=longnames, filemap=mapfile)
        else:
            self.loggy.warning('Only EC23 climatology is supported for comparison with CMIP6 data.')

        self.toc('Plotting')

    @staticmethod
    def pi_worker(util, piclim, face, diag, field_3d, varstat, varlist):
        """
        Main parallel diagnostic worker for performance indices.

        Args:
            util (Supporter): Utility dictionary for remapping and masks.
            piclim (dict): Climatology dictionary.
            face (dict): Interface dictionary.
            diag (Diagnostic): Diagnostic instance.
            field_3d (list): List of 3D fields.
            varstat (dict): Dictionary to store variable statistics.
            varlist (list): List of variables to process.
        """
        loggy = logging.getLogger(__name__)

        for var in varlist:
            # store NaN in dict (can't use defaultdict due to multiprocessing)
            result = init_mydict(diag.seasons, diag.regions)

            if check_var_interface(var, face):
                # get domain
                domain = get_domain(var, face)

                # get masks
                domain_mask = getattr(util, domain + 'mask')

                # check if required variables are there: use interface file
                # check into first file, and load also model variable units
                infile = make_input_filename(var, face, diag)

                # check if var is available
                isavail, varunit = var_is_there(infile, var, face)

                # if the variable is available
                if isavail:
                    # perform the unit conversion extracting offset and factor
                    offset, factor = UnitsHandler(var, org_units=varunit,
                                                  clim=piclim, face=face).offset_factor()

                    # open file: chunking on time only, might be improved
                    if not isinstance(infile, (xr.DataArray, xr.Dataset)):
                        xfield = xr.open_mfdataset(infile, preprocess=xr_preproc, chunks={'time': 12})
                    else:
                        xfield = infile

                    # in case of big files with multi year, be sure of having opened the right records
                    xfield = xfield.sel(time=xfield.time.dt.year.isin(diag.years_joined))

                    # check time axis
                    check_time_axis(xfield.time, diag.years_joined)

                    # get the data-array field for the required var
                    outfield = formula_wrapper(var, face, xfield)

                    # mean over time and fixing of the units
                    for season in diag.seasons:
                        loggy.debug(season)

                        # copy of the full field
                        tmean = outfield.copy(deep=True)

                        # get filenames for climatology
                        clim, vvvv = get_clim_files(piclim, var, diag, season)

                        # open climatology files, fix their metadata
                        cfield = adjust_clim_file(xr.open_mfdataset(clim, preprocess=xr_preproc))
                        vfield = adjust_clim_file(xr.open_mfdataset(vvvv, preprocess=xr_preproc), remove_zero=True)

                        # season selection
                        if season != 'ALL':
                            tmean = tmean.sel(time=tmean.time.dt.season.isin(season))
                            cfield = cfield.sel(time=cfield.time.dt.season.isin(season))
                            vfield = vfield.sel(time=vfield.time.dt.season.isin(season))

                        # averaging
                        tmean = tmean.mean(dim='time')
                        cfield = cfield.load()  # this ensure no crash with multiprocessing
                        vfield = vfield.load()

                        # safe check for old RK08 which has a different format
                        if diag.climatology != 'RK08':
                            cfield = cfield.mean(dim='time')
                            vfield = vfield.mean(dim='time')

                        tmean = tmean * factor + offset

                        # this computation is required to ensure that file access is done in a
                        # correct way. resource errors found if disabled in some specific
                        # configuration
                        tmean = tmean.compute()

                        # apply interpolation, if fixer is available and with different grids
                        fix = getattr(util, f'{domain}fix')
                        remap = getattr(util, f'{domain}remap')

                        if fix:
                            tmean = fix(tmean, keep_attrs=True)
                        try:
                            final = remap(tmean, keep_attrs=True)
                        except ValueError:
                            loggy.error('Cannot interpolate %s with the current weights...', var)
                            continue

                        # vertical interpolation
                        if var in field_3d:
                            # xarray interpolation on plev, forcing to be in Pascal
                            final = final.metpy.convert_coordinate_units('plev', 'Pa')
                            if set(final['plev'].data) != set(cfield['plev'].data):
                                loggy.warning('%s: Need to interpolate vertical levels...', var)
                                final = final.interp(plev=cfield['plev'].data, method='linear')

                                # safety check for missing values
                                sample = final.isel(lon=0, lat=0)
                                if np.sum(np.isnan(sample)) != 0:
                                    loggy.warning(
                                        '%s: You have NaN after the interpolation, this will affect your PIs...', var)
                                    levnan = cfield['plev'].where(np.isnan(sample))
                                    loggy.warning(levnan[~np.isnan(levnan)].data)

                            # zonal mean
                            final = final.mean(dim='lon')

                            # compute PI
                            complete = (final - cfield)**2 / vfield

                            # compute vertical bounds as weights
                            bounds_lev = guess_bounds(complete['plev'], name='plev')
                            bounds = abs(bounds_lev[:, 0] - bounds_lev[:, 1])
                            www = xr.DataArray(bounds, coords=[complete['plev']], dims=['plev'])

                            # vertical mean
                            outarray = complete.weighted(www).mean(dim='plev')

                        # horizontal averaging with land-sea mask
                        else:
                            complete = (final - cfield)**2 / vfield
                            outarray = mask_field(xfield=complete, mask_type=piclim[var]['mask'], dom=domain, mask=domain_mask)

                        # loop on different regions
                        for region in diag.regions:
                            slicearray = select_region(outarray, region)

                            # latitude-based averaging
                            weights = np.cos(np.deg2rad(slicearray.lat))
                            out = slicearray.weighted(weights).mean().data
                            # store the PI
                            result[season][region] = round(float(out), 3)

                            # diagnostic
                            if region == 'Global':
                                loggy.info('PI for %s %s %s %s', region, season, var, result[season][region])

            # nested dictionary, to be redefined as a dict to remove lambdas
            varstat[var] = result


def pi_entry_point():
    """
    Command line interface to run the performance indices calculation.
    """
    # read arguments from command line
    args = parse_arguments(sys.argv[1:], script='pi')

    performance_indices(exp=args.exp, year1=args.year1, year2=args.year2,
                        numproc=args.numproc, loglevel=args.loglevel,
                        climatology=args.climatology,
                        interface=args.interface, config=args.config,
                        model=args.model, ensemble=args.ensemble, outputdir=args.outputdir)


def performance_indices(exp, year1, year2, config='config.yml', loglevel='WARNING',
                        numproc=1, climatology='EC23', interface=None, model=None,
                        ensemble='r1i1p1f1', silent=None, xdataset=None, outputdir=None):
    """
    Wrapper function to compute the performance indices for a given experiment and years.
    """
    pi = PerformanceIndices(exp=exp, year1=year1, year2=year2, config=config,
                            loglevel=loglevel, numproc=numproc, climatology=climatology,
                            interface=interface, model=model, ensemble=ensemble, silent=silent,
                            xdataset=xdataset, outputdir=outputdir)
    pi.prepare()
    pi.run()
    pi.store()
    pi.plot()
