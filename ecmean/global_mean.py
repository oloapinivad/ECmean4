#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    python3 version of ECmean global mean tool.
    Using a reference file from yaml and Xarray

    @author Paolo Davini (p.davini@isac.cnr.it), Sep 2022.
    @author Jost von Hardenberg (jost.hardenberg@polito.it), Sep 2022
"""

__author__ = "Paolo Davini (p.davini@isac.cnr.it), Sep 2022."

import os
import sys
import logging
import argparse
from multiprocessing import Process, Manager
from time import time
from tabulate import tabulate
import numpy as np
import xarray as xr
import yaml
import dask

from ecmean import Diagnostic, Supporter, UnitsHandler
from ecmean.libs.general import weight_split, write_tuning_table, get_domain, \
     check_time_axis, init_mydict, \
          check_var_interface, check_var_climatology, set_multiprocessing_start_method
from ecmean.libs.files import var_is_there, get_inifiles, load_yaml, make_input_filename
from ecmean.libs.formula import formula_wrapper
from ecmean.libs.masks import masked_meansum, select_region
from ecmean.libs.units import units_extra_definition
from ecmean.libs.ncfixers import xr_preproc
from ecmean.libs.parser import parse_arguments
from ecmean.libs.plotting import heatmap_comparison_gm, prepare_clim_dictionaries_gm
from ecmean.libs.loggy import setup_logger

dask.config.set(scheduler="synchronous")


class GlobalMean:
    """
    Attributes:
        exp (str): Experiment name.
        year1 (int): Start year of the experiment.
        year2 (int): End year of the experiment.
        config (str): Path to the configuration file. Default is 'config.yml'.
        loglevel (str): Logging level. Default is 'WARNING'.
        numproc (int): Number of processes to use. Default is 1.
        interface (str): Path to the interface file. Default is None.
        model (str): Model name. Default is None.
        ensemble (str): Ensemble identifier. Default is 'r1i1p1f1'.
        addnan (bool): Whether to add NaNs. Default is False.
        silent (bool): Whether to suppress output. Default is None.
        trend (bool): Whether to compute trends. Default is None.
        line (str): Line identifier. Default is None.
        outputdir (str): Output directory. Default is None.
        xdataset (str): Path to the xdataset. Default is None.
        loggy (logging.Logger): Logger instance.
        diag (Diagnostic): Diagnostic instance.
        face (dict): Interface dictionary.
        ref (dict): Reference dictionary.
        util_dictionary (Supporter): Supporter instance.
        varmean (dict): Dictionary to store variable means.
        vartrend (dict): Dictionary to store variable trends.
        funcname (str): Name of the class.
        start_time (float): Start time for the timer.
    Methods:
        toc(message):
            Update the timer and log the elapsed time.
        prepare():
            Prepare the necessary components for the global mean computation.
        run():
            Run the global mean computation using multiprocessing.
        store():
            Store the computed global mean values in a table and YAML file.
        plot(mapfile=None, figformat='pdf'):
        gm_worker(util, ref, face, diag, varmean, vartrend, varlist):
    """
    def __init__(self, exp, year1, year2, config='config.yml', loglevel='WARNING', numproc=1,
                      interface=None, model=None, ensemble='r1i1p1f1', addnan=False, silent=None,
                      trend=None, line=None, outputdir=None, xdataset=None):
        self.exp = exp
        self.year1 = year1
        self.year2 = year2
        self.config = config
        self.loglevel = loglevel
        self.numproc = numproc
        self.interface = interface
        self.model = model
        self.ensemble = ensemble
        self.addnan = addnan
        self.silent = silent
        self.trend = trend
        self.line = line
        self.outputdir = outputdir
        self.xdataset = xdataset
        self.loggy = setup_logger(level=self.loglevel)
        self.diag = None
        self.face = None
        self.ref = None
        self.util_dictionary = None
        self.varmean = None
        self.vartrend = None
        self.funcname = self.__class__.__name__
        self.start_time = time()
    
    def toc(self, message):
        """Update the timer and log the elapsed time."""
        elapsed_time = time() - self.start_time
        self.start_time = time()
        self.loggy.info('%s time: %.2f seconds', message, elapsed_time)

    def prepare(self):
        """Prepare the necessary components for the global mean computation."""
        plat, mprocmethod = set_multiprocessing_start_method()
        self.loggy.info('Running on %s and multiprocessing method set as "%s"', plat, mprocmethod)

        self.diag = Diagnostic(argparse.Namespace(**self.__dict__))
        self.face = load_yaml(self.diag.interface)
        self.ref = load_yaml(self.diag.reffile)

        check_var_climatology(self.diag.var_all, self.ref.keys())

        os.makedirs(self.diag.tabdir, exist_ok=True)
        os.makedirs(self.diag.figdir, exist_ok=True)

        comp = self.face['model']['component']
        inifiles = get_inifiles(self.face, self.diag)

        units_extra_definition()

        self.util_dictionary = Supporter(comp, inifiles['atm'], inifiles['oce'], areas=True, remap=False)
        self.toc('Preparation')

    def run(self):
        """Run the global mean computaacross all variables on using multiprocessing."""
        mgr = Manager()
        self.varmean = mgr.dict()
        self.vartrend = mgr.dict()
        processes = []

        for varlist in weight_split(self.diag.var_all, self.diag.numproc):
            core = Process(target=self.gm_worker, args=(self.util_dictionary, self.ref, self.face, self.diag,
                                                                        self.varmean, self.vartrend, varlist))
            core.start()
            processes.append(core)

        for proc in processes:
            proc.join()
        self.toc('Computation')

    def yamlfile(self):
        """Define the output YAML filename"""

        yamlfile= self.diag.tabdir / f'global_mean_{self.diag.expname}_{self.diag.modelname}_{self.diag.ensemble}_{self.diag.year1}_{self.diag.year2}.yml'
        return yamlfile

    def store(self, yamlfile=None):
        """Rearrange the data and save the yaml file and the table."""
        global_table = []

        # reorder the data to be stored
        for var in self.diag.var_all:
            gamma = self.ref[var]
            if isinstance(gamma['obs'], dict):
                tabval = gamma['obs']['ALL']['Global']
                outval = str(tabval['mean']) + '\u00B1' + str(tabval['std'])
            else:
                outval = gamma['obs']

            if 'year1' in gamma.keys():
                years = str(gamma['year1']) + '-' + str(gamma['year2'])
            else:
                raise ValueError('Year1 and Year2 are not defined in the reference file')

            out_sequence = [var, gamma['longname'], gamma['units'], self.varmean[var]['ALL']['Global']]
            if self.diag.ftrend:
                out_sequence = out_sequence + [self.vartrend[var]['ALL']['Global']]
            out_sequence = out_sequence + [outval, gamma.get('dataset', ''), years]
            global_table.append(out_sequence)

        head = ['Variable', 'Longname', 'Units', self.diag.modelname]
        if self.diag.ftrend:
            head = head + ['Trend']
        head = head + ['Obs.', 'Dataset', 'Years']

        # save table
        tablefile = os.path.join(self.diag.tabdir,
                                        f'global_mean_{self.diag.expname}_{self.diag.modelname}_{self.diag.ensemble}_{self.diag.year1}_{self.diag.year2}.txt')
        self.loggy.info('Table file is: %s', tablefile)
        with open(tablefile, 'w', encoding='utf-8') as out:
            out.write(tabulate(global_table, headers=head, stralign='center', tablefmt='orgtbl'))

        # reaorder
        self.varmean = {var: self.varmean[var] for var in self.diag.var_all}

        # save yaml fil
        if yamlfile is None:
            yamlfile = self.yamlfile()
        
        self.loggy.info('YAML file is: %s', tablefile)
        with open(yamlfile, 'w', encoding='utf-8') as file:
            yaml.safe_dump(self.varmean, file, default_flow_style=False, sort_keys=False)

    def plot(self, mapfile=None, figformat='pdf'):
        """"
        Plot the global mean values.
        Args:
            mapfile: Path to the output file. If None, it will be defined automatically following ECmean syntax
            figformat: Format of the output file.
        """

        # load yaml file if is missing
        if not self.varmean:
            yamlfile = self.yamlfile()
            self.loggy.info('Loading the stored data from the yaml file %s', yamlfile)
            if os.path.isfile(yamlfile):
                with open(yamlfile, 'r', encoding='utf-8') as file:
                    self.varmean = yaml.safe_load(file)
            else:
                raise FileNotFoundError(f'YAML file {yamlfile} not found')

        # prepare the dictionaries for the plotting
        obsmean, obsstd, data2plot, units_list = prepare_clim_dictionaries_gm(self.varmean, self.ref,
                                                                              self.diag.var_all, self.diag.seasons,
                                                                              self.diag.regions)
        if mapfile is None:
            mapfile = os.path.join(self.diag.figdir,
                                    f'global_mean_{self.diag.expname}_{self.diag.modelname}_r1i1p1f1_{self.diag.year1}_{self.diag.year2}.{figformat}')
        self.loggy.info('Figure file is: %s', mapfile)

        # call the heatmap for plottinh
        heatmap_comparison_gm(data_dict=data2plot, mean_dict=obsmean, std_dict=obsstd,
                                    diag=self.diag, units_list=units_list,
                                    filemap=mapfile, addnan=self.diag.addnan)

        if self.diag.ftable:
            self.loggy.info('Line file is: %s', self.diag.linefile)
            write_tuning_table(self.diag.linefile, self.varmean, self.diag.var_table, self.diag, self.ref)
        self.toc('Plotting')

    @staticmethod
    def gm_worker(util, ref, face, diag, varmean, vartrend, varlist):
        """"
        Workhorse for the global mean computation.

        """
        loggy = logging.getLogger(__name__)

        for var in varlist:
            result = init_mydict(diag.seasons, diag.regions)
            trend = init_mydict(diag.seasons, diag.regions)

            if check_var_interface(var, face):
                domain = get_domain(var, face)
                weights = getattr(util, domain + 'area')
                domain_mask = getattr(util, domain + 'mask')
                infile = make_input_filename(var, face, diag)
                isavail, varunit = var_is_there(infile, var, face)

                if isavail:
                    units_handler = UnitsHandler(var, org_units=varunit, clim=ref, face=face)
                    offset, factor = units_handler.offset, units_handler.factor

                    if not isinstance(infile, (xr.DataArray, xr.Dataset)):
                        xfield = xr.open_mfdataset(infile, preprocess=xr_preproc, chunks={'time': 12})
                    else:
                        xfield = infile

                    xfield = xfield.sel(time=xfield.time.dt.year.isin(diag.years_joined))
                    check_time_axis(xfield.time, diag.years_joined)
                    cfield = formula_wrapper(var, face, xfield).compute()

                    for season in diag.seasons:
                        tfield = cfield.copy(deep=True)
                        if season != 'ALL':
                            tfield = tfield.sel(time=cfield.time.dt.season.isin(season))

                        if diag.ftrend:
                            tfield = tfield.groupby('time.year').mean('time')
                        else:
                            tfield = tfield.mean(dim='time')

                        for region in diag.regions:
                            slicefield = select_region(tfield, region)
                            sliceweights = select_region(weights, region)
                            if isinstance(domain_mask, xr.DataArray):
                                slicemask = select_region(domain_mask, region)
                            else:
                                slicemask = 0.

                            avg = masked_meansum(
                                xfield=slicefield, weights=sliceweights, mask=slicemask,
                                operation=ref[var].get('operation', 'mean'),
                                mask_type=ref[var].get('mask', 'global'),
                                domain=domain)

                            if isinstance(avg, dask.array.core.Array):
                                avg = avg.compute()

                            result[season][region] = float((np.nanmean(avg) + offset) * factor)

                            if diag.ftrend:
                                if len(avg) == len(diag.years_joined):
                                    trend[season][region] = np.polyfit(diag.years_joined, avg, 1)[0]
                            if season == 'ALL' and region == 'Global':
                                loggy.info('Average: %s %s %s %s', var, season, region, result[season][region])

            varmean[var] = result
            vartrend[var] = trend


def gm_entry_point():
    args = parse_arguments(sys.argv[1:], script='gm')
    global_mean(exp=args.exp, year1=args.year1, year2=args.year2, numproc=args.numproc,
                          trend=args.trend, line=args.line, loglevel=args.loglevel,
                          interface=args.interface, config=args.config, model=args.model,
                          ensemble=args.ensemble, addnan=args.addnan, outputdir=args.outputdir)
    print('ECmean4 Global Mean successfully computed!')

def global_mean(exp, year1, year2, config='config.yml', loglevel='WARNING', numproc=1,
                interface=None, model=None, ensemble='r1i1p1f1', addnan=False, silent=None,
                trend=None, line=None, outputdir=None, xdataset=None):
    gm = GlobalMean(exp, year1, year2, config, loglevel, numproc, interface, model, ensemble, addnan, silent, trend, line, outputdir, xdataset)
    gm.prepare()
    gm.run()
    gm.store()
    gm.plot()
    

