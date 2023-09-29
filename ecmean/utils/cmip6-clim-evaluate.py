#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""
    Global mean evaluation of CMIP6 climatology
    It runs and then reads the output of the cmip6 models to provide an assessment of the
    climatology for CMIP6 models to be stored
"""

__author__ = "Paolo Davini (p.davini@isac.cnr.it), Sep 2022."

import os
import yaml
import warnings
import glob
import numpy as np
from ecmean.performance_indices import performance_indices
from ecmean.libs.files import load_yaml
warnings.simplefilter("ignore")

# the 30-year climatological window to be used as a baseline
year1 = 1981
year2 = 2010
expname = 'historical'
refclim = 'EC23'
nprocs = 2
do_compute = True
do_create_clim = True
do_definitive = False
config_file = '../config_CMIP6_PD.yml'
climdir = '../climatology/'

models = ['EC-Earth3', 'IPSL-CM6A-LR', 'FGOALS-g3', 'TaiESM1', 'CanESM5', 'CESM2',
          'MIROC6', 'MPI-ESM1-2-HR', 'AWI-CM-1-1-MR', 'CMCC-CM2-SR5', 'NorESM2-MM', 'GFDL-CM4']

# models with issue in the grid shape
# models= ['ACCESS-CM2']

# models which have not all the data
# models=['UKESM1-0-LL', 'CNRM-CM6-1']


# call the loop of global mean on all the models
if do_compute:
    for model in models:
        print(model)

        if model in ['CNRM-CM6-1', 'UKESM1-0-LL']:
            ensemble = "r1i1p1f2"
        else:
            ensemble = "r1i1p1f1"

        performance_indices(expname, year1, year2, config=config_file, model=model,
                            ensemble=ensemble, numproc=nprocs, climatology=refclim, 
                            loglevel='WARNING')

if do_create_clim:

    cfg = load_yaml(config_file)

    # dictionary with all elements
    full = {}
    for model in models:
        print(model)
        filein = glob.glob(os.path.join(cfg['dirs']['tab'], 'PI4_' + refclim + '_' + expname +
                                        '_' + model + '_r1i1p1f*_' + str(year1) + '_' + str(year2) + '.yml'))
        full[model] = load_yaml(filein[0])


    # idiot averaging
    m0 = models[0]
    out = {}
    for var in full[m0].keys():
        out[var] = {}
        for season in full[m0][var].keys():
            out[var][season] = {}
            for region in full[m0][var][season].keys():
                element = []
                for model in full.keys():
                    if var in full[model]:
                        if not np.isnan(full[model][var][season][region]):
                            element.append(full[model][var][season][region])
                out[var][season][region] = float(round(np.mean(element), 3))

    # clumsy way to get the models for each var
    mout = {}
    for var in full[m0].keys():
        melement = []
        for model in full.keys():
            if var in full[model]:
                if not np.isnan(full[model][var]['ALL']['Global']):
                    melement.append(model)
        mout[var] = melement

    # clim files
    pifile = os.path.join(climdir, refclim, 'pi_climatology_' + refclim + '.yml')
    if not do_definitive:
        update_pifile = os.path.join(climdir, refclim, 'pi_climatology_' + refclim + '_test.yml')
    else:
        update_pifile = os.path.join(climdir, refclim, 'pi_climatology_' + refclim + '.yml')
    piclim = load_yaml(pifile)

    # Update the climatology
    for var in out.keys():
        piclim[var]['cmip6'] = {}
        for season, region_data in out[var].items():
            piclim[var]['cmip6'][season] = {}
            for region, value in region_data.items():
                piclim[var]['cmip6'][season][region] = float(value)
        piclim[var]['cmip6']['models'] = mout[var]
        piclim[var]['cmip6']['year1'] = year1
        piclim[var]['cmip6']['year2'] = year2


    # # update the climatology
    # for var in out.keys():
    #     piclim[var]['cmip6'] = {}
    #     for season in out[var].keys():
    #         piclim[var]['cmip6'][season] = {}
    #         for region in out[var][season].keys():
    #             piclim[var]['cmip6'][season][region] = float(out[var][season][region])
    #     piclim[var]['cmip6']['models'] = mout[var]
    #     piclim[var]['cmip6']['year1'] = year1
    #     piclim[var]['cmip6']['year2'] = year2

    # dump the new file
    with open(update_pifile, 'w', encoding='utf8') as file:
        yaml.safe_dump(piclim, file, sort_keys=False)
