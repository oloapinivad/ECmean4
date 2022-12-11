#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""
    Global mean evaluation of CMIP6 climatology
    It runs and then reads the output of the cmip6 models to provide an assessment of the 
    climatology for CMIP6 models to be stored 
"""

__author__ = "Paolo Davini (p.davini@isac.cnr.it), Sep 2022."

from performance_indices import pi_main, parse_arguments
from ecmean import load_yaml
import os
import sys
import yaml
import warnings
import glob
import numpy as np

warnings.simplefilter("ignore") 

# the 30-year climatological window to be used as a baseline
year1 = 1981
year2 = 2010
expname = 'historical'
refclim = 'EC23'
nprocs = 4
do_compute = False
do_create_clim = True

# models on which we can build the clim
if do_create_clim : 
    models = ['EC-Earth3', 'IPSL-CM6A-LR', 'FGOALS-g3', 'TaiESM1','CanESM5', 'CNRM-CM6-1',
        'MIROC6', 'MPI-ESM1-2-HR', 'AWI-CM-1-1-MR', 'CMCC-CM2-SR5', 'NorESM2-MM']


else : 
    # models working completely
    models = ['EC-Earth3', 'IPSL-CM6A-LR', 'FGOALS-g3', 'TaiESM1','CanESM5', 'CNRM-CM6-1', 
        'MIROC6', 'MPI-ESM1-2-HR', 'AWI-CM-1-1-MR']

    # models with issue in the grid shape for siconc
    #models= ['CMCC-CM2-SR5', 'NorESM2-MM', 'ACCESS-CM2']

    # models which have not all the data
    #models=['UKESM1-0-LL']

    # model whhich does not work 
    #models=['GFDL-CM4', 'CESM2']

# call the loop of global mean on all the models
if do_compute : 
    for model in models : 
        print(model)

        #sys.argv = [expname, str(year1), str(year2), '--config', 'config_CMIP6_PD.yml', '--model', model, '-j', '8']
        #gm_main(sys.argv)
        if model in ['CNRM-CM6-1', 'UKESM1-0-LL'] :
            ensemble = "r1i1p1f2"
        else :
            ensemble = "r1i1p1f1"

        sys.argv = [expname, str(year1), str(year2), '--config', 'config_CMIP6_PD.yml', 
            '--model', model, '-j', str(nprocs), '-k', refclim, '-e', ensemble]
        pi_main(sys.argv)

if do_create_clim : 

    sys.argv = ['historical', str(year1), str(year2), '--config', 'config_CMIP6_PD.yml', '--model', models[0], '-k', refclim]
    args = parse_arguments(sys.argv)
    cfg = load_yaml(args.config)

    # dictionary with all elements
    full = {}
    for model in models:
        print(model)
        filein = glob.glob(os.path.join(cfg['dirs']['tab'], 'PI4_' + refclim + '_' + expname + 
            '_' + model + '_r1i1p1f*_' + str(year1) + '_' + str(year2) + '.yml'))
        full[model] = load_yaml(filein[0])

    # alternative to do sum on nested dictionary
    #full = {key:{key2:{key3:val1+data[key][key2][key3] for key3,val1 in subdic2.items()} for key2,subdic2 in subdic.items()} for key,subdic in full.items()}


    # idiot averaging
    m0 = models[0]
    out = {}
    for var in full[m0].keys() :
        out[var] = {}
        for season in full[m0][var].keys() :
            out[var][season] = {}
            for region in full[m0][var][season].keys() :
                element = []
                for model in full.keys() :
                    if var in full[model] :
                        element.append(full[model][var][season][region])
                out[var][season][region] =  np.nanmean(element)

    # clumsy way to get the models for each var
    mout = {}
    for var in full[m0].keys() :
        melement = []
        for model in full.keys() :
            if var in full[model] :
                melement.append(model)
        mout[var] = melement
    
    # clim files
    pifile = os.path.join(cfg['dirs']['clm'], refclim, 'pi_climatology_' + refclim + '.yml')
    #update_pifile = os.path.join(cfg['dirs']['clm'], refclim, 'pi_climatology_' + refclim + '_update.yml')
    update_pifile = os.path.join(cfg['dirs']['clm'], refclim, 'pi_climatology_' + refclim + '.yml')
    piclim = load_yaml(pifile)

    #update the climatology
    for var in out.keys() :  
        piclim[var]['cmip6'] = {}
        for season in out[var].keys() :
            piclim[var]['cmip6'][season] = {}
            for region in out[var][season].keys() : 
                piclim[var]['cmip6'][season][region] = float(out[var][season][region])
        piclim[var]['cmip6']['models'] = mout[var]
        piclim[var]['cmip6']['year1'] = year1
        piclim[var]['cmip6']['year2'] = year2
       
    # dump the new file
    with open(update_pifile, 'w') as file:
        yaml.safe_dump(piclim, file, sort_keys=False)