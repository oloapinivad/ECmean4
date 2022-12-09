#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""
    global mean evaluation of CMIP6 climatology 
"""

__author__ = "Paolo Davini (p.davini@isac.cnr.it), Sep 2022."

from global_mean import gm_main
from performance_indices import pi_main, parse_arguments
from ecmean import load_yaml
import os
import sys
import pandas as pd
import yaml
import warnings
import glob
from collections import Counter

warnings.simplefilter("ignore") 

# the 30-year climatological window to be used as a baseline
year1 = 1981
year2 = 2010
expname = 'historical'
refclim = 'EC23'
nprocs = 2
do_compute = True
do_create_clim = False

# models on which we can build the clim
if do_create_clim : 
    models = ['EC-Earth3', 'IPSL-CM6A-LR', 'FGOALS-g3', 'TaiESM1','CanESM5', 'CNRM-CM6-1', 'CESM2',
        'MIROC6', 'MPI-ESM1-2-HR', 'AWI-CM-1-1-MR', 'CMCC-CM2-SR5', 'NorESM2-MM']


else : 
    # models working completely
    models = ['EC-Earth3', 'IPSL-CM6A-LR', 'FGOALS-g3', 'TaiESM1','CanESM5', 'CNRM-CM6-1', 'CESM2',
        'MIROC6', 'MPI-ESM1-2-HR', 'AWI-CM-1-1-MR']

    models = ['FGOALS-g3']
    # models with issue in the grid shape for siconc
    #models= ['CMCC-CM2-SR5', 'NorESM2-MM', 'ACCESS-CM2']

    # models which have not all the data
    #models=['UKESM1-0-LL']

    # model whhich does not work 
    #models=['GFDL-CM4']

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

    full = {}
    for model in models:
        filein = glob.glob(os.path.join(cfg['dirs']['tab'], 'PI4_' + refclim + '_' + expname + 
            '_' + model + '_r1i1p1f*_' + str(year1) + '_' + str(year2) + '.yml'))

        # clumsly dictionary merging
        if model in models[0] :
            full = load_yaml(filein[0])
        else : 
            data = load_yaml(filein[0])
            # crazy sum on each element of the dictionary
            full = {key:{key2:{key3:val1+data[key][key2][key3] for key3,val1 in subdic2.items()} for key2,subdic2 in subdic.items()} for key,subdic in full.items()}

    # idiot averaging
    for var in full.keys() :
        for season in full[var].keys() :
            for region in full[var][season].keys() :
                full[var][season][region] = round(full[var][season][region]/len(models),2)

    print(full)

    pifile = os.path.join(cfg['dirs']['clm'], refclim, 'pi_climatology_' + refclim + '.yml')
    #update_pifile = os.path.join(cfg['dirs']['clm'], refclim, 'pi_climatology_' + refclim + '_update.yml')
    update_pifile = os.path.join(cfg['dirs']['clm'], refclim, 'pi_climatology_' + refclim + '.yml')
    piclim = load_yaml(pifile)

    for var in full.keys() :  
        piclim[var]['cmip6'] = {}
        for season in full[var].keys() :
            piclim[var]['cmip6'][season] = {}
            for region in full[var][season].keys() : 
                piclim[var]['cmip6'][season][region] = float(full[var][season][region])
        piclim[var]['cmip6']['models'] = models
        piclim[var]['cmip6']['year1'] = year1
        piclim[var]['cmip6']['year2'] = year2
       
    # dump the new file
    with open(update_pifile, 'w') as file:
        yaml.safe_dump(piclim, file, sort_keys=False)
