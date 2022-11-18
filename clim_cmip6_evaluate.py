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
warnings.simplefilter("ignore") 

# the 30-year climatological window to be used as a baseline
year1 = 1981
year2 = 2010
expname = 'historical'
refclim = 'EC22'
do_compute = False
do_create_clim = True

# models on which we can build the clim
if do_create_clim : 
    models = ['EC-Earth3', 'IPSL-CM6A-LR', 'FGOALS-g3', 'TaiESM1','CanESM5', 'CNRM-CM6-1', 'CESM2',
        'MIROC6', 'MPI-ESM1-2-HR', 'AWI-CM-1-1-MR', 'CMCC-CM2-SR5', 'NorESM2-MM']

else : 
    # models working completely
    models = ['EC-Earth3', 'IPSL-CM6A-LR', 'FGOALS-g3', 'TaiESM1','CanESM5', 'CNRM-CM6-1', 'CESM2',
        'MIROC6', 'MPI-ESM1-2-HR', 'AWI-CM-1-1-MR']

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
            '--model', model, '-j', '4', '-k', refclim, '-e', ensemble]
        pi_main(sys.argv)

if do_create_clim : 
    final = pd.DataFrame()
    for model in models: 

        sys.argv = ['historical', str(year1), str(year2), '--config', 'config_CMIP6_PD.yml', '--model', model, '-k', refclim]
        args = parse_arguments(sys.argv)
        cfg = load_yaml(args.config)
        regions = cfg['PI']['regions']
        fileins = glob.glob(os.path.join(cfg['dirs']['tab'], 'PI4_' + refclim + '_' + expname + 
            '_' + model + '_r1i1p1f*_' + str(year1) + '_' + str(year2) + '.txt'))

        # this is clumsy since it needs to convert the formatted table to pandas
        for filein in fileins : 
            data = pd.read_table(filein, sep='|', skiprows=[1], skipfooter=3, 
                engine='python')
            data.columns = data.columns.str.strip()
            data['Variable'] = data['Variable'].str.strip()
            clean = data[regions + ['Variable']]
            clean['Model'] = [model]*clean.shape[0]
            final = pd.concat([final, pd.DataFrame(clean)])

    varlist = final['Variable'].unique()

    # create a dictionary 
    dict = {}
    for var in varlist :
        dict[var] = {} 
        select = final.loc[final['Variable'] == var]
        for region in regions : 
            select[region] = select[region].astype(float) 
            dict[var][region] = round(select[region].mean(), 2)
        cmip6_models = select[select['Global'].notnull()]['Model']
        dict[var]['models'] = list(cmip6_models)

    pifile = os.path.join(cfg['dirs']['clm'], refclim, 'pi_climatology_' + refclim + '.yml')
    #update_pifile = os.path.join(cfg['dirs']['clm'], refclim, 'pi_climatology_' + refclim + '_update.yml')
    update_pifile = os.path.join(cfg['dirs']['clm'], refclim, 'pi_climatology_' + refclim + '.yml')
    #updated_pifile = pifile
    piclim = load_yaml(pifile)

    #print(dict)

    for var in varlist :  
        piclim[var]['cmip6'] = {}
        for region in regions : 
            piclim[var]['cmip6'][region] = float(dict[var][region])
        piclim[var]['cmip6']['models'] = dict[var]['models']
        piclim[var]['cmip6']['year1'] = year1
        piclim[var]['cmip6']['year2'] = year2
       

    print(piclim)

    with open(update_pifile, 'w') as file:
        yaml.safe_dump(piclim, file, sort_keys=False)
