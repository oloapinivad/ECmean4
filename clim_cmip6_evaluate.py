#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""
    global mean evaluation of CMIP6 climatology 
"""

__author__ = "Paolo Davini (p.davini@isac.cnr.it), Sep 2022."

#from global_mean import gm_main
from performance_indices import pi_main, parse_arguments
from ecmean import load_yaml
import os
import sys
import pandas as pd
import yaml

# the 30-year climatological window to be used as a baseline
year1 = 1981
year2 = 2010
expname = 'historical'
do_compute = True

# an initial list of ten working models
#'AWI-CM-1-1-MR': crash in interpolation
#'CESM2': conflicting values for variable 'lat_bnds' on objects to be combined. You can skip this check by specifying compat='override'.
#'GFDL-CM4': Resulting object does not have monotonic global indexes along dimension lon
#'ACCESS-CM2':ValueError: The horizontal shape of input data is (145, 192), different from that of the regridder (144, 192)!
models = ['EC-Earth3', 'FGOALS-g3', 'IPSL-CM6A-LR', 'TaiESM1','CanESM5',
    'CMCC-CM2-SR5', 'MIROC6', 'MPI-ESM1-2-HR', 'NorESM2-MM']

# call the loop of global mean on all the models
if do_compute : 
    for model in models : 
        print(model)
        sys.argv = [expname, str(year1), str(year2), '--config', 'config_CMIP6_PD.yml', '--model', model]
        #gm_main(sys.argv)
        pi_main(sys.argv)

final = pd.DataFrame()
for model in models: 
    #print(model)
    sys.argv = ['historical', str(year1), str(year2), '--config', 'config_CMIP6_PD.yml', '--model', model]
    args = parse_arguments(sys.argv)
    cfg = load_yaml(args.config)
    filein = os.path.join(cfg['dirs']['tab'], 'PI4_RK08_' + expname +'_' + model + '_r1i1p1f1_1981_2010.txt')
    data = pd.read_table(filein, sep='|', skiprows=[1], skipfooter=3, 
        engine='python')
    clean = data.iloc[:, 1:3 ]
    clean.columns = clean.columns.str.replace(' ', '')
    clean['Var'] = clean['Var'].str.replace(" ", "")
    clean['model'] = [model]*clean.shape[0]
    final = pd.concat([final, pd.DataFrame(clean)])

varlist = final['Var'].unique()

dict = {}
for var in varlist : 
    select = final.loc[final['Var'] == var]
    floater = select['PI'].astype(float)
    dict[var] = floater.mean()

pifile = os.path.join(cfg['dirs']['clm'], 'RK08', 'pi_climatology_RK08.yml')
update_pifile = os.path.join(cfg['dirs']['clm'], 'RK08', 'pi_climatology_RK08_update.yml')
piclim = load_yaml(pifile)


for var in varlist : 
    piclim[var]['cmip6'] = float(dict[var])

print(piclim)

with open(update_pifile, 'w') as file:
    yaml.dump(piclim, file, sort_keys=False)
