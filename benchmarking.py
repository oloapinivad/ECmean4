#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""
    Provide benchmarking running on a set of cores and set of years EC-Earth3 (or any other model) 
    data on a specific machine for both performance indices and global mean
    produce a figure to be automatically placed within the documentation and used as benchmark
    in different phases of the ECmean4 project evolution.
"""

from performance_indices import pi_main
from global_mean import gm_main
import sys
import timeit
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt



expname = 'historical'
refclim = 'EC23'
#models = ['EC-Earth3', 'IPSL-CM6A-LR', 'FGOALS-g3', 'TaiESM1', 'CanESM5', 'CNRM-CM6-1',
#              'MIROC6', 'MPI-ESM1-2-HR', 'AWI-CM-1-1-MR']
models = ['EC-Earth3']
ensemble = "r1i1p1f1"
scripts = ['Global Mean', 'Performance Indices']

# number of processors and years on which looop
nprocs= [1, 2, 4, 6, 8]
nyears = [1, 2, 5, 10, 30]

# file for configuration to be used as a reference
benchconfig = 'docs/config_benchmark.yml'

# number of repetition of each run
nrepeat = 3

# number of years and processors used for the respective figure
nyears_fixed = 10
nprocs_fixed = 4

# flag to define if the saved figure will go into the docs:
# to be done once in a while, or before big merge
do_definitive = False

# loop on models: so far not used
for model in models:

    # define the dataframe
    howmuch = pd.DataFrame(columns = ['script', 'time', 'nprocs', 'nyears'])    

    # loop on the scripts and on the procs
    for script in scripts : 
        for nproc in nprocs :
            if nproc == nprocs_fixed : 
                nn = nyears
            else :
                nn = [nyears_fixed]

            for n in nn : 
                print('Running ' + script + ' with nprocs =', str(nproc) + 'for nyears = ' + str(n))
                print(model)

                year1 = 1980
                year2 = year1 + n - 1
                # separated calls for model (perhaps a specific config file should be added)
                if script == 'Performance Indices' :
                    sys.argv = [expname, str(year1), str(year2), '--config', benchconfig,
                                '--model', model, '-j', str(nproc), '-k', refclim, '-e', ensemble, '-s', '-v', 'warning']
                    single = timeit.timeit(lambda: pi_main(sys.argv), number=nrepeat)
                elif script == 'Global Mean' : 
                    sys.argv = [expname, str(year1), str(year2), '--config', benchconfig,
                                '--model', model, '-j', str(nproc), '-e', ensemble, '-s',  '-v', 'warning']
                    single = timeit.timeit(lambda: gm_main(sys.argv), number=nrepeat)

                # concatenate
                d = pd.DataFrame({'script': [script], 'time': [round(single/nrepeat, 1)], 'nprocs': [nproc], 'nyears': [n]})
                print(d)
                howmuch = pd.concat([howmuch, d])
 
print(howmuch)

# produce the two-panel figure
fig, axs = plt.subplots(2, 1, sharey=True, tight_layout=True, figsize=(10, 10))
palette = ['teal', 'gold']

# first plot: scaling on cores
how1 = howmuch[howmuch["nyears"] == nyears_fixed]
chart1 = sns.barplot(data=how1, x='nprocs', y='time', hue='script', palette = palette,  ax=axs[0])

axs[0].set_title(f' ECmean4 execution time for CMIP6 {model} ({nyears_fixed} years)', fontsize=15)
for i in chart1.containers:
    chart1.bar_label(i,)
axs[0].set_xlabel('Number of cores', fontsize=15)
axs[0].set_ylabel('Execution time (sec)', fontsize=15)
h,l = axs[0].get_legend_handles_labels()
axs[0].legend(h,l)

# second plot: scaling on years
how2 = howmuch[howmuch["nprocs"] == nprocs_fixed]
chart2 = sns.barplot(data=how2, x='nyears', y='time', hue='script', palette = palette, ax=axs[1])

axs[1].set_title(f' ECmean4 execution time for CMIP6 {model} ({nprocs_fixed} nprocs)', fontsize=15)
for i in chart2.containers:
    chart2.bar_label(i,)
axs[1].set_xlabel('Number of years', fontsize=15)
axs[1].set_ylabel('Execution time (sec)', fontsize=15)
h,l = axs[1].get_legend_handles_labels()
axs[1].legend(h,l)

if do_definitive :
    figurepath = '/home/paolo/ECmean4/docs/sphinx/source/_static/benchmark.png'
else :
    figurepath = '/home/paolo/work/figures/ECmean4/benchmark.png'

fig.savefig(figurepath)
