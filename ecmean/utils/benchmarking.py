#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""
    Provide benchmarking running on a set of cores and set of years EC-Earth3 (or any other model)
    data on a specific machine for both performance indices and global mean
    produce a figure to be automatically placed within the documentation and used as benchmark
    in different phases of the ECmean4 project evolution.
    This is made to run on CNR-ISAC machina 'Wilma' only
"""

import shutil
import os
import timeit
from datetime import date
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from ecmean.performance_indices import performance_indices
from ecmean.global_mean import global_mean

expname = 'historical'
refclim = 'EC23'
# models = ['EC-Earth3', 'IPSL-CM6A-LR', 'FGOALS-g3', 'TaiESM1', 'CanESM5', 'CNRM-CM6-1',
#              'MIROC6', 'MPI-ESM1-2-HR', 'AWI-CM-1-1-MR']
models = ['EC-Earth3']
ensemble = "r1i1p1f1"
scripts = ['Global Mean', 'Performance Indices']

# number of processors and years on which looop
nprocs = [1, 2, 4, 6, 8]
nyears = [1, 2, 5, 10, 30]
# nprocs = [2, 4]
# nyears = [1, 10]

# file for configuration to be used as a reference
benchconfig = 'config_benchmark.yml'

# benchmark directory
benchdir = '/home/paolo/work/figures/ECmean4/benchmark'

# number of repetition of each run
nrepeat = 5

# number of years and processors used for the respective figure
nyears_fixed = 10
nprocs_fixed = 4

# flag to define if the saved figure will go into the docs:
# to be done once in a while, or before big merge
do_definitive = True

# loop on models: so far not used
for model in models:

    # define the dataframe
    howmuch = pd.DataFrame(columns=['script', 'time', 'nprocs', 'nyears'])

    # loop on the scripts and on the procs
    for script in scripts:
        for nproc in nprocs:
            if nproc == nprocs_fixed:
                nn = nyears
            else:
                nn = [nyears_fixed]

            for nyear in nn:
                print('Running ' + script + ' with nprocs =', str(nproc) + ' for nyears = ' + str(nyear))
                print(model)

                year1 = 1980
                year2 = year1 + nyear - 1
                # separated calls for model (perhaps a specific config file should be added)
                if script == 'Performance Indices':
                    single = timeit.timeit(lambda: performance_indices(expname, year1, year2, config=benchconfig,
                                                                       model=model, numproc=nproc,
                                                                       climatology=refclim, ensemble=ensemble,
                                                                       loglevel='info'), number=nrepeat)
                    # sys.argv = [expname, str(year1), str(year2), '--config', benchconfig,
                    #            '--model', model, '-j', str(nproc), '-k', refclim, '-e', ensemble, '-v', 'warning']
                    # single = timeit.timeit(lambda: pi_main(sys.argv), number=nrepeat)
                elif script == 'Global Mean':
                    single = timeit.timeit(lambda: global_mean(expname, year1, year2, config=benchconfig,
                                                               model=model, numproc=nproc,
                                                               ensemble=ensemble,
                                                               loglevel='info'), number=nrepeat)
                    # sys.argv = [expname, str(year1), str(year2), '--config', benchconfig,
                    #            '--model', model, '-j', str(nproc), '-e', ensemble, '-v', 'warning']
                    # single = timeit.timeit(lambda: gm_main(sys.argv), number=nrepeat)

                # concatenate
                d = pd.DataFrame({'script': [script], 'time': [round(single / nrepeat, 1)],
                                 'nprocs': [nproc], 'nyears': [nyear]})
                print(d)
                howmuch = pd.concat([howmuch, d])

# create folder
if not os.path.exists(benchdir):
    os.makedirs(benchdir)

# save to file
today = str(date.today())
csvname = os.path.join(benchdir, 'csv-benchmark-' + today + '.csv')
howmuch.to_csv(csvname, sep='\t')
print(howmuch)

# produce the two-panel figure
fig, axs = plt.subplots(2, 1, sharey=True, tight_layout=True, figsize=(10, 10))
palette = ['teal', 'gold']

# first plot: scaling on cores
how1 = howmuch[howmuch["nyears"] == nyears_fixed]
chart1 = sns.barplot(data=how1, x='nprocs', y='time', hue='script', palette=palette, ax=axs[0])

axs[0].set_title(f' ECmean4 execution time for CMIP6 {model} ({nyears_fixed} years)', fontsize=15)
for i in chart1.containers:
    chart1.bar_label(i,)
axs[0].set_xlabel('Number of cores', fontsize=15)
axs[0].set_ylabel('Execution time (sec)', fontsize=15)
axs[0].set(ylim=(0, 300))
hh, ll = axs[0].get_legend_handles_labels()
axs[0].legend(hh, ll)

# second plot: scaling on years
how2 = howmuch[howmuch["nprocs"] == nprocs_fixed]
chart2 = sns.barplot(data=how2, x='nyears', y='time', hue='script', palette=palette, ax=axs[1])

axs[1].set_title(f' ECmean4 execution time for CMIP6 {model} ({nprocs_fixed} nprocs)', fontsize=15)
for i in chart2.containers:
    chart2.bar_label(i,)
axs[1].set_xlabel('Number of years', fontsize=15)
axs[1].set(ylim=(0, 300))
axs[1].set_ylabel('Execution time (sec)', fontsize=15)
hh, ll = axs[1].get_legend_handles_labels()
axs[1].legend(hh, ll)

# save figure for benchmarks
figurepath = os.path.join(benchdir, 'benchmark-' + today + '.png')
fig.savefig(figurepath)

localdir = os.path.dirname(os.path.abspath(__file__))

if do_definitive:
    shutil.copy(figurepath, os.path.join(localdir, '../../docs/sphinx/source/_static/benchmark.png'))
