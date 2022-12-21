#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""
    Provide benchmarking running on a set of cores EC-Earth3 (or any other model) 
    data on a specific machine for both performance indices and global mean
    produce a figure to be automatically placed within the documentation and used as benchmark
    in different phases of the project evolution.
"""

from performance_indices import pi_main
from global_mean import gm_main
import sys
import timeit
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# the 10-year time window to automatically produce a plot
year1 = 1990
year2 = 1999
expname = 'historical'
refclim = 'EC23'
models = ['EC-Earth3']
ensemble = "r1i1p1f1"
scripts = ['Global Mean', 'Performance Indices']
nprocs= [4, 8]
benchconfig = 'docs/config_benchmark.yml'

# number of repetition of each run
nrepeats = 3

# flag to define if the saved figure will go into the docs:
# to be done once in a while, or before big merge
do_definitive = False

# loop on models: so far not used
for model in models:

    # define the dataframe
    howmuch = pd.DataFrame(columns = ['script', 'time', 'nprocs'])    

    # loop on the scripts and on the procs
    for script in scripts : 
        for nproc in nprocs : 
            print(model)
            for nrepeat in range(nrepeats) : 
                # separated calls for model (perhaps a specific config file should be added)
                if script == 'Performance Indices' :
                    sys.argv = [expname, str(year1), str(year2), '--config', benchconfig,
                                '--model', model, '-j', str(nproc), '-k', refclim, '-e', ensemble, '-s']
                    single = timeit.timeit(lambda: pi_main(sys.argv), number=1)
                elif script == 'Global Mean' : 
                    sys.argv = [expname, str(year1), str(year2), '--config', benchconfig,
                                '--model', model, '-j', str(nproc), '-e', ensemble, '-s']
                    single = timeit.timeit(lambda: gm_main(sys.argv), number=1)

                # concatenate
                d = pd.DataFrame({'script': [script], 'time': [round(single, 1)], 'nprocs': [nproc]})
                print(d)
                howmuch = pd.concat([howmuch, d])
 
print(howmuch)

# produce the new figure
fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True, figsize=(10, 5))
palette = ['teal', 'gold']
chart = sns.barplot(data=howmuch, x='nprocs', y='time', hue='script', palette = palette)
nyears = year2 - year1 + 1
axs.set_title(f' ECmean4 execution time for CMIP6 {model} ({nyears} years)', fontsize=15)
for i in chart.containers:
    chart.bar_label(i,)
axs.set_xlabel('Number of cores', fontsize=15)
axs.set_ylabel('Execution time (sec)', fontsize=15)
h,l = axs.get_legend_handles_labels()
axs.legend(h,l)

if do_definitive :
    figurepath = '/home/paolo/ECmean4/docs/sphinx/source/_static/benchmark.png'
else :
    figurepath = '/home/paolo/work/figures/ECmean4/prova.png'

fig.savefig(figurepath)

