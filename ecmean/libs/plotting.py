#!/usr/bin/env python3
'''
Shared functions for XArray ECmean4
'''

##################
# PLOT FUNCTIONS #
##################

import textwrap
import logging
from matplotlib.colors import TwoSlopeNorm
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
import numpy as np
from ecmean.libs.general import dict_to_dataframe, init_mydict

loggy = logging.getLogger(__name__)


def heatmap_comparison_pi(data_dict, cmip6_dict, diag, longnames, filemap: str = None, size_model=14, **kwargs):
    """
    Function to produce a heatmap - seaborn based - for Performance Indices
    based on CMIP6 ratio

    Args:
        data_dict (dict): dictionary of absolute performance indices
        cmip6_dict (dict): dictionary of CMIP6 performance indices
        diag (object): Diagnostic object
        units_list (list): list of units
        filemap (str): path to save the plot
        size_model (int): size of the PIs in the plot

    Keyword Args:
        title (str): title of the plot, overrides default title
    """

    # convert output dictionary to pandas dataframe
    data_table = dict_to_dataframe(data_dict)
    loggy.debug("Data table")
    loggy.debug(data_table)

    # relative pi with re-ordering of rows
    cmip6_table = dict_to_dataframe(cmip6_dict).reindex(longnames)
    relative_table = data_table.div(cmip6_table)

    # compute the total PI mean
    relative_table.loc['Total PI'] = relative_table.mean()

    # reordering columns if info is available
    lll = [(x, y) for x in diag.seasons for y in diag.regions]
    relative_table = relative_table[lll]
    loggy.debug("Relative table")
    loggy.debug(relative_table)

    # defining plot size
    myfield = relative_table
    xfig = len(myfield.columns)
    yfig = len(myfield.index)

    # real plot
    _, axs = plt.subplots(1, 1, sharey=True, tight_layout=True, figsize=(xfig + 5, yfig + 2))

    thr = [0, 1, 5]
    tictoc = [0, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5]

    if 'title' in kwargs:
        title = kwargs['title']
    else:
        title = 'CMIP6 RELATIVE PI'
        title += f" {diag.modelname} {diag.expname} {diag.year1} {diag.year2}"

    tot = len(myfield.columns)
    sss = len(set([tup[1] for tup in myfield.columns]))
    divnorm = TwoSlopeNorm(vmin=thr[0], vcenter=thr[1], vmax=thr[2])
    pal = sns.color_palette("Spectral_r", as_cmap=True)
    chart = sns.heatmap(myfield, norm=divnorm, cmap=pal,
                        cbar_kws={"ticks": tictoc, 'label': title},
                        ax=axs, annot=True, linewidth=0.5, fmt='.2f',
                        annot_kws={'fontsize': size_model, 'fontweight': 'bold'})

    chart = chart.set_facecolor('whitesmoke')
    axs.set_title(title, fontsize=25)
    axs.vlines(list(range(sss, tot + sss, sss)), ymin=-1, ymax=len(myfield.index), colors='k')
    axs.hlines(len(myfield.index) - 1, xmin=-1, xmax=len(myfield.columns), colors='purple', lw=2)
    names = [' '.join(x) for x in myfield.columns]
    axs.set_xticks([x + .5 for x in range(len(names))], names, rotation=45, ha='right', fontsize=16)
    axs.set_yticks([x + .5 for x in range(len(myfield.index))], myfield.index, rotation=0, fontsize=16)
    axs.figure.axes[-1].tick_params(labelsize=15)
    axs.figure.axes[-1].yaxis.label.set_size(15)
    axs.set(xlabel=None)

    if filemap is None:
        filemap = 'PI4_heatmap.pdf'

    # save and close
    plt.savefig(filemap)
    plt.cla()
    plt.close()


def heatmap_comparison_gm(data_dict, mean_dict, std_dict, diag, units_list, filemap=None,
                          addnan=True, size_model=14, size_obs=8, **kwargs):
    """
    Function to produce a heatmap - seaborn based - for Global Mean
    based on season-averaged standard deviation ratio

    Args:
        data_dict (dict): table of model data
        mean_dict (dict): table of observations
        std_dict (dict): table of standard deviation
        diag (dict): diagnostic object
        units_list (list): list of units
        filemap (str): path to save the plot
        addnan (bool): add to the final plots also fields which cannot be compared against observations
        size_model (int): size of the model values in the plot
        size_obs (int): size of the observation values in the plot

    Keyword Args:
        title (str): title of the plot, overrides default title
    """

    # convert the three dictionary to pandas and then add units
    data_table = dict_to_dataframe(data_dict)
    mean_table = dict_to_dataframe(mean_dict)
    std_table = dict_to_dataframe(std_dict)
    for table in [data_table, mean_table, std_table]:
        table.index = table.index + ' [' + units_list + ']'

    loggy.debug("Data table")
    loggy.debug(data_table)

    # define array
    ratio = (data_table - mean_table) / std_table
    if addnan:
        mask = data_table[('ALL', 'Global')].notna()
    else:
        mask = ratio[('ALL', 'Global')].notna()
    clean = ratio[mask]

    # for dimension of plots
    xfig = len(clean.columns)
    yfig = len(clean.index)
    _, axs = plt.subplots(1, 1, sharey=True, tight_layout=True, figsize=(xfig + 5, yfig + 2))

    if 'title' in kwargs:
        title = kwargs['title']
    else:
        title = 'GLOBAL MEAN'
        title += f" {diag.modelname} {diag.expname} {diag.year1} {diag.year2}"

    # set color range and palette
    thr = 10
    tictoc = range(-thr, thr + 1)
    pal = ListedColormap(sns.color_palette("vlag", n_colors=21))
    tot = len(clean.columns)
    sss = len(set([tup[1] for tup in clean.columns]))

    chart = sns.heatmap(clean, annot=data_table[mask], vmin=-thr - 0.5, vmax=thr + 0.5, center=0,
                        annot_kws={'va': 'bottom', 'fontsize': size_model},
                        cbar_kws={'ticks': tictoc, "shrink": .5,
                                  'label': 'Model Bias \n (standard deviation of interannual variability from observations)'},
                        fmt='.2f', cmap=pal)
    if addnan:
        empty = np.where(clean.isna(), 0, np.nan)
        empty = np.where(data_table[mask] == 0, np.nan, empty)
        chart = sns.heatmap(empty, annot=data_table[mask], fmt='.2f',
                            vmin=-thr - 0.5, vmax=thr + 0.5, center=0,
                            annot_kws={'va': 'bottom', 'fontsize': size_model, 'color': 'dimgrey'}, cbar=False,
                            cmap=ListedColormap(['none']))
    chart = sns.heatmap(clean, annot=mean_table[mask], vmin=-thr - 0.5, vmax=thr + 0.5, center=0,
                        annot_kws={'va': 'top', 'fontsize': size_obs, 'fontstyle': 'italic'},
                        fmt='.2f', cmap=pal, cbar=False)
    if addnan:
        empty = np.where(clean.isna(), 0, np.nan)
        empty = np.where(mean_table[mask].isna(), np.nan, empty)
        chart = sns.heatmap(empty, annot=mean_table[mask], vmin=-thr - 0.5, vmax=thr + 0.5, center=0,
                            annot_kws={'va': 'top', 'fontsize': size_obs, 'fontstyle': 'italic', 'color': 'dimgrey'},
                            fmt='.2f', cmap=ListedColormap(['none']), cbar=False)

    chart = chart.set_facecolor('whitesmoke')
    axs.set_title(title, fontsize=25)
    axs.vlines(list(range(sss, tot + sss, sss)), ymin=-1, ymax=len(clean.index), colors='k')
    names = [' '.join(x) for x in clean.columns]
    axs.set_xticks([x + .5 for x in range(len(names))], names, rotation=45, ha='right', fontsize=16)
    axs.set_yticks([x + .5 for x in range(len(clean.index))], clean.index, rotation=0, fontsize=16)
    axs.set_yticklabels(textwrap.fill(y.get_text(), 28) for y in axs.get_yticklabels())
    axs.figure.axes[-1].tick_params(labelsize=15)
    axs.figure.axes[-1].yaxis.label.set_size(15)
    axs.set(xlabel=None)

    if filemap is None:
        filemap = 'Global_Mean_Heatmap.pdf'

    # save and close
    plt.savefig(filemap)
    plt.cla()
    plt.close()


def prepare_clim_dictionaries_pi(data, clim, shortnames):
    """
    Prepare dictionaries for plotting
    Args:
        data: dictionary with data
        clim: dictionary with climatology
        shortnames: list of shortnames
    Returns:
        data2plot: dictionary with data for plotting
        cmip6: dictionary with CMIP6 data
        longnames: list of longnames
    """

    # uniform dictionaries
    filt_piclim = {}
    for k in clim.keys():
        filt_piclim[k] = clim[k]['cmip6']
        for f in ['models', 'year1', 'year2']:
            del filt_piclim[k][f]

    # set longname, reorganize the dictionaries
    data2plot = {clim[var]['longname']: data[var] for var in shortnames}
    cmip6 = {clim[var]['longname']: filt_piclim[var] for var in shortnames}
    longnames = [clim[var]['longname'] for var in shortnames]

    return data2plot, cmip6, longnames


def prepare_clim_dictionaries_gm(data, clim, shortnames, seasons, regions):
    """
    Prepare dictionaries for global mean plotting

    Args:
        data: dictionary with the data
        clim: dictionary with the climatology
        shortnames: list of shortnames
        seasons: list of seasons
        regions: list of regions

    Returns:
        obsmean: dictionary with the mean
        obsstd: dictionary with the standard deviation
        data2plot: dictionary with the data to plot
        units_list: list of units
    """

    # loop on the variables to create obsmean and obsstd
    obsmean = {}
    obsstd = {}
    for var in shortnames:
        gamma = clim[var]
        obs = gamma['obs']

        # extract from yaml table for obs mean and standard deviation
        mmm = init_mydict(seasons, regions)
        sss = init_mydict(seasons, regions)
        # if we have all the obs/std available
        if isinstance(gamma['obs'], dict):
            for season in seasons:
                for region in regions:
                    mmm[season][region] = obs[season][region]['mean']
                    sss[season][region] = obs[season][region]['std']
        # if only global observation is available
        else:
            mmm['ALL']['Global'] = gamma['obs']

        # Assign to obsmean and obsstd using longname as the key
        obsmean[gamma['longname']] = mmm
        obsstd[gamma['longname']] = sss

    # set longname, get units
    data2plot = {clim[var]['longname']: data[var] for var in shortnames}
    units_list = [clim[var]['units'] for var in shortnames]

    return obsmean, obsstd, data2plot, units_list
