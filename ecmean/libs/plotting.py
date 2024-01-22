#!/usr/bin/env python3
'''
Shared functions for XArray ECmean4
'''

##################
# PLOT FUNCTIONS #
##################

import textwrap
from matplotlib.colors import TwoSlopeNorm
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
import numpy as np



def heatmap_comparison_pi(relative_table, diag, filemap, size_model = 14):
    """Function to produce a heatmap - seaborn based - for Performance Indices
    based on CMIP6 ratio"""

    # defining plot size
    myfield = relative_table
    xfig = len(myfield.columns)
    yfig = len(myfield.index)

    # real plot
    _, axs = plt.subplots(1, 1, sharey=True, tight_layout=True, figsize=(xfig+5, yfig+2))

    thr = [0, 1, 5]
    tictoc = [0, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5]
    title = 'CMIP6 RELATIVE PI'


    # axs.subplots_adjust(bottom=0.2)
    # pal = sns.diverging_palette(h_neg=130, h_pos=10, s=99, l=55, sep=3, as_cmap=True)
    tot = len(myfield.columns)
    sss = len(set([tup[1] for tup in myfield.columns]))
    divnorm = TwoSlopeNorm(vmin=thr[0], vcenter=thr[1], vmax=thr[2])
    pal = sns.color_palette("Spectral_r", as_cmap=True)
    # pal = sns.diverging_palette(220, 20, as_cmap=True)
    chart = sns.heatmap(myfield, norm=divnorm, cmap=pal,
                        cbar_kws={"ticks": tictoc, 'label': title},
                        ax=axs, annot=True, linewidth=0.5, fmt='.2f',
                        annot_kws={'fontsize': size_model, 'fontweight': 'bold'})
    
    chart = chart.set_facecolor('whitesmoke')
    axs.set_title(f'{title} {diag.modelname} {diag.expname} {diag.year1} {diag.year2}', fontsize=25)
    axs.vlines(list(range(sss, tot + sss, sss)), ymin=-1, ymax=len(myfield.index), colors='k')
    axs.hlines(len(myfield.index) - 1, xmin=-1, xmax=len(myfield.columns), colors='purple', lw=2)
    names = [' '.join(x) for x in myfield.columns]
    axs.set_xticks([x + .5 for x in range(len(names))], names, rotation=45, ha='right', fontsize=16)
    axs.set_yticks([x + .5 for x in range(len(myfield.index))], myfield.index, rotation=0, fontsize=16)
    axs.figure.axes[-1].tick_params(labelsize=15)
    axs.figure.axes[-1].yaxis.label.set_size(15)
    axs.set(xlabel=None)

    # save and close
    plt.savefig(filemap)
    plt.cla()
    plt.close()


def heatmap_comparison_gm(data_table, mean_table, std_table, diag, filemap, addnan = True,
                          size_model = 14, size_obs = 8):
    """Function to produce a heatmap - seaborn based - for Global Mean
    based on season-averaged standard deviation ratio"""

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
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True, figsize=(xfig+5, yfig+2))

    title = 'GLOBAL MEAN'

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
        empty = np.where(data_table[mask]==0, np.nan, empty)
        chart = sns.heatmap(empty, annot=data_table[mask], fmt='.2f',
                            vmin=-thr - 0.5, vmax=thr + 0.5, center=0,
                            annot_kws={'va': 'bottom', 'fontsize': size_model, 'color':'dimgrey'}, cbar=False,
                            cmap=ListedColormap(['none']))
    chart = sns.heatmap(clean, annot=mean_table[mask], vmin=-thr - 0.5, vmax=thr + 0.5, center=0,
                        annot_kws={'va': 'top', 'fontsize': size_obs, 'fontstyle': 'italic'},
                        fmt='.2f', cmap=pal, cbar=False)
    if addnan:
        empty = np.where(clean.isna(), 0, np.nan)
        empty = np.where(mean_table[mask].isna(), np.nan, empty)
        chart = sns.heatmap(empty, annot=mean_table[mask], vmin=-thr - 0.5, vmax=thr + 0.5, center=0,
                        annot_kws={'va': 'top', 'fontsize': size_obs, 'fontstyle': 'italic', 'color':'dimgrey'},
                        fmt='.2f', cmap=ListedColormap(['none']), cbar=False)

    chart = chart.set_facecolor('whitesmoke')
    axs.set_title(f'{title} {diag.modelname} {diag.expname} {diag.year1} {diag.year2}', fontsize=25)
    axs.vlines(list(range(sss, tot + sss, sss)), ymin=-1, ymax=len(clean.index), colors='k')
    names = [' '.join(x) for x in clean.columns]
    axs.set_xticks([x + .5 for x in range(len(names))], names, rotation=45, ha='right', fontsize=16)
    axs.set_yticks([x + .5 for x in range(len(clean.index))], clean.index, rotation=0, fontsize=16)
    axs.set_yticklabels(textwrap.fill(y.get_text(), 28) for y in axs.get_yticklabels())
    axs.figure.axes[-1].tick_params(labelsize=15)
    axs.figure.axes[-1].yaxis.label.set_size(15)
    axs.set(xlabel=None)

    # save and close
    plt.savefig(filemap)
    plt.cla()
    plt.close()
