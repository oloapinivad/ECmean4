#!/usr/bin/env python3
'''
Shared functions for XArray ECmean4
'''

##################
# PLOT FUNCTIONS #
##################

from matplotlib.colors import TwoSlopeNorm
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
import textwrap


def heatmap_comparison_pi(relative_table, diag, filemap):
    """Function to produce a heatmap - seaborn based - for Performance Indices
    based on CMIP6 ratio"""

    # real plot
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True, figsize=(15, 8))

    thr = [0, 1, 5]
    tictoc = [0, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5]
    title = 'CMIP6 RELATIVE PI'
    myfield = relative_table

    # axs.subplots_adjust(bottom=0.2)
    # pal = sns.diverging_palette(h_neg=130, h_pos=10, s=99, l=55, sep=3, as_cmap=True)
    tot = (len(myfield.columns))
    sss = (len(set([tup[1] for tup in myfield.columns])))
    divnorm = TwoSlopeNorm(vmin=thr[0], vcenter=thr[1], vmax=thr[2])
    pal = sns.color_palette("Spectral_r", as_cmap=True)
    # pal = sns.diverging_palette(220, 20, as_cmap=True)
    chart = sns.heatmap(myfield, norm=divnorm, cmap=pal,
                        cbar_kws={"ticks": tictoc, 'label': title},
                        ax=axs, annot=True, linewidth=0.5,
                        annot_kws={'fontsize': 11, 'fontweight': 'bold'})
    chart = chart.set_facecolor('whitesmoke')
    axs.set_title(f'{title} {diag.modelname} {diag.expname} {diag.year1} {diag.year2}', fontsize=25)
    axs.vlines(list(range(sss, tot + sss, sss)), ymin=-1, ymax=len(myfield.index), colors='k')
    axs.hlines(len(myfield.index) - 1, xmin=-1, xmax=len(myfield.columns), colors='purple', lw=2)
    names = [' '.join(x) for x in myfield.columns]
    axs.set_xticks([x + .5 for x in range(len(names))], names, rotation=45, ha='right', fontsize=15)
    axs.set_yticks([x + .5 for x in range(len(myfield.index))], myfield.index, rotation=0, fontsize=18)
    axs.set(xlabel=None)

    # save and close
    plt.savefig(filemap)
    plt.cla()
    plt.close()


def heatmap_comparison_gm(data_table, mean_table, std_table, diag, filemap):
    """Function to produce a heatmap - seaborn based - for Global Mean
    based on season-averaged standard deviation ratio"""

    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True, figsize=(18, 14))

    ratio = (data_table - mean_table)/std_table


    title = 'GLOBAL MEAN'
    mask = ratio[('ALL','Global')].notna()
    clean = ratio[mask]
    thr = 10
    tictoc = range(-thr, thr+1)
    pal = ListedColormap(sns.color_palette("seismic", n_colors=21))
    tot = (len(clean.columns))
    sss = (len(set([tup[1] for tup in clean.columns])))

    chart = sns.heatmap(clean, annot=data_table[mask], vmin=-thr-0.5, vmax=thr+0.5, center=0,
                        annot_kws={'va':'bottom', 'fontsize': 12}, 
                        cbar_kws={'ticks': tictoc, 'label': 'Model Bias (as standard deviation of interannual variability from observations)'},
                        fmt = '.2f', cmap = pal)
    chart = sns.heatmap(clean, annot=mean_table[mask], vmin=-thr-0.5, vmax=thr+0.5, center=0,
                        annot_kws={'va':'top', 'fontsize': 8, 'fontstyle': 'italic'}, 
                        fmt = '.2f', cmap = pal, cbar = False)
  
    chart = chart.set_facecolor('whitesmoke')
    axs.set_title(f'{title} {diag.modelname} {diag.expname} {diag.year1} {diag.year2}', fontsize=25)
    axs.vlines(list(range(sss, tot + sss, sss)), ymin=-1, ymax=len(clean.index), colors='k')
    names = [' '.join(x) for x in clean.columns]
    axs.set_xticks([x + .5 for x in range(len(names))], names, rotation=45, ha='right', fontsize=16)
    axs.set_yticks([x + .5 for x in range(len(clean.index))], clean.index, rotation=0, fontsize=16)
    axs.set_yticklabels(textwrap.fill(y.get_text(), 20) for y in axs.get_yticklabels())
    axs.figure.axes[-1].tick_params(labelsize=15)
    axs.figure.axes[-1].yaxis.label.set_size(15)
    axs.set(xlabel=None)

    # save and close
    plt.savefig(filemap)
    plt.cla()
    plt.close()