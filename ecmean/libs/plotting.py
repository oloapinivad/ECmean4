#!/usr/bin/env python3
'''
Shared functions for XArray ECmean4
'''

##################
# PLOT FUNCTIONS #
##################

from matplotlib.colors import TwoSlopeNorm
import matplotlib.pyplot as plt
import seaborn as sns


def heatmap_comparison(relative_table, diag, filemap):
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
    axs.set_title(f'{title} {diag.modelname} {diag.year1} {diag.year2}', fontsize=25)
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
