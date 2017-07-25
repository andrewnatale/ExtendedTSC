import sys, os, math, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.ticker import MaxNLocator
from core.TimeseriesCore import TimeseriesCore

dat_root_dir = '/Users/anatale/UCSF/Grabe_lab/data/traak_data/etsc_out_20170721'

report_name = 'traakTM4_npt'

figures_outdir = os.path.join(dat_root_dir, 'figures')

reader = TimeseriesCore(verbose=False)

# load data sets
datasets = {}
datasets['down'] = reader._data_reader(os.path.join(dat_root_dir, '%s_rmsd_tm_helices_to_trek2down_all_frames.dat' % report_name))
datasets['up'] = reader._data_reader(os.path.join(dat_root_dir, '%s_rmsd_tm_helices_to_trek2up_all_frames.dat' % report_name))

datadicts = {}
for key in datasets:
    datadicts[key] = datasets[key].dataset_to_dict()

# setup figure
f,ax = plt.subplots()

stdtime = datadicts['down']['time'].flatten() / 1000.0
convolution = np.ones((50,))/50.0
conv_time = np.convolve(stdtime,convolution,mode='valid')

def formatter(ax, data1, data2, color1a, color2a, color1b, color2b, title=None, xlabel=None, ylabel=None, ylimits=None, hide_xtick=False):
    ax.scatter(stdtime, data1, color=color1a, s=5, edgecolor='none', alpha=0.5, zorder=1)
    ax.plot(conv_time, np.convolve(data1,convolution, mode='valid'), color=color1b, zorder=5)
    ax.scatter(stdtime, data2, color=color2a, s=5, edgecolor='none', alpha=0.5, zorder=1)
    ax.plot(conv_time, np.convolve(data2,convolution, mode='valid'), color=color2b, zorder=5)
    diff = data1-data2
    ax.scatter(stdtime, diff, color='black', s=5, edgecolor='none', alpha=0.1, zorder=1)
    ax.plot(conv_time, np.convolve(diff,convolution, mode='valid'), color='black', zorder=5)
    ax.set_xlim(stdtime[0],stdtime[-1])
    if ylimits:
        ax.set_ylim(ylimits)
    if hide_xtick:
        plt.setp(ax.get_xticklabels(), visible=False)
    locator = MaxNLocator(nbins=6, prune='lower', interger=True)
    ax.yaxis.set_major_locator(locator)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, fontsize=12)
    # gridlines = ax.get_xgridlines()
    # for line in gridlines:
    #     line.set_linestyle('-.')

formatter(
ax,
datadicts['down']['rmsdseries'][0,:],
datadicts['up']['rmsdseries'][0,:],
'plum',
'lightgreen',
'darkorchid',
'green',
title='TRAAK WT + TM4',
xlabel=None,
ylabel='Ca RMSD',
hide_xtick=False,
ylimits=(0,6),
)
ax.text(0.99,0.11, 'to TREK2 down (4XDJ)', color='darkorchid', transform=ax.transAxes,verticalalignment='bottom',horizontalalignment='right')
ax.text(0.99,0.01, 'to TREK2 up (4BW5)', color='green', transform=ax.transAxes,verticalalignment='bottom',horizontalalignment='right')

# show or save plot
plt.show()
#plt.savefig(os.path.join(figures_outdir, '%s_rmsd_and_distance.png' % report_name), bb_inches='tight')
