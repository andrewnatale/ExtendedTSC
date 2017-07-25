import sys, os, math, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.ticker import MaxNLocator
from core.TimeseriesCore import TimeseriesCore

# choose report to plot
rep_idx = int(sys.argv[1])

traj_names = [
'traakWT_full_npt.sim1',
'traakWT_full_npt.sim2',
'traakWT_s1s3_npt.sim1',
'traakWT_s1s3_npt.sim2',
'traakWT_s2s4_npt.sim1',
'traakWT_s2s4_npt.sim2',
'traakG124I_full_npt.sim1',
'traakG124I_full_npt.sim2',
'traakG124I_s1s3_npt.sim1',
'traakG124I_s1s3_npt.sim2',
'traakG124I_s2s4_npt.sim1',
'traakG124I_s2s4_npt.sim2',
'traakTM4_npt',
'traakExtPip2'
]

dat_root_dir = '/Users/anatale/UCSF/Grabe_lab/data/traak_data/etsc_out_20170718'

report_name = traj_names[rep_idx]

figures_outdir = os.path.join(dat_root_dir, 'figures')

reader = TimeseriesCore(verbose=False)

# load data sets
datasets = {}
datasets['basic'] = reader._data_reader(os.path.join(dat_root_dir, '%s_basic_features_all_frames.dat' % report_name))
datasets['rmsd_tm_helices_to_wt'] = reader._data_reader(os.path.join(dat_root_dir, '%s_rmsd_tm_helices_to_wt_all_frames.dat' % report_name))
datasets['rmsd_tm_helices_to_mut'] = reader._data_reader(os.path.join(dat_root_dir, '%s_rmsd_tm_helices_to_mut_all_frames.dat' % report_name))
datasets['rmsd_pore_helices_to_wt'] = reader._data_reader(os.path.join(dat_root_dir, '%s_rmsd_pore_helices_to_wt_all_frames.dat' % report_name))
datasets['rmsd_pore_helices_to_mut'] = reader._data_reader(os.path.join(dat_root_dir, '%s_rmsd_pore_helices_to_mut_all_frames.dat' % report_name))
datasets['rmsd_selectivity_filter_to_wt'] = reader._data_reader(os.path.join(dat_root_dir, '%s_rmsd_selectivity_filter_to_wt_all_frames.dat' % report_name))
datasets['rmsd_selectivity_filter_to_mut'] = reader._data_reader(os.path.join(dat_root_dir, '%s_rmsd_selectivity_filter_to_mut_all_frames.dat' % report_name))

datadicts = {}
for key in datasets:
    datadicts[key] = datasets[key].dataset_to_dict()

# setup figure
plt.figure(figsize=(16,9))
gs0 = GridSpec(1, 2)
gs00 = GridSpecFromSubplotSpec(4, 1, subplot_spec=gs0[0], hspace=0.15)
gs01 = GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[1], hspace=0.15)
#gs0.update(hspace=0.1)
ax1 = plt.subplot(gs00[0,0])
ax2 = plt.subplot(gs00[1,0])
ax3 = plt.subplot(gs00[2,0])
ax4 = plt.subplot(gs00[3,0])
ax5 = plt.subplot(gs01[0,0])
ax6 = plt.subplot(gs01[1,0])
ax7 = plt.subplot(gs01[2,0])

stdtime = datadicts['basic']['time'].flatten() / 1000.0
convolution = np.ones((50,))/50.0
conv_time = np.convolve(stdtime,convolution,mode='valid')

def formatter(ax, data1, data2, color1a, color2a, color1b, color2b, title=None, xlabel=None, ylabel=None, ylimits=None, hide_xtick=False):
    ax.scatter(stdtime, data1, color=color1a, s=5, edgecolor='none', alpha=0.5, zorder=1)
    ax.plot(conv_time, np.convolve(data1,convolution, mode='valid'), color=color1b, zorder=5)
    ax.scatter(stdtime, data2, color=color2a, s=5, edgecolor='none', alpha=0.5, zorder=1)
    ax.plot(conv_time, np.convolve(data2,convolution, mode='valid'), color=color2b, zorder=5)
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
ax1,
datadicts['basic']['fenA'][0,:],
datadicts['basic']['fenB'][0,:],
'skyblue',
'lightcoral',
'blue',
'darkred',
title='M4-M2\' fenestration distance',
xlabel=None,
ylabel='Angstroms',
ylimits=(4,18),
hide_xtick=True
)

ax1.text(0.01,0.99, 'Subunit A', color='blue', transform=ax1.transAxes,verticalalignment='top',horizontalalignment='left')
ax1.text(0.01,0.89, 'Subunit B', color='darkred', transform=ax1.transAxes,verticalalignment='top',horizontalalignment='left')

formatter(
ax2,
datadicts['basic']['trp_gapA'][0,:],
datadicts['basic']['trp_gapB'][0,:],
'skyblue',
'lightcoral',
'blue',
'darkred',
title='P1-M4\' separation distance',
xlabel=None,
ylabel='Angstroms',
ylimits=(4,16),
hide_xtick=True
)

wingA = np.sqrt(np.square(datadicts['basic']['P190_A'][0,:]) + np.square(datadicts['basic']['P190_A'][1,:]))
wingB = np.sqrt(np.square(datadicts['basic']['P190_B'][0,:]) + np.square(datadicts['basic']['P190_B'][1,:]))

formatter(
ax3,
wingA,
wingB,
'skyblue',
'lightcoral',
'blue',
'darkred',
title='M2/M3 bundle extension distance',
xlabel=None,
ylabel='Angstroms',
ylimits=(18,38),
hide_xtick=True
)

vecWA = datadicts['basic']['W262_A_CH'] - datadicts['basic']['W262_A_CB']
vecWA = vecWA / np.linalg.norm(vecWA, axis=0)
orientWA = np.rad2deg(np.arccos(vecWA[2,:]*-1.0))
vecWB = datadicts['basic']['W262_B_CH'] - datadicts['basic']['W262_B_CB']
vecWB = vecWB / np.linalg.norm(vecWB, axis=0)
orientWB = np.rad2deg(np.arccos(vecWB[2,:]*-1.0))

formatter(
ax4,
orientWA,
orientWB,
'skyblue',
'lightcoral',
'blue',
'darkred',
title='W262 orientation',
xlabel='Time (ns)',
ylabel='Degrees',
ylimits=(0,180),
hide_xtick=False
)

formatter(
ax5,
datadicts['rmsd_tm_helices_to_wt']['rmsdseries'][0,:],
datadicts['rmsd_tm_helices_to_mut']['rmsdseries'][0,:],
'plum',
'lightgreen',
'darkorchid',
'green',
title='Transmembrane helices only',
xlabel=None,
ylabel='Ca RMSD',
hide_xtick=True,
ylimits=(0,4),
)
ax5.text(0.99,0.11, 'to 4I9W (WT)', color='darkorchid', transform=ax5.transAxes,verticalalignment='bottom',horizontalalignment='right')
ax5.text(0.99,0.01, 'to 4RUE (G124I)', color='green', transform=ax5.transAxes,verticalalignment='bottom',horizontalalignment='right')

formatter(
ax6,
datadicts['rmsd_pore_helices_to_wt']['rmsdseries'][0,:],
datadicts['rmsd_pore_helices_to_mut']['rmsdseries'][0,:],
'plum',
'lightgreen',
'darkorchid',
'green',
title='Pore helices only',
xlabel=None,
ylabel='Ca RMSD',
hide_xtick=True,
ylimits=(0,1.5)
)

formatter(
ax7,
datadicts['rmsd_selectivity_filter_to_wt']['rmsdseries'][0,:],
datadicts['rmsd_selectivity_filter_to_mut']['rmsdseries'][0,:],
'plum',
'lightgreen',
'darkorchid',
'green',
title='Selectivity filter',
xlabel='Time (ns)',
ylabel='Backbone RMSD',
hide_xtick=False,
ylimits=(0,1.5)
)

plt.suptitle('%s rmsd and distance measurements' % report_name, fontsize=20)

# show or save plot
#plt.show()
plt.savefig(os.path.join(figures_outdir, '%s_rmsd_and_distance.png' % report_name), bb_inches='tight')
