from __future__ import print_function
import sys, os, math, re
import numpy as np
from itertools import cycle
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
from matplotlib.gridspec import GridSpecFromSubplotSpec
from core.TimeseriesDataSet import TimeseriesDataSet as tsds

def dihedral(p1,p2,p3,p4):
    segA = p2 - p1
    segB = p3 - p2
    segC = p4 - p3
    AxB = np.cross(segA, segB, axisa=0, axisb=0, axisc=0)
    BxC = np.cross(segB, segC, axisa=0, axisb=0, axisc=0)
    n1 = AxB / np.linalg.norm(AxB, axis=0)
    n2 = BxC / np.linalg.norm(BxC, axis=0)
    x = np.sum(n1 * n2, axis=0)
    y = np.sum(np.cross(n1, segB / np.linalg.norm(segB, axis=0), axisa=0, axisb=0, axisc=0) * n2, axis=0)
    dihedrals = np.arctan2(y,x) * -1.0
    dihedrals[dihedrals<0] += 2*np.pi
    return dihedrals

def angle(p1,p2,p3):
    v1 = p1 - p2
    v2 = p3 - p2
    n1 = np.linalg.norm(v1, axis=0)
    n2 = np.linalg.norm(v2, axis=0)
    angles = np.arccos(np.sum(v1*v2, axis=0)/(n1*n2))
    return angles

figsize = (16,8)
# if you need more than 12 colors, you can add them here
colors = cycle(['darkred', 'blue', 'forestgreen', 'darkgoldenrod',
                'darkorchid', 'turquoise', 'crimson', 'deepskyblue',
                'darkorange', 'limegreen', 'teal', 'magenta'])

# convolution array to generate a running average
window = 10 # frames
convolution = np.ones((window,))/float(window)

# minor tick spacing for plots (ns)
minor_tick_space = 10

def compare_timeseries_plot(datalist, ax, y_lim=None, y_label=None, plt_label=None):
    # iterate over input data elements
    for idx,elem in enumerate(datalist):
        # pick a new color
        color_now = colors.next()
        # unpack the data tuple
        dataseries, timesteps, datalabel = elem
        # plot the original data as transparent points
        ax.scatter(timesteps, dataseries, color=color_now, s=5, edgecolor='none', alpha=0.5, zorder=1)
        # plot the running average as a darker line on top of the points
        ax.plot(
          np.convolve(timesteps, convolution, mode='valid'),
          np.convolve(dataseries, convolution, mode='valid'),
          color=color_now,
          zorder=5
          )
        # add a basic legend
        ax.text(0.01, 0.99-idx*0.08,
          datalabel, color=color_now,
          transform=ax.transAxes, verticalalignment='top', horizontalalignment='left')
    # set reasonable axis limits
    first_time = np.amin(np.concatenate([i[1] for i in datalist], axis=0))
    last_time = np.amax(np.concatenate([i[1] for i in datalist], axis=0))
    ax.set_xlim(first_time, last_time)
    if y_lim is not None:
        y_upper_lim = y_lim[1]
        y_lower_lim = y_lim[0]
    else:
        y_upper_lim = np.amax(np.concatenate([i[0] for i in datalist], axis=0))*1.2
        y_lower_lim = np.amin(np.concatenate([i[0] for i in datalist], axis=0))*0.8
        if y_lower_lim <= 2.0:
            y_lower_lim = 0.0
    ax.set_ylim(y_lower_lim, y_upper_lim)
    # create axis labels
    ax.set_xlabel('Time (nanoseconds)')
    # add more ticks
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(minor_tick_space))
    if y_label is not None:
        ax.set_ylabel(y_label)
    if plt_label is not None:
        ax.set_title(plt_label)

dat_dir =  '/Users/anatale/UCSF/Grabe_lab/scratch/tm_helix_features'

traj_names = \
[
'traakWT_full_npt.sim1',
'traakWT_s1s3_npt.sim1',
'traakWT_s2s4_npt.sim1',
'traakG124I_full_npt.sim1',
'traakG124I_s1s3_npt.sim1',
'traakG124I_s2s4_npt.sim1',
'traakTM4_npt'
]

# 'traakWT_full_npt.sim2',
# 'traakWT_s1s3_npt.sim2',
# 'traakWT_s2s4_npt.sim2',
# 'traakG124I_full_npt.sim2',
# 'traakG124I_s1s3_npt.sim2',
# 'traakG124I_s2s4_npt.sim2',

def compare_traj():
    # setup plot
    plt.figure(figsize=(12,7))
    gs0 = GridSpec(3, 3)
    gs0.update(hspace=0.4)
    gs0.update(wspace=0.2)
    ax1 = plt.subplot(gs0[0,0])
    ax2 = plt.subplot(gs0[1,0])
    ax3 = plt.subplot(gs0[2,0])
    ax4 = plt.subplot(gs0[0,1])
    ax5 = plt.subplot(gs0[1,1])
    ax6 = plt.subplot(gs0[2,1])
    ax7 = plt.subplot(gs0[0,2])
    ax8 = plt.subplot(gs0[1,2])
    ax9 = plt.subplot(gs0[2,2])

    axlist = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9]

    datasets = {}
    for i in traj_names:
        datasets[i] = (tsds(os.path.join(dat_dir, i+'.m2m3_features.dat')))

    for i,elem in enumerate(traj_names):

        data = [
        (datasets[elem]['m4a_m2b'].series[0], datasets[elem]['time']/1000.0, 'm4a_m2b'),
        (datasets[elem]['m4b_m2a'].series[0], datasets[elem]['time']/1000.0, 'm4b_m2a')
        ]
        compare_timeseries_plot(data, axlist[i], y_lim=(0.0,15.0))

        # data = [
        # (datasets[elem]['A_271_201'].series[0], datasets[elem]['time']/1000.0, 'A_271_201'),
        # (datasets[elem]['B_271_201'].series[0], datasets[elem]['time']/1000.0, 'B_271_201')
        # ]
        # compare_timeseries_plot(data, axlist[i], y_lim=(0.0,15.0))

        # data = [
        # (datasets[elem]['A_271_209'].series[0], datasets[elem]['time']/1000.0, 'A_271_209'),
        # (datasets[elem]['B_271_209'].series[0], datasets[elem]['time']/1000.0, 'B_271_209')
        # ]
        # compare_timeseries_plot(data, axlist[i], y_lim=(0.0,15.0))

        # angleA = angle(datasets[elem]['pivot_A_5'].series,datasets[elem]['pivot_A_6'].series,datasets[elem]['pivot_A_7'].series)
        # angleB = angle(datasets[elem]['pivot_B_5'].series,datasets[elem]['pivot_B_6'].series,datasets[elem]['pivot_B_7'].series)
        # data = [
        # (angleA, datasets[elem]['time']/1000.0, 'm3_A'),
        # (angleB, datasets[elem]['time']/1000.0, 'm3_B')
        # ]
        # compare_timeseries_plot(data, axlist[i], y_lim=(2.0,3.0))

        # angleA = angle(datasets[elem]['pivot_A_1'].series,datasets[elem]['pivot_A_2'].series,datasets[elem]['pivot_A_3'].series)
        # angleB = angle(datasets[elem]['pivot_B_1'].series,datasets[elem]['pivot_B_2'].series,datasets[elem]['pivot_B_3'].series)
        # data = [
        # (angleA, datasets[elem]['time']/1000.0, 'm2top_A'),
        # (angleB, datasets[elem]['time']/1000.0, 'm2top_B')
        # ]
        # compare_timeseries_plot(data, axlist[i], y_lim=(2.0,3.0))

        # angleA = angle(datasets[elem]['pivot_A_2'].series,datasets[elem]['pivot_A_3'].series,datasets[elem]['pivot_A_4'].series)
        # angleB = angle(datasets[elem]['pivot_B_2'].series,datasets[elem]['pivot_B_3'].series,datasets[elem]['pivot_B_4'].series)
        # data = [
        # (angleA, datasets[elem]['time']/1000.0, 'm2bot_A'),
        # (angleB, datasets[elem]['time']/1000.0, 'm2bot_B')
        # ]
        # compare_timeseries_plot(data, axlist[i], y_lim=(2.4,3.4))

        # dihedA = dihedral(datasets[elem]['pivot_A_1'].series,datasets[elem]['pivot_A_2'].series,datasets[elem]['pivot_A_3'].series,datasets[elem]['pivot_A_4'].series)
        # dihedB = dihedral(datasets[elem]['pivot_B_1'].series,datasets[elem]['pivot_B_2'].series,datasets[elem]['pivot_B_3'].series,datasets[elem]['pivot_B_4'].series)
        # data = [
        # (dihedA, datasets[elem]['time']/1000.0, 'm2dih_A'),
        # (dihedB, datasets[elem]['time']/1000.0, 'm2dih_B')
        # ]
        # compare_timeseries_plot(data, axlist[i], y_lim=(0,2*np.pi))

    plt.show()

def compare_features(traj):
    # setup plot
    plt.figure(figsize=(12,7))
    gs0 = GridSpec(3, 1)
    gs0.update(hspace=0.4)
    gs0.update(wspace=0.2)
    ax1 = plt.subplot(gs0[0,0])
    ax2 = plt.subplot(gs0[1,0])
    ax3 = plt.subplot(gs0[2,0])

    dataset = tsds(os.path.join(dat_dir, traj+'.m2m3_features.dat'))

    data1 = [
    (dataset['A_271_201'].series[0], dataset['time']/1000.0, 'A_271_201'),
    (dataset['B_271_201'].series[0], dataset['time']/1000.0, 'B_271_201')
    ]
    compare_timeseries_plot(data1, ax1, y_lim=(0.0,15.0))

    angleA = angle(dataset['pivot_A_5'].series,dataset['pivot_A_6'].series,dataset['pivot_A_7'].series)
    angleB = angle(dataset['pivot_B_5'].series,dataset['pivot_B_6'].series,dataset['pivot_B_7'].series)
    data2 = [
    (angleA, dataset['time']/1000.0, 'm3_A'),
    (angleB, dataset['time']/1000.0, 'm3_B')
    ]
    compare_timeseries_plot(data2, ax2, y_lim=(2.0,3.0))

    angleC = angle(dataset['pivot_A_2'].series,dataset['pivot_A_3'].series,dataset['pivot_A_4'].series)
    angleD = angle(dataset['pivot_B_2'].series,dataset['pivot_B_3'].series,dataset['pivot_B_4'].series)
    data3 = [
    (angleC, dataset['time']/1000.0, 'm2bot_A'),
    (angleD, dataset['time']/1000.0, 'm2bot_B')
    ]
    compare_timeseries_plot(data3, ax3, y_lim=(2.4,3.4))
    plt.show()

compare_features(traj_names[6])
