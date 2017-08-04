import sys, os, math
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from core.DataSet import DataSet

datfile_dir = '/Users/anatale/UCSF/Grabe_lab/data/traak_data/etsc_out_20170718'

# traj_names = \
# [
# 'traakWT_full_npt.sim1',
# 'traakWT_full_npt.sim2',
# 'traakWT_s1s3_npt.sim1',
# 'traakWT_s1s3_npt.sim2',
# 'traakWT_s2s4_npt.sim1',
# 'traakWT_s2s4_npt.sim2',
# 'traakG124I_full_npt.sim1',
# 'traakG124I_full_npt.sim2',
# 'traakG124I_s1s3_npt.sim1',
# 'traakG124I_s1s3_npt.sim2',
# 'traakG124I_s2s4_npt.sim1',
# 'traakG124I_s2s4_npt.sim2',
# 'traakTM4_npt'
# ]

traj_names = \
[
'traakWT_s1s3_npt.sim1',
'traakWT_s1s3_npt.sim2',
'traakWT_s2s4_npt.sim1',
'traakWT_s2s4_npt.sim2',
'traakG124I_s1s3_npt.sim1',
'traakG124I_s1s3_npt.sim2',
'traakG124I_s2s4_npt.sim1',
'traakG124I_s2s4_npt.sim2'
]

feature_names = \
[
'sf_dihed_vector',
]
datfile_dir = datfile_dir + '/state_vectors'

figures_outdir = '/Users/anatale/UCSF/Grabe_lab/data/traak_data/etsc_out_20170718/figures'

try:
    os.makedirs(figures_outdir)
except OSError:
    if not os.path.isdir(figures_outdir):
        raise

# load data
datasets = {}
for elem in traj_names:
    filenames = [os.path.join(datfile_dir, '%s_%s.dat' % (elem, i)) for i in feature_names]
    datasets[elem] = [DataSet(infilename=i).dataset_to_dict() for i in filenames]

layers = \
{
'S0gly':['133A','133B','242A','242B'],
'S1top':['132A','132B','241A','241B'],
'S2S1':['131A','131B','240A','240B'],
'S3S2':['130A','130B','239A','239B'],
'S4S3':['129A','129B','238A','238B'],
'S4bottom':['129A','129B','238A','238B']
}

plt.figure(figsize=(12,7))
gs = GridSpec(2,2)
gs.update(hspace=0.3)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])
#ax3 = plt.subplot(gs[2,0])
ax4 = plt.subplot(gs[0,1])
ax5 = plt.subplot(gs[1,1])
#ax6 = plt.subplot(gs[2,1])

data1wt = np.zeros((6,1,3),dtype=np.float64)
data2wt = np.zeros((4,1,3),dtype=np.float64)
data3wt = np.zeros((6,4,3),dtype=np.float64)

data4mut = np.zeros((6,1,3),dtype=np.float64)
data5mut = np.zeros((4,1,3),dtype=np.float64)
data6mut = np.zeros((6,4,3),dtype=np.float64)

# flips per layer
for i,key1 in enumerate(['S0gly','S1top','S2S1','S3S2','S4S3','S4bottom']):
    print key1
    if key1 == 'S4bottom':
        angle = 'chi1'
    else:
        angle = 'psi'
    wt_counts = np.zeros((3,), dtype=np.float32)
    mut_counts = np.zeros((3,), dtype=np.float32)
    for elem in datasets:
        print elem
        arr_list = []
        for key2 in layers[key1]:
            key3 = '%s_%s_flip' % (key2, angle)
            arr_list.append(datasets[elem][0][key3])
        this_layer = np.stack(arr_list, axis=1)
        this_layer = np.sum(this_layer, axis=1)
        if elem.startswith('traakWT'):
            wt_counts[0] += np.shape(this_layer[this_layer==0])[0]
            wt_counts[1] += np.shape(this_layer[this_layer==1])[0]
            wt_counts[2] += np.shape(this_layer[this_layer>1])[0]
        elif elem.startswith('traakG124I'):
            mut_counts[0] += np.shape(this_layer[this_layer==0])[0]
            mut_counts[1] += np.shape(this_layer[this_layer==1])[0]
            mut_counts[2] += np.shape(this_layer[this_layer>1])[0]
    wt_percentage = wt_counts / np.sum(wt_counts)
    mut_percentage = mut_counts / np.sum(mut_counts)
    data1wt[i,0,:] = wt_percentage
    data4mut[i,0,:] = mut_percentage

# flips per strand
for i,expr in enumerate([re.compile('1..A'),re.compile('2..A'),re.compile('1..B'),re.compile('2..B')]):
    wt_counts = np.zeros((3,), dtype=np.float32)
    mut_counts = np.zeros((3,), dtype=np.float32)
    for elem in datasets:
        print elem
        arr_list = []
        for key in datasets[elem][0]:
            if expr.search(key):
                print 'found', key
                arr_list.append(datasets[elem][0][key])
        this_strand = np.stack(arr_list, axis=1)
        this_strand = np.sum(this_strand, axis=1)
        print this_strand
        if elem.startswith('traakWT'):
            wt_counts[0] += np.shape(this_strand[this_strand==0])[0]
            wt_counts[1] += np.shape(this_strand[this_strand==1])[0]
            wt_counts[2] += np.shape(this_strand[this_strand>1])[0]
        elif elem.startswith('traakG124I'):
            mut_counts[0] += np.shape(this_strand[this_strand==0])[0]
            mut_counts[1] += np.shape(this_strand[this_strand==1])[0]
            mut_counts[2] += np.shape(this_strand[this_strand>1])[0]
    wt_percentage = wt_counts / np.sum(wt_counts)
    mut_percentage = mut_counts / np.sum(mut_counts)
    data2wt[i,0,:] = wt_percentage
    data5mut[i,0,:] = mut_percentage

print data2wt

def complicated_bar_chart(ax, data, labels, title):
    n_groups = np.shape(data)[0]
    group_idx = np.arange(n_groups)
    n_bars = np.shape(data)[1]
    group_width = 0.8
    bar_width = group_width/n_bars
    xtick_list = []
    for i in range(n_bars):
        bar_positions = group_idx+i*bar_width
        ax.bar(bar_positions,data[:,i,0],bottom=0.0,width=bar_width,color='black',edgecolor='black',linewidth=0.5)
        ax.bar(bar_positions,data[:,i,1],bottom=data[:,i,0],width=bar_width,color='darksalmon',edgecolor='black',linewidth=0.5)
        ax.bar(bar_positions,data[:,i,2],bottom=np.sum([data[:,i,0],data[:,i,1]],axis=0),width=bar_width,color='red',edgecolor='black',linewidth=0.5)
        xtick_list.append(bar_positions)
    ax.set_xticks(np.concatenate(xtick_list,axis=0))
    if labels:
        ax.set_xticklabels(labels, rotation=90)
    else:
        ax.set_xticklabels(np.concatenate(xtick_list,axis=0), rotation=90)
    if title:
        ax.set_title(title)
    ax.set_ylim((0.0,1.0))

# data must be NxMx3 array
# N=number of groups
# M=bars per group
complicated_bar_chart(ax1, data1wt, ['S0t','S0S1','S1S2','S2S3','S3S4','S4b'], 'WT')
ax1.set_ylabel('P(flip) by layer')
complicated_bar_chart(ax4, data4mut, ['S0t','S0S1','S1S2','S2S3','S3S4','S4b'], 'G124I')
complicated_bar_chart(ax2, data2wt, ['A1','A2','B1','B2'], None)
ax2.set_ylabel('P(flip) by strand')
complicated_bar_chart(ax5, data5mut, ['A1','A2','B1','B2'], None)
plt.show()
