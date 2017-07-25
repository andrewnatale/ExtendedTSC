import sys, os, math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from core.TimeseriesCore import TimeseriesCore as tsc

datfile_dir = '/Users/anatale/UCSF/Grabe_lab/data/traak_data/etsc_out_20170718'

traj_names = \
[
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
'traakTM4_npt'
]

# traj_names = ['traakWT_full_npt.sim1',]

feature_names = \
[
'filter_carbonyls',
'filter_potassium',
'filter_water'
]

figures_outdir = '/Users/anatale/UCSF/Grabe_lab/data/traak_data/etsc_out_20170718/figures'

try:
    os.makedirs(figures_outdir)
except OSError:
    if not os.path.isdir(figures_outdir):
        raise

reader = tsc()

# load data
datasets = {}
for elem in traj_names:
    filenames = [os.path.join(datfile_dir, '%s_%s_all_frames.dat' % (elem, i)) for i in feature_names]
    tmplist = [reader._data_reader(i) for i in filenames]
    datasets[elem] = [i.dataset_to_dict() for i in tmplist]

carbonyl_keys = \
[
'S0gly',
'S1top',
'S2S1',
'S3S2',
'S4S3',
'S4bottom'
]

compute_dict = {}
for key in traj_names:
    tgt = datasets[key][0]
    # compute reference z coord for alignment
    carbonylZ = np.stack([tgt[carbonyl_key][2,:] for carbonyl_key in carbonyl_keys], axis=0)
    zref = np.sum(carbonylZ, axis=0) / 6.0
    n_frames = np.shape(carbonylZ)[1]
    # align carbonyl coords
    for i in range(np.shape(carbonylZ)[0]):
        carbonylZ[i,:] = carbonylZ[i,:] - zref
    carbonylAvgs = np.sum(carbonylZ, axis=1) / n_frames
    # align potassium and water coords
    potassium = datasets[key][1]['zsearchlist']
    waters = datasets[key][2]['zsearchlist']
    for i in range(np.shape(waters)[0]):
        waters[i,:] = waters[i,:] - zref
    for i in range(np.shape(potassium)[0]):
        potassium[i,:] = potassium[i,:] - zref
    # remove nan and compute histograms
    histK, edgesK = np.histogram(potassium[~np.isnan(potassium)], bins=100, range=(-11.0,9.0), density=True)
    centersK = (edgesK[:-1] + edgesK[1:]) / 2
    histW, edgesW = np.histogram(waters[~np.isnan(waters)], bins=100, range=(-11.0,9.0), density=True)
    centersW = (edgesW[:-1] + edgesW[1:]) / 2
    # return everything as a list
    compute_dict[key] = [n_frames,carbonylAvgs,histK,centersK,histW,centersW]

# total_framesWT = 0
# potassium_listWT = []
# water_listWT = []
# carbonyl_listWT = []
# for key in traj_names:
#     if key.startswith('traakWT'):
#         tgt = datasets[key][0]
#         # compute reference z coord for alignment
#         carbonylZ = np.stack([tgt[carbonyl_key][2,:] for carbonyl_key in carbonyl_keys], axis=0)
#         zref = np.sum(carbonylZ, axis=0) / 6.0
#         total_framesWT += np.shape(carbonylZ)[1]
#         # align carbonyl coords
#         for i in range(np.shape(carbonylZ)[0]):
#             carbonylZ[i,:] = carbonylZ[i,:] - zref
#         # align potassium and water coords
#         potassium = datasets[key][1]['zsearchlist']
#         waters = datasets[key][2]['zsearchlist']
#         for i in range(np.shape(waters)[0]):
#             waters[i,:] = waters[i,:] - zref
#         for i in range(np.shape(potassium)[0]):
#             potassium[i,:] = potassium[i,:] - zref
#         potassium = potassium[~np.isnan(potassium)]
#         waters = waters[~np.isnan(waters)]
#         potassium_listWT.append(potassium)
#         water_listWT.append(waters)
#         carbonyl_listWT.append(carbonylZ)
#
# total_framesMUT = 0
# potassium_listMUT = []
# water_listMUT = []
# carbonyl_listMUT = []
# for key in traj_names:
#     if key.startswith('traakG124I'):
#         tgt = datasets[key][0]
#         # compute reference z coord for alignment
#         carbonylZ = np.stack([tgt[carbonyl_key][2,:] for carbonyl_key in carbonyl_keys], axis=0)
#         zref = np.sum(carbonylZ, axis=0) / 6.0
#         total_framesMUT += np.shape(carbonylZ)[1]
#         # align carbonyl coords
#         for i in range(np.shape(carbonylZ)[0]):
#             carbonylZ[i,:] = carbonylZ[i,:] - zref
#         # align potassium and water coords
#         potassium = datasets[key][1]['zsearchlist']
#         waters = datasets[key][2]['zsearchlist']
#         for i in range(np.shape(waters)[0]):
#             waters[i,:] = waters[i,:] - zref
#         for i in range(np.shape(potassium)[0]):
#             potassium[i,:] = potassium[i,:] - zref
#         potassium = potassium[~np.isnan(potassium)]
#         waters = waters[~np.isnan(waters)]
#         potassium_listMUT.append(potassium)
#         water_listMUT.append(waters)
#         carbonyl_listMUT.append(carbonylZ)
#
# histK_allWT, edgesK_allWT = np.histogram(np.concatenate(potassium_listWT, axis=0), bins=100, range=(-11.0,9.0), density=True)
# centersK_allWT = (edgesK_allWT[:-1] + edgesK_allWT[1:]) / 2
# histW_allWT, edgesW_allWT = np.histogram(np.concatenate(water_listWT, axis=0), bins=100, range=(-11.0,9.0), density=True)
# centersW_allWT = (edgesW_allWT[:-1] + edgesW_allWT[1:]) / 2
# carbonylAvgs_allWT = np.sum(np.concatenate(carbonyl_listWT, axis=1), axis=1) / total_framesWT
#
# histK_allMUT, edgesK_allMUT = np.histogram(np.concatenate(potassium_listMUT, axis=0), bins=100, range=(-11.0,9.0), density=True)
# centersK_allMUT = (edgesK_allMUT[:-1] + edgesK_allMUT[1:]) / 2
# histW_allMUT, edgesW_allMUT = np.histogram(np.concatenate(water_listMUT, axis=0), bins=100, range=(-11.0,9.0), density=True)
# centersW_allMUT = (edgesW_allMUT[:-1] + edgesW_allMUT[1:]) / 2
# carbonylAvgs_allMUT = np.sum(np.concatenate(carbonyl_listMUT, axis=1), axis=1) / total_framesMUT

def generate_average(searchterm):
    count = 0
    for key in traj_names:
        if key.startswith(searchterm):
            if count == 0:
                aggregate_histK = compute_dict[key][2]
                aggregate_histW = compute_dict[key][4]
                aggregate_carbonyls = compute_dict[key][1]
                centers = compute_dict[key][3]
            else:
                aggregate_histK += compute_dict[key][2]
                aggregate_histW += compute_dict[key][4]
                aggregate_carbonyls += compute_dict[key][1]
            count += 1
    aggregate_histK = aggregate_histK / float(count)
    aggregate_histW = aggregate_histW / float(count)
    aggregate_carbonyls = aggregate_carbonyls / float(count)
    return centers,aggregate_histK,aggregate_histW,aggregate_carbonyls

def populate_axis(ax,data,labels):
    centers,histK,histW,carbonyls = data
    ax.plot(centers,histK,c='green',linewidth=3.0, zorder=5)
    ax.plot(centers,histW,c='skyblue',linewidth=3.0, zorder=3)
    ax.set_xlim((-11.0,9.0))
    ax.set_ylim((0.0,0.75))
    ax.set_ylabel('Probability')
    for elem in carbonyls:
        ax.axvline(x=elem,c='red',linewidth=1.0,linestyle='dashed')
    ax.set_title(labels)

plt.figure(figsize=(16,9))
gs = GridSpec(3,1)
gs.update(hspace=0.3)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])
ax3 = plt.subplot(gs[2,0])

populate_axis(ax1,generate_average('traakWT'),'TRAAK WT - average of 6 simulations')
populate_axis(ax2,generate_average('traakG124I'),'TRAAK G124I - average of 6 simulations')
populate_axis(ax3,generate_average('traakTM4'),'TRAAK WT + extended TM4 - 1 simulation, fully loaded')
ax3.set_xlabel('z-coordinate')
ax1.text(-6.1,.7, 'S4', color='red',verticalalignment='top',horizontalalignment='center')
ax1.text(-3.0,.7, 'S3', color='red',verticalalignment='top',horizontalalignment='center')
ax1.text(0.0,.7, 'S2', color='red',verticalalignment='top',horizontalalignment='center')
ax1.text(3.0,.7, 'S1', color='red',verticalalignment='top',horizontalalignment='center')
ax1.text(6.1,.7, 'S0', color='red',verticalalignment='top',horizontalalignment='center')
ax1.text(6.1,.7, 'S0', color='red',verticalalignment='top',horizontalalignment='center')
ax1.text(0.05,.95, 'potassium', color='green',transform=ax1.transAxes,verticalalignment='center',horizontalalignment='center')
ax1.text(0.05,.85, 'water', color='skyblue',transform=ax1.transAxes,verticalalignment='center',horizontalalignment='center')
# #populate_axis(ax2,(centersW_allWT,histK_allWT,histW_allWT,carbonylAvgs_allWT),'wt_preavg')
# populate_axis(ax2,(centersW_allMUT,histK_allMUT,histW_allMUT,carbonylAvgs_allMUT),'mut_preavg')
plt.show()
