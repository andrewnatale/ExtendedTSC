from __future__ import print_function
import sys, os, math, re
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
from sklearn.cluster import DBSCAN
from sklearn import metrics
#from sklearn import manifold
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.gridspec import GridSpec
from matplotlib.gridspec import GridSpecFromSubplotSpec
from analysis.SimpleFeatures import SimpleFeatures
from core.TimeseriesDataSet import TimeseriesDataSet as tsds
from core.MergeDS import merge_along_features

# loop_resids = [(112,113,114,115,116,117,118,119,120,121,122,123,124,125,126),
#                (142,143,144,145,146,147),
#                (148,149,150,151),
#                (251,252,253,254,255,256),
#                (257,258,259,260,261,262,263,264,265,266,267,268,269)]

loop_resids = [(142,143,144,145,146,147,148,149),
               (251,252,253,254,255,256,257,258)]

outdir = '/Users/anatale/UCSF/Grabe_lab/scratch/trek1_clustering'

def featurize_dihedral():
    topo = os.path.abspath(sys.argv[1])
    try:
        traj = os.path.abspath(sys.argv[2])
        segids = ['TK1A','TK1B']
        _,basename  = os.path.split(traj)
        basename = basename.split('.')[0]
    except:
        traj = None
        segids = ['A','B']
        _,basename  = os.path.split(topo)
        basename = basename.split('.')[0]

    for segid in segids:
        for loop_idx in range(len(loop_resids)):
            dha_selections = []
            for resid in loop_resids[loop_idx]:
                phiselectext = 'segid %s and ((resid %d and name C) or (resid %d and name N CA C))' % (segid,resid-1,resid)
                psiselectext = 'segid %s and ((resid %d and name N CA C) or (resid %d and name N))' % (segid,resid,resid+1)
                dha_selections.append(('%d%s_phi' % (resid,segid[-1]), 'dihedral', phiselectext))
                dha_selections.append(('%d%s_psi' % (resid,segid[-1]), 'dihedral', psiselectext))
            f = SimpleFeatures()
            if traj is not None:
                f.load_traj(topo,traj,500)
            else:
                f.load_pdb(topo)
            f.run(dha_selections)
            f.write_data(os.path.join(outdir, '%s_loop_%d-%d_%s_dihedrals' % (basename, loop_resids[loop_idx][0], loop_resids[loop_idx][-1], segid)))

def featurize_cartesian():
    topo = os.path.abspath(sys.argv[1])
    try:
        traj = os.path.abspath(sys.argv[2])
        segids = ['TK1A','TK1B']
        _,basename  = os.path.split(traj)
        basename = basename.split('.')[0]
    except:
        traj = None
        segids = ['A','B']
        _,basename  = os.path.split(topo)
        basename = basename.split('.')[0]

    for segid in segids:
        for loop_idx in range(len(loop_resids)):
            ca_selections = []
            for resid in loop_resids[loop_idx]:
                caselectext = 'segid %s and resid %d and name CA' % (segid,resid)
                ca_selections.append(('%d%s_ca' % (resid,segid[-1]), 'atom', caselectext))
            f = SimpleFeatures()
            if traj is not None:
                f.load_traj(topo,traj,500)
            else:
                f.load_pdb(topo)
            f.run(ca_selections)
            f.write_data(os.path.join(outdir, '%s_loop_%d-%d_%s_cas' % (basename, loop_resids[loop_idx][0], loop_resids[loop_idx][-1], segid)))

def load_and_merge():

    trajnames = [
    'trek1_lowK_apoS3_npt',
    'trek1_lowK_apoS4_npt',
    'trek1_lowK_holoQM_npt',
    'trek1_lowK_holoPC_npt'
    ]

    datasets = []
    for elem in trajnames:
        filelist = []
        for filename in os.listdir(outdir):
            if filename.startswith(elem):
                filelist.append(os.path.join(outdir, filename))
        datasets.append(merge_along_features(filelist))

    return datasets

datasets = load_and_merge()

search1A = '1..A....'
search1B = '1..B....'
search2A = '2..A....'
search2B = '2..B....'

all_strand_1_arrays = []
all_strand_2_arrays = []
for i in datasets:
    newds1A = i.return_trimmed(search1A)
    newds1B = i.return_trimmed(search1B)
    all_strand_1_arrays.append(newds1A.data)
    all_strand_1_arrays.append(newds1B.data)

    newds2A = i.return_trimmed(search2A)
    newds2B = i.return_trimmed(search2B)
    all_strand_2_arrays.append(newds2A.data)
    all_strand_2_arrays.append(newds2B.data)

whichloop = 1
tgt_loop = loop_resids[whichloop]
if whichloop == 0:
    tgt_strand_data = all_strand_1_arrays
elif whichloop == 1:
    tgt_strand_data = all_strand_2_arrays

transformed = []
for arr in tgt_strand_data:
    tmparr = np.empty((arr.shape[0]*2, arr.shape[1]), dtype=np.float64)
    for i in range(arr.shape[0]):
        q1 = np.cos(arr[i,:])
        q2 = np.sin(arr[i,:])
        tmparr[i*2,:] = q1
        tmparr[(i*2)+1,:] = q2
    transformed.append(tmparr)

vectors = np.concatenate(transformed, axis=1)
print(vectors.shape)
print(np.amin(vectors), np.amax(vectors))

def reduce_by_pca():
    pca = PCA(n_components=4)

    pca.fit(vectors.T)
    var = pca.explained_variance_ratio_
    comp = pca.components_
    print(var)
    print(comp.shape)

    plt.figure(figsize=(10,7))
    gs0 = GridSpec(1, 2)
    gs00 = GridSpecFromSubplotSpec(4, 1, subplot_spec=gs0[0], hspace=0.1)
    gs01 = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[1], hspace=0.3)
    ax1 = plt.subplot(gs00[0,0])
    ax2 = plt.subplot(gs00[1,0])
    ax3 = plt.subplot(gs00[2,0])
    ax4 = plt.subplot(gs00[3,0])
    ax5 = plt.subplot(gs01[0,0])
    ax6 = plt.subplot(gs01[1,0])
    axes = [ax1,ax2,ax3,ax4,ax5,ax6]

    for j in range(comp.shape[0]):
        influence = np.empty((comp.shape[1]/2,), dtype=np.float64)
        tgt_comp = comp[j,:]
        tgt_ax = axes[j]
        for i in range(tgt_comp.shape[0]/2):
            influence[i] = np.square(tgt_comp[i*2]) + np.square(tgt_comp[i*2+1])
        xvals = np.arange(comp.shape[1]/4)
        tgt_ax.bar(xvals, influence[0::2], width=0.35, color='darkorange')
        tgt_ax.bar(xvals+0.35, influence[1::2], width=0.35, color='mediumslateblue')
        tgt_ax.set_xticks(xvals+0.35/2)
        tgt_ax.set_ylim((0,0.5))
    ax1.tick_params(labelbottom='off')
    ax2.tick_params(labelbottom='off')
    ax3.tick_params(labelbottom='off')
    ax4.set_xticklabels([str(i) for i in tgt_loop])
    ax1.text(0.01,0.99, 'PC1 - %03.2f' % var[0], color='black', transform=ax1.transAxes,verticalalignment='top',horizontalalignment='left')
    ax1.text(0.99,0.99, 'phi', color='darkorange', transform=ax1.transAxes,verticalalignment='top',horizontalalignment='right')
    ax1.text(0.99,0.89, 'psi', color='mediumslateblue', transform=ax1.transAxes,verticalalignment='top',horizontalalignment='right')
    ax1.set_title('Angle contributions to PCs')
    ax2.text(0.01,0.99, 'PC2 - %03.2f' % var[1], color='black', transform=ax2.transAxes,verticalalignment='top',horizontalalignment='left')
    ax3.text(0.01,0.99, 'PC3 - %03.2f' % var[2], color='black', transform=ax3.transAxes,verticalalignment='top',horizontalalignment='left')
    ax4.text(0.01,0.99, 'PC4 - %03.2f' % var[3], color='black', transform=ax4.transAxes,verticalalignment='top',horizontalalignment='left')
    ax4.set_xlabel('Residue number')

    colors = ['firebrick', 'red', 'darksalmon', 'sienna', 'cadetblue', 'skyblue', 'dodgerblue', 'blue']

    separate_groups = [pca.transform(i.T) for i in transformed]
    # pdb_tform = [pca.transform(i.T) for i in pdb_vectorized]
    # for arr in separate_groups[0:4]:
    #     ax5.scatter(arr[:,0],arr[:,1], s=5, edgecolor='none', color='deepskyblue', alpha=0.7, zorder=1)
    #     ax5.scatter(arr[0,0],arr[0,1], s=50, color='blue', edgecolor='black', marker='*', zorder=5)
    # for arr in pdb_tform[0:2]:
    #     ax5.scatter(arr[0,0],arr[0,1], s=50, color='blue', edgecolor='black', marker='v', zorder=5)
    # for arr in separate_groups[4:]:
    #     ax5.scatter(arr[:,0],arr[:,1], s=5, edgecolor='none', color='lightcoral', alpha=0.7, zorder=1)
    #     ax5.scatter(arr[0,0],arr[0,1], s=50, color='red', edgecolor='black', marker='*', zorder=5)
    # for arr in pdb_tform[2:]:
    #     ax5.scatter(arr[0,0],arr[0,1], s=50, color='red', edgecolor='black', marker='v', zorder=5)
    for idx,arr in enumerate(separate_groups):
        ax5.scatter(arr[:,0],arr[:,1], s=5, color=colors[idx], edgecolor='none', alpha=0.7, zorder=1)
        ax5.scatter(arr[0,0],arr[0,1], s=50, color=colors[idx], edgecolor='black', marker='*', zorder=5)
    ax5.set_xlabel('PC1')
    ax5.set_ylabel('PC2')

    # for arr in separate_groups[0:4]:
    #     ax6.scatter(arr[:,2],arr[:,3], s=5, edgecolor='none', color='deepskyblue', alpha=0.7, zorder=1)
    #     ax6.scatter(arr[0,2],arr[0,3], s=50, color='blue', edgecolor='black', marker='*', zorder=5)
    # for arr in pdb_tform[0:2]:
    #     ax6.scatter(arr[0,2],arr[0,3], s=50, color='blue', edgecolor='black', marker='v', zorder=5)
    # for arr in separate_groups[4:]:
    #     ax6.scatter(arr[:,2],arr[:,3], s=5, edgecolor='none', color='lightcoral', alpha=0.7, zorder=1)
    #     ax6.scatter(arr[0,2],arr[0,3], s=50, color='red', edgecolor='black', marker='*', zorder=5)
    # for arr in pdb_tform[2:]:
    #     ax6.scatter(arr[0,2],arr[0,3], s=50, color='red', edgecolor='black', marker='v', zorder=5)
    for idx,arr in enumerate(separate_groups):
        ax6.scatter(arr[:,2],arr[:,3], s=5, color=colors[idx], edgecolor='none', alpha=0.7, zorder=1)
        ax6.scatter(arr[0,3],arr[0,3], s=50, color=colors[idx], edgecolor='black', marker='*', zorder=5)
    ax6.set_xlabel('PC3')
    ax6.set_ylabel('PC4')
    plt.suptitle('dPCA residues %d - %d' % (tgt_loop[0], tgt_loop[-1]))
    plt.show()
    # #plt.savefig('trek1_dpca_%d-%d.png' % (tgt_loop[0], tgt_loop[-1]), bb_inches='tight')

    # fig3d = plt.figure()
    # ax = fig3d.add_subplot(111, projection='3d')
    # # for arr in separate_groups[0:4]:
    # #     ax.scatter(arr[:,0],arr[:,1], zs=arr[:,2], color='deepskyblue')
    # # for arr in separate_groups[4:]:
    # #     ax.scatter(arr[:,0],arr[:,1], zs=arr[:,2], color='lightcoral')
    # for idx,arr in enumerate(separate_groups):
    #     ax.scatter(arr[:,0],arr[:,1], zs=arr[:,2], color=colors[idx])
    # plt.show()

    return separate_groups

dimreduce = reduce_by_pca()
for elem in dimreduce:
    print(elem.shape)

breakdown = {}
breakdown['apo'] = np.concatenate(dimreduce[0:4], axis=0)
breakdown['holo'] = np.concatenate(dimreduce[4:], axis=0)

breakdown['apoS3'] = np.concatenate(dimreduce[0:2], axis=0)
breakdown['apoS4'] = np.concatenate(dimreduce[2:4], axis=0)

breakdown['holoQM'] = np.concatenate(dimreduce[4:6], axis=0)
breakdown['holoPC'] = np.concatenate(dimreduce[6:], axis=0)

tgt_name = 'apoS4'
tgt = breakdown[tgt_name]

db = DBSCAN(eps=0.3, min_samples=10).fit(tgt)
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

print('Estimated number of clusters for %s: %d' % (tgt_name, n_clusters_))

unique_labels = set(labels)
colors = [plt.cm.Spectral(each)
          for each in np.linspace(0, 1, len(unique_labels))]

gs1 = GridSpec(1, 2)
ax11 = plt.subplot(gs1[0,0])
ax12 = plt.subplot(gs1[0,1])

for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = [0, 0, 0, 1]

    class_member_mask = (labels == k)

    xy = tgt[class_member_mask & core_samples_mask]
    ax11.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=10)

    xy = tgt[class_member_mask & ~core_samples_mask]
    ax11.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=4)

    xy = tgt[class_member_mask & core_samples_mask]
    ax12.plot(xy[:, 2], xy[:, 3], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=10)

    xy = tgt[class_member_mask & ~core_samples_mask]
    ax12.plot(xy[:, 2], xy[:, 3], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=4)

ax11.set_xlabel('PC1')
ax11.set_ylabel('PC2')
ax12.set_xlabel('PC3')
ax12.set_ylabel('PC4')

plt.suptitle('Estimated number of clusters for %s: %d' % (tgt_name,n_clusters_))
plt.show()

# plt.figure(figsize=(10,7))
# gs = GridSpec(2, 1)
# ax10 = plt.subplot(gs[0,0])
# ax11 = plt.subplot(gs[1,0])
# #ax2 = ax1.twinx()
# ax10.plot(datasets[1]['time']/1000.0, separate_groups[2][:,0], color='red')
# ax10.plot(datasets[3]['time']/1000.0, separate_groups[3][:,0], color='black')
# ax11.plot(datasets[0]['time']/1000.0, separate_groups[0][:,1], color='red', linestyle='--')
# ax11.plot(datasets[1]['time']/1000.0, separate_groups[1][:,1], color='black', linestyle='--')
# plt.show()
