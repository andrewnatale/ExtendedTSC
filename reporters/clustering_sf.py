from __future__ import print_function
import sys, os, math, re
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.gridspec import GridSpecFromSubplotSpec
from analysis.SimpleFeatures import SimpleFeatures
from core.TimeseriesDataSet import TimeseriesDataSet as tsds

loop_resids = [(112,113,114,115,116,117,118,119,120,121,122,123,124,125,126),
               (142,143,144,145,146,147),
               (148,149,150,151),
               (251,252,253,254,255,256),
               (257,258,259,260,261,262,263,264,265,266,267,268,269)]

def featurize():
    topo = os.path.abspath(sys.argv[1])
    try:
        traj = os.path.abspath(sys.argv[2])
        segids = ['TK1A','TK1B']
        outdir,basename  = os.path.split(traj)
        basename = basename.split('.')[0]
    except:
        traj = None
        segids = ['A','B']
        outdir,basename  = os.path.split(topo)
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
            f.write_data(os.path.join('clustering_dihed', '%s_loop_%d-%d_%s_dihedrals' % (basename, loop_resids[loop_idx][0], loop_resids[loop_idx][-1], segid)))

def load_and_transform(loop_idx):

    tgt_loop = loop_resids[loop_idx]
    searcher = re.compile('%d-%d_TK1' % (tgt_loop[0], tgt_loop[-1]))
    tracking_info = {}
    counter = 0
    datasets = []
    for filename in os.listdir('clustering_dihed'):
        if searcher.search(filename):
            ds = tsds(infilename=os.path.join('clustering_dihed', filename))
            tracking_info[filename] = (counter, ds.get_length())
            datasets.append(ds)
            counter += 1

    pdb_searcher = re.compile('xtal_clean_loop_%d-%d_' % (tgt_loop[0], tgt_loop[-1]))
    pdb_tracking_info = {}
    pdb_counter = 0
    pdb_datasets = []
    for filename in os.listdir('clustering_dihed'):
        if pdb_searcher.search(filename):
            ds = tsds(infilename=os.path.join('clustering_dihed', filename))
            pdb_tracking_info[filename] = (counter, ds.get_length())
            pdb_datasets.append(ds)
            pdb_counter += 1

    vectorized = []
    for ds in datasets:
        tmparr = np.empty((ds.data.shape[0]*2, ds.data.shape[1]), dtype=np.float64)
        for i in range(ds.data.shape[0]):
            q1 = np.cos(ds.data[i,:])
            q2 = np.sin(ds.data[i,:])
            tmparr[i*2,:] = q1
            tmparr[(i*2)+1,:] = q2
        vectorized.append(tmparr)

    pdb_vectorized = []
    for ds in pdb_datasets:
        tmparr = np.empty((ds.data.shape[0]*2, ds.data.shape[1]), dtype=np.float64)
        for i in range(ds.data.shape[0]):
            q1 = np.cos(ds.data[i,:])
            q2 = np.sin(ds.data[i,:])
            tmparr[i*2,:] = q1
            tmparr[(i*2)+1,:] = q2
        pdb_vectorized.append(tmparr)

    vectors = np.concatenate(vectorized, axis=1)
    print(vectors.shape)
    print(np.amin(vectors), np.amax(vectors))

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

    separate_groups = [pca.transform(i.T) for i in vectorized]
    pdb_tform = [pca.transform(i.T) for i in pdb_vectorized]
    for arr in separate_groups[0:4]:
        ax5.scatter(arr[:,0],arr[:,1], s=5, edgecolor='none', color='deepskyblue', alpha=0.7, zorder=1)
        ax5.scatter(arr[0,0],arr[0,1], s=50, color='blue', edgecolor='black', marker='*', zorder=5)
    for arr in pdb_tform[0:2]:
        ax5.scatter(arr[0,0],arr[0,1], s=50, color='blue', edgecolor='black', marker='v', zorder=5)
    for arr in separate_groups[4:]:
        ax5.scatter(arr[:,0],arr[:,1], s=5, edgecolor='none', color='lightcoral', alpha=0.7, zorder=1)
        ax5.scatter(arr[0,0],arr[0,1], s=50, color='red', edgecolor='black', marker='*', zorder=5)
    for arr in pdb_tform[2:]:
        ax5.scatter(arr[0,0],arr[0,1], s=50, color='red', edgecolor='black', marker='v', zorder=5)
    ax5.set_xlabel('PC1')
    ax5.set_ylabel('PC2')

    for arr in separate_groups[0:4]:
        ax6.scatter(arr[:,2],arr[:,3], s=5, edgecolor='none', color='deepskyblue', alpha=0.7, zorder=1)
        ax6.scatter(arr[0,2],arr[0,3], s=50, color='blue', edgecolor='black', marker='*', zorder=5)
    for arr in pdb_tform[0:2]:
        ax6.scatter(arr[0,2],arr[0,3], s=50, color='blue', edgecolor='black', marker='v', zorder=5)
    for arr in separate_groups[4:]:
        ax6.scatter(arr[:,2],arr[:,3], s=5, edgecolor='none', color='lightcoral', alpha=0.7, zorder=1)
        ax6.scatter(arr[0,2],arr[0,3], s=50, color='red', edgecolor='black', marker='*', zorder=5)
    for arr in pdb_tform[2:]:
        ax6.scatter(arr[0,2],arr[0,3], s=50, color='red', edgecolor='black', marker='v', zorder=5)
    ax6.set_xlabel('PC3')
    ax6.set_ylabel('PC4')
    plt.suptitle('dPCA residues %d - %d' % (tgt_loop[0], tgt_loop[-1]))
    plt.show()
    #plt.savefig('trek1_dpca_%d-%d.png' % (tgt_loop[0], tgt_loop[-1]), bb_inches='tight')

    plt.figure(figsize=(10,7))
    gs = GridSpec(2, 1)
    ax10 = plt.subplot(gs[0,0])
    ax11 = plt.subplot(gs[1,0])
    #ax2 = ax1.twinx()
    ax10.plot(datasets[]['time']/1000.0, separate_groups[2][:,0], color='red')
    ax10.plot(datasets[3]['time']/1000.0, separate_groups[3][:,0], color='black')
    ax11.plot(datasets[0]['time']/1000.0, separate_groups[0][:,1], color='red', linestyle='--')
    ax11.plot(datasets[1]['time']/1000.0, separate_groups[1][:,1], color='black', linestyle='--')
    plt.show()
#featurize()

load_and_transform(3)
