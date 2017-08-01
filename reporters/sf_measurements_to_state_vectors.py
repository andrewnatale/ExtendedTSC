import sys, os, math
import numpy as np
import matplotlib.pyplot as plt
from core.DataSet import DataSet

datfile_dir = '/Users/anatale/UCSF/Grabe_lab/data/traak_data/etsc_out_20170718'

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

feature_names = ['filter_dihedrals',
                 'filter_carbonyls',
                 'filter_potassium',
                 'filter_water']

figures_outdir = '/Users/anatale/UCSF/Grabe_lab/data/traak_data/etsc_out_20170718/figures'

try:
    os.makedirs(figures_outdir)
except OSError:
    if not os.path.isdir(figures_outdir):
        raise

# load data
def load_raw_data():
    datasets = {}
    for elem in traj_names:
        filenames = [os.path.join(datfile_dir, '%s_%s_all_frames.dat' % (elem, i)) for i in feature_names]
        datasets[elem] = [DataSet(infilename=i) for i in filenames]
    return datasets

def dihedral_vis(save=False):
    layers = \
    {
    'S4':['129A','129B','238A','238B'],
    'S3':['130A','130B','239A','239B'],
    'S2':['131A','131B','240A','240B'],
    'S1':['132A','132B','241A','241B'],
    'S0':['133A','133B','242A','242B']
    }
    colors = {'S4':'red','S3':'blue','S2':'green','S1':'gold','S0':'violet'}
    labels = \
    {
    'S4':'S4 - T129,T238',
    'S3':'S3 - I130,V239',
    'S2':'S2 - G131,G240',
    'S1':'S1 - Y132,F241',
    'S0':'S0 - G133,G242'
    }
    for key1 in sorted(datasets.keys()):
        f, ax = plt.subplots()
        target = datasets[key1][0].primaryDS.dataset_to_dict()
        print target
        txt_offset = 0
        for key2 in sorted(layers.keys(),reverse=True):
        #for key2 in ['S1','S2','S3']:
            for res in layers[key2]:
                ax.scatter(target['%s_phi' % res], target['%s_psi' % res], s=2, c=colors[key2], edgecolor='none', alpha=0.3)
            ax.text(0.05,0.05+0.03*txt_offset, labels[key2], color=colors[key2],
                    transform=ax.transAxes,verticalalignment='bottom',horizontalalignment='left', fontsize=8)
            txt_offset += 1
        ax.set_xlim([math.pi * -1.0, math.pi])
        ax.set_ylim([math.pi * -1.0, math.pi])
        ax.set_xlabel('Phi (radians)')
        ax.set_ylabel('Psi (radians)')
        ax.set_title(key1)
        if save:
            f.savefig(os.path.join(figures_outdir,'%s_dihedrals.png' % key1), bb_inches='tight')
        else:
            plt.show()

tsc_writer = tsc()

def vectorize(key):
    target = datasets[key]

    dihedrals = target[0].primaryDS.dataset_to_dict()
    carbonyls = target[1].primaryDS.dataset_to_dict()
    potassium = target[2].primaryDS.dataset_to_dict()
    water = target[3].primaryDS.dataset_to_dict()

    results = np.zeros((9,np.shape(potassium['zsearchlist'])[1]), dtype=int)

    # potassium in the pore
    pp = potassium['zsearchlist']<carbonyls['S4bottom'][2,:]
    pp = np.sum(pp, axis=0)

    # potassium/water in S4
    p4 = np.logical_and(potassium['zsearchlist']>carbonyls['S4bottom'][2,:], potassium['zsearchlist']<carbonyls['S4S3'][2,:])
    p4 = np.sum(p4, axis=0)
    w4 = np.logical_and(water['zsearchlist']>carbonyls['S4bottom'][2,:], water['zsearchlist']<carbonyls['S4S3'][2,:])
    w4 = np.sum(w4, axis=0)

    # potassium/water in S3
    p3 = np.logical_and(potassium['zsearchlist']>carbonyls['S4S3'][2,:], potassium['zsearchlist']<carbonyls['S3S2'][2,:])
    p3 = np.sum(p3, axis=0)
    w3 = np.logical_and(water['zsearchlist']>carbonyls['S4S3'][2,:], water['zsearchlist']<carbonyls['S3S2'][2,:])
    w3 = np.sum(w3, axis=0)

    # potassium/water in S2
    p2 = np.logical_and(potassium['zsearchlist']>carbonyls['S3S2'][2,:], potassium['zsearchlist']<carbonyls['S2S1'][2,:])
    p2 = np.sum(p2, axis=0)
    w2 = np.logical_and(water['zsearchlist']>carbonyls['S3S2'][2,:], water['zsearchlist']<carbonyls['S2S1'][2,:])
    w2 = np.sum(w2, axis=0)

    # potassium/water in S1
    p1 = np.logical_and(potassium['zsearchlist']>carbonyls['S2S1'][2,:], potassium['zsearchlist']<carbonyls['S1top'][2,:])
    p1 = np.sum(p1, axis=0)
    w1 = np.logical_and(water['zsearchlist']>carbonyls['S2S1'][2,:], water['zsearchlist']<carbonyls['S1top'][2,:])
    w1 = np.sum(w1, axis=0)

    # potassium extracellular (S0)
    p0 = potassium['zsearchlist']>carbonyls['S1top'][2,:]
    p0 = np.sum(p0, axis=0)

    # pinched carbonyls at s3
    s3 = np.concatenate([dihedrals['130A_psi'],dihedrals['130B_psi'],dihedrals['239A_psi'],dihedrals['239B_psi']], axis=0) > 0
    s3 = np.sum(s3, axis=0)
    # at s2
    s2 = np.concatenate([dihedrals['131A_psi'],dihedrals['131B_psi'],dihedrals['240A_psi'],dihedrals['240B_psi']], axis=0) < 0
    s2 = np.sum(s2, axis=0)
    # at s1
    s1 = np.concatenate([dihedrals['132A_psi'],dihedrals['132B_psi'],dihedrals['241A_psi'],dihedrals['241B_psi']], axis=0) > 0
    s1 = np.sum(s1, axis=0)

    vectorDS = tsc_writer._custom_dataset()
    vectorDS.copy_metadata(target[0].primaryDS)
    vectorDS.is_mask = True
    vectorDS.add_measurement(('pore_K','count','K+ in site below S4'), width=1)
    vectorDS.add_measurement(('s4_K','count','K+ in S4'), width=1)
    vectorDS.add_measurement(('s3_K','count','K+ in S3'), width=1)
    vectorDS.add_measurement(('s2_K','count','K+ in S2'), width=1)
    vectorDS.add_measurement(('s1_K','count','K+ in S1'), width=1)
    vectorDS.add_measurement(('s0_K','count','K+ in S0'), width=1)
    vectorDS.add_measurement(('s4_W','count','H2Os in S4'), width=1)
    vectorDS.add_measurement(('s3_W','count','H2Os in S3'), width=1)
    vectorDS.add_measurement(('s2_W','count','H2Os in S2'), width=1)
    vectorDS.add_measurement(('s1_W','count','H2Os in S1'), width=1)
    vectorDS.add_measurement(('s3_pinch','count','pinched carbonyls at S3'), width=1)
    vectorDS.add_measurement(('s2_pinch','count','pinched carbonyls at S2'), width=1)
    vectorDS.add_measurement(('s1_pinch','count','pinched carbonyls at S1'), width=1)

    vectorDS.add_collection(np.stack([pp,p4,p3,p2,p1,p0,w4,w3,w2,w1,s3,s2,s1], axis=0))
    vectorDS.setup_timesteps()
    tsc_writer._data_writer(vectorDS,outfilename=os.path.join(datfile_dir, 'state_vectors', '%s_vectorized.dat' % key))

datasets = load_raw_data()
#dihedral_vis(save=True)
try:
    os.makedirs(os.path.join(datfile_dir, 'state_vectors'))
except OSError:
    if not os.path.isdir(os.path.join(datfile_dir, 'state_vectors')):
        raise

for traj in traj_names:
    vectorize(traj)
