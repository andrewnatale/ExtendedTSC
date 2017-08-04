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
datasets = {}
for elem in traj_names:
    filenames = [os.path.join(datfile_dir, '%s_%s_all_frames.dat' % (elem, i)) for i in feature_names]
    datasets[elem] = [DataSet(infilename=i) for i in filenames]

def vectorize_occupancy(key):
    target = datasets[key]

    carbonyls = target[1].dataset_to_dict()
    potassium = target[2].dataset_to_dict()
    water = target[3].dataset_to_dict()

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

    vectorDS = DataSet()
    vectorDS.copy_metadata(target[0])
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

    vectorDS.add_collection(np.stack([pp,p4,p3,p2,p1,p0,w4,w3,w2,w1], axis=0))
    vectorDS.setup_timesteps()
    vectorDS.write_dat(outfilename=os.path.join(datfile_dir, 'state_vectors', '%s_sf_occ_vector.dat' % key))

def generate_dihedral_cutoffs():
    # doesn't automatically generate cutoffs, but creates histograms
    # of the angle of interest that can be examined to extract cutoffs

    layers = \
    {
    'S0top':['133A','133B','242A','242B'],
    'S0S1':['132A','132B','241A','241B'],
    'S1S2':['131A','131B','240A','240B'],
    'S2S3':['130A','130B','239A','239B'],
    'S3S4':['129A','129B','238A','238B'],
    'S4bottom':['129A','129B','238A','238B']
    }

    for key1 in ['S0top','S0S1','S1S2','S2S3','S3S4','S4bottom']:
        if key1 == 'S4bottom':
            angle = 'chi1'
        else:
            angle = 'psi'
        angles = np.zeros((0,))
        for key2 in datasets:
            tmpdict = datasets[key2][0].dataset_to_dict()
            for res in layers[key1]:
                angles = np.concatenate([angles, tmpdict['%s_%s' % (res, angle)][0,:]], axis=0)
        # wrap into (-pi,+pi):
        angles[angles>np.pi] -= np.pi*2
        angles[angles<np.pi*-1] += np.pi*2
        hist, edges = np.histogram(angles, bins=100, range=(-1.0*np.pi,np.pi), density=True)
        centers = (edges[:-1] + edges[1:]) / 2
        # plot
        f, ax = plt.subplots()
        ax.plot(centers,hist,linewidth=2.0)
        ax.set_xlim([math.pi * -1.0, math.pi])
        ax.set_ylim([0,0.02])
        ax.set_xlabel('%s (radians)' % angle)
        ax.set_title(key1)
        plt.show()

def vectorize_angles(key):
    target = datasets[key]
    dihedrals = target[0].dataset_to_dict()

    layers = \
    {
    'S0top':(['133A','133B','242A','242B'], (-1.40,1.75)),
    'S0S1':(['132A','132B','241A','241B'], (-2.5,0.38)),
    'S1S2':(['131A','131B','240A','240B'], (-0.25,3.14)),
    'S2S3':(['130A','130B','239A','239B'], (-2.50,0.16)),
    'S3S4':(['129A','129B','238A','238B'], (-1.0,2.0)),
    'S4bottom':(['129A','129B','238A','238B'], (-0.20,2.50))
    }

    features = []
    arrays = []
    for key1 in ['S0top','S0S1','S1S2','S2S3','S3S4','S4bottom']:
        if key1 == 'S4bottom':
            angle = 'chi1'
        else:
            angle = 'psi'
        low, high = layers[key1][1]
        for res in layers[key1][0]:
            angles = dihedrals['%s_%s' % (res, angle)][0,:]
            # wrap into (-pi,+pi):
            angles[angles>np.pi] -= np.pi*2
            angles[angles<np.pi*-1] += np.pi*2

            true_if_flip = ~np.logical_and(angles>low, angles<high)
            true_if_flip = true_if_flip.astype(np.int8)
            features.append('%s_%s_flip' % (res, angle))
            arrays.append(true_if_flip)

    vectorDS = DataSet()
    vectorDS.copy_metadata(target[0])
    vectorDS.is_mask = True
    for elem in features:
        vectorDS.add_measurement((elem, 'boolean', 'True if flip'), width=1)
    vectorDS.add_collection(np.stack(arrays, axis=0))
    vectorDS.setup_timesteps()
    vectorDS.write_dat(outfilename=os.path.join(datfile_dir, 'state_vectors', '%s_sf_dihed_vector.dat' % key))

try:
    os.makedirs(os.path.join(datfile_dir, 'state_vectors'))
except OSError:
    if not os.path.isdir(os.path.join(datfile_dir, 'state_vectors')):
        raise

#generate_dihedral_cutoffs()

for traj in traj_names:
    vectorize_angles(traj)
    vectorize_occupancy(traj)
