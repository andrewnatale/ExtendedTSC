from __future__ import print_function
import sys, os, shutil
import numpy as np
import MDAnalysis as mda
from analysis.SimpleFeatures import SimpleFeatures
from analysis.VolumeTracker import VolumeTracker
from analysis.ZSearch import ZSearch
from analysis.RMSDseries import RMSDseries

outdir = '/Users/anatale/UCSF/Grabe_lab/scratch/trek_filter_rmsds'
try:
    os.makedirs(outdir)
except OSError:
    if not os.path.isdir(outdir):
        raise

datadir = '/Users/anatale/UCSF/Grabe_lab/data/trek1_ReducedTrj'

traj_names = \
[
'trek1_lowK_apo_S3_npt',
'trek1_lowK_apo_S4_npt',
'trek1_lowK_holo_QM_npt',
'trek1_lowK_holo_PC_npt'
]

pdb_names = ['apo_xtal_clean.pdb','holo_xtal_clean.pdb']

refidx = [0,0,1,1]

filepairs = []
for idx,traj_name in enumerate(traj_names):
    for filename in os.listdir(datadir):
        if filename.startswith('topo.'+traj_name) and filename.endswith('.psf'):
            topo = os.path.abspath(os.path.join(datadir,filename))
        if filename.startswith(traj_name) and filename.endswith('.dcd'):
            traj = os.path.abspath(os.path.join(datadir,filename))
    ref = os.path.abspath(os.path.join(datadir,pdb_names[refidx[idx]]))
    filepairs.append((traj_name,topo,traj,ref))

descriptorlist = [
('m4a_m2b',      'distance',      '((segid TRKA and resid 280) or (segid TRKB and resid 155)) and name CA'),
('m4b_m2a',      'distance',      '((segid TRKB and resid 280) or (segid TRKA and resid 155)) and name CA'),
('A_271_201',    'distance',      'segid TRKA and resid 271 201 and name CG'),
('A_271_209',    'distance',      'segid TRKA and resid 271 209 and name CG'),
('B_271_201',    'distance',      'segid TRKB and resid 271 201 and name CG'),
('B_271_209',    'distance',      'segid TRKB and resid 271 209 and name CG')
]
descriptorlist_tm4 = [
('m4a_m2b',      'distance',      'protein and (resid 253 or resid 410) and name CA'),
('m4b_m2a',      'distance',      'protein and (resid 535 or resid 128) and name CA'),
('A_271_201',    'distance',      'protein and resid 244 174 and name CG'),
('A_271_209',    'distance',      'protein and resid 244 182 and name CG'),
('B_271_201',    'distance',      'protein and resid 526 456 and name CG'),
('B_271_209',    'distance',      'protein and resid 526 464 and name CG')
]

pivot_points = [
(141,142,143,144),
(156,157,158,159),
(174,175,176,177),
(183,184,185,186),
(190,191,192,193),
(207,208,209,210),
(217,218,219,220)
]

for idx,pp in enumerate(pivot_points):
    for segid in ['TRKA', 'TRKB']:
        descriptorlist.append(('pivot_%s_%d' % (segid[-1],idx+1), 'COG', 'segid %s and name CA and resid %d %d %d %d' % (segid, pp[0], pp[1], pp[2], pp[3])))
        if segid == 'TRKA':
            tm4_pp = (pp[0]-27, pp[1]-27, pp[2]-27, pp[3]-27)
        elif segid == 'TRKB':
            tm4_pp = (pp[0]+255, pp[1]+255, pp[2]+255, pp[3]+255)
        descriptorlist_tm4.append(('pivot_%s_%d' % (segid[-1],idx+1), 'COG', 'protein and name CA and resid %d %d %d %d' % (tm4_pp[0], tm4_pp[1], tm4_pp[2], tm4_pp[3])))

os.chdir(outdir)
for pair in filepairs:
    if pair[0].startswith('traakTM4_npt'):
        tgtlist = descriptorlist_tm4
    else:
        tgtlist = descriptorlist
    a = SimpleFeatures()
    a.load_traj(pair[1], pair[2], 50)
    a.run(tgtlist)
    a.write_data(pair[0]+'.m2m3_features')
