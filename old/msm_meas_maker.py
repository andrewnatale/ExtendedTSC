import sys, os, math
import numpy as np
from ExtendedTSC import ExtendedTSC
from ZSearch import ZSearch

# topology and trajectory files
try:
    psffile = sys.argv[1]
    dcdfile = sys.argv[2]
except IndexError:
    psffile = None
    dcdfile = None

input_prefix = '/Users/anatale/school/UCSF/Grabe_lab/data'

outtag = 'msm'
# make output dir in cwd
outdir = os.path.join(os.getcwd(), 'msm_data')
try:
    os.mkdir(outdir)
except OSError:
    if not os.path.isdir(outdir):
        raise

# if I'm not batch processing, I'm probably just testing on a trajectory with the non TM4 topology
if psffile and dcdfile:
    print 'Using inputs specified from command line'
    filepairs = [(os.path.abspath(psffile), os.path.abspath(dcdfile))]
    filepairsTM4 = None
else:
    # hardcode filenames for batch process
    # filepairs = [
    # (os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_full_trimmed/topo.traakG124I_full_npt.psf'),
    #  os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_full_trimmed/traakG124I_full.new_align.500ps.dcd')),
    # (os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_S1S3_trimmed/topo.traakG124I_S1S3_npt.psf'),
    #  os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_S1S3_trimmed/traakG124I_S1S3.new_align.500ps.dcd')),
    # (os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_S2S4_trimmed/topo.traakG124I_S2S4_npt.psf'),
    #  os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_S2S4_trimmed/traakG124I_S2S4.new_align.500ps.dcd')),
    # (os.path.join(input_prefix, 'mackinnon_traak/traakWT_full_trimmed/topo.traakWT_full_npt.psf'),
    #  os.path.join(input_prefix, 'mackinnon_traak/traakWT_full_trimmed/traakWT_full.new_align.500ps.dcd')),
    # (os.path.join(input_prefix, 'mackinnon_traak/traakWT_S1S3_trimmed/topo.traakWT_S1S3_npt.psf'),
    #  os.path.join(input_prefix, 'mackinnon_traak/traakWT_S1S3_trimmed/traakWT_S1S3.new_align.500ps.dcd')),
    # (os.path.join(input_prefix, 'mackinnon_traak/traakWT_S2S4_trimmed/topo.traakWT_S2S4_npt.psf'),
    #  os.path.join(input_prefix, 'mackinnon_traak/traakWT_S2S4_trimmed/traakWT_S2S4.new_align.500ps.dcd'))
    # ]
    filepairs = [
    (os.path.join(input_prefix, 'traakG124I_full_trimmed/topo.traakG124I_full_npt.psf'),
     os.path.join(input_prefix, 'traakG124I_full_trimmed/traakG124I_full_npt.all.200ps.filter_alignZ.dcd')),
    (os.path.join(input_prefix, 'traakWT_full_trimmed/topo.traakWT_full_npt.psf'),
     os.path.join(input_prefix, 'traakWT_full_trimmed/traakWT_full_npt.all.200ps.filter_alignZ.dcd')),
    ]
    filepairsTM4 = [
    (os.path.join(input_prefix, 'traakTM4_trimmed/topo.traakTM4_npt.psf'),
     os.path.join(input_prefix, 'traakTM4_trimmed/traakTM4_npt.all.200ps.filter_alignZ.dcd')),
    ]

# static selections
selections = [
('a_129og1',  'atom',   'segid TRKA and resid 129 and name OG1'),
('a_129o',    'atom',   'segid TRKA and resid 129 and name O'),
('a_130o',    'atom',   'segid TRKA and resid 130 and name O'),
('a_131o',    'atom',   'segid TRKA and resid 131 and name O'),
('a_132o',    'atom',   'segid TRKA and resid 132 and name O'),
('a_238og1',  'atom',   'segid TRKA and resid 238 and name OG1'),
('a_238o',    'atom',   'segid TRKA and resid 238 and name O'),
('a_239o',    'atom',   'segid TRKA and resid 239 and name O'),
('a_240o',    'atom',   'segid TRKA and resid 240 and name O'),
('a_241o',    'atom',   'segid TRKA and resid 241 and name O'),
('b_129og1',  'atom',   'segid TRKB and resid 129 and name OG1'),
('b_129o',    'atom',   'segid TRKB and resid 129 and name O'),
('b_130o',    'atom',   'segid TRKB and resid 130 and name O'),
('b_131o',    'atom',   'segid TRKB and resid 131 and name O'),
('b_132o',    'atom',   'segid TRKB and resid 132 and name O'),
('b_238og1',  'atom',   'segid TRKB and resid 238 and name OG1'),
('b_238o',    'atom',   'segid TRKB and resid 238 and name O'),
('b_239o',    'atom',   'segid TRKB and resid 239 and name O'),
('b_240o',    'atom',   'segid TRKB and resid 240 and name O'),
('b_241o',    'atom',   'segid TRKB and resid 241 and name O')
]

# to make selections in the traakTM4 trajectory, numbers must be adjusted
# for res in chain A, subtract 27 from actual numbers
# for res in chain B, add 255
selections_traakTM4 = [
('a_129og1',  'atom',   'protein and resid 102 and name OG1'),
('a_129o',    'atom',   'protein and resid 102 and name O'),
('a_130o',    'atom',   'protein and resid 103 and name O'),
('a_131o',    'atom',   'protein and resid 104 and name O'),
('a_132o',    'atom',   'protein and resid 105 and name O'),
('a_238og1',  'atom',   'protein and resid 211 and name OG1'),
('a_238o',    'atom',   'protein and resid 211 and name O'),
('a_239o',    'atom',   'protein and resid 212 and name O'),
('a_240o',    'atom',   'protein and resid 213 and name O'),
('a_241o',    'atom',   'protein and resid 214 and name O'),
('b_129og1',  'atom',   'protein and resid 384 and name OG1'),
('b_129o',    'atom',   'protein and resid 384 and name O'),
('b_130o',    'atom',   'protein and resid 385 and name O'),
('b_131o',    'atom',   'protein and resid 386 and name O'),
('b_132o',    'atom',   'protein and resid 387 and name O'),
('b_238og1',  'atom',   'protein and resid 493 and name OG1'),
('b_238o',    'atom',   'protein and resid 493 and name O'),
('b_239o',    'atom',   'protein and resid 494 and name O'),
('b_240o',    'atom',   'protein and resid 495 and name O'),
('b_241o',    'atom',   'protein and resid 496 and name O')
]

# those selections are nice, but all we really want is:
names = {
'S4bottom':['a_129og1','a_238og1','b_129og1','b_238og1'],
'S4S3':['a_129o','a_238o','b_129o','b_238o'],
'S3S2':['a_130o','a_239o','b_130o','b_239o'],
'S2S1':['a_131o','a_240o','b_131o','b_240o'],
'S1top':['a_132o','a_241o','b_132o','b_241o']
}
mod_selections = []
mod_selections_traakTM4 = []
for key in names:
    # rearrange expressions
    mod_selections.append((key, 'COG', '(%s)' % ') or ('.join([i[2] for i in selections if i[0] in names[key]])))
    mod_selections_traakTM4.append((key, 'COG', '(%s)' % ') or ('.join([i[2] for i in selections_traakTM4 if i[0] in names[key]])))
# print mod_selections
# print mod_selections_traakTM4

# dynamic selections
vol_selecttext = 'cyzone 3.0 8.0 -8.0 ((%s))' % ') or ('.join([i[2] for i in selections])
vol_selecttext_traakTM4 = 'cyzone 3.0 8.0 -8.0 ((%s))' % ') or ('.join([i[2] for i in selections_traakTM4])

search_K = 'name POT'
search_HOH = 'name OH2 and resname TIP3'

# print vol_selecttext
# print vol_selecttext_traakTM4

os.chdir(outdir)

if filepairs:
    for pair in filepairs:
        # static measurements
        a = ExtendedTSC()
        a.load_dcd(pair[0],pair[1], traj_stepsize=200)
        a.measures_from_list(mod_selections)
        a.run()
        output = pair[1].split('/')[-1].split('.')[0] + '.' + outtag
        a.write_data(output)
        # search for potassium
        b = ZSearch()
        b.load_dcd(pair[0],pair[1], traj_stepsize=200)
        b.run(vol_selecttext,search_K)
        output = pair[1].split('/')[-1].split('.')[0] + '.' + outtag + '.' + 'potassium'
        b.write_data(output)
        # search for water
        c = ZSearch()
        c.load_dcd(pair[0],pair[1], traj_stepsize=200)
        c.run(vol_selecttext,search_HOH)
        output = pair[1].split('/')[-1].split('.')[0] + '.' + outtag + '.' + 'water'
        c.write_data(output)

if filepairsTM4:
    for pair in filepairsTM4:
        # static measurements
        a = ExtendedTSC()
        a.load_dcd(pair[0],pair[1], traj_stepsize=200)
        a.measures_from_list(mod_selections_traakTM4)
        a.run()
        output = pair[1].split('/')[-1].split('.')[0] + '.' + outtag
        a.write_data(output)
        # search for potassium
        b = ZSearch()
        b.load_dcd(pair[0],pair[1], traj_stepsize=200)
        b.run(vol_selecttext_traakTM4,search_K)
        output = pair[1].split('/')[-1].split('.')[0] + '.' + outtag + '.' + 'potassium'
        b.write_data(output)
        # search for water
        c = ZSearch()
        c.load_dcd(pair[0],pair[1], traj_stepsize=200)
        c.run(vol_selecttext_traakTM4,search_HOH)
        output = pair[1].split('/')[-1].split('.')[0] + '.' + outtag + '.' + 'water'
        c.write_data(output)
