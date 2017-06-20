import sys, os, math
import numpy as np
import ExtendedTSC

# topology and trajectory files
try:
    psffile = sys.argv[1]
    dcdfile = sys.argv[2]
except IndexError:
    psffile = None
    dcdfile = None

input_prefix = '/Volumes/data/andrew'

outtag = 'pore_geo'
# make output dir in cwd
outdir = os.path.join(os.getcwd(), outtag)
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
    filepairs = [
    (os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_full_trimmed/topo.traakG124I_full_npt.psf'),
     os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_full_trimmed/traakG124I_full.new_align.500ps.dcd')),
    (os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_S1S3_trimmed/topo.traakG124I_S1S3_npt.psf'),
     os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_S1S3_trimmed/traakG124I_S1S3.new_align.500ps.dcd')),
    (os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_S2S4_trimmed/topo.traakG124I_S2S4_npt.psf'),
     os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_S2S4_trimmed/traakG124I_S2S4.new_align.500ps.dcd')),
    (os.path.join(input_prefix, 'mackinnon_traak/traakWT_full_trimmed/topo.traakWT_full_npt.psf'),
     os.path.join(input_prefix, 'mackinnon_traak/traakWT_full_trimmed/traakWT_full.new_align.500ps.dcd')),
    (os.path.join(input_prefix, 'mackinnon_traak/traakWT_S1S3_trimmed/topo.traakWT_S1S3_npt.psf'),
     os.path.join(input_prefix, 'mackinnon_traak/traakWT_S1S3_trimmed/traakWT_S1S3.new_align.500ps.dcd')),
    (os.path.join(input_prefix, 'mackinnon_traak/traakWT_S2S4_trimmed/topo.traakWT_S2S4_npt.psf'),
     os.path.join(input_prefix, 'mackinnon_traak/traakWT_S2S4_trimmed/traakWT_S2S4.new_align.500ps.dcd'))
    ]
    filepairsTM4 = [
    (os.path.join(input_prefix, 'traakTM4_trimmed/topo.traakTM4_npt.psf'),
     os.path.join(input_prefix, 'traakTM4_trimmed/traakTM4.new_align.500ps.dcd')),
    ]

# selections for measurements
selections = [
('a_129og1',  'atom',   'segid TRKA and resid 129 and name OG1'),
('a_129o',    'atom',   'segid TRKA and resid 129 and name O'),
('a_130o',    'atom',   'segid TRKA and resid 130 and name O'),
('a_131o',    'atom',   'segid TRKA and resid 131 and name O'),
('a_132o',    'atom',   'segid TRKA and resid 132 and name O'),
('a_133o',    'atom',   'segid TRKA and resid 133 and name O'),
('a_238og1',  'atom',   'segid TRKA and resid 238 and name OG1'),
('a_238o',    'atom',   'segid TRKA and resid 238 and name O'),
('a_239o',    'atom',   'segid TRKA and resid 239 and name O'),
('a_240o',    'atom',   'segid TRKA and resid 240 and name O'),
('a_241o',    'atom',   'segid TRKA and resid 241 and name O'),
('a_242o',    'atom',   'segid TRKA and resid 242 and name O'),
('b_129og1',  'atom',   'segid TRKB and resid 129 and name OG1'),
('b_129o',    'atom',   'segid TRKB and resid 129 and name O'),
('b_130o',    'atom',   'segid TRKB and resid 130 and name O'),
('b_131o',    'atom',   'segid TRKB and resid 131 and name O'),
('b_132o',    'atom',   'segid TRKB and resid 132 and name O'),
('b_133o',    'atom',   'segid TRKB and resid 133 and name O'),
('b_238og1',  'atom',   'segid TRKB and resid 238 and name OG1'),
('b_238o',    'atom',   'segid TRKB and resid 238 and name O'),
('b_239o',    'atom',   'segid TRKB and resid 239 and name O'),
('b_240o',    'atom',   'segid TRKB and resid 240 and name O'),
('b_241o',    'atom',   'segid TRKB and resid 241 and name O'),
('b_242o',    'atom',   'segid TRKB and resid 242 and name O')
]

# to make selections in the traakTM4 trajectory, numbers must be adjusted
# chains A and B do not have different segids
# for res in chain A, subtract 27 from actual numbers
# for res in chain B, add 255
selections_traakTM4 = None
selections_traakTM4 = [
('a_129og1',  'atom',   'protein and resid 102 and name OG1'),
('a_129o',    'atom',   'protein and resid 102 and name O'),
('a_130o',    'atom',   'protein and resid 103 and name O'),
('a_131o',    'atom',   'protein and resid 104 and name O'),
('a_132o',    'atom',   'protein and resid 105 and name O'),
('a_133o',    'atom',   'protein and resid 106 and name O'),
('a_238og1',  'atom',   'protein and resid 211 and name OG1'),
('a_238o',    'atom',   'protein and resid 211 and name O'),
('a_239o',    'atom',   'protein and resid 212 and name O'),
('a_240o',    'atom',   'protein and resid 213 and name O'),
('a_241o',    'atom',   'protein and resid 214 and name O'),
('a_242o',    'atom',   'protein and resid 215 and name O'),
('b_129og1',  'atom',   'protein and resid 384 and name OG1'),
('b_129o',    'atom',   'protein and resid 384 and name O'),
('b_130o',    'atom',   'protein and resid 385 and name O'),
('b_131o',    'atom',   'protein and resid 386 and name O'),
('b_132o',    'atom',   'protein and resid 387 and name O'),
('b_133o',    'atom',   'protein and resid 388 and name O'),
('b_238og1',  'atom',   'protein and resid 493 and name OG1'),
('b_238o',    'atom',   'protein and resid 493 and name O'),
('b_239o',    'atom',   'protein and resid 494 and name O'),
('b_240o',    'atom',   'protein and resid 495 and name O'),
('b_241o',    'atom',   'protein and resid 496 and name O'),
('b_242o',    'atom',   'protein and resid 497 and name O')
]

os.chdir(outdir)

for pair in filepairs:
    a = ExtendedTSC.ExtendedTSC()
    a.load_dcd(pair[0],pair[1])
    a.measures_from_list(selections)
    a.run()
    output = pair[1].split('/')[-1].split('.')[0] + '.' + outtag
    a.write_data(output)

if filepairsTM4:
    for pair in filepairsTM4:
        a = ExtendedTSC.ExtendedTSC()
        a.load_dcd(pair[0],pair[1])
        a.measures_from_list(selections_traakTM4)
        a.run()
        output = pair[1].split('/')[-1].split('.')[0] + '.' + outtag
        a.write_data(output)
