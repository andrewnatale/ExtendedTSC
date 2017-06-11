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

outtag = 'basic'
# make output dir in cwd
outdir = os.path.join(os.getcwd(), 'basic_meas')
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
('S1top',      'COG',      'protein and (resid 132 or resid 241) and name O'),
('S4bottom',   'COG',      'protein and (resid 129 or resid 238) and name OG1'),
('P190_A',     'atom',     'segid TRKA and resid 190 and name CA'),
('P190_B',     'atom',     'segid TRKB and resid 190 and name CA'),
('trp_gapA',   'distance', 'segid TRKA and (resid 117 or resid 259) and name CA'),
('trp_gapB',   'distance', 'segid TRKB and (resid 117 or resid 259) and name CA'),
('fenA',       'distance', '(segid TRKA and resid 280) or (segid TRKB and resid 155) and name CA'),
('fenB',       'distance', '(segid TRKB and resid 280) or (segid TRKA and resid 155) and name CA'),
('expA',       'distance', 'segid TRKA and (resid 169 or resid 278) and name CA'),
('expB',       'distance', 'segid TRKB and (resid 169 or resid 278) and name CA'),
('W262_chi1A', 'dihedral', 'segid TRKA and resid 262 and (name N or name CA or name CB or name CG)'),
('W262_chi2A', 'dihedral', 'segid TRKA and resid 262 and (name CA or name CB or name CG or name CD1)'),
('W262_chi1B', 'dihedral', 'segid TRKB and resid 262 and (name N or name CA or name CB or name CG)'),
('W262_chi2B', 'dihedral', 'segid TRKB and resid 262 and (name CA or name CB or name CG or name CD1)'),
('W262_A_CB',  'atom',     'segid TRKA and resid 262 and name CB'),
('W262_A_CH',  'atom',     'segid TRKA and resid 262 and name CH2'),
('W262_B_CB',  'atom',     'segid TRKB and resid 262 and name CB'),
('W262_B_CH',  'atom',     'segid TRKB and resid 262 and name CH2'),
('pore_width', 'distance', 'protein and resid 277 and name CB')
]

# to make selections in the traakTM4 trajectory, numbers must be adjusted
# chains A and B do not have different segids
# for res in chain A, subtract 27 from actual numbers
# for res in chain B, add 255
selections_traakTM4 = [
('S1top',      'COG',      'protein and (resid 105 or resid 387 or resid 214 or resid 496) and name O'),
('S4bottom',   'COG',      'protein and (resid 102 or resid 384 or resid 211 or resid 493) and name OG1'),
('P190_A',     'atom',     'protein and resid 163 and name CA'),
('P190_B',     'atom',     'protein and resid 445 and name CA'),
('trp_gapA',   'distance', 'protein and (resid 90 or resid 232) and name CA'),
('trp_gapB',   'distance', 'protein and (resid 372 or resid 514) and name CA'),
('fenA',       'distance', "protein and (resid 253 or resid 410) and name CA"),
('fenB',       'distance', "protein and (resid 535 or resid 128) and name CA"),
('expA',       'distance', "protein and (resid 142 or resid 251) and name CA"),
('expB',       'distance', "protein and (resid 424 or resid 533) and name CA"),
('W262_chi1A', 'dihedral', "protein and resid 235 and (name N or name CA or name CB or name CG)"),
('W262_chi2A', 'dihedral', "protein and resid 235 and (name CA or name CB or name CG or name CD1)"),
('W262_chi1B', 'dihedral', "protein and resid 517 and (name N or name CA or name CB or name CG)"),
('W262_chi2B', 'dihedral', "protein and resid 517 and (name CA or name CB or name CG or name CD1)"),
('W262_A_CB',  'atom',     'protein and resid 235 and name CB'),
('W262_A_CH',  'atom',     'protein and resid 235 and name CH2'),
('W262_B_CB',  'atom',     'protein and resid 517 and name CB'),
('W262_B_CH',  'atom',     'protein and resid 517 and name CH2'),
('pore_width', 'distance', 'protein and (resid 250 or resid 532) and name CB')
]

os.chdir(outdir)

for pair in filepairs:
    a = ExtendedTSC.ExtendedTSC()
    a.load_traj(pair[0],pair[1])
    a.measures_from_list(selections)
    a.generate_timeseries()
    output = pair[1].split('/')[-1].split('.')[0] + '.' + outtag
    a.write_data(output)

if filepairsTM4:
    for pair in filepairsTM4:
        a = ExtendedTSC.ExtendedTSC()
        a.load_traj(pair[0],pair[1])
        a.measures_from_list(selections_traakTM4)
        a.generate_timeseries()
        output = pair[1].split('/')[-1].split('.')[0] + '.' + outtag
        a.write_data(output)
