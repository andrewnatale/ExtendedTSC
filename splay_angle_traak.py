import sys, os, math
import numpy as np
import ExtendedTSC
# import mda_plots

# topology and trajectory files
#psffile = sys.argv[1]
#dcdfile = sys.argv[2]
# stepsize in ps
stepsize = 500

# hardcode filenames for batch process
filepairs = [
('/Volumes/data/andrew/lolicato_mutant_traak/traakG124I_full_trimmed/trkG124full.ion.psf',
 '/Volumes/data/andrew/lolicato_mutant_traak/traakG124I_full_trimmed/traakG124I_full_npt.wrap_alignSF.all.500ps.dcd'),
('/Volumes/data/andrew/lolicato_mutant_traak/traakG124I_S1S3_trimmed/trkG124_S1S3.ion.psf',
 '/Volumes/data/andrew/lolicato_mutant_traak/traakG124I_S1S3_trimmed/traakG124I_S1S3_npt.wrap_alignSF.all.500ps.dcd'),
('/Volumes/data/andrew/lolicato_mutant_traak/traakG124I_S2S4_trimmed/trkG124S2S4.ion.psf',
 '/Volumes/data/andrew/lolicato_mutant_traak/traakG124I_S2S4_trimmed/traakG124I_S2S4_npt.wrap_alignSF.all.500ps.dcd'),
('/Volumes/data/andrew/mackinnon_traak/traakWT_full_trimmed/trkWT.popc200.ion.psf',
 '/Volumes/data/andrew/mackinnon_traak/traakWT_full_trimmed/traakWT_full_npt.wrap_alignSF.all.500ps.dcd'),
('/Volumes/data/andrew/mackinnon_traak/traakWT_S1S3_trimmed/trkWTS1S3.ion.psf',
 '/Volumes/data/andrew/mackinnon_traak/traakWT_S1S3_trimmed/traakWT_S1S3_npt.wrap_alignSF.all.500ps.dcd'),
('/Volumes/data/andrew/mackinnon_traak/traakWT_S2S4_trimmed/trkWTS2S4.ion.psf',
 '/Volumes/data/andrew/mackinnon_traak/traakWT_S2S4_trimmed/traakWT_S2S4_npt.wrap_alignSF.all.500ps.dcd')
]

filepairsTM4 = [
('/Volumes/data/andrew/traakTM4_trimmed/traakTM4.popc.ion.psf',
 '/Volumes/data/andrew/traakTM4_trimmed/traakTM4_npt.wrap_alignSF.all.500ps.dcd'),
]

# selections for measurements
selections = [
('s1cog',    'COG',    'protein and (resid 132 or resid 131 or resid 241 or resid 240) and name O'),
('s4cog',    'COG',    'protein and (resid 129 or resid 238) and (name O or name OG1)'),
('W186_A',   'atom',   'segid TRKA and resid 186 and name CA'),
('W186_B',   'atom',   'segid TRKB and resid 186 and name CA'),
('P190_A',   'atom',   'segid TRKA and resid 190 and name CA'),
('P190_B',   'atom',   'segid TRKA and resid 190 and name CA'),
('trp_gapA', 'distance', 'segid TRKA and (resid 117 or resid 259) and name CA'),
('trp_gapB', 'distance', 'segid TRKB and (resid 117 or resid 259) and name CA')
]

# to make selections in the traakTM4 trajectory, numbers must be adjusted
# chains A and B do not have different segids
# for res in chain A, subtract 27 from actual numbers
# for res in chain B, add 255
selections_traakTM4 = [
('s1cog',    'COG',    'protein and (resid 125 or resid 124 or resid 387 or resid 386 or resid 214 or resid 213 or resid 496 or resid 495) and name O'),
('s4cog',    'COG',    'protein and (resid 102 or resid 384 or resid 211 or resid 493) and (name O or name OG1)'),
('W186_A',   'atom',   'protein and resid 159 and name CA'),
('W186_B',   'atom',   'protein and resid 441 and name CA'),
('P190_A',   'atom',   'protein and resid 163 and name CA'),
('P190_B',   'atom',   'protein and resid 445 and name CA'),
('trp_gapA', 'distance', 'protein and (resid 90 or resid 232) and name CA'),
('trp_gapB', 'distance', 'protein and (resid 372 or resid 514) and name CA')
]

os.mkdir('out')

for pair in filepairs:
    a = ExtendedTSC.ExtendedTSC()
    a.measures_from_dcd(selections,pair[0],pair[1],stepsize)
    output = pair[1].split('/')[-1].split('.')[0] + '.splay_meas.dat'
    a.write_dat(os.path.join('out', output))

for pair in filepairsTM4:
    a = ExtendedTSC.ExtendedTSC()
    a.measures_from_dcd(selections_traakTM4,pair[0],pair[1],stepsize)
    output = pair[1].split('/')[-1].split('.')[0] + '.splay_meas.dat'
    a.write_dat(os.path.join('out', output))
