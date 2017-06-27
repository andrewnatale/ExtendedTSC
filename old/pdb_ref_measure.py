import sys, os, math
import numpy as np
from ExtendedTSC import ExtendedTSC

input_prefix = '/Users/anatale/school/UCSF/Grabe_Lab/data'

outtag = 'pdb_ref'
# make output dir in cwd
outdir = os.path.join(os.getcwd(), 'pdb_ref')
try:
    os.mkdir(outdir)
except OSError:
    if not os.path.isdir(outdir):
        raise

traak_measures = [
('S1top',      'COG',      'protein and (resid 132 or resid 241) and name O'),
('S4bottom',   'COG',      'protein and (resid 129 or resid 238) and name OG1'),
('m3tip_A',     'atom',     'segid A and resid 190 and name CA'),
('m3tip_B',     'atom',     'segid B and resid 190 and name CA'),
('trp_gapA',   'distance', 'segid A and (resid 117 or resid 259) and name CA'),
('trp_gapB',   'distance', 'segid B and (resid 117 or resid 259) and name CA'),
('fenA',       'distance', '(segid A and resid 280) or (segid B and resid 155) and name CA'),
('fenB',       'distance', '(segid B and resid 280) or (segid A and resid 155) and name CA'),
('expA',       'distance', 'segid A and (resid 169 or resid 278) and name CA'),
('expB',       'distance', 'segid B and (resid 169 or resid 278) and name CA'),
('m4trp_chi1A', 'dihedral', 'segid A and resid 262 and (name N or name CA or name CB or name CG)'),
('m4trp_chi2A', 'dihedral', 'segid B and resid 262 and (name CA or name CB or name CG or name CD1)'),
('m4trp_chi1B', 'dihedral', 'segid A and resid 262 and (name N or name CA or name CB or name CG)'),
('m4trp_chi2B', 'dihedral', 'segid B and resid 262 and (name CA or name CB or name CG or name CD1)'),
('m4trp_A_CB',  'atom',     'segid A and resid 262 and name CB'),
('m4trp_A_CH',  'atom',     'segid A and resid 262 and name CH2'),
('m4trp_B_CB',  'atom',     'segid B and resid 262 and name CB'),
('m4trp_B_CH',  'atom',     'segid B and resid 262 and name CH2'),
('pore_width', 'distance', 'protein and resid 277 and name CB')
]

trek2_measuresAB = [
('S1top',      'COG',      'protein and (resid 175 or resid 284) and name O'),
('S4bottom',   'COG',      'protein and (resid 172 or resid 281) and name OG1'),
('m3tip_A',     'atom',     'segid A and resid 233 and name CA'),
('m3tip_B',     'atom',     'segid B and resid 233 and name CA'),
('trp_gapA',   'distance', 'segid A and (resid 160 or resid 303) and name CA'),
('trp_gapB',   'distance', 'segid B and (resid 160 or resid 303) and name CA'),
('fenA',       'distance', '(segid A and resid 324) or (segid B and resid 198) and name CA'),
('fenB',       'distance', '(segid B and resid 324) or (segid A and resid 198) and name CA'),
('expA',       'distance', 'segid A and (resid 212 or resid 322) and name CA'),
('expB',       'distance', 'segid B and (resid 212 or resid 322) and name CA'),
('m4trp_chi1A', 'dihedral', 'segid A and resid 306 and (name N or name CA or name CB or name CG)'),
('m4trp_chi2A', 'dihedral', 'segid A and resid 306 and (name CA or name CB or name CG or name CD1)'),
('m4trp_chi1B', 'dihedral', 'segid B and resid 306 and (name N or name CA or name CB or name CG)'),
('m4trp_chi2B', 'dihedral', 'segid B and resid 306 and (name CA or name CB or name CG or name CD1)'),
('m4trp_A_CB',  'atom',     'segid A and resid 306 and name CB'),
('m4trp_A_CH',  'atom',     'segid A and resid 306 and name CH2'),
('m4trp_B_CB',  'atom',     'segid B and resid 306 and name CB'),
('m4trp_B_CH',  'atom',     'segid B and resid 306 and name CH2'),
('pore_width', 'distance', 'protein and resid 321 and name CB')
]

trek2_measuresCD = [
('S1top',      'COG',      'protein and (resid 175 or resid 284) and name O'),
('S4bottom',   'COG',      'protein and (resid 172 or resid 281) and name OG1'),
('m3tip_A',     'atom',     'segid C and resid 233 and name CA'),
('m3tip_B',     'atom',     'segid D and resid 233 and name CA'),
('trp_gapA',   'distance', 'segid C and (resid 160 or resid 303) and name CA'),
('trp_gapB',   'distance', 'segid D and (resid 160 or resid 303) and name CA'),
('fenA',       'distance', '(segid C and resid 324) or (segid D and resid 198) and name CA'),
('fenB',       'distance', '(segid D and resid 324) or (segid C and resid 198) and name CA'),
('expA',       'distance', 'segid C and (resid 212 or resid 322) and name CA'),
('expB',       'distance', 'segid D and (resid 212 or resid 322) and name CA'),
('m4trp_chi1A', 'dihedral', 'segid C and resid 306 and (name N or name CA or name CB or name CG)'),
('m4trp_chi2A', 'dihedral', 'segid C and resid 306 and (name CA or name CB or name CG or name CD1)'),
('m4trp_chi1B', 'dihedral', 'segid D and resid 306 and (name N or name CA or name CB or name CG)'),
('m4trp_chi2B', 'dihedral', 'segid D and resid 306 and (name CA or name CB or name CG or name CD1)'),
('m4trp_A_CB',  'atom',     'segid C and resid 306 and name CB'),
('m4trp_A_CH',  'atom',     'segid C and resid 306 and name CH2'),
('m4trp_B_CB',  'atom',     'segid D and resid 306 and name CB'),
('m4trp_B_CH',  'atom',     'segid D and resid 306 and name CH2'),
('pore_width', 'distance', 'protein and resid 321 and name CB')
]

# each pdb is missing different bits, so this is messy
meas_dict = {
os.path.join(input_prefix, 'pdbs/TRAAKwt_4i9w.pdb'): traak_measures,
os.path.join(input_prefix, 'pdbs/G124I_4rue.pdb'): traak_measures,
os.path.join(input_prefix, 'pdbs/W262S_4ruf.pdb'): traak_measures[0:9]+traak_measures[18:],
os.path.join(input_prefix, 'pdbs/TREK2_down_4xdj_AB.pdb'): trek2_measuresAB,
os.path.join(input_prefix, 'pdbs/TREK2_down_4xdj_CD.pdb'): trek2_measuresCD[0:2]+trek2_measuresCD[3:],
os.path.join(input_prefix, 'pdbs/TREK2_up_4bw5_AB.pdb'): trek2_measuresAB[0:2]+trek2_measuresAB[3:12]+trek2_measuresAB[14:16]+trek2_measuresAB[18:],
os.path.join(input_prefix, 'pdbs/TREK2_up_4bw5_CD.pdb'): trek2_measuresCD[0:2]+trek2_measuresCD[3:12]+trek2_measuresCD[14:16]+trek2_measuresCD[18:]
}

os.chdir(outdir)

for key in meas_dict:
    a = ExtendedTSC()
    a.load_pdb(key)
    a.measures_from_list(meas_dict[key])
    a.run()
    output = key.split('/')[-1].split('.')[0] + '.' + outtag
    a.write_data(output)
