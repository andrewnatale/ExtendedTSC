import sys, os, math
import numpy as np
import ExtendedTSC
# import mda_plots

# topology and trajectory files
psffile = sys.argv[1]
dcdfile = sys.argv[2]
# stepsize in ps
stepsize = 500

output = sys.argv[3]

# selections for distance measurements
selections = [
('fenestrationA',    'distance',    '(segid TRKA and resid 280) or (segid TRKB and resid 155) and name CA'),
('fenestrationB',    'distance',    '(segid TRKB and resid 280) or (segid TRKA and resid 155) and name CA'),
('expansionA',       'distance',    'segid TRKA and (resid 169 or resid 278) and name CA'),
('expansionB',       'distance',    'segid TRKB and (resid 169 or resid 278) and name CA'),
('W262_chi1A',       'dihedral',    'segid TRKA and resid 262 and (name N or name CA or name CB or name CG)'),
('W262_chi2A',       'dihedral',    'segid TRKA and resid 262 and (name CA or name CB or name CG or name CD1)'),
('W262_chi1B',       'dihedral',    'segid TRKB and resid 262 and (name N or name CA or name CB or name CG)'),
('W262_chi2B',       'dihedral',    'segid TRKB and resid 262 and (name CA or name CB or name CG or name CD1)')
]

# to make selections in the traakTM4 trajectory, numbers must be adjusted
# for res in chain A, subtract 27 from actual numbers
# for res in chain B, add 255
selections_traakTM4 = [
('fenestrationA', 'distance', "protein and (resid 253 or resid 410) and name CA"),
('fenestrationB', 'distance', "protein and (resid 535 or resid 128) and name CA"),
('expansionA', 'distance', "protein and (resid 142 or resid 251) and name CA"),
('expansionB', 'distance', "protein and (resid 424 or resid 533) and name CA"),
('W262_chi1A', 'dihedral', "protein and resid 235 and (name N or name CA or name CB or name CG)"),
('W262_chi2A', 'dihedral', "protein and resid 235 and (name CA or name CB or name CG or name CD1)"),
('W262_chi1B', 'dihedral', "protein and resid 517 and (name N or name CA or name CB or name CG)"),
('W262_chi2B', 'dihedral', "protein and resid 517 and (name CA or name CB or name CG or name CD1)")
]

a = ExtendedTSC.ExtendedTSC()
a.measures_from_dcd(selections,psffile,dcdfile,stepsize)
a.write_dat(output)
