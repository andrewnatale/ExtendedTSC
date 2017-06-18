import sys, os, math, re
import numpy as np
import ExtendedTSC
import MDAnalysis as md

psf = '/Users/anatale/school/UCSF/Grabe_Lab/data/traakWT_full_trimmed/topo.traakWT_full_npt.psf'
dcd = '/Users/anatale/school/UCSF/Grabe_Lab/data/traakWT_full_trimmed/traakWT_full.new_align.500ps.dcd'
pdb1 = '/Users/anatale/school/UCSF/Grabe_Lab/TRAAKwt_4i9w.pdb'
pdb2 = '/Users/anatale/school/UCSF/Grabe_Lab/G124I_4rue.pdb'

test = ExtendedTSC.ExtendedTSC()

# test.load_dcd(psf,dcd)
test.load_pdb(pdb1)
test.measures_from_list([('test', 'COM', 'protein and (name CA or name N)'),
                         ('test2', 'distance', 'segid A and (resid 117 or resid 259) and name CA'),
                         ('test3', 'dihedral', 'segid A and resid 262 and (name N or name CA or name CB or name CG)')])
test.run()
print test.primaryDS.data
test.write_data('testA')
