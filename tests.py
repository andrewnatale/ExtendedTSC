import sys, os, math, re
import numpy as np
from ExtendedTSC import ExtendedTSC

a = ExtendedTSC()
a.load_pdb('/Users/anatale/school/UCSF/Grabe_Lab/data/pdbs/G124I_4rue.pdb')
test_sel = [
('S1top',      'COG',      'protein and (resid 132 or resid 241) and name O'),
('S4bottom',   'COG',      'protein and (resid 129 or resid 238) and name OG1')
]
a.measures_from_list(test_sel)
a.run()
a._data_writer(a.primaryDS)
