from __future__ import print_function
import sys, os
import numpy as np
from MDAextensions.analysis.MembraneSurface import MembraneSurface

topo = os.path.abspath(sys.argv[1])
traj = os.path.abspath(sys.argv[2])
output_dir = os.path.abspath(sys.argv[3])

basename = ''.join(traj.split('/')[-1].split('.')[:-1])
print(basename)

s = MembraneSurface()
s.load_traj(topo,traj,50,framerange=(1000,-1,1))
s.run(
50,
100,
'resname POPC and name C22 C32',
'resname POPC and name C218 C316',
solute_sel='protein'
)

os.chdir(output_dir)
s.write_aligned_subset(basename)
s.write_vmd_surfaces(basename)
