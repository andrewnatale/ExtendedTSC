import sys, os, math
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.density import density_from_Universe

padding = 25.0
delta = 2.0

input_prefix = '/Users/anatale/UCSF/Grabe_lab/scratch/density_test'

filelist = [
('topo.traakWT_full_npt.psf', 'wtfull1.1ns.dcd'),
('topo.traakWT_full_npt.psf', 'wtfull2.1ns.dcd'),
('topo.traakWT_s1s3_npt.psf', 'wts1s3.1ns.dcd'),
]

def taildensity(universe, outdir, outtag):
    D = density_from_Universe(universe, delta=delta, padding=padding, atomselection="\
      name C22 or name C23 or name C24 or name C25 or name C26 or\
      name C27 or name C28 or name C29 or name C210 or name C211 or\
      name C212 or name C213 or name C214 or name C215 or name C216 or\
      name C217 or name C218 or\
      name C32 or name C33 or name C34 or name C35 or name C36 or\
      name C37 or name C38 or name C39 or name C310 or name C311 or\
      name C312 or name C313 or name C314 or name C315 or name C316")
    #print D.grid
    print np.shape(D.grid)
    print D.parameters
    print D.edges
    return D.grid
    #D.grid = D.grid / np.amax(D.grid)

    #D.export(os.path.join(outdir,'new_taildensity_%s.dx' % outtag))

def phosdensity(universe, outdir, outtag):
    D = density_from_Universe(universe, delta=delta, padding=padding, atomselection="\
      resname POPC and name P")
    #print D.grid
    print np.shape(D.grid)
    print D.parameters
    return D.grid
    #D.grid = D.grid / np.amax(D.grid)

    #D.export(os.path.join(outdir,'new_phosdensity_%s.dx' % outtag))

stuff = []
for filepair in filelist:
    u = mda.Universe(os.path.join(input_prefix,filepair[0]), os.path.join(input_prefix,filepair[1]))
    stuff.append(taildensity(u,None,None))

for elem in stuff:
    print np.shape(elem)
