import sys, os, math
import numpy as np
import MDAnalysis as md
from MDAnalysis.analysis.density import density_from_Universe

script,psf,dcd = sys.argv

u = md.Universe(psf,dcd)

def taildensity():
    D = density_from_Universe(u, delta=1.0, atomselection="\
      name C22 or name C23 or name C24 or name C25 or name C26 or\
      name C27 or name C28 or name C29 or name C210 or name C211 or\
      name C212 or name C213 or name C214 or name C215 or name C216 or\
      name C217 or name C218 or\
      name C32 or name C33 or name C34 or name C35 or name C36 or\
      name C37 or name C38 or name C39 or name C310 or name C311 or\
      name C312 or name C313 or name C314 or name C315 or name C316")
    #print D.grid
    print np.shape(D.grid)
    print D.grid.sum()
    D.grid = D.grid / D.grid.sum()
    D.export('taildensity.dx')

def phosdensity():
    D = density_from_Universe(u, delta=1.0, atomselection="\
      resname POPC and name P")
    #print D.grid
    print np.shape(D.grid)
    print D.grid.sum()
    D.grid = D.grid / D.grid.sum()
    D.export('phosdensity.dx')

taildensity()
phosdensity()
