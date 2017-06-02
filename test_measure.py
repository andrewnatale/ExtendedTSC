import sys, os, math
import numpy as np
#import ExtendedTSC
import MDAnalysis as md
import MDAnalysis.core.Timeseries as tm

pdbfile = os.path.abspath(sys.argv[1])
u = md.Universe(pdbfile)
print u
a = u.select_atoms('segid A and resid 28')
print a.n_atoms
collection = tm.TimeseriesCollection()
collection.addTimeseries(tm.CenterOfGeometry(a))
collection.compute(u.trajectory)
