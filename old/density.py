import sys, os, math
import numpy as np
import MDAnalysis as md
from MDAnalysis.analysis.density import density_from_Universe

padding = 25.0
delta = 2.0

# topology and trajectory files
try:
    psffile = sys.argv[1]
    dcdfile = sys.argv[2]
except IndexError:
    psffile = None
    dcdfile = None

input_prefix = '/Users/anatale/school/UCSF/Grabe_lab/data'

outtag = 'bin'+str(delta)

# if I'm not batch processing, I'm probably just testing on a trajectory with the non TM4 topology
if psffile and dcdfile:
    print 'Using inputs specified from command line'
    filepairs = [(os.path.abspath(psffile), os.path.abspath(dcdfile))]
    filepairsTM4 = []
else:
    # hardcode filenames for batch process
    filepairs = [
    (os.path.join(input_prefix, 'traakG124I_full_trimmed/topo.traakG124I_full_npt.psf'),
     os.path.join(input_prefix, 'traakG124I_full_trimmed/traakG124I_full_npt.all.200ps.filter_alignZ.dcd')),
    (os.path.join(input_prefix, 'traakWT_full_trimmed/topo.traakWT_full_npt.psf'),
     os.path.join(input_prefix, 'traakWT_full_trimmed/traakWT_full_npt.all.200ps.filter_alignZ.dcd')),
    ]
    filepairsTM4 = [
    (os.path.join(input_prefix, 'traakTM4_trimmed/topo.traakTM4_npt.psf'),
     os.path.join(input_prefix, 'traakTM4_trimmed/traakTM4_npt.all.200ps.filter_alignZ.dcd')),
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
    #D.grid = D.grid / np.amax(D.grid)

    D.export(os.path.join(outdir,'new_taildensity_%s.dx' % outtag))

def phosdensity(universe, outdir, outtag):
    D = density_from_Universe(universe, delta=delta, padding=padding, atomselection="\
      resname POPC and name P")
    #print D.grid
    print np.shape(D.grid)
    #D.grid = D.grid / np.amax(D.grid)

    D.export(os.path.join(outdir,'new_phosdensity_%s.dx' % outtag))

for pair in filepairs+filepairsTM4:
    u = md.Universe(pair[0],pair[1])
    outdir = os.path.dirname(pair[0])
    taildensity(u,outdir,outtag)
    phosdensity(u,outdir,outtag)
