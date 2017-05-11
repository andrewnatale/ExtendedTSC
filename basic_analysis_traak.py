import sys, os, math
import numpy as np
import ExtendedTSC
import traak_plots
from sklearn import cluster
import matplotlib.pyplot as plt

#np.set_printoptions(threshold=np.nan)

# stepsize in ps
stepsize = 500
save = False
outdirname = None

# generate list of input dat files:
indir = os.path.abspath('../measurements/20170508_basic')
infiles = []
for infile in os.listdir(indir):
    if infile.endswith('.dat'):
        infiles.append(infile)
infiles = [os.path.join(indir, i) for i in infiles]

datalist = []
for filename in infiles:
    datalist.append(ExtendedTSC.ExtendedTSC(filename))
for elem in datalist:
    elem.simplify_indexing()

if save:
    plotdir = os.path.join(indir,outdirname)
    os.mkdir(plotdir)
    os.chdir(plotdir)

# timeseries plots for fenestration and expansion
for elem in datalist:
    print elem.trajname
    plotname = elem.trajname.split('/')[-1].split('.')[0]
    dataholder = [
    ('Fenestration: P155 to G280\'', elem.access['fenestrationA'][0,:], elem.access['fenestrationB'][0,:]),
    ('Expansion: G169 to T278', elem.access['expansionA'][0,:] ,elem.access['expansionB'][0,:])
    ]
    traak_plots.timeseries_AB_stack(plotname, dataholder, elem.access['time'][0,:]/1000.0, ylimits=[4.,20.])



# # selections for distance measurements
# selections = [
# ('fenestrationA',    'distance',    '(segid TRKA and resid 280) or (segid TRKB and resid 155) and name CA'),
# ('fenestrationB',    'distance',    '(segid TRKB and resid 280) or (segid TRKA and resid 155) and name CA'),
# ('expansionA',       'distance',    'segid TRKA and (resid 169 or resid 278) and name CA'),
# ('expansionB',       'distance',    'segid TRKB and (resid 169 or resid 278) and name CA'),
# ('W262_chi1A',       'dihedral',    'segid TRKA and resid 262 and (name N or name CA or name CB or name CG)'),
# ('W262_chi2A',       'dihedral',    'segid TRKA and resid 262 and (name CA or name CB or name CG or name CD1)'),
# ('W262_chi1B',       'dihedral',    'segid TRKB and resid 262 and (name N or name CA or name CB or name CG)'),
# ('W262_chi2B',       'dihedral',    'segid TRKB and resid 262 and (name CA or name CB or name CG or name CD1)')
# ]
#
# # to make selections in the traakTM4 trajectory, numbers must be adjusted
# # for res in chain A, subtract 27 from actual numbers
# # for res in chain B, add 255
# selections_traakTM4 = [
# ('fenestrationA', 'distance', "protein and (resid 253 or resid 410) and name CA"),
# ('fenestrationB', 'distance', "protein and (resid 535 or resid 128) and name CA"),
# ('expansionA', 'distance', "protein and (resid 142 or resid 251) and name CA"),
# ('expansionB', 'distance', "protein and (resid 424 or resid 533) and name CA"),
# ('W262_chi1A', 'dihedral', "protein and resid 235 and (name N or name CA or name CB or name CG)"),
# ('W262_chi2A', 'dihedral', "protein and resid 235 and (name CA or name CB or name CG or name CD1)"),
# ('W262_chi1B', 'dihedral', "protein and resid 517 and (name N or name CA or name CB or name CG)"),
# ('W262_chi2B', 'dihedral', "protein and resid 517 and (name CA or name CB or name CG or name CD1)")
# ]
