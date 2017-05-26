import sys, os, math
import numpy as np
import ExtendedTSC
import traak_plots
#from sklearn import cluster
import matplotlib.pyplot as plt

#np.set_printoptions(threshold=np.nan)

# stepsize in ps
stepsize = 500
save = False
outdirname = None

# generate list of input dat files:
indir = os.path.abspath('../../measurements/20170512')
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

# timeseries AB stacks for all traj
for elem in datalist:
    # get wingtip distances
    x1 = np.copy(elem.access['s1cog'])
    x2 = np.copy(elem.access['s4cog'])
    xA = np.copy(elem.access['P190_A'])
    xB = np.copy(elem.access['P190_B'])
    dA = np.linalg.norm(np.cross(xA-x1, xA-x2, axisa=0,axisb=0,axisc=0), axis=0) / np.linalg.norm(x1-x2, axis=0)
    dB = np.linalg.norm(np.cross(xB-x1, xB-x2, axisa=0,axisb=0,axisc=0), axis=0) / np.linalg.norm(x1-x2, axis=0)
    print elem.trajname
    plotname = elem.trajname.split('/')[-1].split('.')[0]
    dataholder = [
    ('Fenestration: P155 to G280\'', elem.access['fenA'][0,:], elem.access['fenB'][0,:]),
    #('Expansion: G169 to T278', elem.access['expA'][0,:] ,elem.access['expB'][0,:]),
    # ('M4-P1 gap: P259 to G117', elem.access['trp_gapA'][0,:], elem.access['trp_gapB'][0,:]),
    ('Wingtip distance: center axis to P190', dA, dB)
    ]
    traak_plots.timeseries_AB_stack(plotname, dataholder, elem.access['time'][0,:]/1000.0, ylimits=[4.,36.])
