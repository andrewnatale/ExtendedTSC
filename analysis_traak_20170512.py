import sys, os, math
import numpy as np
import matplotlib.pyplot as plt
import ExtendedTSC
import traak_plots
#from sklearn import cluster

# generate list of input dat files:
indir = os.path.abspath('measurements05252017')
infiles = []
for infile in os.listdir(indir):
    if infile.endswith('.dat'):
        infiles.append(os.path.join(indir,infile))

datalist = []
for filename in infiles:
    datalist.append(ExtendedTSC.ExtendedTSC(filename))
for elem in datalist:
    elem.simplify_indexing()

save = True
outdirname = 'wingtip_vs_time'

if save:
    plotdir = os.path.join(indir,outdirname)
    # os.mkdir(plotdir)
    os.chdir(plotdir)

# timeseries AB stacks for all traj and all relevant measures
for elem in datalist:
    # calc unit vector through the selectivity filter, SFunitV
    x1 = np.copy(elem.access['S1top'])
    x2 = np.copy(elem.access['S4bottom'])
    SFunitV = (x1-x2) / np.linalg.norm(x1-x2, axis=0)
    # calc wingtip distances, dA and dB
    xA = np.copy(elem.access['P190_A'])
    xB = np.copy(elem.access['P190_B'])
    dA = np.linalg.norm(np.cross(xA-x1, xA-x2, axisa=0,axisb=0,axisc=0), axis=0) / np.linalg.norm(x2-x1, axis=0)
    dB = np.linalg.norm(np.cross(xB-x1, xB-x2, axisa=0,axisb=0,axisc=0), axis=0) / np.linalg.norm(x2-x1, axis=0)
    # calc W262 orientation
    etaA = np.copy(elem.access['W262_A_CH'])
    betaA = np.copy(elem.access['W262_A_CB'])
    etaB = np.copy(elem.access['W262_B_CH'])
    betaB = np.copy(elem.access['W262_B_CB'])
    orientA = np.sum((etaA - betaA) * SFunitV, axis=0)
    orientB = np.sum((etaB - betaB) * SFunitV, axis=0)
    print elem.trajname
    plotname = elem.trajname.split('/')[-1].split('.')[0]
    dataholder1 = [
    # ('M4-M2\' fenestration: P155 to G280\'', elem.access['fenA'][0,:], elem.access['fenB'][0,:]),
    # ('M2-M4 bundle expansion: G169 to T278', elem.access['expA'][0,:] ,elem.access['expB'][0,:]),
    # ('M4-P1 gap: P259 to G117', elem.access['trp_gapA'][0,:], elem.access['trp_gapB'][0,:])
    ('Wingtip distance: pore axis to P190', dA, dB)
    # ('W262 orientation', orientA, orientB)
    ]
    dataholder2 = [
    ('W262 alignment, subunit A', orientA),
    ('W262 alignment, subunit B', orientB)
    ]
    traak_plots.timeseries_AB_stack(plotname, dataholder1, elem.access['time'][0,:]/1000.0, ylimits=[18,38], save=save)
    #traak_plots.timeseries_basic(plotname, dataholder2, elem.access['time'][0,:]/1000.0, ylimits=[-8,8], save=save)
