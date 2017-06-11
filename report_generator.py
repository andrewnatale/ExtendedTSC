import sys, os, math, re
import numpy as np
from itertools import cycle
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import ExtendedTSC

rep_idx = int(sys.argv[1])

reports = ['traakWT_full','traakWT_S1S3','traakWT_S2S4','traakG124I_full','traakG124I_S1S3','traakG124I_S2S4','traakTM4']

report_name = reports[rep_idx]

# load data sets
a = ExtendedTSC.ExtendedTSC(datfile=os.path.join('basic_meas', report_name+'.basic.dat'))
b = ExtendedTSC.ExtendedTSC(datfile=os.path.join('core_volume_tracking', report_name+'.core_vol.dat'),maskfile=os.path.join('core_volume_tracking', report_name+'.core_vol.mask.dat'))
bw = ExtendedTSC.ExtendedTSC(datfile=os.path.join('core_volume_tracking', report_name+'.core_vol.tip3.dat'))

# make data sets indexable
basic = a.primaryDS.simplify_indexing(return_dict=True)
klipid = b.primaryDS.simplify_indexing(return_dict=True)
mask = b.maskDS.simplify_indexing(return_dict=True)
waterZ = bw.primaryDS.simplify_indexing(return_dict=True)

# do some vector math to setup
# calc unit vector through the selectivity filter, SFunitV
zref1 = np.copy(basic['S1top'])
zref2 = np.copy(basic['S4bottom'])
SFunitV = (zref1-zref2) / np.linalg.norm(zref1-zref2, axis=0)
# calc wingtip distances, dA and dB
xA = np.copy(basic['P190_A'])
xB = np.copy(basic['P190_B'])
wingA = np.linalg.norm(np.cross(xA-zref1, xA-zref2, axisa=0,axisb=0,axisc=0), axis=0) / np.linalg.norm(zref2-zref1, axis=0)
wingB = np.linalg.norm(np.cross(xB-zref1, xB-zref2, axisa=0,axisb=0,axisc=0), axis=0) / np.linalg.norm(zref2-zref1, axis=0)
# calc W262 orientation, relative to SFunitV
etaA = np.copy(basic['W262_A_CH'])
betaA = np.copy(basic['W262_A_CB'])
etaB = np.copy(basic['W262_B_CH'])
betaB = np.copy(basic['W262_B_CB'])
orientA = np.sum((etaA - betaA) * SFunitV, axis=0)
orientB = np.sum((etaB - betaB) * SFunitV, axis=0)

# setup plot grid
plt.figure()
gs = GridSpec(3,2)
gs.update(hspace=0.3)
ax1 = plt.subplot(gs[0:2,0])
ax2 = plt.subplot(gs[2,0])
ax3 = plt.subplot(gs[0,1])
ax4 = plt.subplot(gs[1,1])
ax5 = plt.subplot(gs[2,1])
# f,ax1 = plt.subplots()

# plotting functions to populate axes objects
def pore_vol_occupancy(tgt_ax,plotname=None,xlabel='time (ns)'):
    # setup colors for different entities
    Kcolors = ['violet','limegreen','silver','tomato','seagreen','blue','blueviolet','palevioletred']
    Kcolor = cycle(Kcolors).next
    watercolor = 'powderblue'
    lipidcolor = 'gold'
    # regex for sorting data
    find_K = re.compile('POT')
    find_POPC = re.compile('POPC')
    # reference points z=0
    zref = np.reshape(np.copy(klipid['S4bottom'][2,:]), (1,-1))
    # begin plotting to tgt_ax
    for key in klipid:
        # potassium ions on top, each with a unique(ish) color
        if find_K.search(key):
            zcoords = np.reshape(np.copy(klipid[key][2,:]), (1,-1))
            zcoords = zcoords - zref
            boolmask = mask[key].astype(bool)
            masked_zcoords = zcoords[boolmask]
            masked_time = klipid['time'][boolmask] / 1000.0
            tgt_ax.plot(masked_time, masked_zcoords, c=Kcolor(), zorder=3, marker='o', mew=0.0, linestyle="None")
        # carbon atoms of lipid tails, all the same color an below potassium
        elif find_POPC.search(key):
            zcoords = np.reshape(np.copy(klipid[key][2,:]), (1,-1))
            zcoords = zcoords - zref
            boolmask = mask[key].astype(bool)
            masked_zcoords = zcoords[boolmask]
            masked_time = klipid['time'][boolmask] / 1000.0
            tgt_ax.plot(masked_time, masked_zcoords, c=lipidcolor, zorder=2, marker='o', mew=0.0, ms=3, linestyle="None")
        else:
            pass
    # water oxygen atoms, all the same color and below all other points
    for idx,time in np.ndenumerate(waterZ['time']):
        waters = waterZ['water_volume'][:,idx[1]]
        waters = waters[~np.isnan(waters)]
        waters = waters - zref[idx]
        timeblob = np.empty_like(waters)
        timeblob.fill(time/1000.0)
        tgt_ax.plot(timeblob, waters, c=watercolor, zorder=1, marker='o', mew=0.0, ms=3, linestyle="None")
    # limits and labels
    tgt_ax.set_xlim([klipid['time'][0,0]/1000.0, klipid['time'][0,-1]/1000.0])
    tgt_ax.set_ylim([-16,16])
    if xlabel:
        tgt_ax.set_xlabel(xlabel)
    tgt_ax.set_ylabel('Location along pore axis (angstroms)')
    if plotname:
        tgt_ax.set_title(plotname)

def timeseries(tgt_ax,datalist,time,datalabels=None,ylabel=None,xlabel='time (ns)',plotname=None,ylimits=None,multicolor=True):
    if multicolor:
        colors = ['red', 'blue', 'green', 'purple', 'silver', 'gold']
    else:
        colors = ['black','black','black','black','black','black',]
    for idx,elem in enumerate(datalist):
        tgt_ax.plot(time/1000.0, elem, c=colors[idx])
        if datalabels:
            tgt_ax.text(0.99,0.99-0.08*idx, datalabels[idx], color=colors[idx],
                        transform=tgt_ax.transAxes,verticalalignment='top',horizontalalignment='right')
    tgt_ax.set_xlim([ time[0]/1000.0, time[-1]/1000.0 ])
    if ylimits:
        tgt_ax.set_ylim(ylimits)
    if plotname:
        tgt_ax.set_title(plotname)
    if xlabel:
        tgt_ax.set_xlabel(xlabel)
    if ylabel:
        tgt_ax.set_ylabel(ylabel)

pore_vol_occupancy(ax1,plotname='Filter and pore occupancy',xlabel=None)

timeseries(ax2,
           datalist=[basic['pore_width'][0,:]],
           time=basic['time'][0,:],
           datalabels=None,
           ylabel='Angstroms',
           plotname='Pore width at T277 (on M4)',
           ylimits=[5,20],
           multicolor=False
           )

timeseries(ax3,
           datalist=[wingA,wingB],
           time=basic['time'][0,:],
           datalabels=['subunit A','subunit B'],
           ylabel='Angstroms',
           xlabel=None,
           plotname='M2-M3 wing extension (at P190) from pore axis',
           ylimits=[20,40]
           )

timeseries(ax4,
           datalist=[basic['fenA'][0,:],basic['fenB'][0,:]],
           time=basic['time'][0,:],
           datalabels=['subunit A','subunit B'],
           ylabel='Angstroms',
           plotname='M4-M2\' fenestration opening',
           ylimits=[4,16]
           )

timeseries(ax5,
           datalist=[orientA,orientB],
           time=basic['time'][0,:],
           datalabels=['subunit A','subunit B'],
           ylabel='C-eta relative z-position',
           xlabel=None,
           plotname='W262 rotamer orientation',
           ylimits=[-7,7]
           )

plt.suptitle(report_name, fontsize=20)
plt.show()
