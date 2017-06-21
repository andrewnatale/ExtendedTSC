import sys, os, math, re
import numpy as np
from itertools import cycle
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import ExtendedTSC_old

# choose report to plot
rep_idx = int(sys.argv[1])
reports = ['traakWT_full','traakWT_S1S3','traakWT_S2S4','traakG124I_full','traakG124I_S1S3','traakG124I_S2S4','traakTM4']
report_name = reports[rep_idx]
# load data sets
a = ExtendedTSC_old.ExtendedTSC(datfile=os.path.join('basic_meas', report_name+'.basic.dat'))
b = ExtendedTSC_old.ExtendedTSC(datfile=os.path.join('core_volume_tracking', report_name+'.core_vol.dat'),maskfile=os.path.join('core_volume_tracking', report_name+'.core_vol.mask.dat'))
bw = ExtendedTSC_old.ExtendedTSC(datfile=os.path.join('core_volume_tracking', report_name+'.core_vol.tip3.dat'))
# make data sets indexable
basic = a.primaryDS.simplify_indexing(return_dict=True)
klipid = b.primaryDS.simplify_indexing(return_dict=True)
mask = b.maskDS.simplify_indexing(return_dict=True)
waterZ = bw.primaryDS.simplify_indexing(return_dict=True)

# choose reference lines
use_ref = ['TRAAK','TREK2']
try:
    ref_idx = int(sys.argv[2])
    ref_name = use_ref[ref_idx]
except:
    ref_name = None
# load pdb reference measurements
references = ['TRAAKwt_4i9w','G124I_4rue','W262S_4ruf','TREK2_down_4xdj_AB','TREK2_down_4xdj_CD','TREK2_up_4bw5_AB','TREK2_up_4bw5_CD']
ref_data = []
for ref in references:
    ref_data.append(ExtendedTSC_old.ExtendedTSC(datfile=os.path.join('pdb_ref',ref+'.pdb_ref.dat')))

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

def timeseries(tgt_ax,datalist,time,datalabels=None,references=None,reference_labels=None,ylabel=None,xlabel='time (ns)',plotname=None,ylimits=None,multicolor=True):
    # arbitrary color choices
    if multicolor:
        colors = ['red', 'blue', 'green', 'purple']
    else:
        colors = ['black','black','black','black']
    ref_colors = ['silver', 'gold']
    # setup timeseries data
    for idx,elem in enumerate(datalist):
        tgt_ax.plot(time/1000.0, elem, c=colors[idx])
        if datalabels:
            tgt_ax.text(0.99,0.99-0.05*idx, datalabels[idx], color=colors[idx],
                        transform=tgt_ax.transAxes,verticalalignment='top',horizontalalignment='right', fontsize=8)
    # optionally add reference lines
    if references:
        for idx,elem in enumerate(references):
            tgt_ax.axhline(y=elem,c=ref_colors[idx],linewidth=1.0,linestyle='dashed',zorder=0)
            if reference_labels:
                tgt_ax.text(0.01,0.99-0.05*idx, reference_labels[idx], color=ref_colors[idx],
                            transform=tgt_ax.transAxes,verticalalignment='top',horizontalalignment='left', fontsize=8)
    tgt_ax.set_xlim([ time[0]/1000.0, time[-1]/1000.0 ])
    # set more optional params and add labels
    if ylimits:
        tgt_ax.set_ylim(ylimits)
    if plotname:
        tgt_ax.set_title(plotname)
    if xlabel:
        tgt_ax.set_xlabel(xlabel)
    if ylabel:
        tgt_ax.set_ylabel(ylabel)

def process_refs(data_dict):
    # pre-calculated measurements that every structure should have, get vals out of arrays
    pore_width = data_dict['pore_width'][0,0]
    fenA = data_dict['fenA'][0,0]
    fenB = data_dict['fenB'][0,0]
    trp_gapA = data_dict['trp_gapA'][0,0]
    trp_gapB = data_dict['trp_gapB'][0,0]
    expA = data_dict['expA'][0,0]
    expB = data_dict['expB'][0,0]
    # calculations on coordinates and handle missing data
    # calc unit vector through the selectivity filter, SFunitV
    zref1 = np.copy(data_dict['S1top'])
    zref2 = np.copy(data_dict['S4bottom'])
    SFunitV = (zref1-zref2) / np.linalg.norm(zref1-zref2, axis=0)
    # calc wingtip distances, dA and dB
    try:
        xA = np.copy(data_dict['m3tip_A'])
        wing_distA = np.linalg.norm(np.cross(xA-zref1, xA-zref2, axisa=0,axisb=0,axisc=0), axis=0) / np.linalg.norm(zref2-zref1, axis=0)
        wingA = wing_distA[0]
    except KeyError:
        wingA = None
    try:
        xB = np.copy(data_dict['m3tip_B'])
        wing_distB = np.linalg.norm(np.cross(xB-zref1, xB-zref2, axisa=0,axisb=0,axisc=0), axis=0) / np.linalg.norm(zref2-zref1, axis=0)
        wingB = wing_distB[0]
    except KeyError:
        wingB = None
    # calc W262 orientation, relative to SFunitV
    try:
        etaA = np.copy(data_dict['m4trp_A_CH'])
        betaA = np.copy(data_dict['m4trp_A_CB'])
        orientA = np.sum((etaA - betaA) * SFunitV, axis=0)
        trp_oriA = orientA[0]
    except KeyError:
        trp_oriA = None
    try:
        etaB = np.copy(data_dict['m4trp_B_CH'])
        betaB = np.copy(data_dict['m4trp_B_CB'])
        orientB = np.sum((etaB - betaB) * SFunitV, axis=0)
        trp_oriB = orientB[0]
    except KeyError:
        trp_oriB = None
    # idx  0          1    2    3        4        5    6    7     8     9        10
    return pore_width,fenA,fenB,trp_gapA,trp_gapB,expA,expB,wingA,wingB,trp_oriA,trp_oriB

# setup reference points
traakWT_ref = process_refs(ref_data[0].primaryDS.simplify_indexing())
downAB = process_refs(ref_data[3].primaryDS.simplify_indexing())
downCD = process_refs(ref_data[4].primaryDS.simplify_indexing())
upAB = process_refs(ref_data[5].primaryDS.simplify_indexing())
upCD = process_refs(ref_data[6].primaryDS.simplify_indexing())

if ref_name == 'TRAAK':
    ref_labels = ('TRAAKwt_A', 'TRAAKwt_B')
    ref_pw = [traakWT_ref[0],]
    ref_pw_labels = ('TRAAKwt_4i9w',)
    ref_wing = [traakWT_ref[7], traakWT_ref[8]]
    ref_wing_labels = ref_labels
    ref_fen = [traakWT_ref[1], traakWT_ref[2]]
    ref_fen_labels = ref_labels
    ref_trp_ori = [traakWT_ref[9], traakWT_ref[10]]
    ref_trp_ori_labels = ref_labels
elif ref_name == 'TREK2':
    ref_labels = ('TREK2_down','TREK2_up')
    ref_pw = [(downAB[0]+downCD[0])/2.0, (upAB[0]+upCD[0])/2.0]
    ref_pw_labels = ref_labels
    ref_wing = [(downAB[7]+downAB[8]+downCD[8])/3.0, (upAB[8]+upCD[8])/2.0]
    ref_wing_labels = ref_labels
    ref_fen = [(downAB[1]+downAB[2]+downCD[2]+downCD[2])/4.0, (upAB[1]+upAB[2]+upCD[1]+upCD[2])/4.0]
    ref_fen_labels = ref_labels
    ref_trp_ori = [(downAB[9]+downAB[10]+downCD[9]+downCD[10])/4.0]
    ref_trp_ori_labels = ('TREK2_down_avg',)
else:
    ref_labels = None
    ref_pw = None
    ref_pw_labels = ref_labels
    ref_wing = None
    ref_wing_labels = ref_labels
    ref_fen = None
    ref_fen_labels = ref_labels
    ref_trp_ori = None
    ref_trp_ori_labels = None

# setup axes
pore_vol_occupancy(ax1,plotname='Filter and pore occupancy',xlabel=None)

timeseries(
           ax2,
           datalist=[basic['pore_width'][0,:]],
           time=basic['time'][0,:],
           datalabels=None,
           references=ref_pw,
           reference_labels=ref_pw_labels,
           ylabel='Angstroms',
           plotname='Pore width at T277 (on M4)',
           ylimits=[5,20],
           multicolor=False
           )

timeseries(
           ax3,
           datalist=[wingA,wingB],
           time=basic['time'][0,:],
           datalabels=['subunit A','subunit B'],
           references=ref_wing,
           reference_labels=ref_wing_labels,
           ylabel='Angstroms',
           xlabel=None,
           plotname='M2-M3 wing extension (at P190) from pore axis',
           ylimits=[18,40]
           )

timeseries(
           ax4,
           datalist=[basic['fenA'][0,:],basic['fenB'][0,:]],
           time=basic['time'][0,:],
           datalabels=['subunit A','subunit B'],
           references=ref_fen,
           reference_labels=ref_fen_labels,
           ylabel='Angstroms',
           xlabel=None,
           plotname='M4-M2\' fenestration opening',
           ylimits=[4,16]
           )

timeseries(
           ax5,
           datalist=[orientA,orientB],
           time=basic['time'][0,:],
           datalabels=['subunit A','subunit B'],
           references=ref_trp_ori,
           reference_labels=ref_trp_ori_labels,
           ylabel='C-eta relative z-position',
           plotname='W262 indole ring orientation',
           ylimits=[-7,7]
           )

plt.suptitle(report_name, fontsize=20)
plt.show()
