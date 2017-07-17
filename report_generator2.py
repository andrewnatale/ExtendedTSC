import sys, os, math, re
import numpy as np
from itertools import cycle
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
#import ExtendedTSC_old
from core.TimeseriesCore import TimeseriesCore

reports = ['traakWT_full','traakWT_S1S3','traakWT_S2S4','traakG124I_full','traakG124I_S1S3','traakG124I_S2S4','traakTM4']

dat_root_dir = '/Users/anatale/school/UCSF/Grabe_lab/data/features_traak/'
tgt_data = '500ps'

# load data sets
a = TimeseriesCore(datfilename=os.path.join(dat_root_dir, tgt_data, 'basic_meas', 'traakWT_S1S3.basic.dat'))
b = TimeseriesCore(datfilename=os.path.join(dat_root_dir, tgt_data, 'basic_meas', 'traakG124I_S1S3.basic.dat'))
c = TimeseriesCore(datfilename=os.path.join(dat_root_dir, tgt_data, 'basic_meas', 'traakTM4.basic.dat'))

# make data sets indexable
wt_basic = a.primaryDS.dataset_to_dict()
mut_basic = b.primaryDS.dataset_to_dict()
tm4_basic = c.primaryDS.dataset_to_dict()

# choose reference lines
use_ref = ['TRAAK','TREK2']
try:
    #ref_idx = int(sys.argv[2])
    ref_idx = 0
    ref_name = use_ref[ref_idx]
except:
    ref_name = None
# load pdb reference measurements
references = ['TRAAKwt_4i9w','G124I_4rue','W262S_4ruf','TREK2_down_4xdj_AB','TREK2_down_4xdj_CD','TREK2_up_4bw5_AB','TREK2_up_4bw5_CD']
ref_data = []
for ref in references:
    ref_data.append(TimeseriesCore(datfilename=os.path.join(dat_root_dir, 'pdb_ref', ref+'.pdb_ref.dat')))

def proc_dataset(target_dict):
    # do some vector math to setup
    # calc unit vector through the selectivity filter, SFunitV
    zref1 = np.copy(target_dict['S1top'])
    zref2 = np.copy(target_dict['S4bottom'])
    SFunitV = (zref1-zref2) / np.linalg.norm(zref1-zref2, axis=0)
    # calc wingtip distances, dA and dB
    xA = np.copy(target_dict['P190_A'])
    xB = np.copy(target_dict['P190_B'])
    wingA = np.linalg.norm(np.cross(xA-zref1, xA-zref2, axisa=0,axisb=0,axisc=0), axis=0) / np.linalg.norm(zref2-zref1, axis=0)
    wingB = np.linalg.norm(np.cross(xB-zref1, xB-zref2, axisa=0,axisb=0,axisc=0), axis=0) / np.linalg.norm(zref2-zref1, axis=0)
    # calc W262 orientation, relative to SFunitV
    etaA = np.copy(target_dict['W262_A_CH'])
    betaA = np.copy(target_dict['W262_A_CB'])
    etaB = np.copy(target_dict['W262_B_CH'])
    betaB = np.copy(target_dict['W262_B_CB'])
    orientA = np.sum((etaA - betaA) * SFunitV, axis=0)
    orientB = np.sum((etaB - betaB) * SFunitV, axis=0)

    fenA = np.copy(target_dict['fenA'])[0,:]
    fenB = np.copy(target_dict['fenB'])[0,:]
    time = np.copy(target_dict['time'])[0,:]
    #      0     1     2    3    4
    print np.shape(wingA)
    print np.shape(time)
    return wingA,wingB,fenA,fenB,time

# setup plot grid
plt.figure()
gs = GridSpec(3,2)
gs.update(hspace=0.2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])
ax3 = plt.subplot(gs[2,0])
ax4 = plt.subplot(gs[0,1])
ax5 = plt.subplot(gs[1,1])
ax6 = plt.subplot(gs[2,1])

def timeseries(tgt_ax,datalist,time,datalabels=None,references=None,reference_labels=None,ylabel=None,xlabel='time (ns)',plotname=None,ylimits=None,multicolor=True):
    # arbitrary color choices
    if multicolor:
        colors = ['red', 'blue', 'green', 'purple']
    else:
        colors = ['black','black','black','black']
    ref_colors = ['silver', 'black']
    # setup timeseries data
    for idx,elem in enumerate(datalist):
        tgt_ax.plot(time/1000.0, elem, c=colors[idx])
        if datalabels:
            tgt_ax.text(0.99,0.99-0.1*idx, datalabels[idx], color=colors[idx],
                        transform=tgt_ax.transAxes,verticalalignment='top',horizontalalignment='right', fontsize=14)
    # optionally add reference lines
    if references:
        for idx,elem in enumerate(references):
            tgt_ax.axhline(y=elem,c=ref_colors[idx],linewidth=1.0,linestyle='dashed',zorder=0)
            if reference_labels:
                tgt_ax.text(0.01,0.99-0.1*idx, reference_labels[idx], color=ref_colors[idx],
                            transform=tgt_ax.transAxes,verticalalignment='top',horizontalalignment='left', fontsize=14)
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
traakWT_ref = process_refs(ref_data[0].primaryDS.dataset_to_dict())
downAB = process_refs(ref_data[3].primaryDS.dataset_to_dict())
downCD = process_refs(ref_data[4].primaryDS.dataset_to_dict())
upAB = process_refs(ref_data[5].primaryDS.dataset_to_dict())
upCD = process_refs(ref_data[6].primaryDS.dataset_to_dict())

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

wts = proc_dataset(wt_basic)
muts = proc_dataset(mut_basic)
tm4 = proc_dataset(tm4_basic)

timeseries(
           ax1,
           datalist=[wts[0],wts[1]],
           time=wts[4],
           datalabels=['subunit A','subunit B'],
           references=ref_wing,
           reference_labels=['4I9W A', '4I9W B'],
           ylabel=None,
           plotname=None,
           ylimits=[15,40],
           xlabel=None
           )

timeseries(
           ax2,
           datalist=[muts[0],muts[1]],
           time=muts[4],
           datalabels=['subunit A','subunit B'],
           references=ref_wing,
           reference_labels=['4I9W A', '4I9W B'],
           ylabel=None,
           plotname=None,
           ylimits=[15,40],
           xlabel=None
           )

timeseries(
           ax3,
           datalist=[tm4[0],tm4[1]],
           time=tm4[4],
           datalabels=['subunit A','subunit B'],
           references=ref_wing,
           reference_labels=['4I9W A', '4I9W B'],
           ylabel=None,
           xlabel=None,
           plotname=None,
           ylimits=[15,40]
           )

timeseries(
           ax4,
           datalist=[wts[2],wts[3]],
           time=wts[4],
           datalabels=['subunit A','subunit B'],
           references=ref_fen,
           reference_labels=['4I9W A', '4I9W B'],
           ylabel=None,
           xlabel=None,
           plotname=None,
           ylimits=[4,16]
           )

timeseries(
           ax5,
           datalist=[muts[2],muts[3]],
           time=muts[4],
           datalabels=['subunit A','subunit B'],
           references=ref_fen,
           reference_labels=['4I9W A', '4I9W B'],
           ylabel=None,
           xlabel=None,
           plotname=None,
           ylimits=[4,16]
           )

timeseries(
           ax6,
           datalist=[tm4[2],tm4[3]],
           time=tm4[4],
           datalabels=['subunit A','subunit B'],
           references=ref_fen,
           reference_labels=['4I9W A', '4I9W B'],
           ylabel=None,
           plotname=None,
           ylimits=[4,16],
           xlabel=None
           )
#plt.suptitle(report_name, fontsize=20)
plt.show()
