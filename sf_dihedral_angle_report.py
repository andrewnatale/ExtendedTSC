import sys, os, math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from core.TimeseriesCore import TimeseriesCore as tsc

# choose report to plot
rep_idx = int(sys.argv[1])

save=False

datfile_dir = '/Users/anatale/UCSF/Grabe_lab/data/traak_data/etsc_out_20170718'

traj_names = [
'traakWT_full_npt.sim1',
'traakWT_full_npt.sim2',
'traakWT_s1s3_npt.sim1',
'traakWT_s1s3_npt.sim2',
'traakWT_s2s4_npt.sim1',
'traakWT_s2s4_npt.sim2',
'traakG124I_full_npt.sim1',
'traakG124I_full_npt.sim2',
'traakG124I_s1s3_npt.sim1',
'traakG124I_s1s3_npt.sim2',
'traakG124I_s2s4_npt.sim1',
'traakG124I_s2s4_npt.sim2',
'traakTM4_npt',
'traakExtPip2'
]

report_name = traj_names[rep_idx]

feature_names = [
'filter_dihedrals',
]

figures_outdir = '/Users/anatale/UCSF/Grabe_lab/data/traak_data/etsc_out_20170718/figures'

try:
    os.makedirs(figures_outdir)
except OSError:
    if not os.path.isdir(figures_outdir):
        raise

reader = tsc()

# load data
filename = os.path.join(datfile_dir, '%s_filter_dihedrals_all_frames.dat' % report_name)
dataset = reader._data_reader(filename)
datadict = dataset.dataset_to_dict()

stdtime = datadict['time'].flatten() / 1000.0

layers = \
{
'S0top':['133A','133B','242A','242B'],
'S0S1':['132A','132B','241A','241B'],
'S1S2':['131A','131B','240A','240B'],
'S2S3':['130A','130B','239A','239B'],
'S3S4':['129A','129B','238A','238B'],
'S4bottom':['129A','129B','238A','238B']
}

colors = {'S4':'red','S3':'blue','S2':'green','S1':'gold','S0':'violet'}
labels = \
{
'S4':'S4 - T129,T238',
'S3':'S3 - I130,V239',
'S2':'S2 - G131,G240',
'S1':'S1 - Y132,F241',
'S0':'S0 - G133,G242'
}

def angle_space():
    f, ax = plt.subplots()
    txt_offset = 0
    for key2 in sorted(layers.keys(),reverse=True):
    #for key2 in ['S1','S2','S3']:
        for res in layers[key2]:
            ax.scatter(datadict['%s_phi' % res], datadict['%s_psi' % res], s=2, c=colors[key2], edgecolor='none', alpha=0.3)
        ax.text(0.05,0.05+0.03*txt_offset, labels[key2], color=colors[key2],
                transform=ax.transAxes,verticalalignment='bottom',horizontalalignment='left', fontsize=8)
        txt_offset += 1
    ax.set_xlim([math.pi * -1.0, math.pi])
    ax.set_ylim([math.pi * -1.0, math.pi])
    ax.set_xlabel('Phi (radians)')
    ax.set_ylabel('Psi (radians)')
    ax.set_title(report_name)
    if save:
        f.savefig(os.path.join(figures_outdir,'%s_dihedrals.png' % report_name), bb_inches='tight')
    else:
        plt.show()

def psi_grid():
    plt.figure(figsize=(16,9))
    gs = GridSpec(6,4)
    gs.update(hspace=0)
    gs.update(wspace=0)
    all_axes = {}
    for i in [0,1,2,3,4,5]:
        for j in [0,1,2,3]:
            all_axes['ax%d%d' % (i,j)] = plt.subplot(gs[i,j])
    for i,key1 in enumerate(['S0top','S0S1','S1S2','S2S3','S3S4','S4bottom']):
        if key1 == 'S4bottom':
            angle = 'chi1'
        else:
            angle = 'psi'
        for j,key2 in enumerate(layers[key1]):
            tgt_ax = all_axes['ax%d%d' % (i,j)]
            tgt_ax.scatter(stdtime, np.unwrap(datadict['%s_%s' % (key2, angle)][0,:]), s=1, edgecolor='none')
            tgt_ax.set_xlim(stdtime[0],stdtime[-1])
            tgt_ax.set_ylim(math.pi * -1.0, math.pi)
            tgt_ax.text(0.01,0.99, '%s_%s' % (key2, angle), color='black', transform=tgt_ax.transAxes,verticalalignment='top',horizontalalignment='left')
            if i in [0,1,2,3,4]:
                plt.setp(tgt_ax.get_xticklabels(), visible=False)
            if j in [1,2,3]:
                plt.setp(tgt_ax.get_yticklabels(), visible=False)
    plt.suptitle(report_name)
    plt.show()
#angle_space()
psi_grid()
