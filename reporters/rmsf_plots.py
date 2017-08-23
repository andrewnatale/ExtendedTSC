import sys, os, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from core.TimeseriesDataSet import TimeseriesDataSet as tsds

# choose report to plot
rep_idx = int(sys.argv[1])

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

dat_root_dir = '/Users/anatale/UCSF/Grabe_lab/data/traak_data/filter_rmsf'

report_name = traj_names[rep_idx]

# figures_outdir = os.path.join(dat_root_dir, 'figures')

# load data sets
ds = tsds(infilename=os.path.join(dat_root_dir, '%s.50ps.alignYZcenter.rmsf_filter_CA.dat' % report_name))

print(np.amin(ds.data),np.amax(ds.data))

if report_name.startswith('traakTM4_npt'):
    layers = \
    [
    ['A_106','B_388','A_215','B_497'],
    ['A_105','B_387','A_214','B_496'],
    ['A_104','B_386','A_213','B_495'],
    ['A_103','B_385','A_212','B_494'],
    ['A_102','B_384','A_211','B_493'],
    ]
else:
    layers = \
    [
    ['A_133','B_133','A_242','B_242'],
    ['A_132','B_132','A_241','B_241'],
    ['A_131','B_131','A_240','B_240'],
    ['A_130','B_130','A_239','B_239'],
    ['A_129','B_129','A_238','B_238'],
    ]

layer_labels = ['S0 - res 133/242', 'S1 - res 132/241', 'S2 - res 131/240', 'S3 - res 130/239', 'S4 - res 129/238']

# layers.reverse()

styles = [('black','-'), ('black','--'), ('red','-'), ('red','--')]

plt.figure(figsize=(12,7))
gs = GridSpec(5,1)
gs.update(hspace=0)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])
ax3 = plt.subplot(gs[2,0])
ax4 = plt.subplot(gs[3,0])
ax5 = plt.subplot(gs[4,0])
axlist = [ax1,ax2,ax3,ax4,ax5]

for idx1,elem1 in enumerate(layers):
    for idx2,elem2 in enumerate(elem1):
        s = re.compile(elem2)
        for key in ds:
            if s.search(key):
                axlist[idx1].plot(ds['time']/1000.0, ds[key][0], color=styles[idx2][0], linestyle=styles[idx2][1], linewidth=1)
        axlist[idx1].set_xlim([ds['time'][0]/1000.0,ds['time'][-1]/1000.0])
        axlist[idx1].set_ylim((0.2,1.0))
        axlist[idx1].text(0.01,0.99, layer_labels[idx1], color='black', transform=axlist[idx1].transAxes,verticalalignment='top',horizontalalignment='left')
for ax in axlist[0:-1]:
    plt.setp(ax.get_xticklabels(), visible=False)

ax5.set_xlabel('time (ns)')
ax3.set_ylabel('C-alpha RMSF (Angstroms) over 2 ns window')
ax1.set_title(report_name)
#plt.show()
plt.savefig(os.path.join(dat_root_dir, 'figures', '%s_filter_rmsf.png' % report_name), bb_inches='tight')
