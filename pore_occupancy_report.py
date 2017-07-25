import sys, os, math, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from core.TimeseriesCore import TimeseriesCore

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

dat_root_dir = '/Users/anatale/UCSF/Grabe_lab/data/traak_data/etsc_out_20170718'

report_name = traj_names[rep_idx]

figures_outdir = os.path.join(dat_root_dir, 'figures')

reader = TimeseriesCore(verbose=False)

# load data sets
datasets = {}
datasets['carbonyls'] = reader._data_reader(os.path.join(dat_root_dir, '%s_filter_carbonyls_all_frames.dat' % report_name))
datasets['waters'] = reader._data_reader(os.path.join(dat_root_dir, '%s_pore_water_all_frames.dat' % report_name))
datasets['lipids'] = reader._data_reader(os.path.join(dat_root_dir, '%s_pore_lipids_all_frames.dat' % report_name))
datasets['basic'] = reader._data_reader(os.path.join(dat_root_dir, '%s_basic_features_all_frames.dat' % report_name))

datadicts = {}
for key in datasets:
    datadicts[key] = datasets[key].dataset_to_dict()

# set reference z=0 for plotting
zref = datadicts['carbonyls']['S4bottom'][2,:]

# setup figure
plt.figure(figsize=(16,9))
gs = GridSpec(2,2)
gs.update(hspace=0.2)
gs.update(wspace=0.2)
ax1 = plt.subplot(gs[:,0])
ax2 = plt.subplot(gs[0,1])
ax3 = plt.subplot(gs[1,1])

stdtime = datadicts['carbonyls']['time'].flatten() / 1000.0

# ax1 - z coords of everything over time
ax1.set_xlim(stdtime[0],stdtime[-1])
ax1.set_ylim(-16,0)
ax1.set_xlabel('time (ns)')
ax1.set_ylabel('z-coordinate relative to selectivity filter (Angstroms)')

# plot water positions
for idx,time in np.ndenumerate(stdtime):
    waters = datadicts['waters']['zsearchlist'][:,idx]
    waters = waters[~np.isnan(waters)]
    waters = waters - zref[idx]
    timeblob = np.empty_like(waters)
    timeblob.fill(time)
    ax1.scatter(timeblob, waters, color='deepskyblue', s=2, edgecolor='none', zorder=2, alpha=0.5)

# plot lipid positions
for idx,time in np.ndenumerate(stdtime):
    lipids = datadicts['lipids']['zsearchlist'][:,idx]
    lipids = lipids[~np.isnan(lipids)]
    lipids = lipids - zref[idx]
    timeblob = np.empty_like(lipids)
    timeblob.fill(time)
    ax1.scatter(timeblob, lipids, color='orange', s=3, edgecolor='none', zorder=3, alpha=0.5)

# ax2 - counts of water and lipid over time
water_count = ~np.isnan(datadicts['waters']['zsearchlist'])
water_count = np.sum(water_count, axis=0)
lipid_count = ~np.isnan(datadicts['lipids']['zsearchlist'])
lipid_count = np.sum(lipid_count, axis=0)

convolution = np.ones((50,))/50.0
conv_time = np.convolve(stdtime,convolution,mode='valid')
avg_lipid = np.convolve(lipid_count,convolution,mode='valid')
avg_water = np.convolve(water_count,convolution,mode='valid')

print np.amax(water_count), np.amin(water_count)
print np.amax(lipid_count), np.amin(lipid_count)

ax2.scatter(stdtime,water_count,color='deepskyblue', s=10, edgecolor='none', zorder=3, alpha=0.3)
ax2.plot(conv_time,avg_water,zorder=5, color='blue')
ax4 = ax2.twinx()
ax4.scatter(stdtime,lipid_count,color='orange',s=10, edgecolor='none', zorder=3, alpha=0.3)
ax4.plot(conv_time,avg_lipid,zorder=5, color='sienna')
ax2.set_xlim(stdtime[0],stdtime[-1])
ax4.set_xlim(stdtime[0],stdtime[-1])
ax2.set_ylim(0,120)
ax4.set_ylim(0,80)
ax2.set_xlabel('time (ns)')
ax2.set_ylabel('# of water molecules in pore')
ax4.set_ylabel('# of lipid carbons in pore')
ax2.text(0.5,0.99, 'Water', color='deepskyblue',
  transform=ax2.transAxes,verticalalignment='top',horizontalalignment='center', fontsize=12)
ax2.text(0.5,0.94, 'Lipid', color='orange',
  transform=ax2.transAxes,verticalalignment='top',horizontalalignment='center', fontsize=12)
# ax3 - pore_width at T277 on M4
width = datadicts['basic']['pore_width'][0,:]
ax3.scatter(stdtime,width,color='grey',s=10, edgecolor='none', alpha=0.3,zorder=1)
avg_width = np.convolve(width,convolution,mode='valid')
ax3.plot(conv_time,avg_width,color='black',zorder=5)
ax3.set_xlim(stdtime[0],stdtime[-1])
ax3.set_ylim(8.0,22.0)
ax3.set_xlabel('time (ns)')
ax3.set_ylabel('Pore width; T277 CA-CA distance (Angstroms)')

plt.suptitle('%s pore occupancy' % report_name, fontsize=20)
# show or save plot
#plt.show()
plt.savefig(os.path.join(figures_outdir, '%s_pore_occupancy.png' % report_name), bb_inches='tight')
