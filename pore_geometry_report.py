import sys, os, math, re
import numpy as np
from itertools import cycle
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
#import ExtendedTSC_old
from core.TimeseriesCore import TimeseriesCore

# choose report to plot
rep_idx = int(sys.argv[1])
reports = ['traakWT_full_npt','traakWT_s1s3_npt','traakWT_s2s4_npt','traakG124I_full_npt','traakG124I_s1s3_npt','traakG124I_s2s4_npt','traakTM4_npt']
dat_root_dir = '/Users/anatale/school/UCSF/Grabe_lab/data/traj_features'
tgt_data = '50ps/pore_features'

report_name = reports[rep_idx]

figures_dir = os.path.join('/Users/anatale/school/UCSF/Grabe_lab/data/traj_features', tgt_data, 'figures')

reader = TimeseriesCore(verbose=False)

# load data sets
datasets = {}
datasets['dihedrals'] = reader._data_reader(os.path.join(dat_root_dir, tgt_data, '%s_pore_bb_dihedrals_all_frames.dat' % report_name))
datasets['carbonyls'] = reader._data_reader(os.path.join(dat_root_dir, tgt_data, '%s_pore_carbonyls_all_frames.dat' % report_name))
datasets['waters'] = reader._data_reader(os.path.join(dat_root_dir, tgt_data, '%s_pore_water_all_frames.dat' % report_name))
datasets['potassium'] = reader._data_reader(os.path.join(dat_root_dir, tgt_data, '%s_pore_potassium_all_frames.dat' % report_name))
datasets['vectors'] = reader._data_reader(os.path.join(dat_root_dir, tgt_data, '%s_vectorized.dat' % report_name))

datadicts = {}
for key in datasets:
    datadicts[key] = datasets[key].dataset_to_dict()

# set reference z=0 for plotting
zref = (datadicts['carbonyls']['S1top'][2,:]
     + datadicts['carbonyls']['S2S1'][2,:]
     + datadicts['carbonyls']['S3S2'][2,:]
     + datadicts['carbonyls']['S4S3'][2,:]
     + datadicts['carbonyls']['S4bottom'][2,:]) / 5.0
#print np.shape(zref)

# setup figure
f = plt.figure(figsize=(16,9))
ax = f.add_subplot(111)

stdtime = datadicts['carbonyls']['time'].flatten() / 1000.0
ax.set_xlim(stdtime[0],stdtime[-1])
ax.set_ylim(-8.5,8.5)
ax.set_xlabel('time (ns)')
ax.set_ylabel('z-coordinate relative to selectivity filter center (Angstroms)')
plt.title('%s selectivity filter timeseries' % report_name)

# formatter to plot pinching of carbonyls at each filter level
def plot_pinch(tgtkey,maskkey):
    if maskkey is None:
        target = datadicts['carbonyls'][tgtkey][2,:]-zref
        ax.scatter(stdtime, target, c='black', s=2, edgecolor='none')
    else:
        target = datadicts['carbonyls'][tgtkey][2,:]-zref
        mask = datadicts['vectors'][maskkey].flatten()
        mask1 = mask<1
        mask2 = np.logical_and(mask>0,mask<2)
        mask3 = mask>1
        unpinch = target[mask1]
        onepinch = target[mask2]
        multipinch = target[mask3]
        unpinched_time = stdtime[mask1]
        onepinch_time = stdtime[mask2]
        multipinch_time = stdtime[mask3]
        ax.scatter(unpinched_time, unpinch, color='black', s=2, edgecolor='none')
        ax.scatter(onepinch_time, onepinch, color='red', s=2, edgecolor='none')
        ax.scatter(multipinch_time, multipinch, color='red', s=2, edgecolor='none')

# put pinching data into formatter
dha_map = {
'S1top':'s1_pinch',
'S2S1':'s2_pinch',
'S3S2':'s3_pinch',
'S4S3':None,
'S4bottom':None}
for subkey in dha_map:
    plot_pinch(subkey,dha_map[subkey])

# plot water positons
for idx,time in np.ndenumerate(stdtime):
    waters = datadicts['waters']['zsearchlist'][:,idx]
    waters = waters[~np.isnan(waters)]
    waters = waters - zref[idx]
    timeblob = np.empty_like(waters)
    timeblob.fill(time)
    ax.scatter(timeblob, waters, color='deepskyblue', s=2, edgecolor='none')

# plot potassium positions
for idx,time in np.ndenumerate(stdtime):
    potassium = datadicts['potassium']['zsearchlist'][:,idx]
    potassium = potassium[~np.isnan(potassium)]
    potassium = potassium - zref[idx]
    timeblob = np.empty_like(potassium)
    timeblob.fill(time)
    ax.scatter(timeblob, potassium, color='orange', s=3, edgecolor='none')

# show or save plot
plt.show()
#f.savefig(os.path.join(figures_dir, '%s_poreTS.png' % report_name), bb_inches='tight')
