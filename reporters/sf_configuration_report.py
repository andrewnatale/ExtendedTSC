import sys, os, math, re
import numpy as np
import matplotlib.pyplot as plt
from core.DataSet import DataSet

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

# load data sets
datasets = {}
datasets['dihedrals'] = DataSet(infilename=os.path.join(dat_root_dir, '%s_filter_dihedrals_all_frames.dat' % report_name))
datasets['carbonyls'] = DataSet(infilename=os.path.join(dat_root_dir, '%s_filter_carbonyls_all_frames.dat' % report_name))
datasets['waters'] = DataSet(infilename=os.path.join(dat_root_dir, '%s_filter_water_all_frames.dat' % report_name))
datasets['potassium'] = DataSet(infilename=os.path.join(dat_root_dir, '%s_filter_potassium_all_frames.dat' % report_name))
datasets['vectors'] = DataSet(infilename=os.path.join(dat_root_dir, 'state_vectors', '%s_vectorized.dat' % report_name))

datadicts = {}
for key in datasets:
    datadicts[key] = datasets[key].dataset_to_dict()

# set reference z=0 for plotting
zref = (datadicts['carbonyls']['S0gly'][2,:]
     + datadicts['carbonyls']['S1top'][2,:]
     + datadicts['carbonyls']['S2S1'][2,:]
     + datadicts['carbonyls']['S3S2'][2,:]
     + datadicts['carbonyls']['S4S3'][2,:]
     + datadicts['carbonyls']['S4bottom'][2,:]) / 6.0
#print np.shape(zref)

# setup figure
f = plt.figure(figsize=(16,9))
ax = f.add_subplot(111)

stdtime = datadicts['carbonyls']['time'].flatten() / 1000.0
ax.set_xlim(stdtime[0],stdtime[-1])
ax.set_ylim(-11,9)
ax.set_xlabel('time (ns)')
ax.set_ylabel('z-coordinate relative to selectivity filter center (Angstroms)')
plt.title('%s selectivity filter timeseries' % report_name)

# formatter to plot pinching of carbonyls at each filter level
def plot_pinch(tgtkey,maskkey):
    if maskkey is None:
        target = datadicts['carbonyls'][tgtkey][2,:]-zref
        #target = datadicts['carbonyls'][tgtkey][2,:]
        ax.scatter(stdtime, target, c='black', s=2, edgecolor='none')
    else:
        target = datadicts['carbonyls'][tgtkey][2,:]-zref
        #target = datadicts['carbonyls'][tgtkey][2,:]
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
'S4bottom':None,
'S0gly':None}
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
    ax.scatter(timeblob, potassium, color='olivedrab', s=3, edgecolor='none')

# show or save plot
#plt.show()
f.savefig(os.path.join(figures_outdir, '%s_poreTS.png' % report_name), bb_inches='tight')
