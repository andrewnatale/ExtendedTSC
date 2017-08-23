import sys, os, math, re
import numpy as np
from itertools import cycle
import matplotlib.pyplot as plt
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

dat_root_dir = '/Users/anatale/UCSF/Grabe_lab/data/traak_data/etsc_out_20170718'

report_name = traj_names[rep_idx]

figures_outdir = os.path.join(dat_root_dir, 'figures')

# load data sets
datasets = {}
datasets['dihedrals'] = tsds(infilename=os.path.join(dat_root_dir, '%s_filter_dihedrals_all_frames.dat' % report_name))
datasets['carbonyls'] = tsds(infilename=os.path.join(dat_root_dir, '%s_filter_carbonyls_all_frames.dat' % report_name))
datasets['waters'] = tsds(infilename=os.path.join(dat_root_dir, '%s_filter_water_all_frames.dat' % report_name))
datasets['potassium'] = tsds(infilename=os.path.join(dat_root_dir, '%s_filter_potassium_vt_all_frames.dat' % report_name))
datasets['potassium_mask'] = tsds(infilename=os.path.join(dat_root_dir, '%s_filter_potassium_vt_all_frames.mask.dat' % report_name))
datasets['dihed_vector'] = tsds(infilename=os.path.join(dat_root_dir, 'state_vectors', '%s_sf_dihed_vector.dat' % report_name))
datasets['rmsf'] = tsds(infilename=os.path.join('/Users/anatale/UCSF/Grabe_lab/data/traak_data/filter_rmsf', '%s.50ps.alignYZcenter.rmsf_filter_CA.dat' % report_name))

# set reference z=0 for plotting
zref = np.zeros_like(datasets['carbonyls']['S0gly'][2])
zref += datasets['carbonyls']['S0gly'][2]
zref += datasets['carbonyls']['S1top'][2]
zref += datasets['carbonyls']['S2S1'][2]
zref += datasets['carbonyls']['S3S2'][2]
zref += datasets['carbonyls']['S4S3'][2]
zref += datasets['carbonyls']['S4bottom'][2]
zref /= 6.0
#print np.shape(zref)

# setup figure
f = plt.figure(figsize=(12,7))
ax = f.add_subplot(111)

stdtime = datasets['carbonyls']['time'].flatten() / 1000.0
ax.set_xlim(stdtime[0],stdtime[-1])
ax.set_ylim(-11,9)
ax.set_xlabel('time (ns)')
ax.set_ylabel('z-coordinate relative to selectivity filter center (Angstroms)')
plt.title('%s selectivity filter timeseries' % report_name)

# formatter to plot pinching of carbonyls at each filter level
def plot_flip(tgtkey):
    if tgtkey == 'S4bottom':
        angle = 'chi1'
    else:
        angle = 'psi'
    # get avg oxygen z-coordinate and align
    target = datasets['carbonyls'][tgtkey][2]-zref
    #target = datadicts['carbonyls'][tgtkey][2,:]

    mask = np.sum( [datasets['dihed_vector']['%s_%s_flip' % (i, angle)].series.flatten() for i in layers[tgtkey]], axis=0 )
    mask1 = mask<1
    mask2 = mask==1
    mask3 = mask>1
    noflip = target[mask1]
    oneflip = target[mask2]
    multiflip = target[mask3]
    noflip_time = stdtime[mask1]
    oneflip_time = stdtime[mask2]
    multiflip_time = stdtime[mask3]
    ax.scatter(noflip_time, noflip, color='black', s=2, edgecolor='none', zorder=2)
    ax.scatter(oneflip_time, oneflip, color='darksalmon', s=2, edgecolor='none', zorder=2)
    ax.scatter(multiflip_time, multiflip, color='red', s=2, edgecolor='none', zorder=2)

# put pinching data into formatter
layers = \
{
'S0gly':['133A','133B','242A','242B'],
'S1top':['132A','132B','241A','241B'],
'S2S1':['131A','131B','240A','240B'],
'S3S2':['130A','130B','239A','239B'],
'S4S3':['129A','129B','238A','238B'],
'S4bottom':['129A','129B','238A','238B']
}
for key in layers:
    plot_flip(key)

# plot water positons
for idx,time in np.ndenumerate(stdtime):
    waters = datasets['waters']['zsearchlist'].series[:,idx]
    waters = waters[~np.isnan(waters)]
    waters = waters - zref[idx]
    timeblob = np.empty_like(waters)
    timeblob.fill(time)
    ax.scatter(timeblob, waters, color='lightsteelblue', s=2, edgecolor='none', zorder=1)

colors = cycle(['green','fuchsia','orange','olivedrab','darkblue','limegreen','mediumpurple','gold','teal','maroon','sienna','orchid'])

for key in datasets['potassium']:
    if key != 'time':
        zcoords = datasets['potassium'][key][2] - zref
        mask = datasets['potassium_mask'][key][0]
        #print np.shape(zcoords), np.shape(mask)
        sel_zcoords = zcoords[mask==1]
        sel_time = stdtime[mask==1]
        #print np.shape(sel_zcoords)
        ax.scatter(sel_time, sel_zcoords, color=colors.next(), s=3, edgecolor='none', zorder=3)

# show or save plot
plt.show()
#f.savefig(os.path.join(figures_outdir, '%s_filter_ts.png' % report_name), bb_inches='tight')
