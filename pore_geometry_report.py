import sys, os, math, re
import numpy as np
from itertools import cycle
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
#import ExtendedTSC_old
from core.TimeseriesCore import TimeseriesCore

# choose report to plot
rep_idx = int(sys.argv[1])
reports = ['traakWT_full','traakWT_S1S3','traakWT_S2S4','traakG124I_full_npt','traakG124I_S1S3','traakG124I_S2S4','traakTM4']

report_name = reports[rep_idx]
dat_root_dir = '/Users/anatale/school/UCSF/Grabe_lab/data/traj_features/'
tgt_data = '500ps'

# load data sets
a = TimeseriesCore(datfilename=os.path.join(dat_root_dir, tgt_data, 'pore_geo', report_name+'.pore_geo.dat'))
#b = ExtendedTSC.ExtendedTSC(datfile=os.path.join('pore_geo', reports[1]+'.pore_geo.dat'))
c = TimeseriesCore(datfilename=os.path.join(dat_root_dir, tgt_data, 'core_volume_tracking', report_name+'.core_vol.dat'),maskfilename=os.path.join(dat_root_dir, tgt_data, 'core_volume_tracking', report_name+'.core_vol.mask.dat'))
d = TimeseriesCore(datfilename=os.path.join(dat_root_dir, tgt_data, 'core_volume_tracking', report_name+'.core_vol.tip3.dat'))

# make data sets indexable
wt = a.primaryDS.dataset_to_dict()
#mut = b.primaryDS.simplify_indexing(return_dict=True)
klipids = c.primaryDS.dataset_to_dict()
kmask = c.maskDS.dataset_to_dict()
tip3z = d.primaryDS.dataset_to_dict()

# for each sf carbonyl ~plane, calculate area and avg z-displacement
names = {
'D':['a_129og1','a_238og1','b_129og1','b_238og1'],
'E':['a_129o','a_238o','b_129o','b_238o'],
'F':['a_130o','a_239o','b_130o','b_239o'],
'G':['a_131o','a_240o','b_131o','b_240o'],
'H':['a_132o','a_241o','b_132o','b_241o'],
'I':['a_133o','a_242o','b_133o','b_242o']
}

# reference z=0
# zref = np.reshape(np.copy(klipids['S4bottom'][2,:]), (1,-1))
zref = None
count = 0
for key in ['D','E','F','G','H']:
    for elem in names[key]:
        if zref is None:
            zref = wt[elem][2,:]
        else:
            zref = zref + wt[elem][2,:]
        count += 1
zref = np.reshape((zref / float(count)), (1,-1))

print zref
zcoords = {}
area = {}
# avg 4-carbonyl z-coord
for key in names:
    print key,names[key]
    zcoords[key] = (wt[names[key][0]][2,:]+wt[names[key][1]][2,:]+wt[names[key][2]][2,:]+wt[names[key][3]][2,:])/4 - zref[0,:]
    #zcoords[key] = (wt[names[key][0]][2,:]+wt[names[key][1]][2,:]+wt[names[key][2]][2,:]+wt[names[key][3]][2,:])/4

    area[key] = (np.linalg.norm(np.cross(wt[names[key][0]]-wt[names[key][1]], wt[names[key][2]]-wt[names[key][1]], axisa=0,axisb=0,axisc=0), axis=0)
                +np.linalg.norm(np.cross(wt[names[key][0]]-wt[names[key][3]], wt[names[key][2]]-wt[names[key][3]], axisa=0,axisb=0,axisc=0), axis=0)
                +np.linalg.norm(np.cross(wt[names[key][1]]-wt[names[key][0]], wt[names[key][3]]-wt[names[key][0]], axisa=0,axisb=0,axisc=0), axis=0)
                +np.linalg.norm(np.cross(wt[names[key][1]]-wt[names[key][2]], wt[names[key][3]]-wt[names[key][2]], axisa=0,axisb=0,axisc=0), axis=0)
                )/2.0
    area[key] = area[key]
    area[key] = area[key]/40.0

for key in names:
    print key
    print np.amin(zcoords[key]),np.amax(zcoords[key])
    print np.amin(area[key]),np.amax(area[key])

howlong = np.shape(klipids['time'])[1]
reshapetime = np.reshape(klipids['time'], (howlong,)) / 1000.0
# print np.shape(reshapetime)
f, ax = plt.subplots()
ax.plot(reshapetime, zcoords['D'], c='black', lw=1)
ax.fill_between(reshapetime, zcoords['D']-area['D'], zcoords['D']+area['D'], facecolor='green', alpha=0.5)
ax.plot(reshapetime, zcoords['E'], c='black', lw=1)
ax.fill_between(reshapetime, zcoords['E']-area['E'], zcoords['E']+area['E'], facecolor='green', alpha=0.5)
ax.plot(reshapetime, zcoords['F'], c='black', lw=1)
ax.fill_between(reshapetime, zcoords['F']-area['F'], zcoords['F']+area['F'], facecolor='green', alpha=0.5)
ax.plot(reshapetime, zcoords['G'], c='black', lw=1)
ax.fill_between(reshapetime, zcoords['G']-area['G'], zcoords['G']+area['G'], facecolor='green', alpha=0.5)
ax.plot(reshapetime, zcoords['H'], c='black', lw=1)
ax.fill_between(reshapetime, zcoords['H']-area['H'], zcoords['H']+area['H'], facecolor='green', alpha=0.5)
ax.plot(reshapetime, zcoords['I'], c='black', lw=1)
ax.fill_between(reshapetime, zcoords['I']-area['I'], zcoords['I']+area['I'], facecolor='palevioletred', alpha=0.3)

# setup colors for different entities
Kcolors = ['violet','silver','tomato','blue','blueviolet','palevioletred', 'orange']
Kcolor = cycle(Kcolors).next
watercolor = 'powderblue'
# regex for sorting data
find_K = re.compile('POT')
# begin plotting to tgt_ax
for key in klipids:
    # potassium ions on top, each with a unique(ish) color
    if find_K.search(key):
        kzcoords = np.reshape(np.copy(klipids[key][2,:]), (1,-1))
        kzcoords = kzcoords - zref
        #kzcoords = kzcoords
        boolmask = kmask[key].astype(bool)
        masked_kzcoords = kzcoords[boolmask]
        masked_time = klipids['time'][boolmask] / 1000.0
        ax.plot(masked_time, masked_kzcoords, c=Kcolor(), zorder=4, marker='o', mew=0.2, linestyle='None')
# water oxygen atoms, all the same color
for idx,time in np.ndenumerate(tip3z['time']):
    waters = tip3z['water_volume'][:,idx[1]]
    waters = waters[~np.isnan(waters)]
    waters = waters - zref[idx]
    #waters = waters
    timeblob = np.empty_like(waters)
    timeblob.fill(time/1000.0)
    ax.plot(timeblob, waters, c=watercolor, zorder=3, marker='o', mew=0.0, ms=3, linestyle="None")
# limits and labels
ax.set_xlim([klipids['time'][0,0]/1000.0, klipids['time'][0,-1]/1000.0])
ax.set_ylim([-9,9])
plt.show()