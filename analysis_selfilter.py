import sys, os, math
import numpy as np
import re
import ExtendedTSC
import traak_plots
#from sklearn import cluster
import matplotlib.pyplot as plt

# load coordinate data
input_dir = 'filter_occupancy05252017'
prefix = 'traakG124I_S2S4_npt'
coordfile = os.path.join(input_dir, prefix+'atom_coordinates.dat')
a = ExtendedTSC.ExtendedTSC(coordfile)
access = a.simplify_indexing(return_dict=True)

# load mask data
maskfile =  os.path.join(input_dir, prefix+'masks.dat')
b = ExtendedTSC.ExtendedTSC(maskfile,mask=True)
mask = b.simplify_indexing(return_dict=True)

# print info
for key,val in access.iteritems():
    print key, np.shape(val)
print ' '
for key,val in mask.iteritems():
    print key, np.shape(val)
    #print val

# calculate unit vectors in the S4 to S1 directions
S1C = np.copy(access['S1top'])
S4C = np.copy(access['S4bottom'])
SFunitV = (S1C - S4C) / np.linalg.norm(S1C - S4C, axis=0)

colors = ['gold','limegreen','violet','silver', 'orange', 'tomato', 'seagreen', 'blue', 'yellow', 'blueviolet', 'palevioletred', 'darkgoldenrod']
watercolor = 'powderblue'
color_idx = 0
find_water = re.compile('water')
find_K = re.compile('pot')

f, ax = plt.subplots()
for key in access:
    if find_water.search(key):
        heights = np.reshape(np.sum((access[key] - S4C) * SFunitV, axis=0), (1,-1))
        boolmask = mask[key].astype(bool)
        masked_heights = heights[boolmask]
        masked_time = access['time'][boolmask] / 1000.0
        ax.scatter(masked_time, masked_heights, c=watercolor)
    elif find_K.search(key):
        heights = np.reshape(np.sum((access[key] - S4C) * SFunitV, axis=0), (1,-1))
        boolmask = mask[key].astype(bool)
        masked_heights = heights[boolmask]
        masked_time = access['time'][boolmask] / 1000.0
        ax.scatter(masked_time, masked_heights, c=colors[color_idx])
        color_idx += 1
    else:
        pass
ax.set_xlim([access['time'][0,0] / 1000.0 ,access['time'][0,-1] / 1000.0])
ax.set_ylim([-2,15])
ax.set_xlabel('time (ns)')
ax.set_ylabel('Location along pore axis (angstroms)')
ax.set_title(prefix)
# plt.show()
os.chdir(input_dir)
plt.savefig(prefix+'.pore_occupancy.png', bbox_inches='tight')
