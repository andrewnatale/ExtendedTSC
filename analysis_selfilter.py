import sys, os, math
import numpy as np
import re
import ExtendedTSC
import traak_plots
#from sklearn import cluster
import matplotlib.pyplot as plt

# load coordinate data
#infile = 'testing1.dat'
coordfile = 'coordinates_test.dat'
a = ExtendedTSC.ExtendedTSC(coordfile)
access = a.simplify_indexing()

# load mask data
maskfile = 'masks_test.dat'
b = ExtendedTSC.ExtendedTSC(maskfile,mask=True)
mask = b.simplify_indexing()

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

colors = ['gold','limegreen','violet','silver']
watercolor = 'blue'
color_idx = 0
find_water = re.compile('water')
find_K = re.compile('pot')

f, ax = plt.subplots()
for key in access:
    if find_water.search(key):
        heights = np.reshape(np.sum((access[key] - S4C) * SFunitV, axis=0), (1,-1))
        boolmask = mask[key].astype(bool)
        masked_heights = heights[boolmask]
        masked_time = access['time'][boolmask]
        ax.scatter(masked_time, masked_heights, c=watercolor)
    elif find_K.search(key):
        heights = np.reshape(np.sum((access[key] - S4C) * SFunitV, axis=0), (1,-1))
        boolmask = mask[key].astype(bool)
        masked_heights = heights[boolmask]
        masked_time = access['time'][boolmask]
        ax.scatter(masked_time, masked_heights, c=colors[color_idx])
        color_idx += 1
    else:
        pass
ax.set_xlim([access['time'][0,0],access['time'][0,-1]])
ax.set_ylim([-2,15])
plt.show()
