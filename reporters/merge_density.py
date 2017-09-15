from __future__ import print_function
import sys, os, re
import numpy as np
from copy import deepcopy
from gridData import Grid

data_dir = '/Users/anatale/UCSF/Grabe_lab/data/traak_data/density'

searchterm = re.compile('h2o_wt_s')

gridlist = []

#fliplist = [0,1,1,0,0,1] # wt
fliplistS = [1,0,0,1] # wtS
#fliplist = [0,1,1,1,0,1] # mut
#fliplistS = [1,1,0,1] # mutS

for filename in os.listdir(data_dir):
    if re.search(searchterm, filename):
        print(filename)
        gridlist.append(Grid(os.path.join(data_dir,filename)))

print(gridlist)

newgrid = deepcopy(gridlist[0])
newgrid.grid.fill(0.0)

for i,g in enumerate(gridlist):
    if fliplistS[1] == 1:
        newgrid.grid = np.add(newgrid.grid, g.grid[::-1,::-1,:])
    else:
        newgrid.grid = np.add(newgrid.grid, g.grid)

newgrid.grid = newgrid.grid / len(gridlist)

newgrid.export(os.path.join(data_dir,'wt_h2o_s_nosym.dx'), file_format='dx')
