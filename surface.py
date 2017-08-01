from __future__ import print_function
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as mcm
from analysis.MembraneSurface import MembraneSurface

datadir = '/Users/anatale/UCSF/Grabe_lab/scratch/'
topofile = os.path.join(datadir,'tm4_nowater0.pdb')
dcdfile = os.path.join(datadir,'tm4_nowater.dcd')

a = MembraneSurface()
a.load_traj(topofile,dcdfile,50)
a.run(100, 100, 'name C22 C32', 'name C218', stype='interp')

a.write_aligned_subset('tm4_memb_align')

# print(s.Up_hist)
# print(np.mean(s.Up_hist[s.Up_hist>0]), np.median(s.Up_hist[s.Up_hist>0]))

# print(np.amax(s.Up[~np.isnan(s.Up)]), np.amin(s.Up[~np.isnan(s.Up)]))
# print(np.amax(s.Um[~np.isnan(s.Um)]), np.amin(s.Um[~np.isnan(s.Um)]))
#
a.write_vmd_surfaces()
# Up_grid = Grid(surfer.Up, edges=(surfer.xedges,surfer.yedges))
# Up_grid.save('Up_test.pkl')
# Um_grid = Grid(surfer.Um, edges=(surfer.xedges,surfer.yedges))
# Um_grid.save('Um_test.pkl')

# Up = Grid('Up_test.pickle')
# Um = Grid('Um_test.pickle')
# print(Um.grid)

# fig = plt.figure()
# ax = Axes3D(fig)
# ax.plot_surface(s.X, s.Y, s.Up,
#   rstride=1, cstride=1, cmap=mcm.coolwarm, linewidth=0, antialiased=False,
#   vmin=np.amin(s.Up[~np.isnan(s.Up)]), vmax=np.amax(s.Up[~np.isnan(s.Up)]))
# # # ax.plot_surface(s.X, s.Y, s.Up_hist,
# # #   rstride=1, cstride=1, cmap=mcm.coolwarm, linewidth=0, antialiased=False,
# # #   vmin=np.amin(s.Up_hist), vmax=np.amax(s.Up_hist))
# plt.show()
#
# fig = plt.figure()
# ax = Axes3D(fig)
# ax.plot_surface(s.X, s.Y, s.Um,
#   rstride=1, cstride=1, cmap=mcm.coolwarm, linewidth=0, antialiased=False,
#   vmin=np.amin(s.Um[~np.isnan(s.Um)]), vmax=np.amax(s.Um[~np.isnan(s.Um)]))
# # # ax.plot_surface(s.X, s.Y, s.Up_hist,
# # #   rstride=1, cstride=1, cmap=mcm.coolwarm, linewidth=0, antialiased=False,
# # #   vmin=np.amin(s.Up_hist), vmax=np.amax(s.Up_hist))
# plt.show()
