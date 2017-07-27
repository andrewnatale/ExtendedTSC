from __future__ import print_function
import sys, os
import numpy as np
#import scipy.interpolate as scin
#from statsmodels.robust.scale import mad
from gridData import Grid
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase

class _getMDdensity(AnalysisBase):

    def __init__(self, grid_dim, grid_len, density_sel, universe,
      grid_center_sel=None, align_sel=None, **kwargs):
        super(_getMDdensity,self).__init__(universe.trajectory,**kwargs)
        self.grid_dim = grid_dim
        self.grid_len = grid_len
        self.density_sel = density_sel
        self.u = universe
        if grid_center_sel:
            # do something
            #self.grid_center = ???
            pass
        else:
            self.grid_center = (0.0,0.0,0.0)
        if align_sel:
            # do something
            pass

    def _prepare(self):
        # setup vars and data structures
        # setup grid
        xmin = self.grid_center[0]-0.5*self.grid_len
        xmax = self.grid_center[0]+0.5*self.grid_len
        ymin = self.grid_center[1]-0.5*self.grid_len
        ymax = self.grid_center[1]+0.5*self.grid_len
        zmin = self.grid_center[2]-0.5*self.grid_len
        zmax = self.grid_center[2]+0.5*self.grid_len
        self.xedges,self.xspace = np.linspace(xmin,xmax,num=self.grid_dim+1,endpoint=True,retstep=True)
        self.yedges,self.yspace = np.linspace(ymin,ymax,num=self.grid_dim+1,endpoint=True,retstep=True)
        self.zedges,self.zspace = np.linspace(zmin,zmax,num=self.grid_dim+1,endpoint=True,retstep=True)
        self.xcenters = np.convolve(self.xedges,np.array([0.5,0.5]),mode='valid')
        self.ycenters = np.convolve(self.yedges,np.array([0.5,0.5]),mode='valid')
        self.zcenters = np.convolve(self.zedges,np.array([0.5,0.5]),mode='valid')
        self.X, self.Y, self.Z = np.meshgrid(self.xcenters,self.ycenters,self.zcenters)
        # setup atomgroups
        self.density_group = self.u.select_atoms(self.density_sel)
        self.coord_list = []
        self.frame_count = 0

    def _single_frame(self):
        # what to do at each frame
        print(self._ts)
        # copy coords for this timestep
        density_coords = self.density_group.positions.copy()
        # add coords to the list
        self.coord_list.append(density_coords)
        self.frame_count += 1

    def _conclude(self):
        # concatenate lists into arrays
        all_coords = np.concatenate(self.coord_list, axis=0)
        # compute histogram
        self.data,self.edges = np.histogramdd(all_coords,[self.xedges,self.yedges,self.zedges])
        self.data = self.data / self.frame_count
        #self.gridobj = Grid(self.data, edges=self.edges)

psffile = os.path.abspath(sys.argv[1])
dcdfile = os.path.abspath(sys.argv[2])
outname = sys.argv[3]

u = mda.Universe(psffile,dcdfile)

atomselection="\
name C22 or name C23 or name C24 or name C25 or name C26 or\
name C27 or name C28 or name C29 or name C210 or name C211 or\
name C212 or name C213 or name C214 or name C215 or name C216 or\
name C217 or name C218 or\
name C32 or name C33 or name C34 or name C35 or name C36 or\
name C37 or name C38 or name C39 or name C310 or name C311 or\
name C312 or name C313 or name C314 or name C315 or name C316"

dens = _getMDdensity(50,100, atomselection,u)
dens.run()

mygrid = Grid(dens.data, edges=dens.edges)
mygrid.export('%s' % outname, file_format='dx')
