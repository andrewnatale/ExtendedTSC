from __future__ import print_function
import sys, os
import numpy as np
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
            self.grid_center = (0.0, 0.0, 0.0)
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
        # normalize to occupancy per frame
        self.data = self.data / self.frame_count
        # place into a grid object for writing to .dx
        self.gridobj = Grid(self.data, edges=self.edges)
