from __future__ import print_function
import sys, os
import numpy as np
from gridData import Grid
from core.TimeseriesCore import TimeseriesCore
from base.getMDsurfs import _getMDsurfs

class MembraneSurface(TimeseriesCore):

    def run(self, grid_dim, grid_len, surface_sel, midplane_sel, stype='bin' ,grid_center_sel=None):
        if self.primaryDS.framerange is None:
            self.surfer = _getMDsurfs(
              grid_dim,
              grid_len,
              surface_sel,
              midplane_sel,
              self.u,
              stype=stype,
              grid_center_sel=grid_center_sel)
        else:
            self.surfer = _getMDsurfs(
              grid_dim,
              grid_len,
              surface_sel,
              midplane_sel,
              self.u,
              stype=stype,
              grid_center_sel=grid_center_sel,
              start=self.primaryDS.framerange[0],
              stop=self.primaryDS.framerange[1],
              step=self.primaryDS.framerange[2])
        self.surfer.run()
        self._postprocess()

    def write_vmd_surfaces(self):
        self._savevmdmesh('Up')
        self._savevmdmesh('Um')

    def _savevmdmesh(self, target):
        x = self.surfer.X
        y = self.surfer.Y
        if target == 'Up':
            z = self.surfer.Up
        elif target == 'Um':
            z = self.surfer.Um
        with open('%s.tcl' % target, 'w') as meshfile:
            for j in range(self.surfer.grid_dim-1):
                for i in range(self.surfer.grid_dim-1):
                    if np.sum(np.isnan((x[i,j],y[i,j],z[i,j],x[i+1,j],y[i+1,j],z[i+1,j],x[i+1,j+1],y[i+1,j+1],z[i+1,j+1]))) < 1:
                        meshfile.write('draw triangle {%f %f %f} {%f %f %f} {%f %f %f}\n' % \
                          (x[i,j],y[i,j],z[i,j],x[i+1,j],y[i+1,j],z[i+1,j],x[i+1,j+1],y[i+1,j+1],z[i+1,j+1]))
            for j in range(self.surfer.grid_dim-1):
                for i in range(self.surfer.grid_dim-1):
                    if np.sum(np.isnan((x[i,j],y[i,j],z[i,j],x[i,j+1],y[i,j+1],z[i,j+1],x[i+1,j+1],y[i+1,j+1],z[i+1,j+1]))) < 1:
                        meshfile.write('draw triangle {%f %f %f} {%f %f %f} {%f %f %f}\n' % \
                          (x[i,j],y[i,j],z[i,j],x[i,j+1],y[i,j+1],z[i,j+1],x[i+1,j+1],y[i+1,j+1],z[i+1,j+1]))

    def _postprocess(self):
        pass

    def _debug(self):
        pass
