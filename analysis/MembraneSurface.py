from __future__ import print_function
import sys, os
import numpy as np
from gridData import Grid
from MDAnalysis.coordinates.DCD import DCDWriter
from core.TimeseriesCore import TimeseriesCore
from base.getMDsurfs import _getMDsurfs

class MembraneSurface(TimeseriesCore):
    """
    Calculate membrane surface geometry from MD data using _getMDsurfs and write the data to a
    variety of output formats. Can also write .dcd files containing the a subset of the
    trajectory data with the membrane alignment.
    """

    def run(self, grid_dim, grid_len, surface_sel, midplane_sel, additional_save_sel='protein'):
        """Setup and call _getMDsurfs.run() to calculate the upper and lower leaflet surfaces."""

        # executive decision - do averaging every 200ps if using stype='interp'
        if self.primaryDS.traj_stepsize >= 200:
            frequency = 1
        else:
            frequency = 200 // self.primaryDS.traj_stepsize
        self.surfer = _getMDsurfs(
          grid_dim,
          grid_len,
          surface_sel,
          midplane_sel,
          self.u,
          stype='interp',
          interp_freq=1,
          additional_save_sel=additional_save_sel,
          verbose=True,
          start=self.primaryDS.framerange[0],
          stop=self.primaryDS.framerange[1],
          step=self.primaryDS.framerange[2])
        self.surfer.run()
        self._postprocess()

    def write_aligned_subset(self, filename):
        """Write out the coordinates that were aligned during surface calculation to pdb/dcd."""

        self.surfer.u2.trajectory.rewind()
        self.surfer.u2.atoms.write('%s.pdb' % filename)
        with DCDWriter('%s.dcd' % filename, self.surfer.u2.atoms.n_atoms) as w:
            for ts in self.surfer.u2.trajectory:
                w.write_next_timestep(ts)

    def write_vmd_surfaces(self):
        """Write both upper and lower leaflet surfaces to vmd scripts using default names."""

        self._savevmdmesh('Up')
        self._savevmdmesh('Um')

    def _savevmdmesh(self, target):
        """Write target leaflet to vmd script; Up is upper, Um is lower."""

        # select our x,y,z values
        x = self.surfer.X
        y = self.surfer.Y
        if target == 'Up':
            z = self.surfer.Up
        elif target == 'Um':
            z = self.surfer.Um
        # just iterate over the grid and write tcl commands to draw triangles
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
