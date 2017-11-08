#!/usr/bin/env python2
from __future__ import print_function,division
import sys, os
import numpy as np
import scipy.interpolate as scin
#from copy import copy,deepcopy
from gridData import Grid
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.coordinates.DCD import DCDWriter
from MDAextensions.datatools.TimeseriesCore import TimeseriesCore
from MDAextensions.datatools.CustomErrors import LoadError,AnalysisRuntimeError

class MembraneSurface(TimeseriesCore):
    """
    Calculate membrane surface geometry from MD data using _getMDsurfs and write the data to a
    variety of output formats. Can also write .dcd files containing the a subset of the
    trajectory data with the membrane alignment.
    """

    def run(self, grid_dim, grid_len, surface_sel, midplane_sel, solute_sel='protein'):
        """Setup and call _getMDsurfs.run() to calculate the upper and lower leaflet surfaces."""

        # # executive decision - do averaging every 200ps if using stype='interp'
        # if self.primaryDS.traj_stepsize >= 200:
        #     frequency = 1
        # else:
        #     frequency = 200 // self.primaryDS.traj_stepsize
        frequency = 1
        self.surfer = _getMDsurfs(
          grid_dim,
          grid_len,
          surface_sel,
          midplane_sel,
          self.u,
          stype='interp',
          interp_freq=frequency,
          solute_sel=solute_sel,
          verbose=True,
          start=self.primaryDS.framerange[0],
          stop=self.primaryDS.framerange[1],
          step=self.primaryDS.framerange[2])
        self.surfer.run()
        self._postprocess()

    def write_aligned_subset(self, basename):
        """Write out the coordinates that were aligned during surface calculation to pdb/dcd."""

        self.surfer.u2.trajectory.rewind()
        self.surfer.u2.atoms.write('%s.mem_align.pdb' % basename)
        with DCDWriter('%s.mem_align.dcd' % basename, self.surfer.u2.atoms.n_atoms) as w:
            for ts in self.surfer.u2.trajectory:
                w.write_next_timestep(ts)

    def write_vmd_surfaces(self, basename):
        """Write both upper and lower leaflet surfaces to vmd scripts using default names."""

        self._savevmdmesh('Up', basename)
        self._savevmdmesh('Um', basename)

    def _savevmdmesh(self, target, basename):
        """Write target leaflet to vmd script; Up is upper, Um is lower."""

        # select our x,y,z values
        x = self.surfer.X
        y = self.surfer.Y
        if target == 'Up':
            z = self.surfer.Up
        elif target == 'Um':
            z = self.surfer.Um
        # just iterate over the grid and write tcl commands to draw triangles
        with open('%s.%s.tcl' % (basename, target), 'w') as meshfile:
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
        # nothing to do here
        pass


class _getMDsurfs(AnalysisBase):
    """
    MDAnalysis class for extracting membrane hydrophobic surface profiles from simulation data.
    The input trajectory should be prepared by wrapping and aligning the solute of interest to
    the center of the simulation box. It will then be slightly 'dealigned' in order to orient the
    membrane perpendicular to the z-axis to get the best surface. The aligned coordinates are saved
    and can be written to a new trajectory file if desired.

    NOTE that if you align and save a large number of coordinates beyond what is needed for surface
    calculation (using the 'additional_save_sel' keyword arg) this class will use A LOT of mermory -
    potentially much more than the size of the trajectory on disk.

    TODO: allow centering of grid to arbitrary coordinates or to the center of a selection
    """

    def __init__(self, grid_dim, grid_len, surface_sel, midplane_sel, universe, stype='interp',
      interp_freq=1, solute_sel='protein', **kwargs):
        """
        Arguments:
        grid_dim - interger; number of grid points along each axis
        grid_len - float; size of the grid along one edge, in angstroms
        surface_sel - string; MDAnalysis atom selection expression to define the membrane surface
        midplane_sel - string; MDAnalysis atom selection expression used to define a z-coord
            in the middle of the membrane for separating upper and lower leaflets
        universe - MDAnalysis Universe object; contains the trajectory data to process

        Keyword arguments:
        stype - string; 'interp' or 'bin'; method for surface calculation; 'interp' is
            strongly reccomended as it is much faster and produces better looking surfaces
        interp_freq - interger; if using stype='interp' this controls how often a surface
            will be calculated and saved, in time steps (1 == every step)
        solute_sel - string; MDAnalysis atom selection expression -
            center on this selection and add it to the aligned universe - default 'protein'
        """

        super(_getMDsurfs,self).__init__(universe.trajectory,**kwargs)
        self.grid_dim = grid_dim
        self.grid_len = grid_len
        self.surface_sel = surface_sel
        self.midplane_sel = midplane_sel
        self.u = universe
        self.solute_sel = solute_sel
        self.grid_center = (0.0,0.0)
        self.stype = stype
        self.interp_freq = interp_freq

    def _prepare(self):
        # setup vars and data structures
        # setup grid
        xmin = self.grid_center[0] - 0.5*self.grid_len
        xmax = self.grid_center[0] + 0.5*self.grid_len
        ymin = self.grid_center[1] - 0.5*self.grid_len
        ymax = self.grid_center[1] + 0.5*self.grid_len
        self.xedges,self.xspace = np.linspace(xmin,xmax,num=self.grid_dim+1,endpoint=True,retstep=True)
        self.yedges,self.yspace = np.linspace(ymin,ymax,num=self.grid_dim+1,endpoint=True,retstep=True)
        self.xcenters = np.convolve(self.xedges,np.array([0.5,0.5]),mode='valid')
        self.ycenters = np.convolve(self.yedges,np.array([0.5,0.5]),mode='valid')
        self.X, self.Y = np.meshgrid(self.xcenters,self.ycenters)
        # lists to collect timestep data
        if self.stype == 'interp':
            self.tmp_upper_coords_list = []
            self.tmp_lower_coords_list = []
            self.upper_surf_list = []
            self.lower_surf_list = []
        self.upper_coords_list = []
        self.lower_coords_list = []
        self.write_coords_list = []
        self.frame_count = 0
        # atomgroup setup
        self.solute_group = self.u.select_atoms(self.solute_sel)
        self.midplane_group = self.u.select_atoms(self.midplane_sel)
        self.surface_group = self.u.select_atoms(self.surface_sel)
        self.write_group = self.u.select_atoms('(%s)' % ') or ('.join([self.surface_sel, self.midplane_sel, self.solute_sel]))

    def _single_frame(self):
        # increment the counter
        self.frame_count += 1
        # get solute_sel center of geometry for translational alignment
        solute_cog = np.mean(self.solute_group.positions, axis=0)
        # # deepcopy these arrays
        # midplane_coords = deepcopy(self.midplane_group.positions)
        # surface_coords = deepcopy(self.surface_group.positions)
        # write_coords = deepcopy(self.write_group.positions)
        midplane_coords = self.midplane_group.positions - solute_cog
        surface_coords = self.surface_group.positions - solute_cog
        write_coords = self.write_group.positions - solute_cog
        # init rotator on just the midplane
        r = rotator(midplane_coords)
        # rotate
        rot_midplane_coords = r.rotate(midplane_coords)
        rot_surface_coords = r.rotate(surface_coords)
        rot_write_coords = r.rotate(write_coords)
        # now do surface calculation
        # get midplane_ref average z value
        avg_midplane_z = np.mean(rot_midplane_coords[:,2])
        # adust all coords in z to center the membrane
        rot_surface_coords[:,2] = rot_surface_coords[:,2] - avg_midplane_z
        rot_write_coords[:,2] = rot_write_coords[:,2] - avg_midplane_z
        # sort surface atom coords into upper and lower leaflets
        upper_mask = rot_surface_coords[:,2] > 0
        lower_mask = rot_surface_coords[:,2] < 0
        upper_coords = np.stack([
          rot_surface_coords[:,0][upper_mask],
          rot_surface_coords[:,1][upper_mask],
          rot_surface_coords[:,2][upper_mask] ],
          axis=1
          )
        lower_coords = np.stack([
          rot_surface_coords[:,0][lower_mask],
          rot_surface_coords[:,1][lower_mask],
          rot_surface_coords[:,2][lower_mask] ],
          axis=1
          )
        if self.stype == 'interp':
            # add coords to the tmp lists
            self.tmp_upper_coords_list.append(upper_coords)
            self.tmp_lower_coords_list.append(lower_coords)
            # the interp_freq var controls how often we calculate surfs and reset lists
            if self.frame_count % self.interp_freq == 0:
                upp = np.concatenate(self.tmp_upper_coords_list, axis=0)
                low = np.concatenate(self.tmp_lower_coords_list, axis=0)
                # interpolate a surface for the frame and save it
                tUp_surf = scin.griddata((upp[:,0],upp[:,1]), upp[:,2], (self.X,self.Y), method='cubic')
                tUm_surf = scin.griddata((low[:,0],low[:,1]), low[:,2], (self.X,self.Y), method='cubic')
                # save a boolean array indicating defined surface values
                tUp_ct = ~np.isnan(tUp_surf)
                tUm_ct = ~np.isnan(tUm_surf)
                # turn nan into zero so that the surface averages cleanly
                tUp_surf[np.isnan(tUp_surf)] = 0.0
                tUm_surf[np.isnan(tUm_surf)] = 0.0
                # save surface and counts
                self.upper_surf_list.append((tUp_surf, tUp_ct))
                self.lower_surf_list.append((tUm_surf, tUm_ct))
                # reset tmp lists
                self.tmp_upper_coords_list = []
                self.tmp_lower_coords_list = []
        elif self.stype == 'bin':
            # no need to do anything here, coords will be processed in bulk in _conclude
            pass
        # add coords to the list
        self.upper_coords_list.append(upper_coords)
        self.lower_coords_list.append(lower_coords)
        self.write_coords_list.append(rot_write_coords)
        del r, midplane_coords, surface_coords, write_coords

    def _conclude(self):
        # concatenate coordinate lists into arrays
        Up_arr = np.concatenate(self.upper_coords_list, axis=0)
        Um_arr = np.concatenate(self.lower_coords_list, axis=0)
        # compute histograms (perhaps for filtering out low occupancy grid points or other QC)
        self.Up_hist,_ = np.histogramdd(Up_arr[:,0:2],[self.xedges,self.yedges])
        self.Um_hist,_ = np.histogramdd(Um_arr[:,0:2],[self.xedges,self.yedges])
        # finalize surfaces
        if self.stype == 'interp':
            # if there are unprocessed coords, process those first
            if self.frame_count % self.interp_freq != 0:
                upp = np.concatenate(self.tmp_upper_coords_list, axis=0)
                low = np.concatenate(self.tmp_lower_coords_list, axis=0)
                # interpolate a surface for the frame and save it
                tUp_surf = scin.griddata((upp[:,0],upp[:,1]), upp[:,2], (self.X,self.Y), method='cubic')
                tUm_surf = scin.griddata((low[:,0],low[:,1]), low[:,2], (self.X,self.Y), method='cubic')
                # save a boolean array indicating defined surface values
                tUp_ct = ~np.isnan(tUp_surf)
                tUm_ct = ~np.isnan(tUm_surf)
                # turn nan into zero so that the surface averages cleanly
                tUp_surf[np.isnan(tUp_surf)] = 0.0
                tUm_surf[np.isnan(tUm_surf)] = 0.0
                # save surface and counts
                self.upper_surf_list.append((tUp_surf, tUp_ct))
                self.lower_surf_list.append((tUm_surf, tUm_ct))
            # average interpolated surfaces from each timestep
            Up_surf = np.sum([i[0] for i in self.upper_surf_list], axis=0) / np.sum([i[1] for i in self.upper_surf_list], axis=0)
            Um_surf = np.sum([i[0] for i in self.lower_surf_list], axis=0) / np.sum([i[1] for i in self.lower_surf_list], axis=0)
            # remove areas where there was never a count based on histogram
            # NOTE: this is problematic - it eliminates parts of the surface that I can clearly
            # see should have counts by watching the trajectory; for now leave it disabled
            # Up_surf[self.Up_hist==0] = np.nan
            # Um_surf[self.Um_hist==0] = np.nan
        elif self.stype == 'bin':
            # sort coordinates into bins and then cacluate the average value
            # NOTE: this is slow and the surfaces are a lot rougher than when using interpolation
            Up_surf = bin_surf(self.grid_dim, self.xedges, self.yedges, Up_arr)
            Um_surf = bin_surf(self.grid_dim, self.xedges, self.yedges, Um_arr)
        # save surfaces for access
        self.Up = Up_surf
        self.Um = Um_surf
        # this doesn't quite work as advertised in the MDA documentation (see Merge and MemoryReader for details)
        self.u.trajectory.rewind()
        self.u2 = mda.Merge(self.write_group)
        self.u2.transfer_to_memory()
        self.u2.load_new(np.stack(self.write_coords_list, axis=0))

# get a surface by constructing bins and taking a simple average
def bin_surf(grid_dim, xedges, yedges, points):
    tmp_surf = np.empty((grid_dim,grid_dim), dtype=float)
    tmp_surf.fill(np.nan)
    for i in range(grid_dim):
        for j in range(grid_dim):
            #print('binning %d %d' % (i,j))
            # select target bin edges
            x_low = xedges[i]
            x_hig = xedges[i+1]
            y_low = yedges[j]
            y_hig = yedges[j+1]
            # sort coordinates
            in_xbin = np.logical_and(points[:,0]>x_low, points[:,0]<x_hig)
            in_ybin = np.logical_and(points[:,1]>y_low, points[:,1]<y_hig)
            in_bin = np.logical_and(in_xbin,in_ybin)
            # calculate mean z for coords in bin
            tmp_surf[i,j] = np.mean(points[:,2][in_bin])
    return tmp_surf

# calculates a rotation on an input set of coordinates, then can apply that rotation to other sets
# specifically, the rotation aligns the lowest magnitude principle component of the input set
# to the z-axis, while not allowing any rotation about z
class rotator(object):

    def __init__(self, coordinates, z_rotation=0.0):
        """
        Arguments:
        coordinates - array; set of coordinates to use to calculate the rotation matrix for the
            membrane alignment
        z_rotation - float; angle in degrees; the rotation matrix will be modified to add this
            rotation about the z axis in addition to the calculated rotations about x and y
        """

        # transpose and get covariance matrix
        cov = np.cov(coordinates.T)
        # get eigenvectors, values
        evals, evecs = np.linalg.eig(cov)
        # print('starting evecs:')
        # print(evecs)
        # find the eigenvector corresponding to the smallest eigenvalue - this will be normal
        # to the plane of the membrane
        sort_indices = np.argsort(evals)
        small_idx = sort_indices[0]
        x_vs1,y_vs1,z_vs1 = evecs[:,small_idx]
        # make sure it points in the positive z direction so we don't flip the box
        if z_vs1 < 0:
            # print('flipping!')
            x_vs1 = x_vs1 * -1.0
            y_vs1 = y_vs1 * -1.0
            z_vs1 = z_vs1 * -1.0
        # setup rotation matrix - rotation about x
        gamma = -1.0 * np.arctan(y_vs1 / z_vs1)
        # print('gamma:', gamma)
        Rx = np.matrix([[1.0,            0.0,           0.0],
                        [0.0,   np.cos(gamma), np.sin(gamma)],
                        [0.0,  -np.sin(gamma), np.cos(gamma)]])
        Rx_evecs = Rx * evecs
        # print('Rx applied to evecs:')
        # print(Rx_evecs)
        x_vs2,y_vs2,z_vs2 = Rx_evecs.A[:,small_idx]
        # rotation about y
        beta = np.arctan(x_vs2 / z_vs2)
        # print('beta:', beta)
        Ry = np.matrix([[np.cos(beta),   0.0,  -np.sin(beta)],
                        [0.0,            1.0,           0.0],
                        [np.sin(beta),   0.0,  np.cos(beta)]])
        Rxy_evecs = Ry * Rx_evecs
        # print('Ry applied to Rx_evecs:')
        # print(Rxy_evecs)
        # specified z-rotation
        alpha = np.deg2rad(z_rotation)
        Rz = np.matrix([[np.cos(alpha),  np.sin(alpha), 0.0],
                        [-np.sin(alpha), np.cos(alpha), 0.0],
                        [0.0,            0.0,           1.0]])
        # save the rotation matix
        self.Rxyz = Rx * Ry * Rz

    def rotate(self, target):
        """Rotate another set of coordinates by the calculated rotation."""
        Rxyz_target = self.Rxyz * target.T
        return Rxyz_target.A.T
