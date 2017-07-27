from __future__ import print_function
import numpy as np
import scipy.interpolate as scin
from statsmodels.robust.scale import mad
from MDAnalysis.analysis.base import AnalysisBase

class _getMDsurfs(AnalysisBase):

    def __init__(self, grid_dim, grid_len, surface_sel, midplane_sel, universe,
      stype='bin', grid_center_sel=None, **kwargs):
        """
        Arguments:
        grid_dim - integer; number of grid points along each axis
        grid_len - float; size of the grid along one edge, in angstroms
        surface_sel - string; MDAnalysis atom selection expression to define the membrane surface
        midplane_sel - string; MDAnalysis atom selection expression used to define a z-coord
            in the middle of the membrane for separating upper and lower leaflets
        universe - MDAnalysis Universe object; trajectory data to process

        Keyword arguments:
        stype - string; 'interp' or 'bin'; method for surface calculation
        grid_center_sel - string; MDAnalysis atom selection expression; each frame will be
        """
        super(_getMDsurfs,self).__init__(universe.trajectory,**kwargs)
        self.grid_dim = grid_dim
        self.grid_len = grid_len
        self.surface_sel = surface_sel
        self.midplane_sel = midplane_sel
        self.u = universe
        if grid_center_sel:
            # do something
            #self.grid_center = ???
            pass
        else:
            self.grid_center = (0.0,0.0)
        self.stype = stype

    def _prepare(self):
        # setup vars and data structures
        # setup grid
        xmin = self.grid_center[0]-0.5*self.grid_len
        xmax = self.grid_center[0]+0.5*self.grid_len
        ymin = self.grid_center[1]-0.5*self.grid_len
        ymax = self.grid_center[1]+0.5*self.grid_len
        self.xedges,self.xspace = np.linspace(xmin,xmax,num=self.grid_dim+1,endpoint=True,retstep=True)
        self.yedges,self.yspace = np.linspace(ymin,ymax,num=self.grid_dim+1,endpoint=True,retstep=True)
        self.xcenters = np.convolve(self.xedges,np.array([0.5,0.5]),mode='valid')
        self.ycenters = np.convolve(self.yedges,np.array([0.5,0.5]),mode='valid')
        self.X, self.Y = np.meshgrid(self.xcenters,self.ycenters)
        # setup atomgroups
        self.midplane_group = self.u.select_atoms(self.midplane_sel)
        self.surface_group = self.u.select_atoms(self.surface_sel)
        # lists to collect timestep data
        if self.stype == 'interp':
            self.upper_surf_list = []
            self.lower_surf_list = []
        self.upper_coords_list = []
        self.lower_coords_list = []
        self.frame_count = 0

    def _single_frame(self):
        # what to do at each frame
        print(self._ts)
        # copy and transpose coord arrays so they're shaped like what numpy cov and linalg.eig need
        midplane_coords = self.midplane_group.positions.copy()
        midplane_coords = midplane_coords.T
        surface_coords = self.surface_group.positions.copy()
        surface_coords = surface_coords.T
        # calculate membrane alignment rotation matrix
        membrane_coords = np.concatenate([midplane_coords,surface_coords], axis=1)
        # # get center coords
        # x_avg = np.mean(membrane_coords[0,:])
        # y_avg = np.mean(membrane_coords[1,:])
        # z_avg = np.mean(membrane_coords[2,:])
        # # shift to center
        # membrane_coords[0,:] = membrane_coords[0,:] - x_avg
        # membrane_coords[1,:] = membrane_coords[1,:] - y_avg
        # membrane_coords[2,:] = membrane_coords[2,:] - z_avg
        # get covariance matrix
        cov = np.cov(membrane_coords)
        # get eigenvectors
        evals, evecs = np.linalg.eig(cov)
        # find the eigenvector corresponding to the smallest eigenvalue
        sort_indices = np.argsort(evals)
        small_idx = sort_indices[0]
        x_vs,y_vs,z_vs = evecs[:,small_idx]
        # make sure it points in the positive z direction so we don't flip the box
        if z_vs < 0:
            x_vs = x_vs * -1.0
            y_vs = y_vs * -1.0
            z_vs = z_vs * -1.0
        # setup rotation matrix; only rotate about x and y, not z
        alpha = 0.0
        beta = np.arccos(z_vs / (np.sqrt(x_vs**2 + z_vs**2)))
        gamma = np.arccos(z_vs / (np.sqrt(y_vs**2 + z_vs**2)))
        Rz = np.matrix([[np.cos(alpha), np.sin(alpha), 0.0],
                        [-np.sin(alpha), np.cos(alpha), 0.0],
                        [0.0, 0.0, 1.0]])
        Ry = np.matrix([[np.cos(beta), 0.0, -np.sin(beta)],
                        [0.0, 1.0, 0.0],
                        [np.sin(beta), 0.0, np.cos(beta)]])
        Rx = np.matrix([[1.0, 0.0, 0.0],
                        [0.0, np.cos(gamma), np.sin(gamma)],
                        [0.0, -np.sin(gamma), np.cos(gamma)]])
        Rzyx = Rz*Ry*Rx
        # shift, rotate, shift back, and transpose back
        # midplane_coords[0,:] = midplane_coords[0,:] - x_avg
        # midplane_coords[1,:] = midplane_coords[1,:] - y_avg
        # midplane_coords[2,:] = midplane_coords[2,:] - z_avg
        # surface_coords[0,:] = surface_coords[0,:] - x_avg
        # surface_coords[1,:] = surface_coords[1,:] - y_avg
        # surface_coords[2,:] = surface_coords[2,:] - z_avg
        rot_midplane_coords = Rzyx * midplane_coords
        rot_surface_coords = Rzyx * surface_coords
        rot_midplane_coords = rot_midplane_coords.A
        rot_surface_coords = rot_surface_coords.A
        # rot_midplane_coords[0,:] = rot_midplane_coords[0,:] + x_avg
        # rot_midplane_coords[1,:] = rot_midplane_coords[1,:] + y_avg
        # rot_midplane_coords[2,:] = rot_midplane_coords[2,:] + z_avg
        # rot_surface_coords[0,:] = rot_surface_coords[0,:] + x_avg
        # rot_surface_coords[1,:] = rot_surface_coords[1,:] + y_avg
        # rot_surface_coords[2,:] = rot_surface_coords[2,:] + z_avg
        rot_midplane_coords = rot_midplane_coords.T
        rot_surface_coords = rot_surface_coords.T
        # get midplane_ref average z value
        avg_midplane_z = np.mean(rot_midplane_coords[:,2])
        # sort surface atom coords into upper and lower leaflets
        # could pre-sort, but in theory, lipids could flip over the course of trajectory
        upper_mask = rot_surface_coords[:,2] > avg_midplane_z
        lower_mask = rot_surface_coords[:,2] < avg_midplane_z
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
            Up_surf = scin.griddata((upper_coords[:,0],upper_coords[:,1]), upper_coords[:,2], (self.X,self.Y), method='cubic')
            Um_surf = scin.griddata((lower_coords[:,0],lower_coords[:,1]), lower_coords[:,2], (self.X,self.Y), method='cubic')
            Up_ct = ~np.isnan(Up_surf)
            Um_ct = ~np.isnan(Um_surf)
            # are there any nans after interpolation?
            print(np.sum(np.isnan(Up_surf))+np.sum(np.isnan(Um_surf)))
            Up_surf[np.isnan(Up_surf)] = 0.0
            Um_surf[np.isnan(Um_surf)] = 0.0
            self.upper_surf_list.append((Up_surf, Up_ct))
            self.lower_surf_list.append((Um_surf, Um_ct))
        elif self.stype == 'bin':
            pass
        # add coords to the list
        self.upper_coords_list.append(upper_coords)
        self.lower_coords_list.append(lower_coords)
        self.frame_count += 1

    def _conclude(self):
        # concatenate lists into arrays
        Up_arr = np.concatenate(self.upper_coords_list, axis=0)
        Um_arr = np.concatenate(self.lower_coords_list, axis=0)
        # compute histograms (perhaps for filtering out low occupancy grid points)
        self.Up_hist,_ = np.histogramdd(Up_arr[:,0:2],[self.xedges,self.yedges])
        self.Um_hist,_ = np.histogramdd(Um_arr[:,0:2],[self.xedges,self.yedges])
        # either do a simple average of everything in each bin, or use interpolation
        if self.stype == 'interp':
            # average interpolated surfaces from each timestep
            Up_surf = np.sum([i[0] for i in self.upper_surf_list], axis=0) / np.sum([i[1] for i in self.upper_surf_list], axis=0)
            Um_surf = np.sum([i[0] for i in self.lower_surf_list], axis=0) / np.sum([i[1] for i in self.lower_surf_list], axis=0)
        elif self.stype == 'bin':
            Up_surf = np.empty((self.grid_dim,self.grid_dim), dtype=float)
            Um_surf = np.empty((self.grid_dim,self.grid_dim), dtype=float)
            Up_surf.fill(np.nan)
            Um_surf.fill(np.nan)
            for i in range(self.grid_dim):
                for j in range(self.grid_dim):
                    print('binning %d %d' % (i,j))
                    # select target bin edges
                    x_low = self.xedges[i]
                    x_hig = self.xedges[i+1]
                    y_low = self.yedges[j]
                    y_hig = self.yedges[j+1]
                    # sort coordinates
                    Up_in_xbin = np.logical_and(Up_arr[:,0]>x_low, Up_arr[:,0]<x_hig)
                    Up_in_ybin = np.logical_and(Up_arr[:,1]>y_low, Up_arr[:,1]<y_hig)
                    Um_in_xbin = np.logical_and(Um_arr[:,0]>x_low, Um_arr[:,0]<x_hig)
                    Um_in_ybin = np.logical_and(Um_arr[:,1]>y_low, Um_arr[:,1]<y_hig)
                    Up_in_bin = np.logical_and(Up_in_xbin,Up_in_ybin)
                    Um_in_bin = np.logical_and(Um_in_xbin,Um_in_ybin)
                    # calculate median z for coords in bin
                    Up_surf[i,j] = np.mean(Up_arr[:,2][Up_in_bin])
                    Um_surf[i,j] = np.mean(Um_arr[:,2][Um_in_bin])
        # remove areas where there was never a count
        Up_surf[self.Up_hist==0] = np.nan
        Um_surf[self.Um_hist==0] = np.nan
        # save surfaces for access
        self.Up = Up_surf
        self.Um = Um_surf
        #return self.Up, self.Um, self.Up_hist, self.Um_hist, self.X, self.Y
