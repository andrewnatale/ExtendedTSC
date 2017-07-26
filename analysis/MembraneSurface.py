from __future__ import print_function
import sys
import numpy as np
import scipy.interpolate as scin
# tested and working with MDAnalysis-0.16.1
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from core.TrajProcessor import TrajProcessor

"""
step 1 - take input to define grid; dimensions and center
"""

class MembraneSurface(TimeseriesCore):

class _getMDsurfs(AnalysisBase):

    def __init__(self,grid_dim,grid_len,grid_center,edge_ref,midplane_ref,universe,**kwargs):
        super(_getMDsurfs,self).__init__(universe.trajectory,**kwargs)
        self.target = target
        self.reference = reference
        self.u = universe
        self.do_center = center
        self.do_superpose = superpose

    def _prepare(self):
        # setup vars and data structures
        self.rmsd_list = []

    def _single_frame(self):
        # what to do at each frame
        print self._ts
        this_frame_tgt = self.target.positions
        self.rmsd_list.append(rmsd(this_frame_tgt,self.reference, center=self.do_center, superposition=self.do_superpose))

    def _conclude(self):
        # convert to array
        self.rmsd = np.array(self.rmsd_list)
        self.rmsd = np.reshape(self.rmsd, (1,-1))
