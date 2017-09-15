#!/usr/bin/env python2
from __future__ import print_function,division
import sys, os
import numpy as np
# tested and working with MDAnalysis-0.16.1
import MDAnalysis as mda
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.analysis.base import AnalysisBase
from MDAextensions.datatools.TimeseriesCore import TimeseriesCore

class RMSDseries(TimeseriesCore):
    """Calculate RMSDs values for each frame of a trajectory against a specified reference,
    with optional alignment."""

    def __init__(self,**kwargs):
        super(RMSDseries,self).__init__(**kwargs)
        self.center = False
        self.superpose = False
        self.primaryDS.set_static()

    def run(self, selecttext, frame=0, pdbname=None, alt_selecttext=None, align=False):
        """Arguments:
        selecttext - string; MDAnalysis format selection expression; rmsd will be calculated using
            this set of atoms

        Keyword arguments:
        frame - int; trajectory frame to use as a reference structure
        pdbname - string; path to PDB file to use as a reference structure; if specified, this
            reference takes precedence over any frame reference
        alt_selecttext - string; MDAnalysis format selection expression; use if alternate
            expressions are needed for the pdb vs the trajectory topology; the pdb will use
            selecttext and the trajectory will use alt_selecttext
        align - boolean; if True, center and rotate the selection to align with the reference
            before calculating RMSD for each frame
        """

        if self.input_type == None:
            sys.exit('No data has been loaded, cannot run! Exiting...')
        if pdbname:
            # load pdb into its own universe
            ref_u = mda.Universe(pdbname, format='pdb')
            self.primaryDS.rmsd_reference = '%s \"%s\"' % (os.path.abspath(pdbname), selecttext)
            # reference atom selection from pdb
            pdbselect = ref_u.select_atoms(selecttext)
            ref = pdbselect.positions.copy()
            # target atom selection in trajectory - use alt_selecttext variable if given
            if alt_selecttext:
                tgt = self.u.select_atoms(alt_selecttext)
                descriptor = ('rmsdseries', 'RMSD', alt_selecttext)
            else:
                tgt = self.u.select_atoms(selecttext)
                descriptor = ('rmsdseries', 'RMSD', selecttext)
            # selections must match in length
            if len(tgt) != len(pdbselect):
                sys.exit('Atom selections for RMSD calculation are not the same size!\n\
                  PDB selection contains %d atoms, \
                  while trajectory selection contains %d.\n Exiting...' % (len(pdbselect),len(tgt)))
        else:
            self.primaryDS.rmsd_reference = 'frame %d \"%s\"' % (frame, selecttext)
            # target atom selection in trajectory
            tgt = self.u.select_atoms(selecttext)
            # set trajectory to reference frame
            self.u.trajectory[frame]
            # copy reference coordinates from that frame
            ref = tgt.positions.copy()
            # reset trajectory to beginning
            self.u.trajectory.rewind()
            descriptor = ('rmsdseries', 'RMSD', selecttext)
        self.primaryDS.add_feature(descriptor, width=1)
        if align:
            self.center = True
            self.superpose = True
        # run rmsd calculation
        if self.primaryDS.framerange is None:
            series = _simpleRMSD(tgt, ref, self.u, self.center, self.superpose, verbose=True)
        else:
            series = _simpleRMSD(tgt, ref, self.u, self.center, self.superpose, verbose=True,
              start=self.primaryDS.framerange[0],stop=self.primaryDS.framerange[1],step=self.primaryDS.framerange[2])
        series.run()
        self.primaryDS.format_data(series.rmsd)

class _simpleRMSD(AnalysisBase):

    def __init__(self,target,reference,universe,center,superpose,**kwargs):
        super(_simpleRMSD,self).__init__(universe.trajectory,**kwargs)
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
        this_frame_tgt = self.target.positions
        self.rmsd_list.append(rmsd(this_frame_tgt, self.reference, center=self.do_center, superposition=self.do_superpose))

    def _conclude(self):
        # convert to array
        self.rmsd = np.array(self.rmsd_list)
        self.rmsd = np.reshape(self.rmsd, (1,-1))
