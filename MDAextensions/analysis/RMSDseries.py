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

    def run(self, selection_list, ref_frame=0, pdbname=None, align=False):
        #   def run(self, selecttext, frame=0, pdbname=None, alt_selecttext=None, align=False, featurename='rmsdseries'):
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
            self.primaryDS.rmsd_reference = '%s' % os.path.abspath(pdbname)
        else:
            self.primaryDS.rmsd_reference = 'trajectory frame %d' % frame
        # set vars to control alignment
        if align is True:
            self.center = True
            self.superpose = True
        # parse selection list and calc rmsds
        series_list = []
        for elem in selection_list:
            name, selecttext, alt_selecttext = elem
            descriptor = (name, 'RMSDseries', str(selecttext)+';'+str(alt_selecttext)+';'+'align=='+str(align))
            self.primaryDS.add_feature(descriptor, width=1)
            # target atom selection in trajectory
            tgt = self.u.select_atoms(selecttext)
            if pdbname:
                # reference atom coordinates from pdb
                pdbselect = ref_u.select_atoms(alt_selecttext)
                ref = pdbselect.positions.copy()
                if len(tgt) != len(pdbselect):
                    sys.exit('Atom selections for RMSD calculation are not the same size!\n\
                      PDB selection contains %d atoms, \
                      while trajectory selection contains %d.\n Exiting...' % (len(pdbselect),len(tgt)))
            else:
                # set trajectory to reference frame
                self.u.trajectory[ref_frame]
                # copy reference coordinates from that frame
                ref = tgt.positions.copy()
                # reset trajectory to beginning
                self.u.trajectory.rewind()
            # run rmsd calculation
            if self.primaryDS.framerange is None:
                series = _simpleRMSD(tgt, ref, self.u, self.center, self.superpose, verbose=True)
            else:
                series = _simpleRMSD(tgt, ref, self.u, self.center, self.superpose, verbose=True,
                  start=self.primaryDS.framerange[0],stop=self.primaryDS.framerange[1],step=self.primaryDS.framerange[2])
            series.run()
            series_list.append(series.rmsd)
        self.primaryDS.format_data(np.stack(series_list, axis=0))

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
        # self.rmsd = np.reshape(self.rmsd, (1,-1))
