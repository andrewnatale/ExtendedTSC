#!/usr/bin/env python2
from __future__ import print_function,division
import sys, os, datetime
import numpy as np
# tested and working with MDAnalysis-0.16.1
import MDAnalysis as mda
from MDAextensions.datatools.TimeseriesDataSet import TimeseriesDataSet as tds
from MDAextensions.datatools.CustomErrors import LoadError

class TimeseriesCore(object):
    """
    Base class for using MDAnalysis to analyze trajectory data and store it in a database-like
    manner.
    """

    # supported input types
    valid_data_types = ['generic_traj','pdb']

    def __init__(self):
        """
        Setup empty default DataSet objects. Subclasses should use these unless there is a good
        reason to customize.
        """

        self.input_type = None # set later to: 'dcd_traj', 'generic_traj', or 'pdb'
        self.primaryDS = tds()
        self.primaryDS.data_datetime = str(datetime.datetime.now())
        self.primaryDS.data_hostname = os.uname()[1]
        self.maskDS = tds()
        self.maskDS.is_mask = True

    def _custom_dataset(self):
        """Return an initialized DataSet object for customization."""

        custom_dataset = tds()
        custom_dataset.data_datetime = str(datetime.datetime.now())
        custom_dataset.data_hostname = os.uname()[1]
        return custom_dataset

    def write_data(self, fileprefix):
        """
        Writes data for populated default data sets to .dat files in the cwd.

        Arguments:
        fileprefix - string; file names without extension (.dat or .mask.dat
            will be appended automatically)
        """

        if self.primaryDS.populated:
            self.primaryDS._write(outfilename=fileprefix+'.dat')
        if self.maskDS.populated:
            self.maskDS._write(outfilename=fileprefix+'.mask.dat')

    def load_universe(self, universe, traj_stepsize, framerange=None, input_type='generic_traj', toponame=None, trajname=None):
        """
        Load a preinitialized MDAnalysis Universe object. This can be useful if you need to
        customize Universe generation.

        Arguments:
        universe - MDAnalysis Universe object
        traj_stepsize - int; simulation time between frames of input trajectory

        Keyword arguments:
        framerange - tuple; 3 integers (start, stop, and step) describing how to slice the
            trajectory for analysis, default None == (0,-1,1)
        input_type - string; default is generic_traj, which will work for most things, but special types
            like 'dcd_traj' can be passed to let ExtendedTSC know to use optimizations.
            NOTE: PDB loading is not allowed through this method!
        toponame - string; path to topology file - not processed here, just for bookkeeping
        trajname - string; path to trajectory file - not processed here, just for bookkeeping
        """

        if self.input_type == None:
            print('Loading preinitialized Universe object.')
            if toponame:
                self.primaryDS.toponame = toponame
            if trajname:
                self.primaryDS.trajname = trajname
            self.primaryDS.traj_stepsize = traj_stepsize
            # parse framerange - this is an ugly kludge, I blame the AnalysisBase API
            # don't forget, if you mess with this, make cooresponding changes to the _write
            # method in TimeseriesDataSet
            if framerange is None:
                self.primaryDS.framerange = (0, -1, 1)
            else:
                self.primaryDS.framerange = framerange
            if self.primaryDS.framerange[1] == -1:
                self.primaryDS.framerange = (self.primaryDS.framerange[0], None, self.primaryDS.framerange[2])
            self.u = universe
            self.input_type = input_type
        else:
            raise LoadError(0)

    def load_traj(self, topofile, trajfile, traj_stepsize, framerange=None):
        """
        Load topology and trajectory files (any MDAnalysis supported types).

        Arguments:
        topofile - string; path to topology file
        trajfile - string; path to trajectory file
        traj_stepsize - int; simulation time between frames of input trajectory in ps

        Keyword arguments:
        framerange - tuple; 3 integers (start, stop, and step) describing how to slice the
            trajectory for analysis, default None == (0,-1,1)
        """

        if self.input_type == None:
            print('Attempting to load trajectory using files:\n%s\n%s' % (topofile,trajfile))
            self.primaryDS.toponame = os.path.abspath(topofile)
            self.primaryDS.trajname = os.path.abspath(trajfile)
            self.primaryDS.traj_stepsize = traj_stepsize
            # parse framerange - this is an ugly kludge, I blame the AnalysisBase API
            # don't forget, if you mess with this, make cooresponding changes to the .dat
            # writer method in core.TimeseriesDataSet
            if framerange is None:
                self.primaryDS.framerange = (0, -1, 1)
            else:
                self.primaryDS.framerange = framerange
            if self.primaryDS.framerange[1] == -1:
                self.primaryDS.framerange = (self.primaryDS.framerange[0], None, self.primaryDS.framerange[2])
            # init MDA universe
            self.u = mda.Universe(self.primaryDS.toponame,self.primaryDS.trajname)
            self.input_type = 'generic_traj'
        else:
            raise LoadError(0)

    def load_pdb(self, pdbfile):
        """
        Load a structure from a PDB file. PDBs are treated as single-frame trajectories.

        Arguments:
        pdbfile - string; path to pdb format file
        """

        if self.input_type == None:
            print('Attempting to load trajectory using files:\n%s' % pdbfile)
            self.primaryDS.pdbname = os.path.abspath(pdbfile)
            # init universe
            self.u = mda.Universe(self.primaryDS.pdbname, format='pdb')
            self.input_type = 'pdb'
        else:
            raise LoadError(0)
