import sys, os
import numpy as np
# tested and working with MDAnalysis-0.16.1
import MDAnalysis as mda
from TimeseriesCore import TimeseriesCore

class TrajProcessor(TimeseriesCore):
    """Main class for trajectory reading modules that use MDAnalysis as a backend."""

    # supported input types
    valid_data_types = ['dcd_traj','generic_traj','pdb']

    def load_universe(self, universe, traj_stepsize, framerange=None, input_type='generic_traj', toponame=None, trajname=None):
        """Load a preinitialized MDAnalysis Universe object. This can be useful if you need to
    customize Universe generation or want to pass the same Universe to multiple TrajProcessor
    instances for parallel processing.

    WARNING: This needs extra testing to ensure MDA's trajectory API can handle parallel access
    to the same trajectory reader object without exploding!

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
            self.logger.msg('Loading preinitialized Universe object.')
            if toponame:
                self.primaryDS.toponame = toponame
            if trajname:
                self.primaryDS.trajname = trajname
            self.primaryDS.traj_stepsize = traj_stepsize
            self.primaryDS.framerange = framerange
            self.u = universe
            if input_type != 'pdb' and input_type in self.valid_data_types:
                self.input_type = input_type
            else:
                self.logger.exit('Invalid input type for Universe loading! Exiting...')
        else:
            self.logger.err('Can only handle one trajectory or structure per instance! Exiting...')

    def load_dcd(self, psffile, dcdfile, traj_stepsize, framerange=None):
        """Load psf topology and dcd trajectory files. For non dcd format trajectories, use
    load_traj() method instead.

    Arguments:
    psffile - string; path to psf format topology file
    dcdfile_list - string; path to dcd format trajectory file
    traj_stepsize - int; simulation time between frames of input trajectory

    Keyword arguments:
    framerange - tuple; 3 integers (start, stop, and step) describing how to slice the
        trajectory for analysis, default None == (0,-1,1)
    """

        if self.input_type == None:
            self.logger.msg('Attempting to load trajectory using files:\n%s\n%s' % (psffile,dcdfile))
            self.primaryDS.toponame = os.path.abspath(psffile)
            self.primaryDS.trajname = os.path.abspath(dcdfile)
            self.primaryDS.traj_stepsize = traj_stepsize
            self.primaryDS.framerange = framerange
            # init MDA universe
            self.u = mda.Universe(self.primaryDS.toponame,self.primaryDS.trajname,format=u'DCD')
            self.input_type = 'dcd_traj'
        else:
            self.logger.err('Can only handle one trajectory or structure per instance! Exiting...')

    def load_traj(self, topofile, trajfile, traj_stepsize, framerange=None):
        """Load topology and trajectory files (any MDAnalysis supported types).

    Arguments:
    topofile - string; path to topology file
    trajfile - string; path to trajectory file
    traj_stepsize - int; simulation time between frames of input trajectory in ps

    Keyword arguments:
    framerange - tuple; 3 integers (start, stop, and step) describing how to slice the
        trajectory for analysis, default None == (0,-1,1)
    """

        if self.input_type == None:
            self.logger.msg('Attempting to load trajectory using files:\n%s\n%s' % (topofile,trajfile))
            self.primaryDS.toponame = os.path.abspath(topofile)
            self.primaryDS.trajname = os.path.abspath(trajfile)
            self.primaryDS.traj_stepsize = traj_stepsize
            self.primaryDS.framerange = framerange
            # init MDA universe
            self.u = mda.Universe(self.primaryDS.toponame,self.primaryDS.trajname)
            self.input_type = 'generic_traj'
        else:
            self.logger.err('Can only handle one trajectory or structure per instance! Exiting...')

    def load_pdb(self, pdbfile):
        """Load a structure from a PDB file. PDBs are treated as single-frame trajectories.

    Arguments:
    pdbfile - string; path to pdb format file
    """

        if self.input_type == None:
            self.logger.msg('Attempting to load trajectory using files:\n%s' % pdbfile)
            self.primaryDS.pdbname = os.path.abspath(pdbfile)
            # init universe
            self.u = mda.Universe(self.primaryDS.pdbname, format='pdb')
            self.input_type = 'pdb'
        else:
            self.logger.err('Can only handle one trajectory or structure per instance! Exiting...')