import sys, os
from TimeseriesCore import TimeseriesCore
# tested and working with MDAnalysis-0.16.1
import MDAnalysis as md

class TrajInit(TimeseriesCore):
    """Adds trajectory reading modules to TimeseriesCore using MDAnalysis as a backend."""

    # supported input types
    valid_data_types = ['dcd_traj','generic_traj','pdb']

    def load_dcd(self,psffile,dcdfile,traj_stepsize=500,framerange=None):
        """Load psf topology and dcd trajectory files. For non dcd format trajectories, use
    load_traj() method instead.

    Arguments:
    psffile - string; path to psf format topology file
    dcdfile - string; path to dcd format trajectory file

    Keyword arguments:
    traj_stepsize - int; simulation time between frames of input trajectory
    framerange - tuple; 3 integers (start, stop, and step) describing how to slice the
        trajectory for analysis, default None == (0,-1,1)
    """

        if self.input_type == None:
            self.toponame = os.path.abspath(psffile)
            self.trajname = os.path.abspath(dcdfile)
            self.traj_stepsize = traj_stepsize
            self.framerange = framerange
            # init MDA universe
            self.u = md.Universe(self.toponame,self.trajname,format=u'DCD')
            self.input_type = 'dcd_traj'
        else:
            print 'Can only handle one trajectory or structure per instance! Exiting...'
            sys.exit(1)

    def load_traj(self,topofile,trajfile,traj_stepsize=500,framerange=None):
        """Load topology and trajectory files (any MDAnalysis supported types).

    Arguments:
    topofile - string; path to topology file
    trajfile - string; path to trajectory file

    Keyword arguments:
    traj_stepsize - int; simulation time between frames of input trajectory in ps
    framerange - tuple; 3 integers (start, stop, and step) describing how to slice the
        trajectory for analysis, default None == (0,-1,1)
    """

        if self.input_type == None:
            self.toponame = os.path.abspath(topofile)
            self.trajname = os.path.abspath(trajfile)
            self.traj_stepsize = traj_stepsize
            self.framerange = framerange
            # init MDA universe
            self.u = md.Universe(self.toponame,self.trajname)
            self.input_type = 'generic_traj'
        else:
            print 'Can only handle one trajectory or structure per instance! Exiting...'
            sys.exit(1)

    def load_pdb(self,pdbfile):
        """Load a structure from a PDB file. PDBs are treated as single-frame trajectories.

    Arguments:
    pdbfile - string; path to pdb format file
    """

        if self.input_type == None:
            self.pdbname = os.path.abspath(pdbfile)
            # init universe
            self.u = md.Universe(pdbfile, format='pdb')
            self.input_type = 'pdb'
        else:
            print 'Can only handle one trajectory or structure per instance! Exiting...'
            sys.exit(1)
