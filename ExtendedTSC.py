# ExtendedTSC.py - a class to wrap and extend MDAnalysis' TimeseriesCollection module.
# written by Andrew Natale
# probably (definitely) only works on MacOS/Linux/Unix, python 2.7.x, MDAnalysis-0.15.0
# TODO:
#   - allow loading of multiple traj files which form one continuous trajectory
#   - squish bugs
#   - move to MDAnalysis-0.16.1
#   - implement an equivalent of TimeseriesCollection for non psf/dcd format trajectories
#   - modularize for different filetypes
#   - raise python errors instead of calling sys.exit()

import sys, os, math, datetime
import numpy as np
# tested and working with MDAnalysis-0.16.1
import MDAnalysis as md
import MDAnalysis.core.Timeseries as tm
from MDAnalysis.analysis.base import AnalysisBase

class ExtendedTSC(object):
    """Class to wrap and extend MDAnalysis' TimeseriesCollection module."""

    # version info
    version = '0.6.0'

    def __init__(self,datfile=None,maskfile=None):
        """Sets up an ExtendedTSC object, either empty for making measurements
        (takes no init args) or with the names of datfiles to be read."""

        # collect some system data for keeping track of dat files
        self.data_datetime = datetime.datetime.now()
        self.data_hostname = os.uname()[1]
        # setup empty vars
        self.input_type = None # 'dat' 'traj' 'pdb'
        self.toponame = None
        self.trajname = None
        self.traj_stepsize = None # in picoseconds
        self.pdbname = None
        self.framerange = None # all frames, no skipping
        self.primaryDS = _DataSet()
        self.maskDS = _DataSet()
        self.read_mask = False
        # loading dat file(s) - requires no other setup
        if datfile:
            self._data_reader(datfile)
            self.input_type = 'dat'
        if maskfile:
            self._data_reader(maskfile)
            self.input_type = 'dat'

    def load_dcd(self,psffile,dcdfile,traj_stepsize=500,framerange=None):
        """Loads psf topology and dcd trajectory files and initializes MDAnalysis."""

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
        """Loads topology and trajectory files and initializes MDAnalysis.
        This is a generic method which should work with any trajectory types
        supported by MDAnalysis without extra options."""

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

    def load_PBD(self,pdbfile):
        """Loads a structure from a PDB file and initializes MDAnalysis.
        PDBs are treated as though they are a trajectory with a single frame."""

        if self.input_type == None:
            self.pdbname = os.path.abspath(pdbfile)
            # init universe
            self.u = md.Universe(pdbfile, format='pdb')
            self.input_type = 'pdb'
        else:
            print 'Can only handle one trajectory or structure per instance! Exiting...'
            sys.exit(1)

    def measures_from_list(self,selection_list):
        """Loads selection definitions as a list of 3-element tuples: (name, type, selectext).
        This method can be used alone or alongside measures_from_volumesearch to add reference
        measurements to be made alongside the measurements defined by the volume search."""

        # setup primaryDS using selections from a list
        for descriptor in selection_list:
            self.primaryDS.add_measurement(descriptor)

    def measures_from_volumesearch(self,vol_selecttext,search_selecttext,mode='res'):
        """Uses MDAnalysis to do an inital step through the trajectory to identify
        all the atoms/residues that enter a defined volume region in any frame,
        and build a list to measure their coordinates later.
        Records at each timestep whether the atom/reside is present in the volume and saves this
        info into a secondary mask _DataSet object which can be written to a file later.
        """
        # TODO: input type check, shouldn't process pdbfiles
        if self.framerange is None:
            searcher = _VolumeSearch(vol_selecttext,search_selecttext,mode,self.u)
        else:
            searcher = _VolumeSearch(vol_selecttext,search_selecttext,mode,self.u,start=self.framerange[0],stop=self.framerange[1],step=self.framerange[2])
        print searcher.start,searcher.stop,searcher.step,searcher.n_frames
        searcher.run()
        # setup data sets using search results
        for descriptor in searcher.selections:
            self.primaryDS.add_measurement(descriptor)
        for descriptor in searcher.selections_mask:
            # occupancy measure types must have width=1
            self.maskDS.add_measurement(descriptor,width=1)
        # load the masking data from the search
        self.maskDS.add_collection(searcher.mask)

    def run(self):
        # replace generate_timeseries with an input format agnostic method that calls on different timeseries creators
        # and only needs to hand the output arrays
        # step 1: init checks - on DataSet and measurements
        # step 2: query traj format and choose method
        # step 3: pass measuremnt list to a timeseries method and run
        # step 4: format returned array into primaryDS

        # check that the primary DataSet object has been properly initialized
        if not self.primaryDS:
            print 'DataSet not properly initialized! Exiting...'
            sys.exit(1)
        elif self.primaryDS.count_measurements() == 0:
            print 'No measurement descriptors found in DataSet. Exiting...'
            sys.exit(1)
        elif self.primaryDS.populated is True:
            print 'DataSet object already contains an array, cannot generate another! Exiting...'
            sys.exit(1)
        # check that measurement objects behave as expected on the input topology, and set measure widths
        for meas in self.primaryDS.measurements:
            if meas.type == 'atom':
                tmpselect = self.u.select_atoms(meas.selecttext)
                # strictly speaking, TimeseriesCollection can handle having more than one atom
                # in an 'atom' selection - it just measures coordinates for all of them
                # however to keep things simpler downstream, enforce one atom selections here
                if tmpselect.n_atoms != 1:
                    print 'Selection %s \"%s\" found %d atoms instead of 1!\nFix the selection and try again. Exiting now...' \
                      % (mas.type, meas.selecttext, tmpselect.n_atoms)
                    sys.exit(1)
                meas.set_width(3)
            elif meas.type == 'bond':
                # not implemented
                print 'measurement type %s not implemented, exiting' % meas.type
                sys.exit(1)
            elif meas.type == 'angle':
                # not implemented
                print 'measurement type %s not implemented, exiting' % meas.type
                sys.exit(1)
            elif meas.type == 'dihedral':
                tmpselect = self.u.select_atoms(meas.selecttext)
                if tmpselect.n_atoms != 4:
                    print 'Selection %s \"%s\" found %d atoms instead of 4!\nFix the selection and try again. Exiting now...' \
                      % (mas.type, meas.selecttext, tmpselect.n_atoms)
                    sys.exit(1)
                meas.set_width(1)
            elif meas.type == 'distance':
                tmpselect = self.u.select_atoms(meas.selecttext)
                if tmpselect.n_atoms != 2:
                    print 'Selection %s \"%s\" found %d atoms instead of 2!\nFix the selection and try again. Exiting now...' \
                      % (mas.type, meas.selecttext, tmpselect.n_atoms)
                    sys.exit(1)
                meas.set_width(1)
            elif meas.type == 'COG':
                tmpselect = self.u.select_atoms(meas.selecttext)
                if tmpselect.n_atoms == 0:
                    print 'Selection %s \"%s\" found 0 atoms!\nFix the selection and try again. Exiting now...'\
                      % (meas.type, meas.selecttext)
                    sys.exit(1)
                meas.set_width(3)
            elif meas.type == 'COM':
                tmpselect = self.u.select_atoms(meas.selecttext)
                if tmpselect.n_atoms == 0:
                    print 'Selection %s \"%s\" found 0 atoms!\nFix the selection and try again. Exiting now...'\
                      % (meas.type, meas.selecttext)
                    sys.exit(1)
                meas.set_width(3)
            elif meas.type == 'water_dipole':
                # not implemented
                print 'measurement type %s not implemented, exiting' % meas.type
                sys.exit(1)
            else:
                print 'unrecognized measure type %s! Exiting...' % meas.type
                sys.exit(1)
        # check input data format and call the appropriate method to generate the array
        if self.input_type == 'dcd_traj':
            tmp_array = self._DCD_timeseries()
        elif self.input_type == 'pdb':
            tmp_array = self._generic_timeseries()
        elif self.input_type == 'generic_traj':
            tmp_array = self._PDB_measure()
        else:
            sys.exit(1)
        # process the array into a DataSet object
        self.primaryDS.add_collection(tmp_array)
        # create an array of time values (in ps) for plotting and load it into datasets
        starttime = self.framerange[0] * self.traj_stepsize
        endtime = self.framerange[0] * self.traj_stepsize + self.traj_stepsize * self.framerange[2] * (np.shape(tmp_array)[1] - 1)
        self.primaryDS.add_timesteps(np.linspace(float(starttime), float(endtime), num=np.shape(tmp_array)[1]))
        if self.maskDS.populated:
            self.maskDS.add_timesteps(np.linspace(float(starttime), float(endtime), num=np.shape(tmp_array)[1]))
            # compare the maskDS data array to the time array and trim if needed (see NOTE above)
            print 'primary',np.shape(self.primaryDS.data)[1], np.shape(self.primaryDS.time)[0]
            print 'mask',np.shape(self.maskDS.data)[1], np.shape(self.maskDS.time)[0]
            if (np.shape(self.maskDS.data)[1] - np.shape(self.maskDS.time)[0]) == 1:
                print 'trimming...'
                self.maskDS.data = np.delete(self.maskDS.data, -1, axis=1)
            elif ((np.shape(self.maskDS.data)[1] - np.shape(self.maskDS.time)[0]) > 1) or \
                 ((np.shape(self.maskDS.data)[1] - np.shape(self.maskDS.time)[0]) < 0):
                print 'Something is wrong, DataSet arrays do not match (probably a bug). Exiting...'
                sys.exit(1)

    def _DCD_timeseries(self):
        collection = tm.TimeseriesCollection()
        # all the error checking should be done, so assume it's good and load it up

        for meas in self.primaryDS.measurements:
            if meas.type == 'atom':
                collection.addTimeseries(tm.Atom('v', self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'bond':
                # not implemented
                pass
            elif meas.type == 'angle':
                # not implemented
                pass
            elif meas.type == 'dihedral':
                collection.addTimeseries(tm.Dihedral(self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'distance':
                collection.addTimeseries(tm.Distance('r', self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'COG':
                collection.addTimeseries(tm.CenterOfGeometry(self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'COM':
                collection.addTimeseries(tm.CenterOfMass(self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'water_dipole':
                # not implemented
                pass
        # compute Timeseries from universe trajectory and return the resulting array
        # NOTE: There's a bug in MDAnalysis' fast DCD reader code that causes it to drop the last
        # frame it should measure when slicing the trajectory in any way (not the same bug as noted above in
        # 'measures_from_volumesearch'). Haven't tracked it down so we just work around it.
        collection.compute(self.u.trajectory, start=self.framerange[0], stop=self.framerange[1], step=self.framerange[2])
        return collection.data

    def _generic_timeseries(self):
        if self.framerange is None:
            collection = _GenericTSC(self.primaryDS.measurements,self.u)
        else:
            collection = _GenericTSC(self.primaryDS.measurements,self.u,start=self.framerange[0],stop=self.framerange[1],step=self.framerange[2])
        collection.run()
        return collection.data

    def _PDB_measure(self):
        collection = _ProcPDB(self.primaryDS.measurements)
        collection.add_universe(self.u.trajectory)
        collection.run()
        return collection.data

    def water_search(self,vol_selecttext):

        if self.framerange is None:
            searcher = _WaterSearchZ(vol_selecttext,self.u)
        else:
            searcher = _WaterSearchZ(vol_selecttext,self.u,start=self.framerange[0],stop=self.framerange[1],step=self.framerange[2])
        print searcher.start,searcher.stop,searcher.step,searcher.n_frames
        searcher.run()
        # setup data set using search results
        for descriptor in searcher.selections:
            self.primaryDS.add_measurement(descriptor)
        self.primaryDS.add_collection(searcher.coordarray)
        starttime = self.framerange[0] * self.traj_stepsize
        endtime = self.framerange[0] * self.traj_stepsize + self.traj_stepsize * self.framerange[2] * (np.shape(searcher.coordarray)[1] - 1)
        self.primaryDS.add_timesteps(np.linspace(float(starttime), float(endtime), num=np.shape(searcher.coordarray)[1]))
        self.primaryDS.measurements[0].set_width(np.shape(searcher.coordarray)[0])

    def write_data(self,fileprefix):
        """Writes data for all populated DataSets to .dat files in the cwd."""

        if self.primaryDS.populated:
            self._data_writer(fileprefix+'.dat',self.primaryDS)
        if self.maskDS.populated:
            self._data_writer(fileprefix+'.mask.dat',self.maskDS,mask=True)

    def _data_writer(self,outfile,dataset,mask=False):
        """Writes measurements to data file with the following format:

        >header # starts the header selection, should be the first non-comment line in file
        >hostname # system where the file was generated
        >timestamp # time of file generation
        >toponame # OPTIONAL topology file used to generate measurements
        >trajname # OPTIONAL trajectory file(s) used to generate measurments
        >stepsize # OPTIONAL how long (in ps) is each timestep in trajectory
        >pdbname # OPTIONAL pdb file used to generate measurements
        >mask # OPTIONAL is data an occupancy mask, needed for properly those datasets
        >measure # name, type, and selection algebra defining a measurement
        >fields # 'time' and measure names/subnames as column headers
        >endheader # ends the header section
        >data # data values at each timestep for each field
        >end # marks end of file

        NOTE: There should be NO extraneous spaces in any names or entries,
        reading these files depends on splitting on whitespace."""

        with open(outfile, 'w') as datfile:
            # always write this metadata
            datfile.write('>header created by ExtendedTSC-%s\n' % self.version)
            datfile.write('>hostname %s\n' % self.data_hostname)
            datfile.write('>timestamp %s\n' % self.data_datetime)
            # optional fields
            if self.toponame:
                datfile.write('>toponame %s\n' % self.toponame)
            if self.trajname:
                datfile.write('>trajname %s\n' % self.trajname)
            if self.traj_stepsize:
                datfile.write('>stepsize(ps) %d\n' % self.traj_stepsize)
            if self.framerange:
                datfile.write('>framerange %d %d %d (start, stop, skip)\n' % (self.framerange[0],self.framerange[1],self.framerange[2]))
            # if self.pdbname:
            #     datfile.write('>pdbname %d\n' % self.pdbname)
            # obligatory fields
            if mask:
                datfile.write('>mask True\n')
            else:
                datfile.write('>mask False\n')
            for meas in dataset.measurements:
                datfile.write('>measure %s %s \"%s\"\n' % (meas.name, meas.type, meas.selecttext))
            # write out column headers
            datfile.write('>fields time ')
            for meas in dataset.measurements:
                if meas.width == 1:
                    datfile.write('%s ' % meas.name)
                # coordinate data measure types
                elif (meas.width%3 == 0) and (meas.type in ['atom', 'COG', 'COM']):
                    track = 1
                    coords = ['x','y','z']
                    for i in range(meas.width):
                        coord = coords[i%3]
                        datfile.write('%s:%s%d ' % (meas.name,coord,track))
                        if coord == 'z':
                            track += 1
                # for other widths, just use numbers to mark sub-measurements
                else:
                    for i in range(1,meas.width+1):
                        datfile.write('%s:%d ' % (meas.name,i))
            datfile.write('\n')
            datfile.write('>endheader\n')
            # eunmerate over time array and write data
            for idx,elem in np.ndenumerate(dataset.time):
                datfile.write('>data ')
                datfile.write('%s ' % str(elem))
                datfile.write('%s ' % ' '.join([str(i) for i in dataset.data[:,idx[0]]]))
                datfile.write('\n')
            # indcate that file was written completely
            datfile.write('>end')

    def _data_reader(self,infile):
        """Rebuilds array from previously generated measurements in a .dat file.
        Should only be called through __init__ and does not use TimeseriesCollection or MDAnalysis at all."""

        # open and read input file
        with open(infile, 'r') as datfile:
            lines = datfile.readlines()
        lines = [i.strip() for i in lines]
        # superficial check to see if input file was written to completion
        # will definitely not catch all types of problems!
        if lines[-1] != '>end':
            print 'Likely broken input file! Last line is not \">end\"! Exiting...'
            sys.exit(1)
        # read file header
        leader = lines[0].split()[0]
        while leader != '>endheader':
            p = lines.pop(0)
            leader = p.split()[0]
            # skip extraneous or marking lines
            if p.startswith('#') or p == '' or leader == '>header' or leader == '>endheader':
                pass
            # bail if data line is found while reading header
            elif leader == '>data':
                print 'data line found in header, check input file and try again!'
                sys.exit(1)
            # fill basic vars from header lines, not all are required
            elif leader == '>hostname':
                self.data_hostname = ' '.join(p.split()[1:])
            elif leader == '>timestamp':
                self.data_datetime = ' '.join(p.split()[1:])
            elif leader == '>toponame':
                self.toponame = ' '.join(p.split()[1:])
            elif leader == '>trajname':
                self.trajname = ' '.join(p.split()[1:])
            elif leader == '>pdbname':
                self.pdbname = ' '.join(p.split()[1:])
            elif leader == '>mask':
                if p.split()[1] == 'True':
                    self.read_mask = True
                    targetDS = self.maskDS
                elif p.split()[1] == 'False':
                    self.read_mask = False
                    targetDS = self.primaryDS
                else:
                    print 'Cannot determine if file %s is data or mask. Exiting...' % infile
                    sys.exit(1)
                if targetDS.populated:
                    print 'Data set init error with file %s, check your inputs. Exiting...' % infile
                    sys.exit(1)
            elif leader == '>stepsize(ps)':
                self.traj_stepsize = int(p.split()[1])
            elif leader == '>framerange':
                self.framerange = (int(p.split()[1]),int(p.split()[2]),int(p.split()[3]))
            # parse the more complicated lines
            # use '>measure' lines to compose descriptions
            elif leader == '>measure':
                # compose descriptor
                name = p.split()[1]
                meas_type = p.split()[2]
                selection = p.split('\"')[1]
                descriptor = (name, meas_type, selection)
                targetDS.add_measurement(descriptor)
            # use '>fields' line to figure out the width of each measurement
            elif leader == '>fields':
                fields = p.split()[1:]
                for elem in fields:
                    name = elem.split(':')[0]
                    for meas in targetDS.measurements:
                        if name == meas.name:
                            meas.incr_width()
        # process data lines
        tmplist = []
        tmptime = []
        leader = lines[0].split()[0]
        while leader != '>end':
            p = lines.pop(0)
            tmptime.append(float(p.split()[1]))
            if self.read_mask:
                tmplist.append([int(i) for i in p.split()[2:]])
            else:
                tmplist.append([float(i) for i in p.split()[2:]])
            # new leader
            leader = lines[0].split()[0]
        # convert lists to DataSet arrays
        targetDS.add_timesteps(np.array(tmptime))
        targetDS.add_collection(np.array(tmplist).T)

# Auxilary classes called by ExtendedTSC
# These hould not be instantiated directly in scripts!

class _Measurement(object):
    """A class to hold a single TimeseriesCollection measurement set and metadata."""

    def __init__(self, descriptor, width=0):
        self.name, self.type, self.selecttext = descriptor
        self.set_width(width)
        self.series = None

    def set_width(self, width):
        self.width = int(width)

    def incr_width(self):
        self.width += 1

    def add_data(self, array):
        self.series = array

class _DataSet(object):
    """A class to organize a set of measurements generated together."""

    def __init__(self):
        self.measurements = []
        self.populated = False

    def add_measurement(self, descriptor, width=0):
        self.measurements.append(_Measurement(descriptor, width=width))

    def count_measurements(self):
        return len(self.measurements)

    def add_collection(self, array):
        if self.populated:
            print 'Cannot add data to a _DataSet object twice! Exiting...'
            sys.exit(1)
        self.data = array
        # can only add data once
        self.populated = True

    def add_timesteps(self, array):
        self.time = array

    def simplify_indexing(self,return_dict=True):
        """Links _Measurement objects to the data array and (optinally) returns a dict for lookup by measurement name.
        Doesn't add or change any data, just makes access a bit simpler."""

        idx = 0
        self.access = {}
        for meas in self.measurements:
            meas.add_data(self.data[idx:idx+meas.width,:])
            self.access[meas.name] = meas.series
            idx += meas.width
        # plotting doesn't work well without doing this
        self.access['time'] = np.reshape(self.time, (1,-1))
        if return_dict:
            return self.access

class _VolumeSearch(AnalysisBase):
    """Class for steping through a trajectory frame by frame and tracking individual atoms/residues."""

    def __init__(self,vol_selecttext,search_selecttext,mode,universe,**kwargs):
        super(_VolumeSearch,self).__init__(universe.trajectory,**kwargs)
        self.vol_selecttext = vol_selecttext
        self.search_selecttext = search_selecttext
        self.u = universe
        # check mode
        valid_modes = ['res', 'atom']
        if mode in valid_modes:
            self.mode = mode
        else:
            print 'Invalid mode selected for volumetric search! Possible modes: %s.\nExiting...' % ' '.join(valid_modes)
            sys.exit(1)

    def _prepare(self):
        # setup vars and data structures
        self.selection_set = set()
        self.masking_list = []

    def _single_frame(self):
        # what to do at each frame
        print self._ts
        tmpoccupancy = []
        tmpatoms = self.u.select_atoms('(%s) and (%s)' % (self.vol_selecttext, self.search_selecttext))
        for atom in tmpatoms:
            # store unique identifiers for found res/atoms
            if self.mode == 'res':
                self.selection_set.add((atom.segid,atom.resid))
                tmpoccupancy.append((atom.segid,atom.resid))
            elif self.mode == 'atom':
                self.selection_set.add((atom.index,atom.segid,atom.resid,atom.name))
                tmpoccupancy.append(atom.index)
        self.masking_list.append(tmpoccupancy)

    def _conclude(self):
        # setup lists to hold descriptors
        self.selections = []
        self.selections_mask = []
        count = 0
        # compose selection descriptors
        for elem in sorted(self.selection_set):
            count += 1
            if self.mode == 'res':
                tmp_name = '%s_%s' % (str(elem[0]),str(elem[1]))
                tmp_type = 'COM'
                tmp_selecttext = 'segid %s and resid %s' % (str(elem[0]),str(elem[1]))
            elif self.mode == 'atom':
                tmp_name = '%s_%s_%s' % (str(elem[1]),str(elem[2]),str(elem[3]))
                tmp_type = 'atom'
                tmp_selecttext = 'segid %s and resid %s and name %s' % (str(elem[1]),str(elem[2]),str(elem[3]))
            self.selections.append((tmp_name,tmp_type,tmp_selecttext))
            self.selections_mask.append((tmp_name,'occupancy',tmp_selecttext))
        # build occupancy mask array
        self.mask = np.zeros((count,len(self.masking_list)), dtype=int)
        # re-iterate over selection set and check occupancy at each timestep
        for idxA, elem in enumerate(sorted(self.selection_set)):
            if self.mode == 'res':
                identity = elem
            elif self.mode == 'atom':
                identity = elem[0]
            for idxB, occupancy in enumerate(self.masking_list):
                if identity in occupancy:
                    try:
                        self.mask[idxA,idxB] = 1
                    except IndexError:
                        pass

class _WaterSearchZ(AnalysisBase):
    """Similar to _VolumeSearch, but specific for water molecules and does not track
    molecules through whole simulation. Instead, just saves the z-coordinate of each
    water found in the search volume at each timestep."""
    # TODO: specify coordinates (xyz) to save so more than z can be used
    # NOTE: selection is hardcoded to work only for the TIP3 model and CHARMM atom names

    def __init__(self,vol_selecttext,universe,**kwargs):
        super(_WaterSearchZ,self).__init__(universe.trajectory,**kwargs)
        self.vol_selecttext = vol_selecttext
        self.search_selecttext = 'name OH2 and resname TIP3'
        self.u = universe

    def _prepare(self):
        # setup vars and data structures
        self.zcoords = []

    def _single_frame(self):
        # what to do at each frame
        print self._ts
        tmpcoordlist = []
        tmpatoms = self.u.select_atoms('(%s) and (%s)' % (self.vol_selecttext, self.search_selecttext))
        for atom in tmpatoms:
            tmpcoordlist.append(atom.position[2])
        self.zcoords.append(tmpcoordlist)

    def _conclude(self):
        # find the timestep with the most water
        counts = np.zeros((len(self.zcoords)))
        for idx,coordlist in enumerate(self.zcoords):
            count = len(coordlist)
            counts[idx] = count
        # use that number to build array
        dimensions = (int(np.amax(counts)),len(self.zcoords))
        self.coordarray = np.empty(dimensions, dtype=float)
        # fill it with nans
        self.coordarray.fill(np.nan)
        # now fill it with coordinates, leaving nan in the gaps
        for idx1,coordlist in enumerate(self.zcoords):
            for idx0,coord in enumerate(coordlist):
                self.coordarray[idx0,idx1] = coord
        # for consistency make this a list, though it will only ever have one item
        self.selections = []
        self.selections.append(('water_volume', 'z-coordinates', '(%s) and (%s)' % (self.vol_selecttext, self.search_selecttext)))

class _GenericTSC(AnalysisBase):

    def __init__(self,measurements,universe,**kwargs):
        super(_GenericTSC,self).__init__(universe.trajectory,**kwargs)
        self.target_measurements = measurements
        self.u = universe

    def _prepare(self):
        total_width = 0
        # get array dimensions
        for meas in self.target_measurements:
            total_width += meas.width
        self.data = np.empty((total_width,self.n_frames), dtype=float)

    def _single_frame(self):
        print self._ts
        start_idx = 0
        for meas in self.target_measurements:
            tmpgroup = self.u.select_atoms(meas.selecttext)
            tmpwidth = meas.width
            # figure out what property is needed based on type
            if meas.type == 'atom':
                tmpselect= self.u.select_atoms(meas.selecttext)
            elif meas.type == 'bond':
                # not implemented
                pass
            elif meas.type == 'angle':
                # not implemented
                pass
            elif meas.type == 'dihedral':
                collection.addTimeseries(tm.Dihedral(self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'distance':
                collection.addTimeseries(tm.Distance('r', self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'COG':
                collection.addTimeseries(tm.CenterOfGeometry(self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'COM':
                collection.addTimeseries(tm.CenterOfMass(self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'water_dipole':
                # not implemented
                pass

    def _conclude(self):
        pass

class _ProcPDB(object):

    def __init__(self):
        pass
