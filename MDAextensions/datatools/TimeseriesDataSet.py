#!/usr/bin/env python2
from __future__ import print_function,division
from collections import Mapping,Sequence
import numpy as np
import re
from copy import copy,deepcopy

class TimeseriesDataSet(Mapping):
    """
    A class to organize a set of features extracted from the same trajectory.
    Also provides methods to write and read features to/from ExtendedTSC format .dat files.
    """

    current_file_verson = '1.2'
    compatible_file_versions = ['1.1','1.2']

    def __init__(self, infilename=None):
        """
        Keyword arguments:
        infilename - string; path to .dat file
        """
        # the list keeps everything in order until we're set up, then once the array is loaded
        # the features dictionary is how the data should be accessed in normal use
        self.feature_list = []
        self.feature_dict = {}
        # flags
        self.populated = False
        self.is_mask = False
        # metadata
        self.data_datetime = None # data generation timestamp
        self.data_hostname = None # data generation computer hostname
        self.toponame = None # path to topology file
        self.trajname = None # path to trajectory file
        self.traj_stepsize = None # in picoseconds
        self.pdbname = None # path to pdb file
        self.framerange = None # None means all frames, no skipping
        self.feature_list_type = None # 'static' or 'dynamic' - helps with post-processing logic
        self.rmsd_reference = None # only used for RMSD measurement type, either a filename or a frame number
        # load from file
        if infilename:
            self._read(infilename)

    # define abstract methods for Mapping:

    # indexing is by feature name, returns the cooresponding _Feature object
    # if no array has been loaded, this won't give you anything
    def __getitem__(self, i):
        if self.populated:
            return self.feature_dict[i]
        else:
            return None

    # iterate over features, using the list, not the dict (to exclude special elements like 'time')
    def __iter__(self):
        # if self.populated:
        #     return iter(self.feature_dict)
        # else:
        #     return iter(self.feature_list)
        return iter([f.name for f in self.feature_list])

    # __len__() gives length of feature list (so it won't include 'time')
    # other methods can give the shape of the array
    def __len__(self):
        return len(self.feature_list)

    def copy_metadata(self, target):
        """
        Copy metadata from another TimeseriesDataSet instance (i.e. from primaryDS to maskDS).

        Arguments:
        target - initialized TimeseriesDataSet instance to copy from
        """

        self.data_datetime = copy(target.data_datetime)
        self.data_hostname = copy(target.data_hostname)
        self.toponame = copy(target.toponame)
        self.trajname = copy(target.trajname)
        self.traj_stepsize = copy(target.traj_stepsize)
        self.pdbname = copy(target.pdbname)
        self.framerange = copy(target.framerange)
        self.feature_list_type = copy(target.feature_list_type)
        self.rmsd_reference = copy(target.rmsd_reference)

    def set_static(self):
        self.feature_list_type = 'static'

    def set_dynamic(self):
        self.feature_list_type = 'dynamic'

    def add_feature(self, descriptor, width=0):
        if not self.populated:
            self.feature_list.append(_Feature(descriptor, width=width))
        else:
            print('WARNING: Features cannot be added once the Data Set has been formatted!')
            print('WARNING: Ignoring \'add_feature\' operation!')

    def get_width(self):
        if self.populated:
            return np.shape(self.data)[0]
        else:
            return None

    def get_length(self):
        if self.populated:
            return np.shape(self.data)[1]
        else:
            return None

    def get_starttime(self):
        if self.populated:
            return self.time[0]
        else:
            return None

    def get_endtime(self):
        if self.populated:
            return self.time[-1]
        else:
            return None

    def format_data(self, array):
        """Format data using the feature list to allow for dict-like indexing by name."""

        self.data = array
        idx = 0
        for feature in self.feature_list:
            feature.add_data(self.data[idx:idx+feature.width,:])
            self.feature_dict[feature.name] = feature
            idx += feature.width
        # can only add data once
        self.populated = True
        # create an array of time values (in ps) for plotting and load it into datasets
        # you can override these automatic values with the set_time() method if needed
        if self.pdbname:
            self.time = np.array([0,])
        else:
            if self.framerange is None:
                firstframe = 0
                framestep = 1
            else:
                firstframe = self.framerange[0]
                framestep = self.framerange[2]
            starttime = firstframe * self.traj_stepsize
            n_steps = np.shape(self.data)[1]
            endtime = firstframe * self.traj_stepsize + self.traj_stepsize * framestep * (n_steps - 1)
            self.time = np.linspace(float(starttime), float(endtime), num=n_steps, endpoint=True)
        # allow time values to be accessed in the same way as data
        self.feature_dict['time'] = self.time

    def set_time(self, time_array):
        """Time values can usually be calculated from metadata and the number of frames, but
        if they need to be overridden, pass the new time values in a 1D array to this method."""

        self.time = time_array
        self.feature_dict['time'] = self.time

    def _write(self, outfilename=None):
        """
        Writes contents of the DataSet to stdout or a file with the following format:

        >header # starts the header selection, should be the first non-comment line in file
        >version # file format version for checking compatibility
        >hostname # system hostname of machine where the file was generated
        >timestamp # time of file generation
        >toponame # OPTIONAL topology file used to generate measurements
        >trajname # OPTIONAL trajectory file(s) used to generate measurements
        >stepsize(ps) # OPTIONAL how long (in ps) is each timestep in trajectory
        >pdbname # OPTIONAL pdb file used to generate measurements
        >rmsd_reference # OPTIONAL describes the reference structure used for rmsd calculations
        >framerange # specifies the trajectory slicing used to generate data
        >mask # True if the data is occupancy/counts (i.e. an interger array), needed to properly load those datasets
        >feature_list_type # static or dynamic - needed when merging dat files
        >feature # name, type, and selection algebra defining a feature (formerly '>measure')
        >fields # 'time' and measure names:subnames as column headers
        >endheader # ends the header section
        >data # data values at each timestep for each field
        >end # marks end of file

        NOTE: There should be NO extraneous spaces in any names or entries,
        reading these files depends on splitting on whitespace.

        Keyword arguments:
        outfilename - string; output file name; if not specified, lines are printed to stdout
        """

        if outfilename:
            print('Writing DataSet to file: %s' % outfilename)
        # build a list of lines to write
        output_lines = []
        # always write this metadata
        output_lines.append('>header\n')
        output_lines.append('>version %s\n' % str(self.current_file_verson))
        output_lines.append('>hostname %s\n' % self.data_hostname)
        output_lines.append('>timestamp %s\n' % self.data_datetime)
        # optional fields
        if self.toponame:
            output_lines.append('>toponame %s\n' % self.toponame)
        if self.trajname:
            output_lines.append('>trajname %s\n' % self.trajname)
        if self.traj_stepsize:
            output_lines.append('>stepsize(ps) %d\n' % self.traj_stepsize)
        if self.pdbname:
            output_lines.append('>pdbname %s\n' % self.pdbname)
        if self.rmsd_reference:
            output_lines.append('>rmsd_reference %s\n' % self.rmsd_reference)
        # obligatory fields
        # framerange requires special care so as to not break file loading
        if self.framerange is None:
            output_lines.append('>framerange 0 -1 1 (start, stop, step)\n')
        elif self.framerange[1] is None:
            output_lines.append('>framerange %d -1 %d (start, stop, step)\n' % (self.framerange[0],self.framerange[2]))
        else:
            output_lines.append('>framerange %d %d %d (start, stop, step)\n' % (self.framerange[0],self.framerange[1],self.framerange[2]))
        if self.is_mask:
            output_lines.append('>mask True\n')
        else:
            output_lines.append('>mask False\n')
        if self.feature_list_type:
            output_lines.append('>feature_list_type %s\n' % self.feature_list_type)
        for feature in self.feature_list:
            output_lines.append('>feature %s %s \"%s\"\n' % (feature.name, feature.type, feature.selecttext))
        # write out column headers
        output_lines.append('>fields time ')
        for feature in self.feature_list:
            if feature.width == 1:
                output_lines.append('%s ' % feature.name)
            # coordinate data measure types
            elif (feature.width%3 == 0) and (feature.type in ['atom', 'COG', 'COM']):
                track = 1
                coords = ['x','y','z']
                for i in range(feature.width):
                    coord = coords[i%3]
                    output_lines.append('%s:%s%d ' % (feature.name,coord,track))
                    if coord == 'z':
                        track += 1
            # for other widths, just use numbers to mark sub-measurements
            else:
                for i in range(1,feature.width+1):
                    output_lines.append('%s:%d ' % (feature.name,i))
        output_lines.append('\n')
        output_lines.append('>endheader\n')
        # eunmerate over time array and write data
        for idx,elem in np.ndenumerate(self.time):
            output_lines.append('>data ')
            output_lines.append('%s ' % str(elem))
            output_lines.append('%s ' % ' '.join([str(i) for i in self.data[:,idx[0]]]))
            output_lines.append('\n')
        # indcate that file was written completely
        output_lines.append('>end')
        # either write the file or just print what would be written
        if outfilename is None:
            print(''.join(output_lines))
        else:
            with open(outfilename, 'w') as datfile:
                datfile.write(''.join(output_lines))
            print('Finished writing output to: %s' % outfilename)

    def _read(self, infilename, enforce_version=False):
        """
        Rebuild DataSet from saved measurements in a .dat file.

        Arguments:
        infile - string; path to .dat file

        Keyword Arguments:
        enforce_version - boolean; if True turn on strict version checking - reading will fail if
            file version is not among those given in 'self.compatible_file_versions'
        """

        print('Reading data file: %s' % infilename)
        tmpversion = None
        # open and read input file
        with open(infilename, 'r') as datfile:
            lines = datfile.readlines()
        lines = [i.strip() for i in lines]
        # superficial check to see if input file was written to completion
        # will definitely not catch all types of problems!
        if lines[-1] != '>end':
            sys.exit('Likely broken input file! Last line is not \">end\"! Exiting...')
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
                sys.exit('Data line found in header, check input file and try again! Exiting...')
            # fill basic vars from header lines, not all are required
            elif leader == '>version':
                tmpversion = p.split()[1]
            elif leader == '>hostname':
                self.data_hostname = ' '.join(p.split()[1:])
            elif leader == '>timestamp':
                self.data_datetime = ' '.join(p.split()[1:])
            elif leader == '>toponame':
                self.toponame = p.split()[1]
            elif leader == '>trajname':
                self.trajname = p.split()[1]
            elif leader == '>pdbname':
                self.pdbname = p.split()[1]
            elif leader == '>rmsd_reference':
                self.rmsd_reference = ' '.join(p.split()[1:])
            elif leader == '>mask':
                if p.split()[1] == 'True':
                    read_mask = True
                elif p.split()[1] == 'False':
                    read_mask = False
                else:
                    sys.exit('Cannot determine if file %s is data or mask. Exiting...' % infilename)
            elif leader == '>stepsize(ps)':
                self.traj_stepsize = int(p.split()[1])
            # despite the mess in every other method that deals with framerange,
            # when loading from .dat files, it should always be 3 intergers
            elif leader == '>framerange':
                self.framerange = (int(p.split()[1]),int(p.split()[2]),int(p.split()[3]))
            elif leader == '>feature_list_type':
                self.feature_list_type = p.split()[1]
            # parse the more complicated lines
            # use '>feature' lines to compose descriptions
            # '>measure' is an old name for the same field
            elif leader == '>measure' or leader == '>feature':
                # compose descriptor
                name = p.split()[1]
                feature_type = p.split()[2]
                selection = p.split('\"')[1]
                descriptor = (name, feature_type, selection)
                self.add_feature(descriptor)
            # use '>fields' line to figure out the width of each feature
            elif leader == '>fields':
                fields = p.split()[1:]
                for elem in fields:
                    name = elem.split(':')[0]
                    for feature in self.feature_list:
                        if name == feature.name:
                            feature.incr_width()
        # version check
        self._check_file_version(tmpversion, enforce_version)
        # process data lines
        tmplist = []
        tmptime = []
        leader = lines[0].split()[0]
        while leader != '>end':
            p = lines.pop(0)
            tmptime.append(float(p.split()[1]))
            if read_mask:
                tmplist.append([int(i) for i in p.split()[2:]])
            else:
                tmplist.append([float(i) for i in p.split()[2:]])
            # new leader
            leader = lines[0].split()[0]
        # convert lists to arrays
        self.format_data(np.array(tmplist).T)
        # use time values from file
        self.feature_dict['time'] = np.array(tmptime)
        if read_mask:
            self.is_mask = True
        print('Finished reading file: %s' % infilename)

    def _check_file_version(self, version_number, enforce):
        # in some cases, we must fail on any mismatch
        if (enforce == True) and not (version_number in self.compatible_file_versions):
            sys.exit('File version %s cannot be processed! Exiting...' % version_number)
        elif version_number is None:
            print('Warning! Cannot detect file version!')
        elif not (version_number in self.compatible_file_versions):
            # in the future, some versions may become obsolete and processing should stop here if detected
            print('Warning! File version is not current! Some operations may not be supported.')

    def trim(self, regex_term, exclude=False):
        """
        Regenerate internal data structures from a subset of _Features, either keeping or rejecting
        _Features based on whether the name variable matches the regex_term. if exclude is False,
        keep only matching feature, if it is True, keep only non matching features.

        This operation is destuctive by design, use with care.
        """

        # match feature names by regex
        searcher = re.compile(regex_term)
        saved = []
        for feature in self.feature_list:
            if exclude is False:
                if searcher.match(feature.name):
                    saved.append(feature)
                else:
                    pass
            elif exclude is True:
                if searcher.match(feature.name):
                    pass
                else:
                    saved.append(feature)
        # merge individual feature series into a new primary array
        newdata = np.concatenate([feature.series for feature in saved], axis=0)
        # overwrite old data structures with the trimmed new ones
        self.feature_list = saved
        self.format_data(newdata)

    def match(self, regex_term, exclude=False):
        """
        Dry run for trim(), just print matching feature names, don't remove anything.
        """

        # match feature names by regex
        searcher = re.compile(regex_term)
        for feature in self.feature_list:
            if exclude is False:
                if searcher.match(feature.name):
                    print(feature.name)
                else:
                    pass
            elif exclude is True:
                if searcher.match(feature.name):
                    pass
                else:
                    print(feature.name)

    def return_trimmed(self, regex_term, exclude=False):
        """
        Return a new trimmed TimeseriesDataSet object without modifying the original. Deep copy
        _Features so that new doesn't reference the original.
        """

        new = TimeseriesDataSet()
        new.copy_metadata(self)
        # match feature names by regex
        searcher = re.compile(regex_term)
        saved = []
        for feature in self.feature_list:
            if exclude is False:
                if searcher.match(feature.name):
                    saved.append(deepcopy(feature))
                else:
                    pass
            elif exclude is True:
                if searcher.match(feature.name):
                    pass
                else:
                    saved.append(deepcopy(feature))
        new.feature_list = saved
        # merge individual feature series into a new primary array
        newdata = np.concatenate([feature.series for feature in saved], axis=0)
        new.format_data(newdata)
        return new

class _Feature(Sequence):
    """A class to hold a single measurement and keep track of its description and data."""

    # there are a lot of these, and the list changes a lot, so no checks;
    # this is for reference only
    valid_feature_types = ['atom', 'dihedral', 'COM', 'COG', 'distance', 'zsearchlist',
      'rmsdseries', 'occupancy', 'boolean', 'count']

    def __init__(self, descriptor, width):
        self.name, self.type, self.selecttext = descriptor
        self.set_width(width)
        self.series = None

    def __getitem__(self, i):
        if self.series is not None:
            return self.series[i,:]
        else:
            return None

    def __len__(self):
        return self.width

    def set_width(self, width):
        self.width = int(width)

    def incr_width(self):
        self.width += 1

    def add_data(self, array_ref):
        self.series = array_ref

    def get_descriptor(self):
        return (self.name, self.type, self.selecttext)
