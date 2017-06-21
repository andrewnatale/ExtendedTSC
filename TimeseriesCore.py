import sys, os, datetime
import numpy as np

class TimeseriesCore(object):
    """ExtendedTSC base class. On its own, this class can be used to read (and write) ExtendedTSC
    format data files. It is subclassed in other modules to provide trajectory processing
    functions with output to .dat files."""

    # version info
    file_format_version = '1.0'

    def __init__(self,datfile=None,maskfile=None):
        """Initialize TimeseriesCore and optionally load dat files.

    Keyword arguments:
    datfile - string; path to ExtendedTSC format dat file (non-mask)
    maskfile - string; path to ExtendedTSC format dat file (mask)
    """

        # collect some system data for keeping track of dat files
        self.data_datetime = None
        self.data_hostname = None
        # setup empty vars
        self.input_type = None # 'dat' 'traj' 'pdb'
        self.toponame = None
        self.trajname = None
        self.traj_stepsize = None # in picoseconds
        self.pdbname = None
        self.framerange = None # all frames, no skipping
        # init empty DataSet objects
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
            self._datasets_to_dicts()
        else:
            self._collect_sys_info()

    def write_data(self,fileprefix):
        """Writes data for all populated DataSets to .dat files in the cwd.

    Arguments:
    fileprefix - string; file names without extension (.dat or .mask.dat
        will be appended automatically)
    """

        if self.primaryDS.populated:
            self._data_writer(fileprefix+'.dat',self.primaryDS)
        if self.maskDS.populated:
            self._data_writer(fileprefix+'.mask.dat',self.maskDS,mask=True)

    def _collect_sys_info(self):
        """Collect some system data for tracking dat files."""

        self.data_datetime = datetime.datetime.now()
        self.data_hostname = os.uname()[1]

    def _setup_time(self):
        """Format time values (in picoseconds) after loading DataSet objects."""

        if (not self.primaryDS.populated) or (self.input_type == 'dat'):
            print 'Cannot setup DataSet time values! No data or wrong input type! Exiting...'
            sys.exit(1)
        # create an array of time values (in ps) for plotting and load it into datasets
        if self.input_type == 'pdb':
            self.primaryDS.add_timesteps(np.array([0,]))
        else:
            if self.framerange is None:
                firstframe = 0
                framestep = 1
            else:
                firstframe = self.framerange[0]
                framestep = self.framerange[2]
            starttime = firstframe * self.traj_stepsize
            n_steps = np.shape(self.primaryDS.data)[1]
            endtime = firstframe * self.traj_stepsize + self.traj_stepsize * framestep * (n_steps - 1)
            self.primaryDS.add_timesteps(np.linspace(float(starttime), float(endtime), num=n_steps))
            if self.maskDS.populated:
                self.maskDS.add_timesteps(np.linspace(float(starttime), float(endtime), num=n_steps))

    def _datasets_to_dicts(self):
        """Call the _simplify_indexing method of each DataSet object to build
    dictionaries based on measurement names"""

        if self.primaryDS.populated:
            self.primary = self.primaryDS._simplify_indexing()
        if self.maskDS.populated:
            self.mask = self.primaryDS._simplify_indexing()

    def _data_writer(self,dataset,outfile=None,mask=False):
        """Writes contents of a DataSet to a file with the following format:

    >header # starts the header selection, should be the first non-comment line in file
    >version # file format version for checking compatibility
    >hostname # system hostname of machine where the file was generated
    >timestamp # time of file generation
    >toponame # OPTIONAL topology file used to generate measurements
    >trajname # OPTIONAL trajectory file(s) used to generate measurments
    >stepsize # OPTIONAL how long (in ps) is each timestep in trajectory
    >pdbname # OPTIONAL pdb file used to generate measurements
    >mask # OPTIONAL is the data an occupancy mask, needed to properly load those datasets
    >measure # name, type, and selection algebra defining a measurement
    >fields # 'time' and measure names/subnames as column headers
    >endheader # ends the header section
    >data # data values at each timestep for each field
    >end # marks end of file

    NOTE: There should be NO extraneous spaces in any names or entries,
    reading these files depends on splitting on whitespace.

    Arguments:
    outfile - string; full output file name
    dataset - DataSet object; target data to process

    Keyword arguments:
    mask - boolean; indicates whether the target data is a boolean mask
    """

        # build a list of lines to write
        output_lines = []
        # always write this metadata
        output_lines.append('>header\n')
        output_lines.append('>version %s\n' % str(self.file_format_version))
        output_lines.append('>hostname %s\n' % self.data_hostname)
        output_lines.append('>timestamp %s\n' % self.data_datetime)
        # optional fields
        if self.toponame:
            output_lines.append('>toponame %s\n' % self.toponame)
        if self.trajname:
            output_lines.append('>trajname %s\n' % self.trajname)
        if self.traj_stepsize:
            output_lines.append('>stepsize(ps) %d\n' % self.traj_stepsize)
        if self.framerange:
            output_lines.append('>framerange %d %d %d (start, stop, skip)\n' % (self.framerange[0],self.framerange[1],self.framerange[2]))
        if self.pdbname:
            output_lines.append('>pdbname %s\n' % self.pdbname)
        # obligatory fields
        if mask:
            output_lines.append('>mask True\n')
        else:
            output_lines.append('>mask False\n')
        for meas in dataset.measurements:
            output_lines.append('>measure %s %s \"%s\"\n' % (meas.name, meas.type, meas.selecttext))
        # write out column headers
        output_lines.append('>fields time ')
        for meas in dataset.measurements:
            if meas.width == 1:
                output_lines.append('%s ' % meas.name)
            # coordinate data measure types
            elif (meas.width%3 == 0) and (meas.type in ['atom', 'COG', 'COM']):
                track = 1
                coords = ['x','y','z']
                for i in range(meas.width):
                    coord = coords[i%3]
                    output_lines.append('%s:%s%d ' % (meas.name,coord,track))
                    if coord == 'z':
                        track += 1
            # for other widths, just use numbers to mark sub-measurements
            else:
                for i in range(1,meas.width+1):
                    output_lines.append('%s:%d ' % (meas.name,i))
        output_lines.append('\n')
        output_lines.append('>endheader\n')
        # eunmerate over time array and write data
        for idx,elem in np.ndenumerate(dataset.time):
            output_lines.append('>data ')
            output_lines.append('%s ' % str(elem))
            output_lines.append('%s ' % ' '.join([str(i) for i in dataset.data[:,idx[0]]]))
            output_lines.append('\n')
        # indcate that file was written completely
        output_lines.append('>end')

        if outfile is None:
            print ''.join(output_lines)
        else:
            with open(outfile, 'w') as datfile:
                datfile.write(''.join(output_lines))

    def _data_reader(self,infile):
        """Rebuild DataSet from saved measurements in a .dat file.

    Arguments:
    infile - string; name of file to read
    """

        tmpversion = None
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
        # version check
        self._check_version(tmpversion)
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

    def _check_version(self,version_number):
        if version_number is None:
            print 'Warning! Cannot detect filetype version.'
        if version_number != self.file_format_version:
            # in the future, some versions may become obsolete and processing should stop here
            pass

# Auxilary core classes used by TimeseriesCore

class _Measurement(object):
    """A class to hold a single measurement and keep track of its description and data."""

    def __init__(self, descriptor, width):
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
    """A class to organize a set of measurements generated together on one pass through a trajectory."""

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
        """Links _Measurement objects to the data array and (optionally) returns a dict for lookup
    by measurement name. Doesn't add or change any data, just makes access a bit simpler."""

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
