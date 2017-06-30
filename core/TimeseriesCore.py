import sys, os, datetime
import numpy as np
import DataSet

class TimeseriesCore(object):
    """ExtendedTSC base class. On its own, this class can be used to read (and write) ExtendedTSC
    format data files. It is subclassed in other modules to provide trajectory processing
    functions with output to .dat files."""

    # version info
    file_format_version = '1.0'

    def __init__(self, datfilename=None, maskfilename=None, verbose=True, log=None):
        """Initialize TimeseriesCore and optionally load dat files.

    Keyword arguments:
    datfile - string; path to ExtendedTSC format dat file (non-mask)
    maskfile - string; path to ExtendedTSC format dat file (mask)
    verbose - boolean; if True, print to stdout while running
    log - boolean; if True, write script output to a log file in the cwd
    """

        # init logger
        self.logger = _TSLogger(verbose, log)
        self.input_type = None # 'dat' 'traj' 'pdb'
        # load a single dat/mask pair through __init__(), otherwise use load()
        if datfilename and not maskfilename:
            self.primaryDS = self._data_reader(datfilename, False)
            self.maskDS = None
            self.input_type = 'dat'
        elif datfilename and maskfilename:
            self.primaryDS = self._data_reader(datfilename, False)
            self.maskDS = self._data_reader(maskfilename, True)
            self.input_type = 'dat'
        elif maskfilename and not datfilename:
            self.logger.err('Cannot load a mask file without a corresponding dat file! Exiting...')
        else:
            # not sure if this should be left to subclasses?
            self._init_datasets()

    def write_data(self, fileprefix):
        """Writes data for all populated DataSets to .dat files in the cwd.

    Arguments:
    fileprefix - string; file names without extension (.dat or .mask.dat
        will be appended automatically)
    """

        if self.primaryDS.populated:
            self._data_writer(self.primaryDS, outfilename=fileprefix+'.dat')
        if self.maskDS.populated:
            self._data_writer(self.maskDS, outfilename=fileprefix+'.mask.dat')

    def _init_datasets(self):
        """Setup empty default DataSet objects. Subclasses should use these unless there is a good
    reason to customize."""

        self.primaryDS = DataSet._DataSet()
        self.primaryDS.data_datetime = str(datetime.datetime.now())
        self.primaryDS.data_hostname = os.uname()[1]
        self.maskDS = DataSet._DataSet()
        self.maskDS.is_mask = True

    # def _datasets_to_dicts(self):
    #     """Call the _simplify_indexing method of each DataSet object to build
    # dictionaries based on measurement names"""
    #
    #     if self.primaryDS.populated:
    #         self.primary = self.primaryDS.simplify_indexing()
    #     if self.maskDS.populated:
    #         self.mask = self.primaryDS.simplify_indexing()

    def _custom_dataset(self):
        """Return an initialized DataSet object for customization."""

        custom_dataset = DataSet._DataSet()
        custom_dataset.data_datetime = str(datetime.datetime.now())
        custom_dataset.data_hostname = os.uname()[1]
        return custom_dataset

    def _data_writer(self, dataset, outfilename=None):
        """Writes contents of a DataSet to a file with the following format:

    >header # starts the header selection, should be the first non-comment line in file
    >version # file format version for checking compatibility
    >hostname # system hostname of machine where the file was generated
    >timestamp # time of file generation
    >toponame # OPTIONAL topology file used to generate measurements
    >trajname # OPTIONAL trajectory file(s) used to generate measurements
    >stepsize # OPTIONAL how long (in ps) is each timestep in trajectory
    >pdbname # OPTIONAL pdb file used to generate measurements
    >mask # OPTIONAL is the data an occupancy mask, needed to properly load those datasets
    >measure # name, type, and selection algebra defining a measurement
    >fields # 'time' and measure names:subnames as column headers
    >endheader # ends the header section
    >data # data values at each timestep for each field
    >end # marks end of file

    NOTE: There should be NO extraneous spaces in any names or entries,
    reading these files depends on splitting on whitespace.

    Arguments:
    dataset - DataSet object; target data to process

    Keyword arguments:
    outfilename - string; output file name; if not specified, file lines are printed to stdout
    """

        if outfilename:
            self.logger.msg('Writing DataSet to file: %s' % outfilename)
        # build a list of lines to write
        output_lines = []
        # always write this metadata
        output_lines.append('>header\n')
        output_lines.append('>version %s\n' % str(self.file_format_version))
        output_lines.append('>hostname %s\n' % dataset.data_hostname)
        output_lines.append('>timestamp %s\n' % dataset.data_datetime)
        # optional fields
        if dataset.toponame:
            output_lines.append('>toponame %s\n' % dataset.toponame)
        if dataset.trajname:
            output_lines.append('>trajname %s\n' % dataset.trajname)
        if dataset.traj_stepsize:
            output_lines.append('>stepsize(ps) %d\n' % dataset.traj_stepsize)
        if dataset.framerange:
            output_lines.append('>framerange %d %d %d (start, stop, skip)\n' % (dataset.framerange[0],dataset.framerange[1],dataset.framerange[2]))
        if dataset.pdbname:
            output_lines.append('>pdbname %s\n' % dataset.pdbname)
        # obligatory fields
        if dataset.is_mask:
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
        # either write the file or just print what would be written
        if outfilename is None:
            print ''.join(output_lines)
        else:
            with open(outfilename, 'w') as datfile:
                datfile.write(''.join(output_lines))
            self.logger.msg('Finished writing output to: %s' % outfilename)

    def _data_reader(self, infilename, mask):
        """Rebuild DataSet from saved measurements in a .dat file.

    Arguments:
    infile - string; name of file to read
    """

        self.logger.msg('Reading data file: %s' % infilename)
        if mask:
            self.logger.msg('Expecting file %s to be of special type \'mask\'' % infilename)
        tmpversion = None
        tmpDataSet = DataSet._DataSet()
        # open and read input file
        with open(infilename, 'r') as datfile:
            lines = datfile.readlines()
        lines = [i.strip() for i in lines]
        # superficial check to see if input file was written to completion
        # will definitely not catch all types of problems!
        if lines[-1] != '>end':
            self.logger.err('Likely broken input file! Last line is not \">end\"! Exiting...')
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
                self.logger.err('Data line found in header, check input file and try again!')
            # fill basic vars from header lines, not all are required
            elif leader == '>version':
                tmpversion = p.split()[1]
            elif leader == '>hostname':
                tmpDataSet.data_hostname = ' '.join(p.split()[1:])
            elif leader == '>timestamp':
                tmpDataSet.data_datetime = ' '.join(p.split()[1:])
            elif leader == '>toponame':
                tmpDataSet.toponame = p.split()[1]
            elif leader == '>trajname':
                tmpDataSet.trajname = p.split()[1]
            elif leader == '>pdbname':
                tmpDataSet.pdbname = p.split()[1]
            elif leader == '>mask':
                if p.split()[1] == 'True':
                    read_mask = True
                elif p.split()[1] == 'False':
                    read_mask = False
                else:
                    self.logger.err('Cannot determine if file %s is data or mask. Exiting...' % infilename)
                if read_mask != mask:
                    self.logger.err('Mismatch between expected filetype (mask or dat) and what was found in %s. Exiting...' % infilename)
            elif leader == '>stepsize(ps)':
                tmpDataSet.traj_stepsize = int(p.split()[1])
            elif leader == '>framerange':
                tmpDataSet.framerange = (int(p.split()[1]),int(p.split()[2]),int(p.split()[3]))
            # parse the more complicated lines
            # use '>measure' lines to compose descriptions
            elif leader == '>measure':
                # compose descriptor
                name = p.split()[1]
                meas_type = p.split()[2]
                selection = p.split('\"')[1]
                descriptor = (name, meas_type, selection)
                tmpDataSet.add_measurement(descriptor)
            # use '>fields' line to figure out the width of each measurement
            elif leader == '>fields':
                fields = p.split()[1:]
                for elem in fields:
                    name = elem.split(':')[0]
                    for meas in tmpDataSet.measurements:
                        if name == meas.name:
                            meas.incr_width()
        # version check
        self._check_file_version(tmpversion)
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
        # convert lists to DataSet arrays
        tmpDataSet.add_timesteps(np.array(tmptime))
        tmpDataSet.add_collection(np.array(tmplist).T)
        if read_mask:
            tmpDataSet.is_mask = True
        self.logger.msg('Finished reading file: %s' % infilename)
        return tmpDataSet

    def _check_file_version(self, version_number):
        if version_number is None:
            self.logger.msg('Warning! Cannot detect filetype version!')
        elif version_number != self.file_format_version:
            # in the future, some versions may become obsolete and processing should stop here
            pass

class _TSLogger(object):
    """A simple class to log and print actions taken by TimeseriesCore and child objects."""

    def __init__(self, verbose, log):
        self.verbose = verbose
        self.log_file = None
        if log:
            self.log_file = open(os.path.join(os.getcwd(), log), 'w')

    def msg(self, msgtext):
        if self.log_file:
            self.log_file.write(msgtext)
            self.log_file.write('\n')
        if self.verbose:
            print msgtext

    def err(self, msgtext):
        if self.log_file:
            self.log_file.write('Error: '+msgtext)
            self.log_file.write('\n')
            self.log_file.close()
        sys.exit('Error: '+msgtext)
