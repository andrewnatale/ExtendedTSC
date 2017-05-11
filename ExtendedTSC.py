# ExtendedTSC.py - a small module to wrap and extend MDAnalysis' TimeseriesCollection module.
# v0.1.1
# TODO:
#   - allow to define a range of frames to analyze
#   - allow loading of multiple dcd files into a single trajectory
#   - implement additional measurement types
#   - implement extension of a TimeseriesCollection with more complicated calculations

import sys, os, math, datetime
import numpy as np
# tested with MDAnalysis-0.15.0
# MDAnalysis-0.16.0 WILL NOT WORK WITH THIS SCRIPT
import MDAnalysis as md
import MDAnalysis.core.Timeseries as tm

# wrap and extend MDAnalysis' TimeseriesCollection module
class ExtendedTSC(object):

    def __init__(self,readfile=None):
        self.populated = False
        self.toponame = None
        self.trajname = None
        self.stepsize = None
        self.data_datetime = None
        self.data_hostname = None
        self.measurements = []
        if readfile:
            self.read_dat(readfile)


    # use MDAnalysis to generate measurement array from simulation data
    def measures_from_dcd(self,selections,psffile,dcdfile,stepsize):
        # check 'populated' flag
        if self.populated:
            print 'data arrays already populated, cannot do so again!'
            sys.exit(1)
        # tracking info
        self.toponame = os.path.abspath(psffile)
        self.trajname = os.path.abspath(dcdfile)
        self.stepsize = stepsize
        self.data_datetime = datetime.datetime.now()
        self.data_hostname = os.uname()[1]
        # init Measurement objects into a list
        for descriptor in selections:
            self.measurements.append(Measurement(descriptor))
        # init mda universe object
        self.u = md.Universe(psffile,dcdfile,topology_format=u'PSF')
        # load measurements
        self.collection = tm.TimeseriesCollection()
        for meas in self.measurements:
            if meas.type == 'atom':
                tmpselect = self.u.select_atoms(meas.selecttext)
                meas.set_width(tmpselect.n_atoms * 3)
                self.collection.addTimeseries(tm.Atom('v', self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'bond':
                # not implemented
                print 'measurement type %s not implemented, exiting' % meas.type
                sys.exit(1)
            elif meas.type == 'angle':
                # not implemented
                print 'measurement type %s not implemented, exiting' % meas.type
                sys.exit(1)
            elif meas.type == 'dihedral':
                meas.set_width(1)
                self.collection.addTimeseries(tm.Dihedral(self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'distance':
                meas.set_width(1)
                self.collection.addTimeseries(tm.Distance('r', self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'COG':
                meas.set_width(3)
                self.collection.addTimeseries(tm.CenterOfGeometry(self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'COM':
                # not implemented
                print 'measurement type %s not implemented, exiting' % meas.type
                sys.exit(1)
            elif meas.type == 'water_dipole':
                # not implemented
                print 'measurement type %s not implemented, exiting' % meas.type
                sys.exit(1)
            else:
                # warn, but continue
                print 'unrecognized measure type!'
        # compute from trajectory
        self.collection.compute(self.u.trajectory, start=0, stop=-1)
        # create an array of time values (in ns) for plotting
        self.timesteps = np.linspace(0,
                                     float(np.shape(self.collection.data)[1]-1)*float(stepsize),
                                     np.shape(self.collection.data)[1]
                                     )
        # toggle 'populated' flag
        self.populated = True

    def read_dat(self,infile):
        # rebuild array from previously generated measurements in a .dat file
        # check 'populated' flag
        if self.populated:
            print 'data arrays already populated, cannot do so again!'
            sys.exit(1)
        with open(infile, 'r') as datfile:
            lines = datfile.readlines()
        lines = [i.strip() for i in lines]
        # superficial check to see if input file was written to completion
        # will definitely not catch all types of errors!
        if lines[-1] != '>end':
            print 'possibly broken input file!'
            print 'last line is not \">end\"!'
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
            # fill basic vars from header lines
            elif leader == '>hostname':
                self.data_hostname = ' '.join(p.split()[1:])
            elif leader == '>timestamp':
                self.data_datetime = ' '.join(p.split()[1:])
            elif leader == '>toponame':
                self.toponame = ' '.join(p.split()[1:])
            elif leader == '>trajname':
                self.trajname = ' '.join(p.split()[1:])
            elif leader == '>stepsize(ps)':
                self.stepsize = int(p.split()[1])
            # parse the more complicated lines
            # use '>measure' lines to create Measurement objects
            elif leader == '>measure':
                # compose descriptor
                name = p.split()[1]
                meas_type = p.split()[2]
                selection = p.split('\"')[1]
                descriptor = (name, meas_type, selection)
                self.measurements.append(Measurement(descriptor))
            # use '>fields' line to figure out the width of each measurement
            elif leader == '>fields':
                fields = p.split()[1:]
                for elem in fields:
                    name = elem.split(':')[0]
                    for meas in self.measurements:
                        if name == meas.name:
                            meas.incr_width()
        # process data lines
        tmplist = []
        tmptime = []
        leader = lines[0].split()[0]
        while leader != '>end':
            p = lines.pop(0)
            tmptime.append(float(p.split()[1]))
            tmplist.append([float(i) for i in p.split()[2:]])
            # new leader
            leader = lines[0].split()[0]
        # convert to arrays
        self.timesteps = np.array(tmptime)
        self.collection = DummyCollection(np.array(tmplist).T)
        # toggle 'populated' flag
        self.populated = True

    def write_dat(self,outfile):
        # write measurements to data file with the following format:
        # # start a line containing a comment with '#'
        # >header # starts the header selection, should be the first non-comment line in file
        # # header sections may contain the following fields:
        # >hostname # computer where the file was generated
        # >timestamp # time of file generation
        # >toponame # psf file used to generate measurements
        # >trajname # dcd file(s) used to generate measurments
        # >stepsize # how long (in ps) is each timestep in trajectory
        # >measure # name, type, and selection algebra defining a measurement
        # >fields # 'time' and measure names as column headers
        # >endheader # ends the header section
        # >data # values at each timestep for each field
        # >end # marks end of file to help catch write errors
        if not self.populated:
            print 'no data to write!'
            sys.exit(1)
        # write header with tracking and selection information
        with open(outfile, 'w') as datfile:
            datfile.write('>header\n')
            if self.data_hostname:
                datfile.write('>hostname %s\n' % self.data_hostname)
            if self.data_datetime:
                datfile.write('>timestamp %s\n' % self.data_datetime)
            if self.toponame:
                datfile.write('>toponame %s\n' % self.toponame)
            if self.trajname:
                datfile.write('>trajname %s\n' % self.trajname)
            if self.stepsize:
                datfile.write('>stepsize(ps) %s\n' % self.stepsize)
            for meas in self.measurements:
                datfile.write('>measure %s %s \"%s\"\n' % (meas.name, meas.type, meas.selecttext))
            # write out column headers
            datfile.write('>fields time ')
            for meas in self.measurements:
                if meas.width == 1:
                    datfile.write('%s ' % meas.name)
                # if width is 3, the data are likely a single x,y,z coordinate
                elif meas.width == 3:
                    datfile.write('%s:x ' % meas.name)
                    datfile.write('%s:y ' % meas.name)
                    datfile.write('%s:z ' % meas.name)
                # if width is otherwise divisible by 3, data is a group of coordinates
                elif meas.width%3 == 0:
                    track = 1
                    coords = ['x','y','z']
                    for i in range(meas.width):
                        coord = coords[i%3]
                        datfile.write('%s:%s%d ' % (meas.name,coord,track))
                        if coord == 'z':
                            track += 1
                # for other widths (should these exist?), just use numbers to mark sub-measurements
                else:
                    for i in range(1,meas.width+1):
                        datfile.write('%s:%d ' % (meas.name,i))
            datfile.write('\n')
            datfile.write('>endheader\n')
            # eunmerate over time array and write data
            for idx,elem in np.ndenumerate(self.timesteps):
                datfile.write('>data ')
                datfile.write('%s ' % str(elem))
                datfile.write('%s ' % ' '.join([str(i) for i in self.collection.data[:,idx[0]]]))
                datfile.write('\n')
            # indcate that file was written completely
            datfile.write('>end')

    def simplify_indexing(self):
        # link Measurement objects to array and return a dict for lookup by name
        # doesn't add or change any data, just makes access a bit simpler
        # this function might be buggy... use with caution!
        # if I were sure the width determination always worked well, this would
        # be allowed to happen automatically when generating or loading data
        idx = 0
        self.access = {}
        for meas in self.measurements:
            try:
                meas.add_data(self.collection.data[idx:idx+meas.width,:])
            except IndexError:
                meas.add_data(self.collection.data[idx:,:])
            self.access[meas.name] = meas.ts
            idx += meas.width
        # plotting doesn't work well without doing this
        self.access['time'] = np.reshape(self.timesteps, (1,-1))
        #return access

# container object for a single Timeseries measurement
class Measurement(object):

    def __init__(self, descriptor):
        self.name, self.type, self.selecttext = descriptor
        self.width = 0
        self.ts = None

    def set_width(self, int_width):
        self.width = int_width

    def incr_width(self):
        self.width += 1

    def add_data(self, array):
        self.ts = array

# spoof an object from MDAnalysis when reading data from a text file
# keeps the code above simpler and makes data manipulations more consistent
class DummyCollection(object):

    def __init__(self, array):
        self.data = array
