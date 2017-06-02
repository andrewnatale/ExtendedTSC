# ExtendedTSC.py - a small module to wrap and extend MDAnalysis' TimeseriesCollection module.
# v0.3.0
# probably only works on MacOS/Linux/Unix
# TODO:
#   - allow to define a range of frames to analyze - DONE v0.2.0
#   - allow loading of multiple dcd files which form one continuous trajectory
#   - implement additional measurement types allowed by TimeseriesCollection
#   - implement extension of a TimeseriesCollection with more complicated calculations - partially DONE v0.3.0 - added interface to write arbitrary dat files
#   - pass a Universe object as arg instead of creating our own?
#   - getting messy, needs a refactor and incorporate volume tracking as a core function
#   - load pdb

import sys, os, math, datetime
import numpy as np
# tested with MDAnalysis-0.15.0
# MDAnalysis-0.16.0 WILL NOT WORK WITH THIS SCRIPT
import MDAnalysis as md
import MDAnalysis.core.Timeseries as tm

# wrap and extend MDAnalysis' TimeseriesCollection module
class ExtendedTSC(object):

    def __init__(self,readfile=None,mask=False):
        self.version = '0.3.0'
        self.populated = False
        self.deriv_populated = False
        self.toponame = None
        self.trajname = None
        self.traj_stepsize = None
        self.startframe = None
        self.stopframe = None
        self.frameskip = None
        self.data_datetime = None
        self.data_hostname = None
        self.measurements = []
        if readfile:
            self.read_dat(readfile,mask=mask)

    # use MDAnalysis to generate measurement array from simulation data
    def measures_from_dcd(self,selections,psffile,dcdfile,traj_stepsize=500,startframe=0,stopframe=-1,frameskip=1):
        # check 'populated' flag
        if self.populated:
            print 'data arrays already populated, cannot do so again!'
            sys.exit(1)
        # tracking info
        self.toponame = os.path.abspath(psffile)
        self.trajname = os.path.abspath(dcdfile)
        self.traj_stepsize = traj_stepsize
        self.startframe = startframe
        self.stopframe = stopframe
        self.frameskip = frameskip
        self.data_datetime = datetime.datetime.now()
        self.data_hostname = os.uname()[1]
        # init Measurement objects into a list
        for descriptor in selections:
            self.measurements.append(Measurement(descriptor))
        # init mda universe object
        self.u = md.Universe(psffile,dcdfile,topology_format=u'PSF')
        # load measurements
        self.collection = tm.TimeseriesCollection()
        # this loop identifies all the possible measurements that can be used with TimeseriesCollection
        # however not all are implemented (really easy to do, I just haven't used some at all)
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
                meas.set_width(3)
                self.collection.addTimeseries(tm.CenterOfMass(self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'water_dipole':
                # not implemented
                print 'measurement type %s not implemented, exiting' % meas.type
                sys.exit(1)
            else:
                print 'unrecognized measure type %s!' % meas.type
                sys.exit(1)
        # compute from trajectory
        # NOTE: I think there's a bug in MDAnalysis' DCD reader that causes it to drop the last
        # frame it should measure when slicing the trajectory in any way - need to confirm this
        # everything should still work mostly as advertised
        self.collection.compute(self.u.trajectory, start=self.startframe, stop=self.stopframe, skip=self.frameskip)
        # create an array of time values (in ps) for plotting
        starttime = self.startframe * self.traj_stepsize
        endtime = self.startframe * self.traj_stepsize + self.traj_stepsize * self.frameskip * (np.shape(self.collection.data)[1] - 1)
        self.timesteps = np.linspace(float(starttime), float(endtime), num=np.shape(self.collection.data)[1])
        # toggle 'populated' flag
        self.populated = True

    def derivative_collection(self,descriptors,widths,array,time=None):
        # basically just a way to interface with the file IO functions with arbitrary data
        # the main way to use it would be to generate an inital collection, from a trajectory or a dat file,
        # then do some calculations to derive new timeseries-like data structures
        # then feed that back into this function to be able to write those new calculations
        # to a dat file in the same way as the inital ones
        # descriptors is a list of 3 string tuples (won't be processed by MDA, so can simply be descriptive)
        # widths is a list of the same length as descriptors, containing the width of each measure
        # array contains all the data that will be linked with the measurements
        # NOTE: no error checking done here, this functions relies on the user to properly format the imput!
        if time is not None:
            self.timesteps = time
        self.deriv_measurements = []
        for desc in descriptors:
            self.deriv_measurements.append(Measurement(desc))
        for idx,elem in enumerate(widths):
            self.deriv_measurements[idx].set_width(elem)
        self.deriv_collection = DummyCollection(array)
        # toggle 'populated' flag
        self.deriv_populated = True

    def read_dat(self,infile,mask=False):
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
                self.traj_stepsize = int(p.split()[1])
            elif leader == '>framerange':
                self.startframe = int(p.split()[1])
                self.stopframe = int(p.split()[2])
                self.frameskip = int(p.split()[3])
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
            if mask is True:
                tmplist.append([int(i) for i in p.split()[2:]])
            else:
                tmplist.append([float(i) for i in p.split()[2:]])
            # new leader
            leader = lines[0].split()[0]
        # convert to arrays
        self.timesteps = np.array(tmptime)
        self.collection = DummyCollection(np.array(tmplist).T)
        # toggle 'populated' flag
        self.populated = True

    def simplify_indexing(self,derivative=False,return_dict=False):
        # link Measurement objects to array and return a dict for lookup by name
        # doesn't add or change any data, just makes access a bit simpler
        # this function might be buggy... use with caution!
        # if I were sure the width determination always worked well, this would
        # be allowed to happen automatically when generating or loading data
        if derivative is True:
            tgt_meas = self.deriv_measurements
            tgt_collection = self.deriv_collection
        else:
            tgt_meas = self.measurements
            tgt_collection = self.collection
        idx = 0
        self.access = {}
        for meas in tgt_meas:
            try:
                meas.add_data(tgt_collection.data[idx:idx+meas.width,:])
            except IndexError:
                meas.add_data(tgt_collection.data[idx:,:])
            self.access[meas.name] = meas.ts
            idx += meas.width
        # plotting doesn't work well without doing this
        self.access['time'] = np.reshape(self.timesteps, (1,-1))
        if return_dict:
            return self.access

    def write_dat(self,outfile,derivative=False):
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
        if derivative is True:
            if not self.deriv_populated:
                print 'no data to write!'
                sys.exit(1)
            tgt_meas = self.deriv_measurements
            tgt_collection = self.deriv_collection
        else:
            if not self.populated:
                print 'no data to write!'
                sys.exit(1)
            tgt_meas = self.measurements
            tgt_collection = self.collection
        # write header with tracking and selection information
        with open(outfile, 'w') as datfile:
            datfile.write('>header created by ExtendedTSC-%s\n' % self.version)
            datfile.write('>hostname %s\n' % self.data_hostname)
            datfile.write('>timestamp %s\n' % self.data_datetime)
            datfile.write('>toponame %s\n' % self.toponame)
            datfile.write('>trajname %s\n' % self.trajname)
            datfile.write('>stepsize(ps) %d\n' % self.traj_stepsize)
            datfile.write('>framerange %d %d %d (start, stop, skip)\n' % (self.startframe, self.stopframe, self.frameskip))
            for meas in tgt_meas:
                datfile.write('>measure %s %s \"%s\"\n' % (meas.name, meas.type, meas.selecttext))
            # write out column headers
            datfile.write('>fields time ')
            for meas in tgt_meas:
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
                datfile.write('%s ' % ' '.join([str(i) for i in tgt_collection.data[:,idx[0]]]))
                datfile.write('\n')
            # indcate that file was written completely
            datfile.write('>end')

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
        self.isdummy = True
