import numpy as np

class _DataSet(object):
    """A class to organize a set of measurements generated together on one pass through a trajectory."""

    def __init__(self):
        self.measurements = []
        self.populated = False
        self.is_mask = False
        # metadata
        self.data_datetime = None # data generation timestamp
        self.data_hostname = None # data generation computer hostname
        self.toponame = None # path to topology file
        self.trajname = None # path to trajectory file
        self.traj_stepsize = None # in picoseconds
        self.pdbname = None # path to pdb file
        self.framerange = None # all frames, no skipping

    def copy_metadata(self, target):
        """Copy metadata from another DataSet instance (i.e. from primaryDS to maskDS)."""

        self.data_datetime = target.data_datetime
        self.data_hostname = target.data_hostname
        self.toponame = target.toponame
        self.trajname = target.trajname
        self.traj_stepsize = target.traj_stepsize
        self.pdbname = target.pdbname
        self.framerange = target.framerange

    def add_measurement(self, descriptor, width=0):
        self.measurements.append(_Measurement(descriptor, width=width))

    def count_measurements(self):
        return len(self.measurements)

    def get_width(self):
        return np.shape(self.data)[0]

    def get_length(self):
        return np.shape(self.data)[1]

    def get_starttime(self):
        return self.time[0]

    def get_endtime(self):
        return self.time[-1]

    def add_collection(self, array):
        self.data = array
        # can only add data once
        self.populated = True

    def add_timesteps(self, array):
        self.time = array

    def dataset_to_dict(self):
        """Links _Measurement objects to the data array and returns a dict for lookup
    by measurement name. Doesn't add or change any data, just makes access a bit simpler."""

        idx = 0
        access = {}
        for meas in self.measurements:
            meas.add_data(self.data[idx:idx+meas.width,:])
            access[meas.name] = meas.series
            idx += meas.width
        # plotting doesn't work well without doing this
        access['time'] = np.reshape(self.time, (1,-1))
        return access

    def setup_timesteps(self):
        """Generate time values (in picoseconds) after being populated."""

        if not self.populated:
            self.logger.err('Cannot setup DataSet time values! No data or wrong input type! Exiting...')
        # create an array of time values (in ps) for plotting and load it into datasets
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
            self.time = np.linspace(float(starttime), float(endtime), num=n_steps)

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
