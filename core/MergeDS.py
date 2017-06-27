import sys, os, datetime
import numpy as np
from TimeseriesCore import TimeseriesCore
from TimeseriesCore import _TSLogger

class MergeDS(TimeseriesCore):
    """Adds methods to merge related DataSet objects."""

    def __init__(self, verbose=True, log=None):
        # init logger
        self.logger = _TSLogger(verbose, log)
        # init datasets
        self.input_type = 'dat'
        self._init_datasets()

    def merge_along_time(self, datlist, strict_checking=False):
        """Load a list of datasets that contain the same features and merge them along the time
    axis to generate one longer timeseries.

    NOTE: Checking needs work, it is very basic!

    WARNING: This cannot presently handle mask DataSets (and is very likely to fail on their
    cooresponding primary DataSets). This is because when the classes that create these sets
    (so far only VolumeTracker) are run in parallel on different chunks of the same trajectory,
    they will find and track different objects, resulting in data sets with very different shapes
    and features.
    It wouldn't be too hard to make this work (all the information that you would get from running
    just one process on the full traj is still there) but I just don't need it yet.

    Arguments:
    datlist - list of strings; paths to dat files to merge

    Keyword arguments:
    strict_checking - boolean; if True do extra checks before merge - NOT IMPLEMENTED
    """

        pad_with_nan = False
        # load specified files into lists of DataSet objects
        dats = []
        #masks = []
        for filename in datlist:
            dats.append(self._data_reader(filename, False))
        # if masklist:
        #     for filename in masklist:
        #         masks.append(self._data_reader(filename, True))
        # extract parameters from the first DataSet in the list and compare to all the others
        test_stride = dats[0].framerange[2]
        test_width = dats[0].get_width()
        test_n_features = dats[0].count_measurements()
        first_frame = dats[0].framerange[0]
        last_frame = dats[0].framerange[1]
        # permissive checks
        for elem in dats[1:]:
            # check stride
            if (elem.framerange[2] != test_stride):
                self.logger.exit('Failed to merge DataSets due to stride mismatch! Exiting...')
            # check number of features
            if (elem.count_measurements() != test_n_features):
                self.logger.exit('Failed to merge DataSets due to different feature counts! Exiting...')
            # check width
            if (elem.get_width() != test_width) and (test_n_features == 1):
                self.logger.msg('Warning, data width mismatch! This is expected for some DataSets as long as the feature count = 1.\nArrays will be padded with \'nan\' as needed to match.')
                pad_with_nan = True
            elif (elem.get_width() != test_width):
                self.logger.exit('Failed to merge DataSets due to different data widths and feature count > 1! Exiting...')
            # check if last frame of preceding file matches first frame of this one
            if (last_frame != elem.framerange[0]):
                self.logger.msg('Warning! Apparent overlap or missing frames during merge!')
            last_frame = elem.framerange[1]
        # if strict checking is requested, do some addtional checks
        if strict_checking:
            # checks to implement:
            # topo/traj name matching
            # matching features
            # more stringent time continuity check
            pass
        # empty default DataSets should have been pre-initialized by __init__()
        # copy the metadata from the first DataSet in the list
        self.primaryDS.copy_metadata(dats[0])
        # get the full framerange
        self.primaryDS.framerange = (first_frame, last_frame, test_stride)
        # copy measurement objects
        for meas in dats[0].measurements:
            self.primaryDS.add_measurement((meas.name, meas.type, meas.selecttext), meas.width)
        # add padding (if needed), concatenate, and load data and time arrays
        if pad_with_nan is True:
            all_widths = [i.get_width() for i in dats]
            maxwidth = np.amax(all_widths)
            for elem in dats:
                width_diff = maxwidth - elem.get_width()
                if width_diff > 0:
                    padding = np.empty((width_diff,elem.get_length()), dtype=float)
                    padding.fill(np.nan)
                    elem.data = np.concatenate([elem.data, padding], axis=0)
                elif width_diff == 0:
                    self.primaryDS.measurements[0].set_width(elem.get_width())
        self.primaryDS.add_collection(np.concatenate([i.data for i in dats], axis=1))
        self.primaryDS.add_timesteps(np.concatenate([i.time for i in dats], axis=0))
        ## rinse and repeat for maskDS
        # if masklist:
        #     self.maskDS.copy_metadata(self.primaryDS)
        #     for meas in masks[0].measurements:
        #         self.maskDS.add_measurement((meas.name, meas.type, meas.selecttext), meas.width)
        #     self.maskDS.add_collection(np.concatenate([i.data for i in masks], axis=1))
        #     self.maskDS.add_timesteps(np.concatenate([i.time for i in masks], axis=0))

    def merge_along_features(self, datlist, strict_checking=False):
        """Load a list of datasets that cover the same time window of the same trajectory and merge
        them along the features axis."""
        pass # not yet implemented
