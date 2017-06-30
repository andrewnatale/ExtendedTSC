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

    def merge_along_time(self, datlist):
        """Load a list of datasets that contain the same features and merge them along the time
    axis to generate one longer timeseries.

    WARNING: Datfiles must be input in the correct time-order! This method has no way to sort
        them itself.

    Arguments:
    datlist - list of strings; paths to dat files to merge
    """

        # control switches
        pad_with_nan = False
        # load specified files into lists of DataSet objects
        dats = []
        for filename in datlist:
            dats.append(self._data_reader(filename, enforce_version=True))
        # extract parameters from the first DataSet in the list and compare to all the others
        test_width = dats[0].get_width()
        test_n_features = dats[0].count_measurements()
        first_frame = dats[0].framerange[0]
        last_frame = dats[0].framerange[1]
        test_stride = dats[0].framerange[2]
        test_is_mask = dats[0].is_mask
        test_feature_list_type = dats[0].feature_list_type
        # checking is fairly permissive, only stops if merging would result in broken or nonsense
        # DataSets, but does not protect against many types of input errors
        # first make sure parameters that absolutely must match actually do
        for elem in dats[1:]:
            if (elem.framerange[2] != test_stride) \
              or (elem.is_mask != test_is_mask) \
              or (elem.test_feature_list_type != elem.feature_list_type):
                self.logger.exit('Failed to merge DataSets due to property mismatch! Exiting...')
        # conditional checks
        for elem in dats[1:]:
            # static feature lists checks
            if test_feature_list_type == 'static':
                if elem.count_measurements() != test_n_features:
                    self.logger.exit('Failed to merge DataSets due to different feature counts! Exiting...')
                elif (elem.get_width() != test_width) and (test_n_features != 1):
                    self.logger.exit('Failed to merge DataSets due to different data widths and feature count > 1! Exiting...')
                elif (elem.get_width() != test_width) and (test_n_features == 1):
                    self.logger.msg('Warning, data width mismatch! This is expected for some DataSets as long as the feature count = 1.\nArrays will be padded with \'nan\' as needed to match.')
                    pad_with_nan = True
            elif test_feature_list_type == 'dynamic':
                pass
            else:
                self.logger.exit('Cannot detect feature_list_type! Exiting...')
            # check if last frame of preceding file matches first frame of this one
            if last_frame != elem.framerange[0]:
                self.logger.msg('Warning! Apparent overlap or missing frames during merge!')
            last_frame = elem.framerange[1]
        # empty default DataSets should have been pre-initialized by __init__()
        # copy the metadata from the first DataSet in the list
        self.primaryDS.copy_metadata(dats[0])
        self.primaryDS.is_mask = test_is_mask
        # get the full framerange
        self.primaryDS.framerange = (first_frame, last_frame, test_stride)
        # at this point have to diverge and use different code depending on 'static' vs 'dynamic'
        if test_feature_list_type == 'static':
            # copy measurement objects
            for meas in dats[0].measurements:
                self.primaryDS.add_measurement((meas.name, meas.type, meas.selecttext), meas.width)
            # add padding (if needed), concatenate, and load data and time arrays
            # have to use np.nan because 0 could be an actual data values in float datasets
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
        # this is much more complicated, and probably more dangerous too
        if test_feature_list_type == 'dynamic':
            # collect all measurement names into a set and figure out total size
            featureset = set()
            total_length = 0
            for elem in dats:
                for meas in elem.measurements:
                    featureset.add((meas.name, meas.type, meas.selecttext, meas.get_width()))
                total_length += elem.get_length()
            # get the final width of the array
            total_width = np.sum([i[3] for i in featureset])
            # sort features by name and assign indices, this will be the final ordering
            feature_dict = {}
            start_idx = 0
            for elem in sorted(featureset, key=lambda x: x[0]):
                feature_dict[elem[0]] = ((elem[0], elem[1], elem[2], start_idx, start_idx+elem[3]))
                start_idx += elem[3]
            # setup final array
            if test_is_mask is True:
                tmparray = np.zeros((total_width,total_length), dtype=int)
            else:
                tmparray = np.empty((total_width,total_length), dtype=float)
                tmparray.fill(np.nan)
            # iterate over datasets to load the array
            start_time_idx = 0
            for elem in dats:
                end_time_idx = start_time_idx + elem.get_length()
                tmp_dict = elem.dataset_to_dict()
                for key in tmp_dict:
                    if key in feature_dict.keys():
                        start_width_idx = feature_dict[key][3]
                        end_width_idx = feature_dict[key][4]
                    tmparray[start_width_idx:end_width_idx,start_time_idx:end_time_idx] = tmp_dict[key]
                start_time_idx = end_time_idx
            # finish up
            for elem in sorted(feature_dict.values(), key=lambda x: x[3]):
                self.primaryDS.add_measurement((elem[0], elem[1], elem[2]), width=elem[4]-elem[3])
            self.primaryDS.add_collection(tmparray)
            self.primaryDS.setup_timesteps()

    def merge_along_features(self, datlist, strict_checking=False):
        """Load a list of datasets that cover the same time window of the same trajectory and merge
        them along the features axis."""
        pass # not yet implemented
