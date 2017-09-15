#!/usr/bin/env python2
from __future__ import print_function,division
import sys, os, datetime
import re
import numpy as np
from MDAextensions.datatools.TimeseriesDataSet import TimeseriesDataSet as tsds

def merge_along_time(datlist):
    """
    Load a list of dataset files that contain the same features and merge them along the time
    axis to generate one longer timeseries. Automatically sorts the input files by starting
    frame.

    NOTE: For specific types of datasets, features do not always have to be exactly the same.
    This function should be able to determine whether a merge operation is permissible based
    on dataset metadata.

    Arguments:
    datlist - list of strings; paths to dat files to merge
    """

    # control flags
    pad_with_nan = False
    # load specified files into lists of DataSet objects
    dats = []
    for filename in datlist:
        print('Loading %s for merge operation...' % filename)
        newDataSet = tsds()
        newDataSet._read(filename, enforce_version=True)
        dats.append(newDataSet)
    # sort by start frame
    sorted_dats = sorted(dats, key=lambda ds: ds.framerange[0])
    dats = sorted_dats
    # extract parameters from the first DataSet in the list and compare to all the others
    test_width = dats[0].get_width()
    test_n_features = len(dats[0])
    first_frame = dats[0].framerange[0]
    last_frame = dats[0].framerange[1]
    test_stride = dats[0].framerange[2]
    test_is_mask = dats[0].is_mask
    test_feature_list_type = dats[0].feature_list_type
    # checking is fairly permissive, only stops if merging would result in totally broken DataSets
    # DOES NOT protect against all types of input errors
    # first make sure parameters that absolutely must match actually do
    for ds in dats[1:]:
        if (ds.framerange[2] != test_stride) \
          or (ds.is_mask != test_is_mask) \
          or (ds.feature_list_type != test_feature_list_type):
            sys.exit('Failed to merge DataSets due to property mismatch! Exiting...')
    # conditional checks
    for ds in dats[1:]:
        # static feature lists checks
        if test_feature_list_type == 'static':
            if len(ds) != test_n_features:
                sys.exit('Failed to merge DataSets due to different feature counts! Exiting...')
            elif (ds.get_width() != test_width) and (test_n_features != 1):
                sys.exit('Failed to merge DataSets due to different data widths and feature count > 1! Exiting...')
            elif (ds.get_width() != test_width) and (test_n_features == 1):
                print('Warning, data width mismatch! This is expected for some DataSets as long as the feature count = 1.\nArrays will be padded with \'nan\' as needed to match.')
                pad_with_nan = True
        elif test_feature_list_type == 'dynamic':
            pass
        else:
            print('Cannot detect feature_list_type! Exiting...')
            sys.exit(1)
        # check if last frame of preceding file matches first frame of this one
        if last_frame != ds.framerange[0]:
            print('Warning! Apparent overlap or missing frames during merge!')
        last_frame = ds.framerange[1]
    # create empty dataset and copy the metadata from the first DataSet in the list
    mergeDS = tsds()
    mergeDS.copy_metadata(dats[0])
    mergeDS.is_mask = test_is_mask
    # get the full framerange
    mergeDS.framerange = (first_frame, last_frame, test_stride)
    # at this point we must diverge and use different code depending on 'static' vs 'dynamic'
    if test_feature_list_type == 'static':
        # copy features
        for feature in dats[0].feature_list:
            mergeDS.add_feature((feature.name, feature.type, feature.selecttext), width=feature.width)
        # add padding (if needed), concatenate, load data, and generate time array
        # have to use np.nan because 0 could be an actual data values in float datasets
        if pad_with_nan is True:
            all_widths = [i.get_width() for i in dats]
            maxwidth = np.amax(all_widths)
            mergeDS.feature_list[0].set_width(maxwidth)
            for ds in dats:
                width_diff = maxwidth - ds.get_width()
                if width_diff > 0:
                    padding = np.empty((width_diff,ds.get_length()), dtype=float)
                    padding.fill(np.nan)
                    ds.data = np.concatenate([ds.data, padding], axis=0)
                elif width_diff == 0:
                    pass
        mergeDS.format_data(np.concatenate([ds.data for ds in dats], axis=1))
    # this is much more complicated, and probably more bug prone, but it seems to work ok
    elif test_feature_list_type == 'dynamic':
        # collect all feature names into a set and figure out total size
        featureset = set()
        total_length = 0
        for ds in dats:
            for feature in ds.feature_list:
                featureset.add((feature.name, feature.type, feature.selecttext, feature.width))
            total_length += ds.get_length()
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
        for ds in dats:
            end_time_idx = start_time_idx + ds.get_length()
            for key in ds.keys():
                if key in feature_dict.keys():
                    start_width_idx = feature_dict[key][3]
                    end_width_idx = feature_dict[key][4]
                    tmparray[start_width_idx:end_width_idx,start_time_idx:end_time_idx] = ds[key][:]
            del start_width_idx
            del end_width_idx
            start_time_idx = end_time_idx
        # finish up
        for elem in sorted(feature_dict.values(), key=lambda x: x[3]):
            mergeDS.add_feature((elem[0], elem[1], elem[2]), width=elem[4]-elem[3])
        mergeDS.format_data(tmparray)
    # return the merged DataSet
    return mergeDS

def merge_along_features(datlist):
    """
    Load a list of datasets that cover the same time window of the same trajectory and merge
    them along the features axis.

    Arguments:
    datlist - list of strings; paths to dat files to merge
    """

    # load specified files into lists of DataSet objects
    dats = []
    for filename in datlist:
        print('Loading %s for merge operation...' % filename)
        newDataSet = tsds()
        newDataSet._read(filename, enforce_version=True)
        dats.append(newDataSet)
    # extract parameters from the first DataSet in the list and compare to all the others
    test_length = dats[0].get_length()
    first_frame = dats[0].framerange[0]
    last_frame = dats[0].framerange[1]
    test_stride = dats[0].framerange[2]
    test_is_mask = dats[0].is_mask
    test_feature_list_type = dats[0].feature_list_type
    # checking is fairly permissive, only stops if merging would result in totally broken DataSets
    # DOES NOT protect against all types of input errors
    # first make sure parameters that absolutely must match actually do
    for ds in dats[1:]:
        if (ds.get_length() != test_length) \
          or (ds.framerange[0] != first_frame) \
          or (ds.framerange[1] != last_frame) \
          or (ds.framerange[2] != test_stride) \
          or (ds.is_mask != test_is_mask) \
          or (ds.feature_list_type != test_feature_list_type):
            sys.exit('Failed to merge DataSets due to property mismatch! Exiting...')
    # create empty dataset and copy the metadata from the first DataSet in the list
    mergeDS = tsds()
    mergeDS.copy_metadata(dats[0])
    mergeDS.is_mask = test_is_mask
    # now merge - much easier than timewise merge
    for ds in dats:
        for feature in ds.feature_list:
            mergeDS.add_feature((feature.name, feature.type, feature.selecttext), width=feature.width)
    mergeDS.format_data(np.concatenate([ds.data for ds in dats], axis=0))
    return mergeDS
