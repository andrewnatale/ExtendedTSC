#!/usr/bin/env python2
from __future__ import print_function,division

class LoadError(Exception):

    errmsgs = {
    0: 'Can only handle one trajectory or structure per instance! Exiting...',
    1: 'No data has been loaded, cannot run! Exiting...',
    }

    def __init__(self, errtype, bad_input=None):
        self.errtype = errtype
        self.input = bad_input
        if self.input is not None:
            msg = '%s %s' % (self.errmsgs[self.errtype], self.input)
        else:
            msg = self.errmsgs[self.errtype]
        super(LoadError, self).__init__(msg)

class DataSetReadError(Exception):

    errmsgs = {
    0: 'Likely broken input file! Last line is not \">end\"! Bad file:',
    1: 'Data line found in file header, check input and try again! Bad line:',
    2: 'Cannot determine if file is float data or integer mask! Bad file:',
    3: 'DataSet file version cannot be processed! Unsupported version:',
    4: 'Failed to merge DataSets due to property mismatch! Exiting...',
    5: 'Cannot detect feature_list_type for DataSet:',
    }

    def __init__(self, errtype, bad_input=None):
        self.errtype = errtype
        self.input = bad_input
        if self.input is not None:
            msg = '%s %s' % (self.errmsgs[self.errtype], self.input)
        else:
            msg = self.errmsgs[self.errtype]
        super(DataSetReadError, self).__init__(msg)

class AnalysisRuntimeError(Exception):

    errmsgs = {
    0: 'No feature descriptors found in DataSet. Exiting...',
    1: 'DataSet object already contains an array, cannot generate another! Exiting...',
    2: 'Bad selection, wrong number of atoms:',
    3: 'Feature not implemented:',
    }

    def __init__(self, errtype, bad_input=None):
        self.errtype = errtype
        self.input = bad_input
        if self.input is not None:
            msg = '%s %s' % (self.errmsgs[self.errtype], self.input)
        else:
            msg = self.errmsgs[self.errtype]
        super(AnalysisRuntimeError, self).__init__(msg)
