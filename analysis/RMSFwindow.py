from __future__ import print_function
import sys, os
import numpy as np
from scipy.signal import convolve as sp_conv
# tested and working with MDAnalysis-0.16.1
from MDAnalysis.analysis.base import AnalysisFromFunction
from core.TimeseriesCore import TimeseriesCore

class RMSFwindow(TimeseriesCore):
    """
    Calculate root-mean-squared fluctuations of specified atoms in a sliding window across a
    trajectory. Input trajectory should be pre-aligned.
    """

    def __init__(self,**kwargs):
        super(RMSFwindow,self).__init__(**kwargs)
        self.primaryDS.set_static()

    def run(self, selecttext, windowsize):

        if self.input_type == None:
            sys.exit('No data has been loaded, cannot run! Exiting...')
        # setup a temporary dataset
        tmpDS = self._custom_dataset()
        tmpDS.copy_metadata(self.primaryDS)
        # create feature descriptors
        self.descriptor_list = []
        atoms = self.u.select_atoms(selecttext)
        for atom in atoms:
            tmpDS.add_feature(('%s_%s_%s' % (atom.segid, atom.resid, atom.name), 'rmsf_over_%d_frame_window' % windowsize, 'segid %s and resid %s and name %s' % (atom.segid, atom.resid, atom.name)), width=3)
        position_data = AnalysisFromFunction(lambda x: x.positions, atoms).run().results
        # hold the coordinates in the temporary dataset
        tmpDS.format_data(position_data.reshape(position_data.shape[0], -1).T)
        # use array convolution to get rolling averages - set up convolution array
        convolution = np.empty((windowsize,), dtype=np.float64)
        convolution.fill(1.0)
        convolution = convolution / float(windowsize)
        # average time values
        time_conv = sp_conv(tmpDS['time'], convolution, mode='valid')
        # average coordinates
        running_avgs = []
        for i in range(tmpDS.data.shape[0]):
            i_conv = sp_conv(tmpDS.data[i,:],convolution, mode='valid')
            running_avgs.append(i_conv)
        running_avg_array = np.stack(running_avgs, axis=0)
        # setup array to hold final values
        final_data = np.empty((running_avg_array.shape[0]/3, running_avg_array.shape[1]))
        # iterate over time windows
        for i in range(running_avg_array.shape[1]):
            # broadcasting doesn't work here, so make a stretched array of avg coords manually
            i_avg_coord = np.repeat(running_avg_array[:,i].reshape(1,-1).T, windowsize, axis=1)
            window_coords = tmpDS.data[:,i:i+windowsize]
            diff_coords = window_coords - i_avg_coord
            # iterate over atoms and calculate rmsf
            for j in range(running_avg_array.shape[0]/3):
                j_diff_coords = diff_coords[j*3:j*3+3,:]
                disp = np.linalg.norm(j_diff_coords, axis=0)
                disp_sq = np.square(disp)
                ij_rmsf = np.sqrt(np.mean(disp_sq))
                final_data[j,i] = ij_rmsf
        # load up the dataset
        for elem in tmpDS.feature_list:
            self.primaryDS.add_feature(elem.get_descriptor(), width=1)
        self.primaryDS.format_data(final_data)
        self.primaryDS.time = time_conv
        self.primaryDS.feature_dict['time'] = self.primaryDS.time
