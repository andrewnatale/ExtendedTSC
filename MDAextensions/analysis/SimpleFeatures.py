#!/usr/bin/env python2
from __future__ import print_function,division
import sys
import numpy as np
# tested and working with MDAnalysis-0.16.1
from MDAnalysis.analysis.base import AnalysisBase
from MDAextensions.datatools.TimeseriesCore import TimeseriesCore
from MDAextensions.datatools.CustomErrors import LoadError,AnalysisRuntimeError

class SimpleFeatures(TimeseriesCore):
    """
    Class for analyzers that need to make certain simple types of features using GenTimeseries
    (works like the deprecated TimeseriesCollection module). Takes as input a list of selection
    descriptors. When subclassed, the run() method can be redefined to generate _Feature objects
    in some way other than loading a list prior to generating a timeseries."""

    def run(self,selection_list):
        """
        Arguments:
        selection_list - list; a list of 3-string tuples of the form:
            (name, type, selectext);
            name is a identifier
            type is a code recognized by the run() method
            selecttext is an MDAnalysis format atom selector expression
        """

        if self.input_type == None:
            raise LoadError(1)
        # setup primaryDS using selections from a list
        for descriptor in selection_list:
            self.primaryDS.add_feature(descriptor)
        self.primaryDS.set_static()
        self._generate_timeseries()

    def _generate_timeseries(self):
        """Featurize a trajectory based on a loaded or generated selection list, then
        poulate a DataSet object with the results."""

        # check that the primary DataSet object has been properly initialized
        if len(self.primaryDS) == 0:
            raise AnalysisRuntimeError(0)
        elif self.primaryDS.populated is True:
            raise AnalysisRuntimeError(1)
        # check that feature selections behave as expected on the input topology, and set widths
        for feature in self.primaryDS.feature_list:
            if feature.type == 'atom':
                tmpselect = self.u.select_atoms(feature.selecttext)
                if tmpselect.n_atoms != 1:
                    raise AnalysisRuntimeError(2, feature.selecttext)
                feature.set_width(3)
            elif feature.type == 'bond':
                # not implemented
                raise AnalysisRuntimeError(3, feature.type)
            elif feature.type == 'angle':
                # not implemented
                raise AnalysisRuntimeError(3, feature.type)
            elif feature.type == 'dihedral':
                tmpselect = self.u.select_atoms(feature.selecttext)
                if tmpselect.n_atoms != 4:
                    raise AnalysisRuntimeError(2, feature.selecttext)
                feature.set_width(1)
            elif feature.type == 'distance':
                tmpselect = self.u.select_atoms(feature.selecttext)
                if tmpselect.n_atoms != 2:
                    raise AnalysisRuntimeError(2, feature.selecttext)
                feature.set_width(1)
            elif feature.type == 'COG':
                tmpselect = self.u.select_atoms(feature.selecttext)
                if tmpselect.n_atoms == 0:
                    raise AnalysisRuntimeError(2, feature.selecttext)
                feature.set_width(3)
            elif feature.type == 'COM':
                tmpselect = self.u.select_atoms(feature.selecttext)
                if tmpselect.n_atoms == 0:
                    raise AnalysisRuntimeError(2, feature.selecttext)
                feature.set_width(3)
            elif feature.type == 'water_dipole':
                # not implemented
                raise AnalysisRuntimeError(3, feature.type)
            else:
                raise AnalysisRuntimeError(3, feature.type)
        # generate the array
        if self.primaryDS.framerange is None:
            collection = GenTimeseries(self.primaryDS.feature_list, self.u, verbose=True)
        else:
            collection = GenTimeseries(self.primaryDS.feature_list, self.u, verbose=True,
              start=self.primaryDS.framerange[0], stop=self.primaryDS.framerange[1], step=self.primaryDS.framerange[2])
        collection.run()
        # process the array into a DataSet object
        self.primaryDS.format_data(collection.data)

class GenTimeseries(AnalysisBase):
    """
    Measures a variety of properties on a single pass through a trajectory. Similar in principle
    to the DCD only (and now deprecated) MDAnalysis.core.Timeseries.TimeseriesCollection module, but
    works with any trajectory type (including pdb files) and uses the MDAnalysis AnalysisBase API.

    NOTE: In order to provide 'TimeseriesCollection'-like functionality (i.e. simultaneous
    generation of different feature types), this class explicitly relies on the
    MDAextensions.datatools.TimeseriesDataSet._Feature API, and so it will not work with MDAnalysis alone.
    """

    def __init__(self,feature_list,universe,**kwargs):
        super(GenTimeseries,self).__init__(universe.trajectory,**kwargs)
        self.target_features = feature_list
        self.u = universe

    def _prepare(self):
        self.atomgroups = {}
        self.coordinates = {}
        total_width = 0
        # setup empty arrays
        for feature in self.target_features:
            total_width += feature.width
            self.atomgroups[feature.name] = self.u.select_atoms(feature.selecttext)
            self.coordinates[feature.name] = np.empty((self.atomgroups[feature.name].n_atoms * 3, self.n_frames), dtype=float)
        self.data = np.empty((total_width,self.n_frames), dtype=float)

    def _single_frame(self):
        for key in self.atomgroups:
            self.coordinates[key][:,self._frame_index] = np.concatenate([atom.position for atom in self.atomgroups[key]], axis=0)

    def _conclude(self):
        start_feature_idx = 0
        for feature in self.target_features:
            end_feature_idx = start_feature_idx + feature.width
            if feature.type == 'atom':
                # just load the coordinates
                self.data[start_feature_idx:end_feature_idx,:] = self.coordinates[feature.name]
            elif feature.type == 'bond':
                # not implemented
                pass
            elif feature.type == 'angle':
                # not implemented
                pass
            elif feature.type == 'dihedral':
                self.data[start_feature_idx:end_feature_idx,:] = self._dihedral(self.coordinates[feature.name])
            elif feature.type == 'distance':
                self.data[start_feature_idx:end_feature_idx,:] = self._distance(self.coordinates[feature.name])
            elif feature.type == 'COG':
                count = self.atomgroups[feature.name].n_atoms
                self.data[start_feature_idx:end_feature_idx,:] = self._center_of_geometry(count, self.coordinates[feature.name])
            elif feature.type == 'COM':
                count = self.atomgroups[feature.name].n_atoms
                weights = [atom.mass for atom in self.atomgroups[feature.name]]
                self.data[start_feature_idx:end_feature_idx,:] = self._center_of_mass(count, self.coordinates[feature.name], weights)
            elif feature.type == 'water_dipole':
                # not implemented
                pass
            start_feature_idx += feature.width

    def _dihedral(self, array):
        # expect to get 4 points
        p1 = array[0:3,:]
        p2 = array[3:6,:]
        p3 = array[6:9,:]
        p4 = array[9:,:]
        segA = p2 - p1
        segB = p3 - p2
        segC = p4 - p3
        AxB = np.cross(segA, segB, axisa=0, axisb=0, axisc=0)
        BxC = np.cross(segB, segC, axisa=0, axisb=0, axisc=0)
        n1 = AxB / np.linalg.norm(AxB, axis=0)
        n2 = BxC / np.linalg.norm(BxC, axis=0)
        x = np.sum(n1 * n2, axis=0)
        y = np.sum(np.cross(n1, segB / np.linalg.norm(segB, axis=0), axisa=0, axisb=0, axisc=0) * n2, axis=0)
        # flip the sign to get the same thing TimeseriesCollection gets ???
        return np.arctan2(y,x) * -1.0

    def _distance(self, array):
        # expect to get 2 points
        p1 = array[0:3,:]
        p2 = array[3:,:]
        return np.linalg.norm(p1-p2, axis=0)

    def _center_of_geometry(self, count, array):
        cog = array[0:3,:]
        for i in range(1,count,1):
            cog += array[i*3:(i*3)+3,:]
        return cog / float(count)

    def _center_of_mass(self, count, array, weights):
        com = array[0:3,:] * weights[0]
        for i in range(1,count,1):
            com += array[i*3:(i*3)+3,:] * weights[i]
        return com / float(np.sum(weights))
