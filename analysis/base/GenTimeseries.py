import numpy as np
from MDAnalysis.analysis.base import AnalysisBase

class GenTimeseries(AnalysisBase):
    """
    Measures a variety of properties on a single pass through a trajectory. Similar in principle
    to the DCD only (and now deprecated) MDAnalysis.core.Timeseries.TimeseriesCollection module, but
    works with any trajectory type (including pdb files) and uses the MDAnalysis AnalysisBase API.

    NOTE: In order to provide 'TimeseriesCollection'-like functionality (i.e. simultaneous
    generation of different feature types), this class explicitly relies on the
    ExtendedTSC.TimeseriesDataSet API, and so it will not work with MDAnalysis alone.
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
