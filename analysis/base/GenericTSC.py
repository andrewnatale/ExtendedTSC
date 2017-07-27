import numpy as np
from MDAnalysis.analysis.base import AnalysisBase

class _GenericTSC(AnalysisBase):
    """Measures a variety of properties on a single pass through a trajectory. Similar in principle
    to MDAnalysis.core.Timeseries.TimeseriesCollection, but works with any trajectory type
    (including pdb files) and uses the AnalysisBase API."""

    def __init__(self,measurements,universe,**kwargs):
        super(_GenericTSC,self).__init__(universe.trajectory,**kwargs)
        self.target_measurements = measurements
        self.u = universe

    def _prepare(self):
        self.atomgroups = {}
        self.coordinates = {}
        total_width = 0
        # setup empty arrays
        for meas in self.target_measurements:
            total_width += meas.width
            self.atomgroups[meas.name] = self.u.select_atoms(meas.selecttext)
            self.coordinates[meas.name] = np.empty((self.atomgroups[meas.name].n_atoms * 3, self.n_frames), dtype=float)
        self.data = np.empty((total_width,self.n_frames), dtype=float)

    def _single_frame(self):
        #print(self._ts)
        for key in self.atomgroups:
            self.coordinates[key][:,self._frame_index] = np.concatenate([atom.position for atom in self.atomgroups[key]], axis=0)

    def _conclude(self):
        start_meas_idx = 0
        for meas in self.target_measurements:
            end_meas_idx = start_meas_idx + meas.width
            if meas.type == 'atom':
                # just load the coordinates
                self.data[start_meas_idx:end_meas_idx,:] = self.coordinates[meas.name]
            elif meas.type == 'bond':
                # not implemented
                pass
            elif meas.type == 'angle':
                # not implemented
                pass
            elif meas.type == 'dihedral':
                self.data[start_meas_idx:end_meas_idx,:] = self._dihedral(self.coordinates[meas.name])
            elif meas.type == 'distance':
                self.data[start_meas_idx:end_meas_idx,:] = self._distance(self.coordinates[meas.name])
            elif meas.type == 'COG':
                count = self.atomgroups[meas.name].n_atoms
                self.data[start_meas_idx:end_meas_idx,:] = self._center_of_geometry(count, self.coordinates[meas.name])
            elif meas.type == 'COM':
                count = self.atomgroups[meas.name].n_atoms
                weights = [atom.mass for atom in self.atomgroups[meas.name]]
                self.data[start_meas_idx:end_meas_idx,:] = self._center_of_mass(count, self.coordinates[meas.name], weights)
            elif meas.type == 'water_dipole':
                # not implemented
                pass
            start_meas_idx += meas.width

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
