import numpy as np
from TrajProcessor import TrajProcessor
# tested and working with MDAnalysis-0.16.1
import MDAnalysis.core.Timeseries as tm
from MDAnalysis.analysis.base import AnalysisBase

class ExtendedTSC(TrajProcessor):
    """Wrap and extend MDAnalysis' TimeseriesCollection module."""

    def measures_from_list(self,selection_list):
        """Load selection definitions. This method can be used alone or alongside
    measures_from_volumesearch() to add reference measurements to be made alongside the
    measurements defined by the volume search.

    Arguments:
    selection_list - list; a list of 3-string tuples of the form:
        (name, type, selectext);
        name is a identifier
        type is a code recognized by the run() method
        selecttext is an MDAnalysis format atom selector expression
    """

        # setup primaryDS using selections from a list
        for descriptor in selection_list:
            self.primaryDS.add_measurement(descriptor)

    def measures_from_volumesearch(self,vol_selecttext,search_selecttext,mode='atom'):
        """Step through a trajectory finding all the atoms/residues within a defined volume and
    create a list for measuring their positions. Populate a mask DataSet object with binary
    values for each found entity at each timestep indicating its presence in search volume.

    Arguments:
    vol_selecttext - string; MDAnalysis geometric selection expression; can be multiple volumes
        chained with and/or
    search_selecttext - string; MDAnalysis atom selection expression; defines the atom types to
        be searched for

    Keyword arguments:
    mode - string; either 'res' or 'atom'; format the resulting measurement list to give
        coordinates for each found atom, or the center of mass of found residues
    """

        # TODO: input type check, shouldn't process pdbfiles
        self.maskDS._copy_metadata(self.primaryDS)
        if self.primaryDS.framerange is None:
            searcher = _VolumeSearch(vol_selecttext,search_selecttext,mode,self.u)
        else:
            searcher = _VolumeSearch(vol_selecttext,search_selecttext,mode,self.u,start=self.primaryDS.framerange[0],stop=self.primaryDS.framerange[1],step=self.primaryDS.framerange[2])
        print searcher.start,searcher.stop,searcher.step,searcher.n_frames
        searcher.run()
        # setup data sets using search results
        for descriptor in searcher.selections:
            self.primaryDS.add_measurement(descriptor)
        for descriptor in searcher.selections_mask:
            # occupancy measure types must have width=1
            self.maskDS.add_measurement(descriptor,width=1)
        # load the masking data from the search
        self.maskDS.add_collection(searcher.mask)
        self.maskDS._setup_timesteps()

    def run(self):
        """Make measurements on a trajectory based on a loaded or generated selection list, then
    poulate a DataSet object with the results."""

        # check that the primary DataSet object has been properly initialized
        if not self.primaryDS:
            self.logger.err('DataSet not properly initialized! Exiting...')
        elif self.primaryDS.count_measurements() == 0:
            self.logger.err('No measurement descriptors found in DataSet. Exiting...')
        elif self.primaryDS.populated is True:
            self.logger.err('DataSet object already contains an array, cannot generate another! Exiting...')
        # check that measurement objects behave as expected on the input topology, and set measure widths
        for meas in self.primaryDS.measurements:
            if meas.type == 'atom':
                tmpselect = self.u.select_atoms(meas.selecttext)
                # strictly speaking, TimeseriesCollection can handle having more than one atom
                # in an 'atom' selection - it just measures coordinates for all of them
                # however to keep things simpler downstream, enforce one atom selections here
                if tmpselect.n_atoms != 1:
                    self.logger.err('Selection %s \"%s\" found %d atoms instead of 1!\nFix the selection and try again. Exiting now...' \
                      % (meas.type, meas.selecttext, tmpselect.n_atoms))
                meas.set_width(3)
            elif meas.type == 'bond':
                # not implemented
                self.logger.err('measurement type %s not implemented, exiting' % meas.type)
            elif meas.type == 'angle':
                # not implemented
                self.logger.err('measurement type %s not implemented, exiting' % meas.type)
            elif meas.type == 'dihedral':
                tmpselect = self.u.select_atoms(meas.selecttext)
                if tmpselect.n_atoms != 4:
                    self.logger.err('Selection %s \"%s\" found %d atoms instead of 4!\nFix the selection and try again. Exiting now...' \
                      % (meas.type, meas.selecttext, tmpselect.n_atoms))
                meas.set_width(1)
            elif meas.type == 'distance':
                tmpselect = self.u.select_atoms(meas.selecttext)
                if tmpselect.n_atoms != 2:
                    self.logger.err('Selection %s \"%s\" found %d atoms instead of 2!\nFix the selection and try again. Exiting now...' \
                      % (meas.type, meas.selecttext, tmpselect.n_atoms))
                meas.set_width(1)
            elif meas.type == 'COG':
                tmpselect = self.u.select_atoms(meas.selecttext)
                if tmpselect.n_atoms == 0:
                    self.logger.err('Selection %s \"%s\" found 0 atoms!\nFix the selection and try again. Exiting now...'\
                      % (meas.type, meas.selecttext))
                meas.set_width(3)
            elif meas.type == 'COM':
                tmpselect = self.u.select_atoms(meas.selecttext)
                if tmpselect.n_atoms == 0:
                    self.logger.err('Selection %s \"%s\" found 0 atoms!\nFix the selection and try again. Exiting now...'\
                      % (meas.type, meas.selecttext))
                meas.set_width(3)
            elif meas.type == 'water_dipole':
                # not implemented
                self.logger.err('measurement type %s not implemented, exiting' % meas.type)
            else:
                self.logger.err('unrecognized measure type %s! Exiting...' % meas.type)
        # check input data format and call the appropriate method to generate the array
        if self.input_type == 'dcd_traj':
            tmp_array = self._DCD_timeseries()
        elif self.input_type == 'generic_traj':
            tmp_array = self._generic_timeseries()
        elif self.input_type == 'pdb':
            tmp_array = self._generic_timeseries()
        else:
            self.logger.err('Cannot determine input data format! Check trajectory initialization!')
        # process the array into a DataSet object
        self.primaryDS.add_collection(tmp_array)
        self.primaryDS._setup_timesteps()

    def _DCD_timeseries(self):
        """Use TimeseriesCollection to make measurements (fast but DCD only)."""

        collection = tm.TimeseriesCollection()
        # all the error checking should be done, so assume it's good and load it up
        for meas in self.primaryDS.measurements:
            if meas.type == 'atom':
                collection.addTimeseries(tm.Atom('v', self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'bond':
                # not implemented
                pass
            elif meas.type == 'angle':
                # not implemented
                pass
            elif meas.type == 'dihedral':
                collection.addTimeseries(tm.Dihedral(self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'distance':
                collection.addTimeseries(tm.Distance('r', self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'COG':
                collection.addTimeseries(tm.CenterOfGeometry(self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'COM':
                collection.addTimeseries(tm.CenterOfMass(self.u.select_atoms(meas.selecttext)))
            elif meas.type == 'water_dipole':
                # not implemented
                pass
        # compute Timeseries from universe trajectory and return the resulting array
        # NOTE: There's a bug in MDAnalysis' fast DCD reader code that causes it to drop the last
        # frame it should measure when slicing the trajectory in any way (not the same bug as noted above in
        # 'measures_from_volumesearch'). Haven't tracked it down so we just work around it - ???still active in 0.16.1???
        if self.primaryDS.framerange is None:
            collection.compute(self.u.trajectory)
        else:
            collection.compute(self.u.trajectory, start=self.primaryDS.framerange[0], stop=self.primaryDS.framerange[1], step=self.primaryDS.framerange[2])
        return collection.data

    def _generic_timeseries(self):
        """Use _GenericTSC to make measurements (any trajectory format)."""

        if self.primaryDS.framerange is None:
            collection = _GenericTSC(self.primaryDS.measurements,self.u)
        else:
            collection = _GenericTSC(self.primaryDS.measurements,self.u,start=self.primaryDS.framerange[0],stop=self.primaryDS.framerange[1],step=self.primaryDS.framerange[2])
        collection.run()
        return collection.data

class _VolumeSearch(AnalysisBase):
    """Class for steping through a trajectory frame by frame and tracking individual
    atoms/residues."""

    def __init__(self,vol_selecttext,search_selecttext,mode,universe,**kwargs):
        super(_VolumeSearch,self).__init__(universe.trajectory,**kwargs)
        self.vol_selecttext = vol_selecttext
        self.search_selecttext = search_selecttext
        self.u = universe
        # check mode
        valid_modes = ['res', 'atom']
        if mode in valid_modes:
            self.mode = mode
        else:
            errmsg = 'Invalid mode selected for volumetric search! Possible modes: %s.\nExiting...' % ' '.join(valid_modes)
            sys.exit(errmsg)

    def _prepare(self):
        # setup vars and data structures
        self.vol_group = self.u.select_atoms('(%s) and (%s)' % (self.vol_selecttext, self.search_selecttext), updating=True)
        self.selection_set = set()
        self.masking_list = []

    def _single_frame(self):
        # what to do at each frame
        print self._ts
        tmpoccupancy = []
        for atom in self.vol_group:
            # store unique identifiers for found res/atoms
            if self.mode == 'res':
                self.selection_set.add((atom.segid,atom.resid))
                tmpoccupancy.append((atom.segid,atom.resid))
            elif self.mode == 'atom':
                self.selection_set.add((atom.index,atom.segid,atom.resid,atom.name))
                tmpoccupancy.append(atom.index)
        self.masking_list.append(tmpoccupancy)

    def _conclude(self):
        # setup lists to hold descriptors
        self.selections = []
        self.selections_mask = []
        count = 0
        # compose selection descriptors
        for elem in sorted(self.selection_set):
            count += 1
            if self.mode == 'res':
                tmp_name = '%s_%s' % (str(elem[0]),str(elem[1]))
                tmp_type = 'COM'
                tmp_selecttext = 'segid %s and resid %s' % (str(elem[0]),str(elem[1]))
            elif self.mode == 'atom':
                tmp_name = '%s_%s_%s' % (str(elem[1]),str(elem[2]),str(elem[3]))
                tmp_type = 'atom'
                tmp_selecttext = 'segid %s and resid %s and name %s' % (str(elem[1]),str(elem[2]),str(elem[3]))
            self.selections.append((tmp_name,tmp_type,tmp_selecttext))
            self.selections_mask.append((tmp_name,'occupancy',tmp_selecttext))
        # build occupancy mask array
        self.mask = np.zeros((count,len(self.masking_list)), dtype=int)
        # re-iterate over selection set and check occupancy at each timestep
        for idxA, elem in enumerate(sorted(self.selection_set)):
            if self.mode == 'res':
                identity = elem
            elif self.mode == 'atom':
                identity = elem[0]
            for idxB, occupancy in enumerate(self.masking_list):
                if identity in occupancy:
                    try:
                        self.mask[idxA,idxB] = 1
                    except IndexError:
                        pass

class _GenericTSC(AnalysisBase):
    """Measures a variety of properties on a single pass through a trajectory. Similar in function
    to TimeseriesCollection, but works with any trajectory type (including pdb files) and uses
    the AnalysisBase API."""

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
        #print self._ts
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
