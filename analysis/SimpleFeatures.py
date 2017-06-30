import numpy as np
# tested and working with MDAnalysis-0.16.1
import MDAnalysis.core.Timeseries as tm
from core.TrajProcessor import TrajProcessor
from base.GenericTSC import _GenericTSC

class SimpleFeatures(TrajProcessor):
    """Reference class for analyzers that need to make certain simple types of measurements using
    either MDAnalysis.core.Timeseries.TimeseriesCollection (fast - for DCD files) or the builtin
    _GenericTSC (slower, but works with anuything MDAnalysis can read). Use a list of selection
    descriptors. When subclassed, the run() method can be redefined to generate measurement objects
    in some way other than loading a list before generating a timeseries."""

    def run(self,selection_list):
        """Arguments:
    selection_list - list; a list of 3-string tuples of the form:
        (name, type, selectext);
        name is a identifier
        type is a code recognized by the run() method
        selecttext is an MDAnalysis format atom selector expression
    """

        # setup primaryDS using selections from a list
        for descriptor in selection_list:
            self.primaryDS.add_measurement(descriptor)
        self.primaryDS.set_static()
        self._generate_timeseries()

    def _generate_timeseries(self):
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
        self.primaryDS.setup_timesteps()

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
