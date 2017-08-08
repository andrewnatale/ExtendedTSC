from __future__ import print_function
import sys
import numpy as np
# tested and working with MDAnalysis-0.16.1
#import MDAnalysis.core.Timeseries as tm # will be deprecated in MDAnalysis 0.17.0
from core.TimeseriesCore import TimeseriesCore
from base.GenTimeseries import GenTimeseries

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
            sys.exit('No data has been loaded, cannot run! Exiting...')
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
            sys.exit('No feature descriptors found in DataSet. Exiting...')
        elif self.primaryDS.populated is True:
            sys.exit('DataSet object already contains an array, cannot generate another! Exiting...')
        # check that feature selections behave as expected on the input topology, and set widths
        for feature in self.primaryDS.feature_list:
            if feature.type == 'atom':
                tmpselect = self.u.select_atoms(feature.selecttext)
                if tmpselect.n_atoms != 1:
                    sys.exit('Selection %s \"%s\" found %d atoms instead of 1!\nFix the selection and try again. Exiting now...' \
                      % (feature.type, feature.selecttext, tmpselect.n_atoms))
                feature.set_width(3)
            elif feature.type == 'bond':
                # not implemented
                sys.exit('feature type %s not implemented, exiting' % feature.type)
            elif feature.type == 'angle':
                # not implemented
                sys.exit('feature type %s not implemented, exiting' % feature.type)
            elif feature.type == 'dihedral':
                tmpselect = self.u.select_atoms(feature.selecttext)
                if tmpselect.n_atoms != 4:
                    sys.exit('Selection %s \"%s\" found %d atoms instead of 4!\nFix the selection and try again. Exiting now...' \
                      % (feature.type, feature.selecttext, tmpselect.n_atoms))
                feature.set_width(1)
            elif feature.type == 'distance':
                tmpselect = self.u.select_atoms(feature.selecttext)
                if tmpselect.n_atoms != 2:
                    sys.exit('Selection %s \"%s\" found %d atoms instead of 2!\nFix the selection and try again. Exiting now...' \
                      % (feature.type, feature.selecttext, tmpselect.n_atoms))
                feature.set_width(1)
            elif feature.type == 'COG':
                tmpselect = self.u.select_atoms(feature.selecttext)
                if tmpselect.n_atoms == 0:
                    sys.exit('Selection %s \"%s\" found 0 atoms!\nFix the selection and try again. Exiting now...'\
                      % (feature.type, feature.selecttext))
                feature.set_width(3)
            elif feature.type == 'COM':
                tmpselect = self.u.select_atoms(feature.selecttext)
                if tmpselect.n_atoms == 0:
                    sys.exit('Selection %s \"%s\" found 0 atoms!\nFix the selection and try again. Exiting now...'\
                      % (feature.type, feature.selecttext))
                feature.set_width(3)
            elif feature.type == 'water_dipole':
                # not implemented
                sys.exit('feature type %s not implemented, exiting' % feature.type)
            else:
                sys.exit('unrecognized measure type %s! Exiting...' % feature.type)
        # generate the array
        if self.primaryDS.framerange is None:
            collection = GenTimeseries(self.primaryDS.feature_list, self.u, verbose=True)
        else:
            collection = GenTimeseries(self.primaryDS.feature_list, self.u, verbose=True,
              start=self.primaryDS.framerange[0], stop=self.primaryDS.framerange[1], step=self.primaryDS.framerange[2])
        collection.run()
        # process the array into a DataSet object
        self.primaryDS.format_data(collection.data)
