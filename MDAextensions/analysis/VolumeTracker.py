#!/usr/bin/env python2
from __future__ import print_function,division
import sys
import numpy as np
# tested and working with MDAnalysis-0.16.1
from MDAnalysis.analysis.base import AnalysisBase
from MDAextensions.analysis.SimpleFeatures import SimpleFeatures

class VolumeTracker(SimpleFeatures):
    """
    Dynamicaly generate a list of atoms or residues using a volume built from MDAnalysis
    geometric selection expressions, then measure their positions at every timestep. Additionaly
    populate a mask DataSet object with interger values for each found entity at each timestep
    indicating its presence (1) or abscence (0) in the search volume.
    """

    def run(self,vol_selecttext,search_selecttext):
        """
        Arguments:
        vol_selecttext - string; MDAnalysis geometric selection expression; can be multiple volumes
            chained with and/or
        search_selecttext - string; MDAnalysis atom selection expression; defines the atom types to
            be searched for
        """

        if self.input_type == None:
            sys.exit('No data has been loaded, cannot run! Exiting...')
        # TODO: input type check, shouldn't process pdbfiles
        self.maskDS.copy_metadata(self.primaryDS)
        self.primaryDS.set_dynamic()
        self.maskDS.set_dynamic()
        if self.primaryDS.framerange is None:
            searcher = _VolumeSearch(vol_selecttext, search_selecttext, self.u, verbose=True)
        else:
            searcher = _VolumeSearch(vol_selecttext, search_selecttext, self.u, verbose=True,
              start=self.primaryDS.framerange[0], stop=self.primaryDS.framerange[1], step=self.primaryDS.framerange[2])
        searcher.run()
        # setup data sets using search results
        for descriptor in searcher.selections:
            self.primaryDS.add_feature(descriptor)
        for descriptor in searcher.selections_mask:
            # occupancy measure types must have width=1
            self.maskDS.add_feature(descriptor,width=1)
        # load the masking data from the search
        self.maskDS.format_data(searcher.mask)
        # now go get the coordinates
        self._generate_timeseries()

class _VolumeSearch(AnalysisBase):
    """
    Class for steping through a trajectory frame by frame and tracking the
    atoms that enter a specified volume.
    """

    def __init__(self,vol_selecttext,search_selecttext,universe,**kwargs):
        super(_VolumeSearch,self).__init__(universe.trajectory,**kwargs)
        self.vol_selecttext = vol_selecttext
        self.search_selecttext = search_selecttext
        self.u = universe

    def _prepare(self):
        # setup vars and data structures
        self.vol_group = self.u.select_atoms('(%s) and (%s)' % (self.vol_selecttext, self.search_selecttext), updating=True)
        self.selection_set = set()
        self.masking_list = []

    def _single_frame(self):
        # what to do at each frame
        tmpoccupancy = []
        for atom in self.vol_group:
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
            #count += 1
            tmp_name = '%s_%s_%s_%s' % (str(elem[0]),str(elem[1]),str(elem[2]),str(elem[3]))
            tmp_type = 'atom'
            # select by index, because segid/resid may not be unique
            tmp_selecttext = 'bynum %s' % str(elem[0])
            self.selections.append((tmp_name,tmp_type,tmp_selecttext))
            self.selections_mask.append((tmp_name,'occupancy',tmp_selecttext))
        # build occupancy mask array
        self.mask = np.zeros((len(self.selection_set),len(self.masking_list)), dtype=int)
        # re-iterate over selection set and check occupancy at each timestep
        for idxA, elem in enumerate(sorted(self.selection_set)):
            identity = elem[0]
            for idxB, occupancy in enumerate(self.masking_list):
                if identity in occupancy:
                    self.mask[idxA,idxB] = 1
