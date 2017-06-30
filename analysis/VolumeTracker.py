import numpy as np
# tested and working with MDAnalysis-0.16.1
from MDAnalysis.analysis.base import AnalysisBase
from SimpleFeatures import SimpleFeatures

class VolumeTracker(SimpleFeatures):
    """Dynamicaly generate a list of atoms or residues using a volume built from MDAnalysis
    geometric selection expressions, then measure their positions at every timestep. Additionaly
    populate a mask DataSet object with boolan values for each found entity at each timestep
    indicating its presence (1) or abscence (0) in the search volume."""

    def run(self,vol_selecttext,search_selecttext,mode='atom'):
        """Arguments:
    vol_selecttext - string; MDAnalysis geometric selection expression; can be multiple volumes
        chained with and/or
    search_selecttext - string; MDAnalysis atom selection expression; defines the atom types to
        be searched for

    Keyword arguments:
    mode - string; either 'res' or 'atom'; format the resulting measurement list to give
        coordinates for each found atom, or the center of mass coordinates of found residues
    """

        # TODO: input type check, shouldn't process pdbfiles
        self.maskDS.copy_metadata(self.primaryDS)
        self.primaryDS.set_dynamic()
        self.maskDS.set_dynamic()
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
        self.maskDS.setup_timesteps()
        self._generate_timeseries()

class _VolumeSearch(AnalysisBase):
    """Class for steping through a trajectory frame by frame and tracking the
    atoms or residues that enter a specified volume."""

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
