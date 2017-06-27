import numpy as np
# tested and working with MDAnalysis-0.16.1
from MDAnalysis.analysis.base import AnalysisBase
from core.TrajProcessor import TrajProcessor

class ZSearch(TrajProcessor):
    """Similar to _VolumeSearch, but saves only the z-coordinate of each atom found in the search
    volume at each timestep."""

    def run(self,vol_selecttext,search_selecttext):
        """Arguments:
    vol_selecttext - string; MDAnalysis geometric selection expression;
        can be multiple volumes chained with and/or
    search_selecttext - string; MDAnalysis atom selection expression;
        defines what to look for in the selected volume
    """

        if self.primaryDS.framerange is None:
            searcher = _SearchZ(vol_selecttext,search_selecttext,self.u)
        else:
            searcher = _SearchZ(vol_selecttext,search_selecttext,self.u,start=self.primaryDS.framerange[0],stop=self.primaryDS.framerange[1],step=self.primaryDS.framerange[2])
        searcher.run()
        # setup data set using search results
        for descriptor in searcher.selections:
            self.primaryDS.add_measurement(descriptor)
        self.primaryDS.add_collection(searcher.coordarray)
        self.primaryDS.setup_timesteps()
        self.primaryDS.measurements[0].set_width(np.shape(searcher.coordarray)[0])

class _SearchZ(AnalysisBase):
    # TODO: specify coordinates (xyz) to save so more than z can be used

    def __init__(self,vol_selecttext,search_selecttext,universe,**kwargs):
        super(_SearchZ,self).__init__(universe.trajectory,**kwargs)
        self.vol_selecttext = vol_selecttext
        self.search_selecttext = search_selecttext
        self.u = universe

    def _prepare(self):
        # setup vars and data structures
        self.zcoords = []
        self.vol_group = self.u.select_atoms('(%s) and (%s)' % (self.vol_selecttext, self.search_selecttext), updating=True)

    def _single_frame(self):
        # what to do at each frame
        print self._ts
        tmpcoordlist = []
        for atom in self.vol_group:
            tmpcoordlist.append(atom.position[2])
        self.zcoords.append(tmpcoordlist)

    def _conclude(self):
        # find the timestep with the most atoms
        counts = np.zeros((len(self.zcoords)))
        for idx,coordlist in enumerate(self.zcoords):
            count = len(coordlist)
            counts[idx] = count
        # use that number to build array
        dimensions = (int(np.amax(counts)),len(self.zcoords))
        self.coordarray = np.empty(dimensions, dtype=float)
        # fill it with nans
        self.coordarray.fill(np.nan)
        # now fill it with coordinates, leaving nan in the gaps
        for idx1,coordlist in enumerate(self.zcoords):
            for idx0,coord in enumerate(coordlist):
                self.coordarray[idx0,idx1] = coord
        # for consistency make this a list, though it will only ever have one item
        self.selections = []
        self.selections.append(('zsearchlist', 'z-coord', '(%s) and (%s)' % (self.vol_selecttext, self.search_selecttext)))
