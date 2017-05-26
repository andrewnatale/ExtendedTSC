import sys, os, math, datetime
import numpy as np
# tested with MDAnalysis-0.15.0
# MDAnalysis-0.16.0 WILL NOT WORK WITH THIS SCRIPT
import MDAnalysis as md
from MDAnalysis.analysis.base import AnalysisBase
import ExtendedTSC

psffile = '/Users/anatale/school/UCSF/Grabe_Lab/k2p/test_measure/trkWT.popc200.ion.psf'
dcdfile = '/Users/anatale/school/UCSF/Grabe_Lab/k2p/test_measure/traakWT_full_npt.wrap_alignSF.all.500ps.dcd'
frameskip = 1

# volume selections for selectivity filter

# select via a cylinder - faster but the trajectory must be aligned such that the z-axisw is parallel to the filter
selecttext = '(cyzone 3.5 15 -2 protein and (resid 129 or resid 238) and name OG1) and (name POT or (name OH2 and resname TIP3))'
# selecttext = '(cyzone 3.5 15 -2 protein and (resid 129 or resid 238) and name OG1) and (name POT)'

# select via 5 spheres centered on the S0-S4 positions, much slower than using a cylinder but requires no trajectory alignment
# selecttext = '(\
# (sphzone 3 protein and (resid 129 or resid 238) and (name OG1 or name O)) or \
# (sphzone 3 protein and (resid 129 or resid 238 or resid 130 or resid 239) and name O) or \
# (sphzone 3 protein and (resid 130 or resid 239 or resid 131 or resid 240) and name O) or \
# (sphzone 3 protein and (resid 131 or resid 240 or resid 132 or resid 241) and name O) or \
# (sphzone 3 protein and (resid 132 or resid 241 or resid 133 or resid 242) and name O)\
# ) and \
# (name POT or (name OH2 and resname TIP3)) \
# '

# build a list of the ions and water molecules that are present in the fiter at any timestep
u = md.Universe(psffile,dcdfile,topology_format=u'PSF')

# use MDAnalysis' AnalysisBase class to build a subclass for locating the atoms we want
class SelFilterMeasure(AnalysisBase):
    def __init__(self):
        # checks
        self.check = True

    def _prepare(self):
        # setup vars and data structures
        self.filter_set = set()
        self.masking_data = []

    def _single_frame(self):
        # what to do at each frame
        print self._ts
        tmpoccupancy = []
        tmpatoms = u.select_atoms(selecttext)
        for atom in tmpatoms:
            # store unique identifier
            self.filter_set.add((atom.index,atom.name,atom.resid,atom.segid))
            tmpoccupancy.append(atom.index)
        self.masking_data.append(tmpoccupancy)

sel_filter = SelFilterMeasure()
# print test.check
# print isinstance(test, AnalysisBase)
# print isinstance(test, SelFilterMeasure)
sel_filter._setup_frames(u.trajectory, step=frameskip)
print sel_filter.start, sel_filter.stop, sel_filter.step, sel_filter.nframes
sel_filter.run()
print 'identified %d atoms:' % len(sel_filter.filter_set)
for elem in sorted(sel_filter.filter_set):
    print elem

print sel_filter.masking_data
# build selection lists for use with ExtendedTSC
selections = []
selections_masks = []
water_idx = 0
pot_idx = 0
for elem in sorted(sel_filter.filter_set):
    if elem[1] == 'OH2':
        water_idx += 1
        tmp_name = 'filter_water%d' % water_idx
        tmp_type = 'COM'
        tmp_selecttext = 'segid %s and resid %d' % (elem[3],elem[2])
    elif elem[1] == 'POT':
        pot_idx += 1
        tmp_name = 'filter_pot%d' % pot_idx
        tmp_type = 'atom'
        tmp_selecttext = 'segid %s and resid %d' % (elem[3],elem[2])
    selections.append((tmp_name,tmp_type,tmp_selecttext))
    selections_masks.append((tmp_name,'mask',tmp_selecttext))
# add some filter reference points in case they come in handy later
selections.append(('S4bottom', 'COG', 'protein and (resid 129 or resid 238) and name OG1'))
selections.append(('S1top', 'COG', 'protein and (resid 132 or resid 241) and name O'))
print selections

# use ExtendedTSC to make measurements (sets up its own internal 'universe' object)
a = ExtendedTSC.ExtendedTSC()
a.measures_from_dcd(selections,psffile,dcdfile,frameskip=frameskip)

# format masking data
masks = np.zeros((water_idx+pot_idx,np.shape(a.timesteps)[0]), dtype=int)
for idxA, elem in enumerate(sorted(sel_filter.filter_set)):
    atom_idx = elem[0]
    for idxB, occupancy in enumerate(sel_filter.masking_data):
        if atom_idx in occupancy:
            try:
                masks[idxA,idxB] = 1
            except IndexError:
                pass

mask_widths = np.ones((water_idx+pot_idx),dtype=int).tolist()
a.derivative_collection(selections_masks,mask_widths,masks)
a.write_dat('coordinates_test.dat')
a.write_dat('masks_test.dat',derivative=True)
