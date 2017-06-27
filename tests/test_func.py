import os
from analysis.SimpleFeatures import SimpleFeatures as sf
from analysis.VolumeTracker import VolumeTracker as vt
from analysis.ZSearch import ZSearch as zs
from core.TimeseriesCore import TimeseriesCore as tc
from core.MergeDS import MergeDS as mds

input_data = '/Users/anatale/school/UCSF/Grabe_lab/data'
psffile = os.path.join(input_data, 'traakG124I_full_trimmed/topo.traakG124I_full_npt.psf')
dcdfile = os.path.join(input_data, 'traakG124I_full_trimmed/traakG124I_full_npt.all.500ps.filter_alignZ.dcd')
sample_dats = '/Users/anatale/school/UCSF/Grabe_lab/code/ExtendedTSC/sample_dats'

def test_sf():
    selectlist = [
    ('S1top',      'COG',      'protein and (resid 132 or resid 241) and name O'),
    ('S4bottom',   'COG',      'protein and (resid 129 or resid 238) and name OG1')
    ]

    a = sf()
    a.load_dcd(psffile,dcdfile,traj_stepsize=500,framerange=(0,-1,100))
    a.run(selectlist)
    a._data_writer(a.primaryDS)

def test_vt():
    vol_selecttext = 'sphzone 8 (protein and (resid 132 or resid 241) and name O)'
    search_selecttext = 'name POT'

    a = vt()
    a.load_dcd(psffile,dcdfile,traj_stepsize=500,framerange=(0,-1,100))
    a.run(vol_selecttext,search_selecttext)
    a._data_writer(a.primaryDS)
    a._data_writer(a.maskDS)

def test_zs():
    vol_selecttext = 'sphzone 8 (protein and (resid 132 or resid 241) and name O)'
    search_selecttext = 'name OH2 and resname TIP3'

    a = zs()
    a.load_dcd(psffile,dcdfile,traj_stepsize=500,framerange=(0,-1,100))
    a.run(vol_selecttext,search_selecttext)
    a._data_writer(a.primaryDS)

def test_core_read():
    datfile = os.path.join(sample_dats, 'test_vtrack_set_frames0to50.dat')
    maskfile = os.path.join(sample_dats, 'test_vtrack_set_frames0to50.mask.dat')

    a = tc(datfilename=datfile, maskfilename=maskfile)
    a._data_writer(a.primaryDS)
    a._data_writer(a.maskDS)

def test_simple_merge():
    datfiles = ['test_simple_set_frames0to50.dat','test_simple_set_frames50to100.dat','test_simple_set_frames100to150.dat']
    datfiles = [os.path.join(sample_dats,i) for i in datfiles]
    c = mds()
    c.merge_along_time(datfiles)
    c._data_writer(c.primaryDS, outfilename='test1.dat')

def test_padded_merge():
    datfiles = ['test_zsearch_set_frames0to50.dat','test_zsearch_set_frames50to100.dat','test_zsearch_set_frames100to150.dat','test_zsearch_set_frames150to200.dat']
    datfiles = [os.path.join(sample_dats,i) for i in datfiles]
    c = mds()
    c.merge_along_time(datfiles)
    c._data_writer(c.primaryDS, outfilename='test2.dat')

def test_mask_merge():
    pass

def test_multiproc():
    pass
