import os
from analysis.SimpleFeatures import SimpleFeatures as sf
from analysis.VolumeTracker import VolumeTracker as vt
from analysis.ZSearch import ZSearch as zs
from core.TimeseriesCore import TimeseriesCore as tc

input_data = '/Users/anatale/school/UCSF/Grabe_lab/data'
psffile = os.path.join(input_data, 'traakG124I_full_trimmed/topo.traakG124I_full_npt.psf')
dcdfile = os.path.join(input_data, 'traakG124I_full_trimmed/traakG124I_full_npt.all.500ps.filter_alignZ.dcd')

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
    datfile = os.path.abspath('tests/inputs/traakTM4_npt.core_vol.dat')
    maskfile = os.path.abspath('tests/inputs/traakTM4_npt.core_vol.mask.dat')

    a = tc(datfilename=datfile, maskfilename=maskfile)
    a._data_writer(a.primaryDS)
    a._data_writer(a.maskDS)

def test_merge():
    datfiles = ['traakTM4_npt.test01.dat','traakTM4_npt.test02.dat','traakTM4_npt.test03.dat']
    datfiles = [os.path.join(input_prefix,i) for i in datfiles]
    c = TimeseriesCore()
    c.merge_along_time(datfiles)
    print c.primaryDS.__dict__
    c._data_writer(c.primaryDS, outfile='test_output/merge_test.dat')

def test_multiproc():
    pass
