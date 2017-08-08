from __future__ import print_function
import os
from tests.data import PDB,DCD

try:
    os.mkdir('out_test')
except OSError:
    pass
os.chdir('out_test')

def test_simple():
    from analysis.SimpleFeatures import SimpleFeatures
    selections = [
    ('sel1', 'COG',         'resid 10 11 12'),
    ('sel2', 'COM',         'resid 74 75 76'),
    ('sel3', 'distance',    'resid 10 76 and name CA'),
    ('sel4', 'dihedral',    '(resid 45 and name N CA C) or (resid 46 and name N)'),
    ('sel5', 'atom',        'resid 15 and name O')
    ]

    for i,chunk in enumerate([(0,33),(33,66),(66,-1)]):
        a = SimpleFeatures()
        a.load_traj(PDB,DCD,traj_stepsize=100,framerange=(chunk[0],chunk[1],1))
        a.run(selections)
        a.write_data('test_simple_%d' % i)

def test_vtrack():
    from analysis.VolumeTracker import VolumeTracker
    vol_selecttext = 'sphzone 10.0 (resid 25 and name CA)'
    search_selecttext = 'name CB'

    for i,chunk in enumerate([(0,33),(33,66),(66,-1)]):
        a = VolumeTracker()
        a.load_traj(PDB,DCD,traj_stepsize=100,framerange=(chunk[0],chunk[1],1))
        a.run(vol_selecttext,search_selecttext)
        a.write_data('test_vtrack_%d' % i)

def test_zsearch():
    from analysis.ZSearch import ZSearch
    vol_selecttext = 'sphzone 10.0 (resid 25 and name CA)'
    search_selecttext = 'name CB'

    for i,chunk in enumerate([(0,33),(33,66),(66,-1)]):
        a = ZSearch()
        a.load_traj(PDB,DCD,traj_stepsize=100,framerange=(chunk[0],chunk[1],1))
        a.run(vol_selecttext,search_selecttext)
        a.write_data('test_zsearch_%d' % i)

def test_rmsd():
    from analysis.RMSDseries import RMSDseries

    # frame
    a = RMSDseries()
    a.load_traj(PDB,DCD,traj_stepsize=100)
    a.run('name CA', frame=50)
    a.write_data('test_rmsd_frameref')

    # pdb
    b = RMSDseries()
    b.load_traj(PDB,DCD,traj_stepsize=100)
    b.run('name N', pdbname=PDB, align=True)
    b.write_data('test_rmsd_pdbref')

def time_merge():
    from core.MergeDS import merge_along_time

    datlist = []
    for i in [0,1,2]:
        datlist.append('test_simple_%d.dat' % i)

    new = merge_along_time(datlist)
    new._write(outfilename='merge.dat')
    # if sorting works, this should be the same as above
    rev = merge_along_time(datlist.__reversed__())
    rev._write(outfilename='merge_rev.dat')

    datlist = []
    for i in [0,1,2]:
        datlist.append('test_zsearch_%d.dat' % i)

    pad = merge_along_time(datlist)
    pad._write(outfilename='merge_pad.dat')

    datlist = []
    for i in [0,1,2]:
        datlist.append('test_vtrack_%d.dat' % i)

    dynamic = merge_along_time(datlist)
    dynamic._write(outfilename='merge_dyn.dat')

    datlist = []
    for i in [0,1,2]:
        datlist.append('test_vtrack_%d.mask.dat' % i)

    dynamic = merge_along_time(datlist)
    dynamic._write(outfilename='merge_dyn.mask.dat')

def feature_merge():
    from core.MergeDS import merge_along_features

    merge = merge_along_features(['test_simple_0.dat', 'test_zsearch_0.dat'])
    merge._write(outfilename='feat_merge.dat')

if __name__ == '__main__':
    test_simple()
    test_vtrack()
    test_zsearch()
    test_rmsd()
    time_merge()
    feature_merge()
