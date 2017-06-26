import sys, os, math, re
import numpy as np
from ExtendedTSC import ExtendedTSC
from TimeseriesCore import TimeseriesCore

def test_etsc():
    input_prefix = '/Users/anatale/school/UCSF/Grabe_lab/data'

    psffile = os.path.join(input_prefix, 'traakG124I_full_trimmed/topo.traakG124I_full_npt.psf')
    dcdfile = os.path.join(input_prefix, 'traakG124I_full_trimmed/traakG124I_full_npt.all.200ps.filter_alignZ.dcd')

    selections = [
    ('S1top',      'COG',      'protein and (resid 132 or resid 241) and name O'),
    ('S4bottom',   'COG',      'protein and (resid 129 or resid 238) and name OG1')
    ]

    a = ExtendedTSC()
    a.load_dcd(psffile,dcdfile,traj_stepsize=200,framerange=(0,-1,100))
    a.measures_from_list(selections)
    a.run()
    a._data_writer(a.primaryDS)
    print a.maskDS.__dict__

def test_core():
    datfile = os.path.abspath('core_volume_tracking/traakTM4_npt.core_vol.dat')
    maskfile = os.path.abspath('core_volume_tracking/traakTM4_npt.core_vol.mask.dat')

    b = TimeseriesCore(datfile=datfile, maskfile=maskfile)
    #b._data_writer(a.primaryDS)
    b._data_writer(b.maskDS)

def test_merge():
    datfiles = ['testfiles/traakTM4_npt.test01.dat','testfiles/traakTM4_npt.test02.dat','testfiles/traakTM4_npt.test03.dat']

    c = TimeseriesCore()
    c.merge_along_time(datfiles)
    print c.primaryDS.__dict__
    c._data_writer(c.primaryDS)

#test_etsc()
#test_core()
test_merge()
