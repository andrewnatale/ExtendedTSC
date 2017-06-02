import sys, os, math
import numpy as np
import ExtendedTSC

# topology and trajectory files
try:
    psffile = sys.argv[1]
    dcdfile = sys.argv[2]
except IndexError:
    psffile = None
    dcdfile = None

#input_prefix = '/Volumes/data/andrew'
input_prefix = '/mnt/hd_scratch/traak_data'

outtag = 'SelFilter_vol'
# make output dir in cwd
outdir = os.getcwd()
#os.mkdir(outdir)

# if I'm not batch processing, I'm probably just testing on a trajectory with the non TM4 topology
if psffile and dcdfile:
    print 'Using inputs specified from command line'
    filepairs = [(os.path.abspath(psffile), os.path.abspath(dcdfile))]
    filepairsTM4 = None
else:
    # hardcode filenames for batch process
    filepairs = [
    (os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_full_trimmed/trkG124full.ion.psf'),
     os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_full_trimmed/traakG124I_full_npt.wrap_alignSF.all.500ps.dcd')),
    (os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_S1S3_trimmed/trkG124_S1S3.ion.psf'),
     os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_S1S3_trimmed/traakG124I_S1S3_npt.wrap_alignSF.all.500ps.dcd')),
    (os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_S2S4_trimmed/trkG124S2S4.ion.psf'),
     os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_S2S4_trimmed/traakG124I_S2S4_npt.wrap_alignSF.all.500ps.dcd')),
    (os.path.join(input_prefix, 'mackinnon_traak/traakWT_full_trimmed/trkWT.popc200.ion.psf'),
     os.path.join(input_prefix, 'mackinnon_traak/traakWT_full_trimmed/traakWT_full_npt.wrap_alignSF.all.500ps.dcd')),
    (os.path.join(input_prefix, 'mackinnon_traak/traakWT_S1S3_trimmed/trkWTS1S3.ion.psf'),
     os.path.join(input_prefix, 'mackinnon_traak/traakWT_S1S3_trimmed/traakWT_S1S3_npt.wrap_alignSF.all.500ps.dcd')),
    (os.path.join(input_prefix, 'mackinnon_traak/traakWT_S2S4_trimmed/trkWTS2S4.ion.psf'),
     os.path.join(input_prefix, 'mackinnon_traak/traakWT_S2S4_trimmed/traakWT_S2S4_npt.wrap_alignSF.all.500ps.dcd'))
    ]
    filepairsTM4 = [
    (os.path.join(input_prefix, 'traakTM4_trimmed/traakTM4.popc.ion.psf'),
     os.path.join(input_prefix, 'traakTM4_trimmed/traakTM4_npt.wrap_alignSF.all.500ps.dcd')),
    ]

# selections for measurements
ref_selections = [
('S1top',      'COG',      'protein and (resid 132 or resid 241) and name O'),
('S4bottom',   'COG',      'protein and (resid 129 or resid 238) and name OG1')
]

#vol_shape = 'cyzone 3.5 15 -2'
vol_shape = 'sphzone 3'
vol_center =  'protein and (resid 129 or resid 238) and name OG1'
#search_selectext = 'name POT or (name OH2 and resname TIP3)'
search_selectext = 'name POT'
# to make selections in the traakTM4 trajectory, numbers must be adjusted
# chains A and B do not have different segids
# for res in chain A, subtract 27 from actual numbers
# for res in chain B, add 255
ref_selections_traakTM4 = [
('S1top',      'COG',      'protein and (resid 105 or resid 387 or resid 214 or resid 496) and name O'),
('S4bottom',   'COG',      'protein and (resid 102 or resid 384 or resid 211 or resid 493) and name OG1')
]

vol_center_traakTM4 = 'protein and (resid 102 or resid 384 or resid 211 or resid 493) and name OG1'
search_selectext_traakTM4 = 'name POT or (name OH2 and resname TIP3)'

for pair in filepairs:
    a = ExtendedTSC.ExtendedTSC()
    a.load_traj(pair[0],pair[1])
    a.measures_from_list(ref_selections)
    a.measures_from_volumesearch(vol_shape,vol_center,search_selectext)
    a.generate_timeseries()
    output = pair[1].split('/')[-1].split('.')[0] + '.' + outtag
    a.write_data(output)

if filepairsTM4:
    for pair in filepairsTM4:
        a = ExtendedTSC.ExtendedTSC()
        a.load_traj(pair[0],pair[1])
        a.measures_from_list(ref_selections_traakTM4)
        a.measures_from_volumesearch(vol_shape,vol_center_traakTM4,search_selectext_traakTM4)
        a.generate_timeseries()
        output = pair[1].split('/')[-1].split('.')[0] + '.' + outtag
        a.write_data(output)
