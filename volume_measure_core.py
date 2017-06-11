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

input_prefix = '/Volumes/data/andrew'

outtag = 'core_vol'
# make output dir in cwd
outdir = os.path.join(os.getcwd(), 'core_volume_tracking')
try:
    os.mkdir(outdir)
except OSError:
    if not os.path.isdir(outdir):
        raise

# if I'm not batch processing, I'm probably just testing on a trajectory with the non TM4 topology
if psffile and dcdfile:
    print 'Using inputs specified from command line'
    filepairs = [(os.path.abspath(psffile), os.path.abspath(dcdfile))]
    filepairsTM4 = None
else:
    # hardcode filenames for batch process
    filepairs = [
    (os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_full_trimmed/topo.traakG124I_full_npt.psf'),
     os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_full_trimmed/traakG124I_full.new_align.500ps.dcd')),
    (os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_S1S3_trimmed/topo.traakG124I_S1S3_npt.psf'),
     os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_S1S3_trimmed/traakG124I_S1S3.new_align.500ps.dcd')),
    (os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_S2S4_trimmed/topo.traakG124I_S2S4_npt.psf'),
     os.path.join(input_prefix, 'lolicato_mutant_traak/traakG124I_S2S4_trimmed/traakG124I_S2S4.new_align.500ps.dcd')),
    (os.path.join(input_prefix, 'mackinnon_traak/traakWT_full_trimmed/topo.traakWT_full_npt.psf'),
     os.path.join(input_prefix, 'mackinnon_traak/traakWT_full_trimmed/traakWT_full.new_align.500ps.dcd')),
    (os.path.join(input_prefix, 'mackinnon_traak/traakWT_S1S3_trimmed/topo.traakWT_S1S3_npt.psf'),
     os.path.join(input_prefix, 'mackinnon_traak/traakWT_S1S3_trimmed/traakWT_S1S3.new_align.500ps.dcd')),
    (os.path.join(input_prefix, 'mackinnon_traak/traakWT_S2S4_trimmed/topo.traakWT_S2S4_npt.psf'),
     os.path.join(input_prefix, 'mackinnon_traak/traakWT_S2S4_trimmed/traakWT_S2S4.new_align.500ps.dcd'))
    ]
    filepairsTM4 = [
    (os.path.join(input_prefix, 'traakTM4_trimmed/topo.traakTM4_npt.psf'),
     os.path.join(input_prefix, 'traakTM4_trimmed/traakTM4.new_align.500ps.dcd')),
    ]

# selections for measurements
ref_selections = [
('S1top',      'COG',      'protein and (resid 132 or resid 241) and name O'),
('S4bottom',   'COG',      'protein and (resid 129 or resid 238) and name OG1')
]

vol_selecttext = '(cyzone 4 16 0 (protein and (resid 129 or resid 238) and name OG1)) or \
                  (cyzone 8 0 -16 (protein and (resid 129 or resid 238) and name OG1))'

search_selectext = 'name POT or\
      name C22 or name C23 or name C24 or name C25 or name C26 or\
      name C27 or name C28 or name C29 or name C210 or name C211 or\
      name C212 or name C213 or name C214 or name C215 or name C216 or\
      name C217 or name C218 or\
      name C32 or name C33 or name C34 or name C35 or name C36 or\
      name C37 or name C38 or name C39 or name C310 or name C311 or\
      name C312 or name C313 or name C314 or name C315 or name C316'

# to make selections in the traakTM4 trajectory, numbers must be adjusted
# chains A and B do not have different segids
# for res in chain A, subtract 27 from actual numbers
# for res in chain B, add 255
ref_selections_traakTM4 = [
('S1top',      'COG',      'protein and (resid 105 or resid 387 or resid 214 or resid 496) and name O'),
('S4bottom',   'COG',      'protein and (resid 102 or resid 384 or resid 211 or resid 493) and name OG1')
]

vol_selecttext_traakTM4 = '(cyzone 4 16 0 (protein and (resid 102 or resid 384 or resid 211 or resid 493) and name OG1)) or \
                           (cyzone 8 0 -16 (protein and (resid 102 or resid 384 or resid 211 or resid 493) and name OG1))'

os.chdir(outdir)

for pair in filepairs:
    a = ExtendedTSC.ExtendedTSC()
    a.load_traj(pair[0],pair[1])
    a.measures_from_list(ref_selections)
    a.measures_from_volumesearch(vol_selecttext,search_selectext,mode='atom')
    a.generate_timeseries()
    output = pair[1].split('/')[-1].split('.')[0] + '.' + outtag
    a.write_data(output)
    b = ExtendedTSC.ExtendedTSC()
    b.load_traj(pair[0],pair[1])
    b.water_search(vol_selecttext)
    output = pair[1].split('/')[-1].split('.')[0] + '.' + outtag + '.' + 'tip3'
    b.write_data(output)

if filepairsTM4:
    for pair in filepairsTM4:
        a = ExtendedTSC.ExtendedTSC()
        a.load_traj(pair[0],pair[1])
        a.measures_from_list(ref_selections_traakTM4)
        a.measures_from_volumesearch(vol_selecttext_traakTM4,search_selectext,mode='atom')
        a.generate_timeseries()
        output = pair[1].split('/')[-1].split('.')[0] + '.' + outtag
        a.write_data(output)
        b = ExtendedTSC.ExtendedTSC()
        b.load_traj(pair[0],pair[1])
        b.water_search(vol_selecttext_traakTM4)
        output = pair[1].split('/')[-1].split('.')[0] + '.' + outtag + '.' + 'tip3'
        b.write_data(output)
