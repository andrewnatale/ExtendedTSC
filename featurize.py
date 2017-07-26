from __future__ import print_function
import sys, os, shutil
import numpy as np
import multiprocessing
import MDAnalysis as mda
from analysis.SimpleFeatures import SimpleFeatures
from analysis.VolumeTracker import VolumeTracker
from analysis.ZSearch import ZSearch
from analysis.RMSDseries import RMSDseries
from core.MergeDS import merge_along_time

# load config file
# provides 3 dicts; 'options', 'universe_recipe', and 'feature_sets'
configfile = sys.argv[1]
execfile(configfile)
#print options,universe_recipe,feature_sets

feature_set_types = [
'simple', # use SimpleFeatures
'vtrack', # use VolumeTracker
'zsearch', # use Zsearch
'rmsd' # use RMSDseries
]

# prep output dir and go there
try:
    os.makedirs(options['output_prefix'])
except OSError:
    if not os.path.isdir(options['output_prefix']):
        raise
os.chdir(options['output_prefix'])

# optionally copy input files to another location (like a ramdisk) before loading
if options['copy_to']:
    try:
        os.makedirs(options['copy_to'])
    except OSError:
        if not os.path.isdir(options['copy_to']):
            raise
    shutil.copy(os.path.join(options['input_prefix'], universe_recipe['toponame']), options['copy_to'])
    shutil.copy(os.path.join(options['input_prefix'], universe_recipe['trajname']), options['copy_to'])
    input_prefix = options['copy_to']
else:
    input_prefix = options['input_prefix']

# init universe, but only to get n_frames
throwaway = mda.Universe(
  os.path.join(input_prefix, universe_recipe['toponame']),
  os.path.join(input_prefix, universe_recipe['trajname']))

n_frames = throwaway.trajectory.n_frames
del throwaway

# figure out how many chunks to break each big job into
# NOTE: the stopframe variable is non-inclusive!
# to get frames 0-49 you must pass a framerange of (0, 50, 1)
startframe = None
stopframe = None
chunks = []
while stopframe < n_frames:
    if startframe is None:
        startframe = 0
        stopframe = options['chunksize']
    else:
        startframe = startframe + options['chunksize']
        stopframe = stopframe + options['chunksize']
    if stopframe >= n_frames-1:
        chunks.append((startframe, -1, 1))
    else:
        chunks.append((startframe, stopframe, 1))

# build job array
job_array = []
for key in feature_sets:
    for idx,framerange in enumerate(chunks):
        task_id = str(idx + 1)
        while len(task_id) < 4:
            task_id = '0'+task_id
        job_array.append((key, feature_sets[key], framerange, task_id))
#print job_array

def job_runner(opts):
    feature_set_name, feature_set_options, framerange, task_id = opts
    fst = feature_set_options['feature_set_type']
    # output tag
    outname = '%s_%s_%s' % (options['job_name'], feature_set_name, task_id)
    # get feature set type and init
    if fst == 'simple':
        a = SimpleFeatures(verbose=True,log=outname+'.log')
    elif fst == 'vtrack':
        a = VolumeTracker(verbose=True,log=outname+'.log')
    elif fst == 'zsearch':
        a = ZSearch(verbose=True,log=outname+'.log')
    elif fst == 'rmsd':
        a = RMSDseries(verbose=True,log=outname+'.log')
    # re-load universe
    a.load_dcd(
      os.path.join(input_prefix, universe_recipe['toponame']),
      os.path.join(input_prefix, universe_recipe['trajname']),
      universe_recipe['stepsize'],
      framerange=framerange
      )
    # choose .run() options and calculate
    if fst == 'simple':
        a.run(feature_set_options['descriptorlist'])
    elif fst == 'vtrack':
        a.run(feature_set_options['volselectext'], feature_set_options['searchselectext'])
    elif fst == 'zsearch':
        a.run(feature_set_options['volselectext'], feature_set_options['searchselectext'])
    elif fst == 'rmsd':
        a.run(
          feature_set_options['selecttext'],
          frame=feature_set_options['ref_frame'],
          pdbname=feature_set_options['pdbname'],
          alt_selecttext=feature_set_options['alt_selecttext'],
          align=feature_set_options['align']
          )
    a.write_data(outname)

try:
    mppool = multiprocessing.Pool(int(options['num_proc']))
    mppool.map(job_runner, job_array)
except:
    print('Something went wrong in the multiprocessing pool!!')
    print('Output files may be missing or contain errors!!')
else:
    # merge data sets
    for key in feature_sets:
        mergelist = []
        maskmergelist = []
        for filename in sorted(os.listdir(os.getcwd())):
            # since this script wrote the files, they should have very predictable filename lengths
            if filename.startswith('%s_%s' % (options['job_name'], key)) \
              and filename.endswith('.dat') \
              and len(filename) == (len(options['job_name']+key)+10):
                mergelist.append(filename)
            if filename.startswith('%s_%s' % (options['job_name'], key)) \
              and filename.endswith('mask.dat') \
              and len(filename) == (len(options['job_name']+key)+15):
                maskmergelist.append(filename)
        b = merge_along_time(sorted(mergelist))
        b.write_dat(outfilename='%s_%s_all_frames' % (options['job_name'], key))
        if len(maskmergelist) > 0:
            c = merge_along_time(sorted(maskmergelist))
            c.write_dat(outfilename='%s_%s_all_frames.mask' % (options['job_name'], key))
finally:
    if options['copy_to']:
        if options['copy_to'] != options['input_prefix']:
            print('Cleaning up...')
            os.remove(os.path.join(input_prefix, universe_recipe['toponame']))
            os.remove(os.path.join(input_prefix, universe_recipe['trajname']))
