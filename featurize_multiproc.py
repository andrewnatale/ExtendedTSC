import sys, os, math, ast
import numpy as np
import multiprocessing
import MDAnalysis as md
from ExtendedTSC import ExtendedTSC
from ZSearch import ZSearch

# load config file
configfile = 'example_multiproc_config.py'
execfile(configfile)

feature_set_types = [
'basic', # use ExtendedTSC
'volumesearch', # use ExtendedTSC
'zsearch' # use Zsearch
]

# init universe
if io_options['ram_fs']:
    #copy to ram
    pass
else:
    input_prefix = io_options['input_prefix']

u = md.Universe(
  os.path.join(input_prefix, universe_recipe['toponame']),
  os.path.join(input_prefix, universe_recipe['trajname']))

startframe = None
stopframe = None
frameranges = []
while stopframe != u.trajectory.n_frames - 1:
    if startframe is None:
        startframe = 0
        stopframe = trajectory_options['chunksize'] - 1
    else:
        startframe = startframe + trajectory_options['chunksize']
        stopframe = stopframe + trajectory_options['chunksize']
    if stopframe >= u.trajectory.n_frames:
        stopframe = u.trajectory.n_frames - 1
    frameranges.append((startframe, stopframe, 1))

job_array = []
for key in feature_sets:
    for framerange in frameranges:
        job_array.append([key, feature_sets[key], framerange])

def job_runner(feature_set_name, feature_set_options, framerange):
    # get feature set type
    if feature_set_options['feature_set_type'] == 'basic':
        a = ExtendedTSC(verbose=False,log=True)
        a.load_universe(
          u,
          universe_recipe['stepsize'],
          framerange=framerange,
          input_type=universe_recipe['input_type'],
          toponame=os.path.join(input_prefix, universe_recipe['toponame']),
          trajname=os.path.join(input_prefix, universe_recipe['trajname']))
        a.measures_from_list(feature_set_options['descriptorlist'])
        a.run()
        a.write_data('%s_%s_%d_%d' % (feature_set_name, feature_set_options['feature_set_type'], framerange[0], framerange[1]))
    elif feature_set_options['feature_set_type'] == 'volumesearch':
        a = ExtendedTSC(verbose=False,log=True)
        a.load_universe(
          u,
          universe_recipe['stepsize'],
          framerange=framerange,
          input_type=universe_recipe['input_type'],
          toponame=os.path.join(input_prefix, universe_recipe['toponame']),
          trajname=os.path.join(input_prefix, universe_recipe['trajname']))
        a.measures_from_volumesearch(feature_set_options['volselectext'], feature_set_options['searchselectext'])
        a.run()
        a.write_data('%s_%s_%d_%d' % (feature_set_name, feature_set_options['feature_set_type'], framerange[0], framerange[1]))
    elif feature_set_options['feature_set_type'] == 'zsearch':
        a = ZSearch(verbose=False,log=True)
        a.load_universe(
          u,
          universe_recipe['stepsize'],
          framerange=framerange,
          input_type=universe_recipe['input_type'],
          toponame=os.path.join(input_prefix, universe_recipe['toponame']),
          trajname=os.path.join(input_prefix, universe_recipe['trajname']))
        a.run(feature_set_options['volselectext'], feature_set_options['searchselectext'])
        a.write_data('%s_%s_%d_%d' % (feature_set_name, feature_set_options['feature_set_type'], framerange[0], framerange[1]))

mppool = multiprocessing.Pool(int(io_options['num_processes']))
mppool.map(job_runner, job_array)
