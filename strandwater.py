options = \
{
'job_name' : 'traak_up',
'input_prefix' : '/Users/anatale/UCSF/Grabe_lab/data',
'copy_to': None,
'output_prefix' : '/Users/anatale/UCSF/Grabe_lab/scratch/filter_water_traak_up',
'num_proc' : 3,
'chunksize' : 100,
'stride' : 1
}

universe_recipe = \
{
'toponame' : 'topo.traak_up_nohydrogen.psf',
'trajname' : 'traak_up_nohydrogen.dcd',
'stepsize' : 50
}

feature_sets = \
{
'strand1A' :
    {
    'feature_set_type' : 'vtrack',
    'volselectext' : '(sphzone 9.0 (protein and resid 103 and name CA)) and (sphzone 9.0 (protein and resid 106 and name CB))',
    'searchselectext' : 'name OW'
    },
'strand1B' :
    {
    'feature_set_type' : 'vtrack',
    'volselectext' : '(sphzone 9.0 (protein and resid 392 and name CA)) and (sphzone 9.0 (protein and resid 395 and name CB))',
    'searchselectext' : 'name OW'
    },
'strand2A' :
    {
    'feature_set_type' : 'vtrack',
    'volselectext' : '(sphzone 9.0 (protein and resid 212 and name CA)) and (sphzone 9.0 (protein and resid 215 and name CB))',
    'searchselectext' : 'name OW'
    },
'strand2B' :
    {
    'feature_set_type' : 'vtrack',
    'volselectext' : '(sphzone 9.0 (protein and resid 501 and name CA)) and (sphzone 9.0 (protein and resid 504 and name CB))',
    'searchselectext' : 'name OW'
    },
}
