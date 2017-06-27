io_options = \
{
'input_prefix' : '/Users/anatale/school/UCSF/Grabe_lab/data/traakTM4_trimmed',
'copy_to': None,
'output_prefix' : None,
'num_processes' : 2
}

universe_recipe = \
{
'toponame' : 'topo.traakTM4_npt.psf',
'trajname' : 'traakTM4_npt.all.500ps.filter_alignZ.dcd',
'stepsize' : 500,
'input_type' : 'dcd_traj',
'extraopts' : None
}

feature_sets = \
{
'basic_set' :
    {
    'feature_set_type' : 'basic',
    'descriptorlist' : ''
    }
,
'vsearch_set' :
    {
    'feature_set_type' : 'volumesearch',
    'volselectext' : '',
    'searchselectext' : ''
    }
,
'zsearch_set' :
    {
    'feature_set_type' : 'zsearch',
    'volselectext' : '',
    'searchselectext' : ''
    }
}

trajectory_options = \
{
'chunksize' : 50,
'stride' : 1
}
