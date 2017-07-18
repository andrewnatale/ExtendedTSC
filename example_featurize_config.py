options = \
{
'job_name' : 'test',
'input_prefix' : '/path/to/indir',
'copy_to': '/path/to/ramdisk',
'output_prefix' : '/path/to/outdir',
'num_proc' : 5,
'chunksize' : 50,
'stride' : 1
}

universe_recipe = \
{
'toponame' : 'topo.file.name.psf',
'trajname' : 'traj.file.name.dcd',
'stepsize' : 500,
'input_type' : 'dcd_traj'
}

feature_sets = \
{
'simple_set' :
    {
    'feature_set_type' : 'simple',
    'descriptorlist' : ''
    }
,
'vtrack_set' :
    {
    'feature_set_type' : 'vtrack',
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
,
'rmsd_set' :
    {
    'feature_set_type' : 'rmsd',
    'selecttext' : '',
    'ref_frame' : 0,
    'pdbname' : None,
    'alt_selecttext' : None,
    'align' : False
    }
}
