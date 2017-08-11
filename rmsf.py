from __future__ import print_function
import sys, os
from analysis.RMSFwindow import RMSFwindow

topo = os.path.abspath(sys.argv[1])
traj = os.path.abspath(sys.argv[2])
output_dir = os.path.abspath(sys.argv[3])
tag = sys.argv[4]

basename = ''.join(traj.split('/')[-1].split('.')[:-1])
print(basename)

a = RMSFwindow()
a.load_traj(topo,traj,50)
# 41 frames == ~2ns time window
# using an odd number means the windows are centered on existing time steps
a.run('protein and (resid 129-133 or resid 238-242) and name CA', 41)

os.chdir(output_dir)
a.write_data('%s.%s' % (basename, tag))
