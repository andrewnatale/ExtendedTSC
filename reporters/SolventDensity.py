from __future__ import print_function
import sys, os
import MDAnalysis as mda
from base.getMDdensity import _getMDdensity

psffile = os.path.abspath(sys.argv[1])
dcdfile = os.path.abspath(sys.argv[2])
outname = sys.argv[3]

u = mda.Universe(psffile,dcdfile)

lipid_hphobic_sel="\
name C22 or name C23 or name C24 or name C25 or name C26 or\
name C27 or name C28 or name C29 or name C210 or name C211 or\
name C212 or name C213 or name C214 or name C215 or name C216 or\
name C217 or name C218 or\
name C32 or name C33 or name C34 or name C35 or name C36 or\
name C37 or name C38 or name C39 or name C310 or name C311 or\
name C312 or name C313 or name C314 or name C315 or name C316"

water_ox_sel = "resname TIP3 and name OH2"

pot_sel = 'resname POT'

# dens_h2o = _getMDdensity(32, 32, water_ox_sel, u, verbose=True)
# dens_h2o.run()
#
# dens_h2o.gridobj.export('h2o_%s' % outname, file_format='dx')
#
# dens_lipid = _getMDdensity(75, 100, lipid_hphobic_sel, u, verbose=True)
# dens_lipid.run()
#
# dens_lipid.gridobj.export('lipid_%s' % outname, file_format='dx')

dens_pot = _getMDdensity(20,40,pot_sel, u, verbose=True)
dens_pot.run()
dens_pot.gridobj.export('pot_%s' % outname, file_format='dx')
