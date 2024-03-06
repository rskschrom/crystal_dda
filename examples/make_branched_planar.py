from crystal_dda.polygons import BranchedPlanar
from crystal_dda.shapes import Crystal
import matplotlib.pyplot as plt
import numpy as np

# create branched planar crystal polygon
a_axis_max = 3.
a_axis = 2.5
a_axis_core = 0.3
frac_tip = 0.2
frac_bcov = 0.3
frac_gap = 0.6
nsub = 5
bp = BranchedPlanar(a_axis, a_axis_max, a_axis_core,
                    frac_tip, frac_bcov, frac_gap, nsub)

# create 3d crystal
c_axis = 0.2
cr = Crystal(bp, c_axis)

# fill with dipoles
dip_len = 0.04
cr.create_dipoles(dip_len)
cr.write_dipoles('crystal.txt')

# convert dipoles to physical dimensions
crx = dip_len*cr.ix
cry = dip_len*cr.iy
crx = crx-np.mean(crx)
cry = cry-np.mean(cry)

# plot dipoles in 2d and associated crystal polygon
plt.scatter(crx, cry, c='k', s=4.)
plt.plot(bp.x, bp.y, 'm-')
ax = plt.gca()
ax.set_aspect(1.)
plt.savefig('crystal.png')
