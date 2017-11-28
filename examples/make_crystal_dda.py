'''
Test code to call DDA input file generating code.

Robert Schrom @ 11/2017
'''
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from crystal_dda import branched_planar_dda

# set values to create branched planar crystal with
a = 0.2
amax = 3.
ac = 0.5

fb = 0.4
ft = 0.2
fg = 0.7
fmb = 0.4

nsb = 5
asp = 20.
ag = amax*fg

fname, afrac = branched_planar_dda(a, asp, amax, ac, ag, ft, fb, fmb, nsb)
print fname, afrac

# plot dipole locations
dda_data = np.genfromtxt(fname, skip_header=3)
x3d = dda_data[:,0]
y3d = dda_data[:,1]
z3d = dda_data[:,2]

# get 2d slice of x y values
zval = z3d[0]
x2d = x3d[z3d==zval]
y2d = y3d[z3d==zval]
plt.scatter(x2d, y2d, c='b', s=1, edgecolor='')
ax = plt.gca()
ax.set_aspect(1.)
plt.savefig('dipole_locations.png')
