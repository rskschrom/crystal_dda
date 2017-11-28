'''
Code to generate DDA files and compare accuracy of analytical area fraction.

Robert Schrom @ 11/2017
'''
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from crystal_dda.crystal_dda import branched_planar_dda
from crystal_dda.polygons import make_branched_planar, make_hexagon
from crystal_dda.geometry import afrac_dda_subregion

# set values to create branched planar crystal with
amax = 3.
ac = 0.05

fb = 0.6
ft = 0.4
fg = 0.3

nsb = 11
asp = 20.
ag = amax*fg
nxp = 300

# determine main branch fraction width (give same width as sub-branches, or 1)
wt = ft/2.*amax
wsb = 1./((nsb-1)/fb+1.)*(amax-ac-wt)
wmb = min(max(wsb/2., ac/2.), min(wsb/2., ac/2.))
fmb = wmb/(ac/2.)

# loop over a axis lengths for given crystal
numa = 21
avals = np.linspace(ac, amax, numa)
afrac_anl_vals = np.empty(numa)
afrac_dda_vals = np.empty(numa)
xdda, ydda = make_branched_planar(amax, ac, ag, ft, fb, fmb, nsb, 0.)

for i, a in enumerate(avals):
    fname, afrac = branched_planar_dda(a, asp, amax, ac, ag, ft, fb, nsb, nxp)

    # compare area fraction of actual DDA particle
    xhex, yhex = make_hexagon(a)
    afrac_dda = afrac_dda_subregion(xhex, yhex, xdda, ydda)
    print 'a:{:.2e}\tanalytical: {:.3f}\tdda: {:.3f}'.format(a, afrac, afrac_dda)
    afrac_anl_vals[i] = afrac
    afrac_dda_vals[i] = afrac_dda

# plot
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
mpl.rc('text', usetex=True)

plt.plot(avals, afrac_anl_vals, 'r--', lw=3., label='analytical')
plt.plot(avals, afrac_dda_vals, 'b--', lw=3., label='dipole shape')
plt.legend()

ax = plt.gca()
ax.set_xlabel('a-axis length (mm)', fontsize=32)
ax.set_ylabel('area fraction', fontsize=32)

ax.tick_params(axis='both', which='major', labelsize=28, pad=20)
ax.set_xticklabels(ax.get_xticks())
ax.set_yticklabels(ax.get_yticks())
ax.grid(color='k', linestyle=(0.5, [2,6]), linewidth=1.)

plt.savefig('afrac_accuracy.png')
