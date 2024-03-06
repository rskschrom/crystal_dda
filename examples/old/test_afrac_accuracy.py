'''
Code to generate DDA files and compare accuracy of analytical area fraction.

Robert Schrom @ 11/2017
'''
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from crystal_dda.geometry import afrac_branched, afrac_dda_subregion
import os

# set values to create branched planar crystal with
amax = 3.
ac = 0.5

fb = 0.4
ft = 0.3
fg = 0.5

nsb = 5
nxp = 300
nzp = 5

# loop over a axis lengths for given crystal
numa = 201
avals = np.linspace(0.1, amax, numa)
afrac_new = np.empty([numa])
afrac_old = np.empty([numa])

for i, a in enumerate(avals):
    afrac_new[i] = afrac_branched(a, amax, ac, ft, fb, fg, nsb)
    if a > ac:
        afrac_old[i] = (ac/a)**2.+(1.-(ac/a)**2.)*0.27
    else:
        afrac_old[i] = 1.

# plot
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
mpl.rc('text', usetex=True)

plt.plot(avals, afrac_new*917., 'r--', lw=3., label='New effective density')
plt.plot(avals, afrac_old*917., 'b--', lw=3., label='Old effective density')
plt.legend(fontsize=26)

ax = plt.gca()
ax.set_xlabel('a-axis length (mm)', fontsize=26)
ax.set_ylabel('Effective density (kg m$^{\sf{-3}}$)', fontsize=26)

ax.tick_params(axis='both', which='major', labelsize=24, pad=20)
ax.set_xticklabels(ax.get_xticks())
ax.set_yticklabels(ax.get_yticks())
ax.grid(color='k', linestyle=(0.5, [2,6]), linewidth=1.)

plt.savefig('afrac.png')
os.system('convert -trim {} {}'.format('afrac.png', 'afrac.png'))
