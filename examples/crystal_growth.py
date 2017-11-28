import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import crystal_dda.polygons as poly
from geom.crystal_dda import rotate, in_polygon
from matplotlib import path
import os

# set branched planar crystal properties
amax = 3.
ac = 0.5
fg = 0.7
ft = 0.2
fb = 0.4
nsb = 11

ag = fg*amax
wt = ft*amax/2.

# determine main branch fraction width (give same width as sub-branches, or 1)
wt = ft/2.*amax
wsb = 1./((nsb-1)/fb+1.)*(amax-ac-wt)
wmb = min(max(wsb/2., ac/2.), min(wsb/2., ac/2.))
fmb = wmb/(ac/2.)
print fmb, wmb, wsb

# test points
numxp = 300
numyp = 300
x2d, y2d = np.meshgrid(np.linspace(-amax, amax, numxp),
                       np.linspace(-amax, amax, numyp), indexing='ij')
xp = x2d.flatten()
yp = y2d.flatten()
indicator = np.empty([numxp*numyp])

# create bounding polygons for two regions
xcore, ycore = poly.make_hexagon(ac)
xg, yg = poly.make_hexagon(ag)
xbound = np.array([-amax/2., -amax/2.+wt, 0., amax/2.-wt, amax/2.])
ybound = np.array([np.sqrt(3.)*amax/2., np.sqrt(3.)*amax/2., np.sqrt(3.)*ag/2.,
                   np.sqrt(3.)*amax/2., np.sqrt(3.)*amax/2.])

# rotations for bounding star shape
xboundr, yboundr = rotate(xbound, ybound, -60.)
xbound = np.concatenate((xbound, xboundr[1:]))
ybound = np.concatenate((ybound, yboundr[1:]))

for i in range(4):
    xboundr, yboundr = rotate(xboundr, yboundr, -60.)
    xbound = np.concatenate((xbound, xboundr[1:]))
    ybound = np.concatenate((ybound, yboundr[1:]))

# create branched planar crystal
x, y = poly.make_branched_planar(amax, ac, ag, ft, fb, fmb, nsb, 0.)
inbranched = in_polygon(x, y, xp, yp)
xp_br = xp[inbranched]
yp_br = yp[inbranched]

# get crystal shapes at various sizes
numa = 11
avals = np.linspace(ac, amax, numa)

for a in avals:
    # subset full-size crystal to a size
    print a
    xhex, yhex = poly.make_hexagon(a)
    inhex = in_polygon(xhex, yhex, xp_br, yp_br)
    xp_sub = xp_br[inhex]
    yp_sub = yp_br[inhex]

    #plot
    plt.scatter(xp_sub, yp_sub, c='b', s=2, edgecolor='')
    plt.plot(xcore, ycore, 'r--', lw=2.)
    plt.plot(xg, yg, 'k--', lw=2.)
    plt.plot(xbound, ybound, 'g--', lw=2.)

    ax = plt.gca()
    ax.set_aspect(1.)
    ax.set_xlim([-amax, amax])
    ax.set_ylim([-amax, amax])

    ax.grid(color='k', linestyle=(0.5, [2,6]), linewidth=1.)
    imgname = 'crystal{:.1f}.png'.format(a)
    plt.savefig(imgname, dpi=150)
    os.system('convert -trim {} {}'.format(imgname, imgname))
