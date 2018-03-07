'''
Create DDA input file for branched planar crystal.

Robert Schrom @ 11/2017
'''
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from polygons import make_hexagon, make_branched_planar
import geometry as geom
import os

# create branched planar crystal
def branched_planar_dda(a, amax, ac, ft, fb, fg, nsb, numxp, numzp):
    # test points
    numyp = numxp
    x2d, y2d = np.meshgrid(np.linspace(-a, a, numxp),
                           np.linspace(-a, a, numyp), indexing='ij')
    xp = x2d.flatten()
    yp = y2d.flatten()

    # make hexagon first
    xhex, yhex = make_hexagon(a)
    inhex = geom.in_polygon(xhex, yhex, xp, yp)
    xp_hex = xp[inhex]
    yp_hex = yp[inhex]

    # determine main branch fraction width (give same width as sub-branches, or 1)
    fmb = geom.frac_main_branch(amax, ac, ft, fb, nsb)

    # determine which hexagon points are in branched planar
    xbr, ybr = make_branched_planar(amax, ac, ft, fb, fg, nsb, 0.)
    inbranched = geom.in_polygon(xbr, ybr, xp_hex, yp_hex)
    xp_br = xp_hex[inbranched]
    yp_br = yp_hex[inbranched]
    afrac = float(len(xp_br))/float(len(xp_hex))

    # rescale to dda domain
    dx = 2.*a/(numxp-1)
    dy = 2.*a/(numyp-1)
    xp_dda = (xp_br+a)/dx
    yp_dda = (yp_br+a)/dy

    # set number of z dipoles
    delta_x = np.max(xp_dda)-np.min(xp_dda)
    delta_y = np.max(yp_dda)-np.min(yp_dda)
    z1d = np.arange(numzp)+1

    # add z dimension 
    x3, z3 = np.meshgrid(xp_dda, z1d, indexing='ij')
    y3, z3 = np.meshgrid(yp_dda, z1d, indexing='ij')
    xdip = x3.flatten()
    ydip = y3.flatten()
    zdip = z3.flatten()
    numdip = len(xdip)

    # write to file
    fname = 'crystal{:.1f}.txt'.format(a)
    f = open(fname, "w")
    f.write('#generated by crystal_dda.py\n')
    f.write('#shape: \'read\'\n')
    f.write('#box size: {:.0f}x{:.0f}x{:.0f}\n'.format(delta_x, delta_y, numzp))

    for i in range(numdip):
        f.write(' {:.0f} {:.0f} {:.0f}\n'.format(xdip[i], ydip[i], zdip[i]))
    f.close()

    return fname, afrac

