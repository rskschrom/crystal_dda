'''
Create DDA input file for branched planar crystal.

Robert Schrom @ 11/2017
'''
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import path
from polygons import make_hexagon, make_branched_planar
import geometry as geom
import os

# rotate x and y points
def rotate(x, y, angle):
    ang_rad = angle/180.*np.pi
    xrot = x*np.cos(ang_rad)-y*np.sin(ang_rad)
    yrot = x*np.sin(ang_rad)+y*np.cos(ang_rad)
    return xrot, yrot

# function to take polygon and points and return indicator
def in_polygon(xpoly, ypoly, xpoints, ypoints):
    p = path.Path([(xpoly[i], ypoly[i]) for i in range(len(xpoly))])
    points = np.empty([len(xpoints), 2])
    points[:,0] = xpoints
    points[:,1] = ypoints
    indicator = p.contains_points(points)
    return indicator

# create branched planar crystal
def branched_planar_dda(a, asp, amax, ac, ag, ft, fb, fmb, nsb):
    # test points
    numxp = 300
    numyp = 300
    x2d, y2d = np.meshgrid(np.linspace(-amax, amax, numxp),
                           np.linspace(-amax, amax, numyp), indexing='ij')
    xp = x2d.flatten()
    yp = y2d.flatten()

    # make crystal and subset to specific size a
    x, y = make_branched_planar(amax, ac, ag, ft, fb, fmb, nsb, 0.)
    inbranched = in_polygon(x, y, xp, yp)
    xp_br = xp[inbranched]
    yp_br = yp[inbranched]

    xhex, yhex = make_hexagon(a)
    inhex = in_polygon(xhex, yhex, xp_br, yp_br)
    xp_sub = xp_br[inhex]
    yp_sub = yp_br[inhex]

    # rescale to dda domain
    dx = 2.*amax/(numxp-1)
    dy = 2.*amax/(numyp-1)
    xp_dda = (xp_br+amax)/dx
    yp_dda = (yp_br+amax)/dy

    # set aspect ratio and thickness
    delta_x = np.max(xp_dda)-np.min(xp_dda)
    delta_y = np.max(yp_dda)-np.min(yp_dda)
    maxd = max(delta_x, delta_y)
    thickness = int(maxd/asp)
    print maxd, thickness, float(maxd)/float(thickness)
    z1d = np.arange(thickness)+1

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
    f.write('#box size: {:.0f}x{:.0f}x{:.0f}\n'.format(delta_x, delta_y, thickness))

    for i in range(numdip):
        f.write(' {:.0f} {:.0f} {:.0f}\n'.format(xdip[i], ydip[i], zdip[i]))
    f.close()

    # get analytical area fraction
    afrac = geom.afrac_branched(a, amax, ac, ag, ft, fb, fmb, nsb)
    return fname, afrac

