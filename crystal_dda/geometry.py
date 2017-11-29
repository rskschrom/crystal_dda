'''
Geometry functions

Robert Schrom @ 11/2017
'''
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import path

# function to take polygon and points and return indicator
def in_polygon(xpoly, ypoly, xpoints, ypoints):
    p = path.Path([(xpoly[i], ypoly[i]) for i in range(len(xpoly))])
    points = np.empty([len(xpoints), 2])
    points[:,0] = xpoints
    points[:,1] = ypoints
    indicator = p.contains_points(points)
    return indicator

# rotation
def rotate(x, y, angle):
    angle = angle*np.pi/180.
    xrot = np.cos(angle)*x-np.sin(angle)*y
    yrot = np.sin(angle)*x+np.cos(angle)*y
    return xrot, yrot

# flip over y-axis
def flipy(x, y):
    xflip = -x
    yflip = y
    return xflip, yflip

# get equation of line from p1 to p2
def points2eqn(p1x, p1y, p2x, p2y):
    m = (p2y-p1y)/(p2x-p1x)
    b = p1y-m*p1x
    return m, b

# get point of intersection of two lines
def intersection(m1, b1, m2, b2):
    xint = (b2-b1)/(m1-m2)
    yint = m1*xint+b1
    return xint, yint

# get area fraction of stellar
def afrac_stellar(fbranch):
    afrac = 1.-(1.-fbranch)**2.
    return afrac

# get area fraction of branched planar
def afrac_branched(a, amax, ac, ag, ft, fb, fmb, nsb):
    # set particle sizes (mm)
    na = 500
    avals = np.linspace(0., a, na)

    # calculate area fraction (deposition)
    afrac_dep = np.ma.masked_all([na])
    a2branch = avals[avals>=ag]
    a2greg = avals[(avals>=ac)&(avals<ag)]
    afrac_dep[avals<ac] = 1.
    afrac_dep[(avals>=ac)&(avals<ag)] = fb+fmb*2.*(1.-fb)/(a2greg/ac+1.)
    afrac_dep[avals>=ag] = fb/(amax-ag)*(ft*amax*(1.-ag/a2branch)+
                           ag*(amax/a2branch-1.))+fmb*2.*(1.-fb)/(a2branch/ac+1.)

    # calculate true area fraction
    darea = np.sqrt(3.)/4*(avals[1:]**3.-avals[0:-1]**3.)
    afrac = np.sum(afrac_dep[0:na-1]*darea[0:na-1])/(np.sqrt(3.)/4.*a**3.)
    return afrac

# formulation like microphysics code
def afrac_branched_alt(a, amax, ac, ag, ft, fb, fmb, nsb):
            if (a.lt.ac) then
               RhoDep = 917.
            endif
            if (a.ge.ac.and.a.lt.ag) then
               RhoDep = 917.*fb+917.*fmb*2.*(1.-fb)/(a/ac+1.)
            endif
            if (acur.ge.ag) then
               print *, ag/acur, amax/a
               RhoDep = fb/(amax-ag)*ag*(amax/a-1.)
               RhoDep = 917.*fb/(amax-ag)*(ft*amax*(1.-ag/a)+ag*(amax/a-1.))+917.*fmb*2.*(1.-fb)/(a/ac+1.)
               RhoDep = 917.*fb
            endif

# calculate area fraction in hexagonal region within dda polygon using random points
def afrac_dda_subregion(xhex, yhex, xdda_poly, ydda_poly):
    # probability method for estimating area fraction
    numtrials = 20000
    a = np.max(xhex)
    rx = 2.*a*(np.random.rand(numtrials)-0.5)
    ry = 2.*a*(np.random.rand(numtrials)-0.5)

    # determine random points in hexagon
    inhex = in_polygon(xhex, yhex, rx, ry)
    rx_hex = rx[inhex]
    ry_hex = ry[inhex]
    numhex = len(rx_hex)

    # determine which hexagon points are in polygon defining full size dda polygon
    inbranched = in_polygon(xdda_poly, ydda_poly, rx_hex, ry_hex)
    rx_br = rx_hex[inbranched]
    ry_br = ry_hex[inbranched]
    numbr = len(rx_br)

    frac = float(numbr)/float(numhex)
    return frac
