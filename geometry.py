'''
Geometry functions

Robert Schrom @ 11/2017
'''
import numpy as np

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
    afrac_dep[avals<ac] = 1.
    afrac_dep[(avals>=ac)&(avals<ag)] = fb
    afrac_dep[avals>=ag] = fb/(amax-ag)*(ft*amax*(1.-ag/a2branch)+ag*(amax/a2branch-1.))

    # calculate true area fraction
    darea = np.sqrt(3.)/4*(avals[1:]**2.-avals[0:-1]**2.)
    afrac = np.sum(afrac_dep[0:na-1]*darea[0:na-1])/(np.sqrt(3.)/4.*a**2.)
    return afrac
