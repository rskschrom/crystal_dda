'''
Create DDA input file for branched planar crystal.

Robert Schrom @ 11/2017
'''
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from . polygons import make_hexagon, make_branched_planar
from . import geometry as geom
import os

# create full-size crystals in x-y plane
def branched_planar_struct(amax=3., ac=0.1, ft=0.4, fb=0.5,
                           fg=0.3, nsb=5, numxp=150):
    # test points
    numyp = numxp
    x2d, y2d = np.meshgrid(np.linspace(-amax, amax, numxp),
                           np.linspace(-amax, amax, numyp), indexing='ij')
    xp = x2d.flatten()
    yp = y2d.flatten()

    # make hexagon first
    xhex, yhex = make_hexagon(amax)
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

    return xp_br, yp_br, afrac

# get dda file of subset full-size branched planar
def subset_branched_dda(a, amax, xp_br, yp_br, afrac, numxp=150, numzp=9):
    # make hexagon first
    xhex, yhex = make_hexagon(a)
    inhex = geom.in_polygon(xhex, yhex, xp_br, yp_br)
    xp_sub = xp_br[inhex]
    yp_sub = yp_br[inhex]

    # rescale to dda domain
    dx = 2.*amax/(numxp-1)
    dy = 2.*amax/(numxp-1)
    xp_dda = xp_sub/dx
    yp_dda = yp_sub/dy
    xp_dda = xp_dda-np.min(xp_dda)+1
    yp_dda = yp_dda-np.min(yp_dda)+1

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
        f.write(' {:.1f} {:.1f} {:.1f}\n'.format(xdip[i], ydip[i], zdip[i]))
    f.close()

    return fname, afrac

# create branched planar crystal
def branched_planar_dda(a=3., amax=3., ac=0.1, ft=0.4, fb=0.5,
                        fg=0.3, nsb=5, numxp=150, numzp=9, outdir='',ind=0):
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
    fname = f'crystal_{ind:04d}_{a:.1f}.txt'
    f = open(f'{outdir}{fname}', "w")
    f.write('#generated by crystal_dda.py\n')
    f.write('#shape: \'read\'\n')
    f.write('#box size: {:.0f}x{:.0f}x{:.0f}\n'.format(delta_x, delta_y, numzp))

    for i in range(numdip):
        f.write(' {:.0f} {:.0f} {:.0f}\n'.format(xdip[i], ydip[i], zdip[i]))
    f.close()

    return fname, afrac

# create branched planar crystal
def branched_planar_dda_chunk(a=3., amax=3., ac=0.1, ft=0.4, fb=0.5,
                              fg=0.3, nsb=5, numxp=150, numzp=9, nchunk=10,
                              chxi=0, chyi=0, fname='crystal.txt'):
    # test chunk of points
    x1d = np.linspace(-a, a, numxp)
    y1d = np.linspace(-a, a, numxp)
    nxch = int(numxp/nchunk)+1
    xchunk = x1d[nxch*chxi:min(numxp,nxch*(chxi+1))]
    ychunk = y1d[nxch*chyi:min(numxp,nxch*(chyi+1))]
    #print(nxch, nxch*chxi, nxch*(chxi+1), numxp)

    x2d, y2d = np.meshgrid(xchunk, ychunk, indexing='ij')
    xp = x2d.flatten()
    yp = y2d.flatten()

    # make hexagon first
    xhex, yhex = make_hexagon(a)
    inhex = geom.in_polygon(xhex, yhex, xp, yp)
    xp_hex = xp[inhex]
    yp_hex = yp[inhex]
    ndiphex = len(xp_hex)

    # determine main branch fraction width (give same width as sub-branches, or 1)
    fmb = geom.frac_main_branch(amax, ac, ft, fb, nsb)

    # determine which hexagon points are in branched planar
    xbr, ybr = make_branched_planar(amax, ac, ft, fb, fg, nsb, 0.)
    inbranched = geom.in_polygon(xbr, ybr, xp_hex, yp_hex)
    xp_br = xp_hex[inbranched]
    yp_br = yp_hex[inbranched]
    ndip = len(xp_br)

    if ndip>0:
        # rescale to dda domain
        dx = 2.*a/(numxp-1)
        dy = 2.*a/(numxp-1)
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
        f = open(fname, "a")
        for i in range(numdip):
            f.write(' {:.0f} {:.0f} {:.0f}\n'.format(xdip[i], ydip[i], zdip[i]))
        f.close()

    return ndip, ndiphex, fname

# create plate crystal
def plate_dda(a, numxp, numzp):
    # test points
    numyp = numxp
    x2d, y2d = np.meshgrid(np.linspace(-a, a, numxp),
                           np.linspace(-a, a, numyp), indexing='ij')
    xp = x2d.flatten()
    yp = y2d.flatten()

    # make hexagon
    xhex, yhex = make_hexagon(a)
    inhex = geom.in_polygon(xhex, yhex, xp, yp)
    xp_hex = xp[inhex]
    yp_hex = yp[inhex]

    # rescale to dda domain
    dx = 2.*a/(numxp-1)
    dy = 2.*a/(numyp-1)
    xp_dda = (xp_hex+a)/dx
    yp_dda = (yp_hex+a)/dy

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
    fname = 'plate.txt'
    f = open(fname, "w")
    f.write('#generated by crystal_dda.py\n')
    f.write('#shape: \'read\'\n')
    f.write('#box size: {:.0f}x{:.0f}x{:.0f}\n'.format(delta_x, delta_y, numzp))

    for i in range(numdip):
        f.write(' {:.0f} {:.0f} {:.0f}\n'.format(xdip[i], ydip[i], zdip[i]))
    f.close()

    return fname

# get area fraction of real branched planar crystal
def branched_planar_afrac(a, amax, ac, ft, fb, fg, nsb, numxp, numzp):
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

    return afrac

# create branched planar crystal
def branched_planar_dipoles(a, amax, ac, ft, fb, fg, nsb, numxp, numzp):
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
    thick = 2.*a*float(numzp)/float(numxp)
    z1d = np.linspace(0., thick, numzp)

    # add z dimension 
    x3, z3 = np.meshgrid(xp_br, z1d, indexing='ij')
    y3, z3 = np.meshgrid(yp_br, z1d, indexing='ij')
    xdip = x3.flatten()
    ydip = y3.flatten()
    zdip = z3.flatten()

    return xdip, ydip, zdip
