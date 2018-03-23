import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import crystal_dda.polygons as poly
import crystal_dda.geometry as geom
from matplotlib import path
import os

# swap dipoles
def swap_dda(xbr, ybr, xfree, yfree, numswap):
    # get random dipole locations in branched crystal
    numbr = len(xbr)
    brind = np.arange(numbr)
    np.random.shuffle(brind)
    rndind1 = brind[:numswap]

    # get random dipole locations in free space of hexagon
    numfree = len(xfree)
    freeind = np.arange(numfree)
    np.random.shuffle(freeind)
    rndind2 = freeind[:numswap]

    # swap branched-planar dipoles with free spaces
    xtmp = xbr[rndind1]
    ytmp = ybr[rndind1]
    xbr[rndind1] = xfree[rndind2]
    ybr[rndind1] = yfree[rndind2]
    xfree[rndind2] = xtmp
    yfree[rndind2] = ytmp

    return xbr, ybr, xfree, yfree

# get number of neighbor dipoles
def get_neighbors(xbr, ybr, diplen):
    numdip = len(xbr)
    numneighbors = np.empty([numdip], dtype=int)

    for i in range(numdip):
        dist = np.sqrt((xbr-xbr[i])**2.+(ybr-ybr[i])**2.)
        numneighbors[i] = len(dist[dist<=1.5*diplen])-1
        #print numneighbors[i]

    return numneighbors

# swap dipoles that are close
def swap_dda_close(xbr, ybr, xfree, yfree, swapfrac, diplen):
    # get indices of dipoles with many neighbors
    numdip = len(xbr)
    numneigh = 15
    numneighbors = get_neighbors(xbr, ybr, diplen*2.)
    print np.max(numneighbors), np.median(numneighbors)

    closeind = np.arange(numdip)[numneighbors>=numneigh]
    numclose = len(closeind)
    numswap = int(numclose*swapfrac)
    np.random.shuffle(closeind)
    rndind1 = closeind[:numswap]
    '''
    # get random dipole locations in branched crystal
    numbr = len(xbr)
    brind = np.arange(numbr)
    np.random.shuffle(brind)
    rndind1 = brind[:numswap]
    '''
    # get random dipole locations in free space of hexagon
    numfree = len(xfree)
    freeind = np.arange(numfree)
    np.random.shuffle(freeind)
    rndind2 = freeind[:numswap]

    # swap branched-planar dipoles with free spaces
    xtmp = xbr[rndind1]
    ytmp = ybr[rndind1]
    xbr[rndind1] = xfree[rndind2]
    ybr[rndind1] = yfree[rndind2]
    xfree[rndind2] = xtmp
    yfree[rndind2] = ytmp

    return xbr, ybr, xfree, yfree, numneighbors

# set branched planar crystal properties
amax = 3.
ac = 0.5

fb = 0.6
ft = 0.4
fg = 0.6

nsb = 5

ag = fg*amax+(1.-fg)*ac
wt = ft*amax/2.

# test points
numxp = 200
numyp = 200
x2d, y2d = np.meshgrid(np.linspace(-amax, amax, numxp),
                       np.linspace(-amax, amax, numyp), indexing='ij')
xp = x2d.flatten()
yp = y2d.flatten()
indicator = np.empty([numxp*numyp])
diplen = float((np.max(xp)-np.min(xp))/numxp)
print diplen

# create bounding polygons for two regions
xcore, ycore = poly.make_hexagon(ac)
xg, yg = poly.make_hexagon(ag)
xbound = np.array([-amax/2., -amax/2.+wt, 0., amax/2.-wt, amax/2.])
ybound = np.array([np.sqrt(3.)*amax/2., np.sqrt(3.)*amax/2., np.sqrt(3.)*ag/2.,
                   np.sqrt(3.)*amax/2., np.sqrt(3.)*amax/2.])

# rotations for bounding star shape
xboundr, yboundr = geom.rotate(xbound, ybound, -60.)
xbound = np.concatenate((xbound, xboundr[1:]))
ybound = np.concatenate((ybound, yboundr[1:]))

for i in range(4):
    xboundr, yboundr = geom.rotate(xboundr, yboundr, -60.)
    xbound = np.concatenate((xbound, xboundr[1:]))
    ybound = np.concatenate((ybound, yboundr[1:]))

# create dda grid of hexagon
xhex, yhex = poly.make_hexagon(amax)
inhex = geom.in_polygon(xhex, yhex, xp, yp)
xp_hex = xp[inhex]
yp_hex = yp[inhex]

# create branched planar crystal
x, y = poly.make_branched_planar(amax, ac, ft, fb, fg, nsb, 0.)
inbranched = geom.in_polygon(x, y, xp_hex, yp_hex)
xp_br = xp_hex[inbranched]
yp_br = yp_hex[inbranched]
numdip = len(yp_br)
numhex = len(yp_hex)
nn = get_neighbors(xp_br, yp_br, 2.*diplen)

# get free dipole loccations
xp_free = xp_hex[inbranched==False]
yp_free = yp_hex[inbranched==False]

# get crystal shapes at various sizes
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
mpl.rc('text', usetex=True)
numswap = 800
numiter = 25

# determine scaling factor for number of swaps after each iteration
sfrac = float(numdip-numswap)/float(numdip)

for i in range(numiter):
    # shuffle dipoles around after initial plot
    print '{:d}/{:d}'.format(i+1, numiter)
    if i>0:
        swfrac = 0.3*(float(i)/float(numiter))**(1.)-0.3*(float(i)/float(numiter))**(3.)
        xp_br, yp_br, xp_free, yp_free, nn = swap_dda_close(xp_br, yp_br, xp_free, yp_free, swfrac, diplen)
        print len(xp_br), numdip

    # plot
    plt.figure(i)
    plt.scatter(xp_br, yp_br, c=nn, cmap='seismic', s=29, edgecolor='', vmin=0, vmax=26)
    #plt.scatter(xp_br, yp_br, c='b', s=29, edgecolor='')
    cb = plt.colorbar()
    cb.set_label('\# of close dipoles', fontsize=28)
    cb_la = [ti.get_text().replace('$', '') for ti in cb.ax.get_yticklabels()]
    cb.ax.set_yticklabels(cb_la, fontsize=22)

    plt.plot(xhex, yhex, 'k--', lw=4.)

    ax = plt.gca()
    ax.set_aspect(1.)
    ax.set_xlim([-amax, amax])
    ax.set_ylim([-amax, amax])

    #ax.set_xlabel('x distance (mm)')
    #ax.set_ylabel('y distance (mm)')

    ax.tick_params(axis='both', which='major', labelsize=28, pad=20)
    ax.set_xticklabels(ax.set_xticks([]))
    ax.set_yticklabels(ax.set_yticks([]))

    ax.grid(color='k', linestyle=(0.5, [2,6]), linewidth=0.)
    imgname = 'smear_{:03d}.png'.format(i)
    plt.savefig(imgname, dpi=60)
    os.system('convert -trim {} {}'.format(imgname, imgname))

    # decrease number of swaps
    #numswap = int(numswap*sfrac)
    print numswap, sfrac
