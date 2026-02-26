import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from crystal_dda.polygons import BranchedPlanar
from crystal_dda.shapes import Crystal

# get points inside the filled crystal
def get_crystal_points(a, a_max, ac, ft, fb, fg, nsb, dip_len, c=0.1):
    bp = bp = BranchedPlanar(a, a_max, ac, ft, fb, fg, nsb)
    
    # create 3d crystal
    cr = Crystal(bp, c)
    cr.create_dipoles(dip_len)

    # convert dipoles to physical dimensions and center them
    crx = dip_len*cr.ix
    cry = dip_len*cr.iy
    crx = crx-np.mean(crx)
    cry = cry-np.mean(cry)
    return crx, cry

# randomly sample crystal parameters
nrnd = 100
rng = np.random.default_rng()
ac = 0.2*rng.random(nrnd)+0.05
fb = 0.8*rng.random(nrnd)+0.1
ft = rng.random(nrnd)
fg = rng.random(nrnd)

# create crystals and plot (color by fb)
nplot = 100
a = 3.
nsb = 4
fig = plt.figure(figsize=(15,15))
cmap = plt.get_cmap('magma_r')

# set dipole length for crystal
dip_len = 0.1

for i in range(nplot):
    crx, cry = get_crystal_points(a, a, ac[i], ft[i], fb[i], fg[i], nsb, dip_len)

    # plot crystal in panel
    ax = fig.add_subplot(10,10,i+1)
    plt.scatter(crx, cry, c=np.expand_dims(cmap(fb[i]), axis=0), s=0.5)
    ax.set_aspect(1.)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.set_xlim([-a, a])
    ax.set_ylim([-a, a])

# add colorbar for fb
plt.scatter(crx+1.e9, cry+1.e9, c=crx, cmap='magma_r', vmin=0., vmax=1.)

fig.subplots_adjust(right=0.8)
cb_ax = fig.add_axes([0.83, 0.1, 0.03, 0.8])
cb = plt.colorbar(cax=cb_ax)
cb.ax.tick_params(labelsize=18)
cb.set_label('$f_b$', fontsize=26)
imgname = 'crystals.png'
plt.savefig(imgname, bbox_inches="tight")
