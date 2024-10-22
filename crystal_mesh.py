import numpy as np
import pyvista as pv
import polygons as pg
import geometry as geom
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

# create crystal shape
a = 1.5
c = 0.1
amax = 3.
ac = 0.2

fb = 0.3
ft = 0.2
fg = 0.7

nsb = 5

xf, yf = pg.branched_planar(amax, a, ac, ft, fb, fg, nsb)
#x, y = pg.branched_planar_components(amax, a, ac, ft, fb, fg, nsb)

# create mesh of polygon surface
npoint = len(xf)
points = np.vstack((xf, yf)).T

segs = np.stack([np.arange(npoint), np.arange(npoint) + 1], axis=1) % npoint

A = dict(vertices=points, segments=segs)
B = tr.triangulate(A,)
vert_tr = B['vertices']
faces_tr = B['triangles']

ntr = faces_tr.shape[0]
faces_tr = np.concatenate((np.full([ntr,1], 3), faces_tr), axis=1)

x_tr0 = vert_tr[:,0]
y_tr0 = vert_tr[:,1]
npoint = len(x_tr0)
z_tr0 = np.full([npoint], -c)

surf = pv.PolyData(np.vstack((x_tr0,y_tr0,z_tr0)).T, faces_tr)
surf.plot(show_edges=True)

