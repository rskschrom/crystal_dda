import numpy as np
import pyvista as pv
import polygons as pg
import geometry as geom
import matplotlib.pyplot as plt
import triangle as tr
import pymeshfix as pmf

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

'''
# triangulate flipped points to complete sector
points_flip = points[:]
points_flip[:,0] = -points_flip[:,0]

A = dict(vertices=points_flip, segments=segs)
B = tr.triangulate(A, 'qpa0.05')
vert_tr_flip = B['vertices']
faces_tr_flip = B['triangles']

x_tr0_flip = vert_tr_flip[:,0]
y_tr0_flip = vert_tr_flip[:,1]

x_tr0 = np.concatenate((x_tr0,x_tr0_flip))
y_tr0 = np.concatenate((y_tr0,y_tr0_flip))

faces_tr_flip = np.concatenate((np.full([faces_tr_flip.shape[0],1], 3), faces_tr_flip), axis=1)

faces_tr0 = np.concatenate((faces_tr,faces_tr_flip))
faces_tr0[ntr:,1:] = faces_tr0[ntr:,1:]+ntr

ntr = faces_tr0.shape[0]
print(ntr, x_tr0.shape)

z_tr0 = np.full([ntr], -c)

points_tr = np.vstack((x_tr0,y_tr0,z_tr0)).T
surf = pv.PolyData(points_tr, faces_tr0)
surf.plot(show_edges=True)
'''
'''
# rotate sector 60 degrees and concatenate with previous sector
x_tr = x_tr0[:]
y_tr = y_tr0[:]
faces_tr = faces_tr0[:]

for i in range(5):
    xrot, yrot = geom.rotate(x_tr0, y_tr0, -60.*(i+1))
    x_tr = np.concatenate((x_tr, xrot))
    y_tr = np.concatenate((y_tr, yrot))

    faces_tr = np.concatenate((faces_tr,faces_tr0))
    faces_tr[(i+1)*ntr:,1:] = faces_tr[(i+1)*ntr:,1:]+(i+1)*ntr
    
z_tr = np.full([6*ntr], c)
points_tr = np.vstack((x_tr,y_tr,z_tr)).T
surf_top = pv.PolyData(points_tr, faces_tr)

mesh = surf_top.extrude((0.,0.,2*c), capping=True)
mesh.flip_normals()
mesh.plot(show_edges=True)
print(mesh.is_manifold)
'''
'''
#surf_top.plot_normals(mag=0.1, faces=True, show_edges=True)

#surf_top.save('top.stl')

# get bottom polygon face
for i in range(6*ntr):
    f1 = faces_tr[i,1]
    f3 = faces_tr[i,3]
    faces_tr[i,1] = f3
    faces_tr[i,3] = f1

points_tr[:,2] = -points_tr[:,2]
surf_bot = pv.PolyData(points_tr, faces_tr)
#surf_bot.plot_normals(mag=0.1, faces=True, show_edges=True)

# create meshes on side of particle
xsp = np.array([])
ysp = np.array([])
zsp = np.array([])
npf = len(xf)

for i in range(npf-1):
    xsp = np.concatenate((xsp,np.array([xf[i],xf[i+1],xf[i+1],xf[i]])))
    ysp = np.concatenate((ysp,np.array([yf[i],yf[i+1],yf[i+1],yf[i]])))
    zsp = np.concatenate((zsp,np.array([c,c,-c,-c])))
  
    if i==0:
        # triangles
        #faces_s = np.array([[4*i,4*i+1,4*i+3,],
        #                    [4*i+1,4*i+2,4*i+3]])
                            
        # quad
        faces_s = np.array([[4*i,4*i+1,4*i+2,4*i+3]])
    else:
        # triangles
        #faces_s = np.concatenate((faces_s,np.array([[4*i,4*i+1,4*i+3,],
        #                                            [4*i+1,4*i+2,4*i+3]])), axis=0)
        
        # quad
        faces_s = np.concatenate((faces_s,np.array([[4*i,4*i+1,4*i+2,4*i+3]])), axis=0)

points_s = np.vstack((xsp,ysp,zsp)).T
faces_s = np.concatenate((np.full([faces_s.shape[0],1], 4), faces_s), axis=1)
surf_s = pv.PolyData(points_s, faces_s)

# remove duplicate faces
nface = faces_s.shape[0]
face_centers = np.empty([nface,3])
face_norms = surf_s.cell_normals

for i in range(nface):
    face_centers[i,:] = np.mean(points_s[faces_s[i,1:],:], axis=0)

tol = a/1.e6
face_data = np.concatenate((face_centers/tol,np.abs(face_norms)/tol), axis=1)
_, ui, uc = np.unique(face_data.astype(int), return_index=True, return_counts=True, axis=0)

# keep only faces that occur exactly once
ui = ui[uc==1]
faces_s = faces_s[ui,:]
surf_s = pv.PolyData(points_s, faces_s)
mesh = surf_s.extrude((0.,0.,1), capping=True)
mesh.plot(show_edges=True)
#surf_s.plot_normals(mag=0.1, faces=True, show_edges=True)



'''
'''
meshfix = pmf.MeshFix(surf_s)
meshfix.repair(verbose=True)
mesh_fix = meshfix.mesh
#mesh_fix.plot(show_edges=True)
'''
'''
# merge side, top, and bottom
#surf_merge = surf_top.merge(surf_s)
#surf_merge = surf_merge.merge(surf_bot)
surf_merge = surf_top+surf_bot+surf_s
surf_merge = surf_merge.subdivide(3)
#surf_merge = surf_merge.clean()
points_sfc = pv.PolyData(surf_merge.cell_centers())

sfc_test = points_sfc.reconstruct_surface()

#surf_merge.plot_normals(mag=0.1, faces=True, show_edges=True)
print(sfc_test.is_manifold)

#voxels = pv.voxelize(surf_merge, density=2*a/150)
edges = surf_merge.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False)
# repair and save mesh

meshfix = pmf.MeshFix(sfc_test)
meshfix.repair(verbose=True)
mesh_fix = meshfix.mesh
print(mesh_fix.is_manifold)

surf_merge.save('crystal.stl')

# plot
pl = pv.Plotter()
#pl.add_mesh(surf, show_edges=True, line_width=1)
#pl.add_mesh(surf_s, show_edges=True, line_width=1)
#pl.add_mesh(surf_merge, show_edges=True, line_width=1)
#pl.add_mesh(edges, color='r', line_width=1.)
pl.add_mesh(mesh_fix)
pl.show()
'''
