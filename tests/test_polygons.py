import numpy as np
from crystal_dda.polygons import Hexagon

def test_hexagon():
    hx = Hexagon(1.)
    x = hx.x
    y = hx.y
    points = np.vstack((x,y))
    points_test = np.array([[1.,0.5,-0.5,-1.,-0.5,0.5,1.],
                           [0.,0.8660254,0.8660254,0.,-0.8660254,-0.8660254,0.]])
    diff = np.sum((points-points_test)**2.)
    assert diff<1.e-6
