import numpy as np
from crystal_dda.polygon_utils import hexagon

def test_hexagon():
    x, y = hexagon(1.)
    points = np.vstack((x,y))
    points_test = np.array([[1.,0.5,-0.5,-1.,-0.5,0.5,1.],
                           [0.,0.8660254,0.8660254,0.,-0.8660254,-0.8660254,0.]])
    diff = np.sum((points-points_test)**2.)
    assert diff<1.e-6
