import numpy as np
from crystal_dda.polygon_utils import hexagon

class Polygon():
    '''
    The basic class for closed polygons.
    
    Parameters
    ----------
    x : ndarray, list
        The set of x points defining the polygon.
    y : ndarray, list
        The set of y points defining the polygon.
        
    Returns
    -------
    Polygon
        The returned `Polygon` class.
    '''
    def __init__(self, x, y):
        self.x = x
        self.y = y
        
class Hexagon(Polygon):
    '''
    The class for hexagons.
    
    Parameters
    ----------
    a : float
        The hexagon side length.
        
    Returns
    -------
    Hexagon
        The returned `Hexagon` class.
    '''
    def __init__(self, a):
        self.a = a
        x, y = hexagon(a)
        super(Hexagon, self).__init__(x, y)
