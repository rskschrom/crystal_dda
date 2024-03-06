import numpy as np
from crystal_dda.polygon_utils import hexagon

class Polygon():
    '''
    The basic class for closed polygons.
    
    Parameters
    ----------
    x : ndarray
        The set of x points defining the polygon.
    y : ndarray
        The set of y points defining the polygon.
        
    Returns
    -------
    Polygon
        The returned `Polygon` class.
    '''
    def __init__(self, x, y):
        self.x = x
        self.y = y
        
    def size(self):
        '''
        Get the maximum dimension for the polygon.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        dmax : float
            The maximum distance between polygon points.
        '''
        npoint = len(self.x)
        
        dmax = 0.
        
        # iterate over points and find largest distances
        for i in range(npoint):
            dmax_test = np.max((self.x[i]-self.x)**2.+(self.y[i]-self.y)**2.)
            dmax = max(dmax, dmax_test)
            
        dmax = np.sqrt(dmax)
        return dmax
        
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
