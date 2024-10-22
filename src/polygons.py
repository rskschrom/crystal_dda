import numpy as np
from crystal_dda.polygon_utils import hexagon, stellar, branched_planar
from crystal_dda.geometry import in_polygon

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

    def within(self, points):
        '''
        Return a boolean index list corresponding to whether each point is within the polygon boundary.

        Parameters
        ----------
        points : ndarray
            The (N,2) array of x-y points that are being evaluated.

        Returns
        -------
        indicator : ndarray
            The (N) length array of boolean indices.
        '''
        indicator = in_polygon(self.x, self.y, points[:,0], points[:,1])
        return indicator
        
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
        
class Stellar(Polygon):
    '''
    The class for stellar crystal polygons.
    
    Parameters
    ----------
    fbranch : float
        The fraction of the stellar crystal branch width relative to the hexagon side length.
    a : float
        The parent hexagon side length.
        
    Returns
    -------
    Stellar
        The returned `Stellar` class.
    '''
    def __init__(self, fbranch, a):
        self.fbranch = fbranch
        self.a = a
        x, y = stellar(fbranch, a)
        super(Stellar, self).__init__(x, y)
        
class BranchedPlanar(Polygon):
    '''
    The class for branched planar crystal polygons.
    
    Parameters
    ----------
    a : float
        The crystal planar axis length.
    amax : float
        The maximum size that characterizes the crystal geometry.
    ac : float
        The hexagonal core axis length.
    ft : float
        The main branch tip width fraction.
    fb : float
        The branch fractional coverage.
    fg : float
        The gap fraction of the crystal.
    nsb : int
        The number of subbranches of the crystal.
        
    Returns
    -------
    BranchedPlanar
        The returned `BranchedPlanar` class.
    '''
    def __init__(self, a, amax, ac, ft, fb, fg, nsb):
        self.a = a
        self.amax = amax
        self.ac = ac
        self.ft = ft
        self.fb = fb
        self.fg = fg
        self.nsb = nsb
        x, y = branched_planar(a, amax, ac, ft, fb, fg, nsb)
        super(BranchedPlanar, self).__init__(x, y)
