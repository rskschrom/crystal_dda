import matplotlib.pyplot as plt
import numpy as np
from crystal_dda.polygons import Polygon, Hexagon
from crystal_dda.geometry import in_polygon

class Crystal:
    '''
    Class to create a crystal shape for ADDA input.
    
    Parameters
    ----------
    poly : Polygon
        The `Polygon` object defining the 2D planar shape of the crystal.
    c : float
        The crystal basal-length dimension (c = h/2, where h is the crystal thickness).
        
    Returns
    -------
    Crystal
        The returned `Crystal` object.
    '''
    def __init__(self, poly, c):
        self.poly = poly
        self.c = c
        self.ix = np.array([])
        self.iy = np.array([])
        self.iz = np.array([])
        self.ndip = 0
        
    def create_dipoles(self, dip_len):
        '''
        Fill the crystal with dimensionless dipole indices based on the crystal physical dimensions and the dipole spacing `dip_len`.
        
        Parameters
        ----------
        dip_len : float
            The dipole spacing.
        
        Returns
        -------
        None
        '''
        # fill in polygon with dipoles on grid spanned by polygon dmax
        dmax = self.poly.size()
        ndipx = int(dmax/dip_len)+1
        xg2, yg2 = np.meshgrid(np.linspace(-dmax/2., dmax/2, ndipx),
                               np.linspace(-dmax/2., dmax/2., ndipx), indexing='ij')
        xgf= xg2.flatten()
        ygf = yg2.flatten()

        inpoly = in_polygon(self.poly.x, self.poly.y, xgf, ygf)
        x_dip_poly = xgf[inpoly]
        y_dip_poly = ygf[inpoly]

        # rescale to dda domain (with min index values of 1)
        ix_dip_poly = x_dip_poly/dip_len
        iy_dip_poly = y_dip_poly/dip_len
        ix_dip_poly = ix_dip_poly-np.min(ix_dip_poly)+1
        iy_dip_poly = iy_dip_poly-np.min(iy_dip_poly)+1

        # extrude 2d dipole indices to 3d along z axis
        ndipz = int(2.*self.c/dip_len)+1
        iz1d = np.arange(ndipz)+1
        
        ix3, iz3 = np.meshgrid(ix_dip_poly, iz1d, indexing='ij')
        iy3, iz3 = np.meshgrid(iy_dip_poly, iz1d, indexing='ij')
        ix_dip = ix3.flatten()
        iy_dip = iy3.flatten()
        iz_dip = iz3.flatten()
        ndip = len(ix_dip)
        
        # set dipole indices for crystal
        self.ix = ix_dip
        self.iy = iy_dip
        self.iz = iz_dip
        self.ndip = ndip
        
        return

    def write_dipoles(self, file_name):
        '''
        Write the dipole indices to a file.
        
        Parameters
        ----------
        file_name : str
            The output file name.
        
        Returns
        -------
        None
        '''
        if len(self.ix)>0:
            # get bounding box of dipoles
            dx = np.max(self.ix)
            dy = np.max(self.iy)
            dz = np.max(self.iz)
        
            header = "#generated by crystal-dda\n"+\
                     "#shape: 'read'\n"+\
                    f"#box size: {dx:.0f}x{dy:.0f}x{dz:.0f}\n"
            np.savetxt(file_name, np.c_[self.ix, self.iy, self.iz],
                   fmt=('%d','%d','%d'), header=header, comments='')
        else:
            # no dipoles have been filled
            print('error: dipoles have not been filled for this crystal.'
                  ' try running create_dipoles() first.')
        return
