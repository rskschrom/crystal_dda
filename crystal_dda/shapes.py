class Crystal:
    '''
    Class to create a crystal shape for ADDA input.
    
    Parameters
    ----------
    a : float
        The crystal a-axis dimension.
    c : float
        The crystal c-axis dimension.
        
    Returns
    -------
    Crystal
        The returned `Crystal` object.
    '''
    def __init__(self, a, c):
        self.a = a
        self.c = c
        self.ix = None
        self.iy = None
        self.iz = None
        
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
        return
