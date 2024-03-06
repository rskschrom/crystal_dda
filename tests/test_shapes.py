import numpy as np
from crystal_dda.polygons import Hexagon
from crystal_dda.shapes import Crystal

# create crystal object
def create_crystal():
    hx = Hexagon(1.)
    cr = Crystal(hx, 0.1)
    return cr
    
def test_crystal():
    cr = create_crystal()
    assert 0 == 0
    
def test_create_dipoles():
    cr = create_crystal()
    cr.create_dipoles(0.04)
    assert cr.ndip == 9702
    
def test_write_dipoles():
    cr = create_crystal()
    cr.create_dipoles(0.04)
    cr.write_dipoles('test.txt')
    assert 0 == 0
