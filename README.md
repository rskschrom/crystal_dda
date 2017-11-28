# crystal_dda
A package for generating branched planar crystal polygons and input text files for DDA code

# Installation
Clone package to local machine. Then cd into directory and install with:

```
python setup.py install --user
```

Numpy and Matplotlib are required dependencies.

# Usage

The main function to create a branched planar dda file is ```branched_planar_dda``` which we import with:

```
from crystal_dda.crystal_dda import branched_planar_dda
```

Next, we set the values defining the branched planar crystal shape appropriately. For example:

```
a = 0.2
amax = 3.
ac = 0.5

fb = 0.4
ft = 0.2
fg = 0.7
fmb = 0.4

nsb = 5
asp = 20.
ag = amax*fg
```
Finaly, we create the DDA input file with:

```
fname, afrac = branched_planar_dda(a, asp, amax, ac, ag, ft, fb, fmb, nsb)
```

where ```fname``` is the name of the input file that we just created and ```afrac``` is the estimated area fraction of the crystal. See the examples directory for more information. We can now run the DDA code with this file as an input to calculate the scattering.
