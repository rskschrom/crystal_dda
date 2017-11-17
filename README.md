# crystal_dda
Modules and script for generating branched planar crystal polygons and input text files for DDA code

# Installation
Download all python files to a a local machine. Nicer Python package installation to come (maybe). Need Numpy and Matplotlib to run these files.

# Usage
Set the variables in ``make_crystal_dda.py`` appropriately and then run with:

```
python make_crystal_dda.py
```
The name of the file you just generated will be printed out as well as the area fraction of the crystal. Then you can run the DDA code with this file as an input.
