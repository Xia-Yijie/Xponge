import XpongeLib as xlib
import os

datapath = os.path.split(xlib.__file__)[0]
print(datapath)
xlib._parmchk2("temp.mol2", "mol2", "temp.frcmod", datapath, 0, 1, 1)
