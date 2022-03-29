# Xponge

## \_\_init\_\_

# XpongeLib

XpongeLib is the C/C++ compiled library for Xponge.

## \_\_init\_\_

### \_parmchk2
``` python
XpongeLib._parmchk2(i, iformat, o, datapath, print_all, print_dihedral_contain_X,  gaffORgaff2)
```
    This function is used to check the force field parameters for gaff and gaff2
#### Input:
#####i
    input file name
##### iformat
    input file format (prepi, prepc, ac, mol2, frcmod, leaplog)
##### o
    output file name
##### datapath
    the path to store the parmchk data
##### print\_all
    print out all force field parameters including those in the parmfile or not. 0 for no and 1 for yes.
##### print\_dihedral\_contain\_X
    print out parameters that matching improper dihedral parameters that contain 'X' in the force field parameter file. 0 for no and 1 for yes.
##### gaffORgaff2
    1 for gaff and 2 for gaff2

#### Output:
    None
