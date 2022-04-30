[TOC]

# Xponge

a **package** for molecular modelling

## \_\_Init\_\_

the basic **module** of Xponge

### Xponge._GlobalSetting

This **class** is used to set the global settings.

#### properties
- `farthest_bonded_force`: an **int**, the atom index distance of farthest bonded force. 4 for 'dihedral' and 5 for 'cmap'.
- `BondedForces`: the **list** to store all of the bonded force type.
- `BondedForcesMap`: the **dict** to store all of the name - bonded force type pairs. 
- `VirtualAtomTypes`: the **dict** to store all of the virtual atom type - dependent atom number pairs.
- `UnitMapping`: the **dict** to store all of the units.
- `PDBResidueNameMap`:  the **dict** to store all of the residue name mapping when loading pdb files.
- `HISMap`: the **dict** to store HIS residue information in different protonated state.
- `nocenter`: a **bool**ean value to ensure whether change the coordinate to the center when building.
- `boxspace`: a **float** to set the periodic box space between the inner molecules when building.
- 
#### instance methods

##### Add_PDB_Residue_Name_Mapping

```python
Add_PDB_Residue_Name_Mapping(self, place, pdb_name, real_name)
```
This **function** is used to add the residue name mapping to the property `PDBResidueNameMap`.

###### Input

- `palce`: a **str**ing, "head" or "tail".

- `pdb_name`: the residual name written in the pdb files.

- `real_name`: the residual name used in Xponge 

###### Output

- `None`
  
###### Example

In `Xponge.forcefield.AMBER.ff14SB`

``` python
residues = "ALA ARG ASN ASP CYS CYX GLN GLU GLY HID HIE HIP ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HIS".split()

for res in residues:
    #...    
    GlobalSetting.Add_PDB_Residue_Name_Mapping("head", res, "N" + res)
    GlobalSetting.Add_PDB_Residue_Name_Mapping("tail", res, "C" + res)
```

##### Set_Invisible_Bonded_Forces

`Set_Invisible_Bonded_Forces(self, types)`

This **function** is used to remove elements from the property `BondedForces`, and disables the corresponding types of bonded forces when building.

###### Input

- `types`: a **list** of **str**ings, the type names of the bonded forces.


###### Output

- `None`

###### Example

In `Xponge.forcefield.CHARMM27.__init__`

```python
GlobalSetting.Set_Invisible_Bonded_Forces(["improper"])
```

The code above disables the periodic improper dihedral ("improper")

##### Set_Visible_Bonded_Forces

`Set_Visible_Bonded_Forces(self, types)`

This **function** is used to replace  the property `BondedForces`, and disables the types of bonded forces except named here when building.

###### Input

- `types`: a **list** of **str**ings, the type names of the bonded forces.

###### Output

- `None`

#### static methods

##### Add_Unit_Transfer_Function

```python
Add_Unit_Transfer_Function(sometype)
```

This function is used to return a function to add a static method  `_unit_transfer` for a class. It is recommended used as a **decorator**. The origin `_unit_transfer`  method will be kept.

###### Input

- `sometype`: the class
- `somefunction`: the method

###### Example

In `Xponge.forcefield.BASE.LJ`

```python
@GlobalSetting.Add_Unit_Transfer_Function(LJType)
def LJ_Unit_Transfer(self):
    if self.A != None and self.B != None:
        if self.B == 0 or self.A == 0:
            self.sigma = 0
            self.epsilon = 0
        else:
            self.sigma = (self.A / self.B) ** (1/6)
            self.epsilon = 0.25 * self.B * self.sigma ** (-6)
        self.A = None
        self.B = None
    if self.sigma != None:
        self.rmin = self.sigma * (4 ** (1/12) / 2)
        self.sigma = None
```

The code above set the unit transfer function to convert`A`,`B` and `sigma` to `epsilon` and `rmin` of the class `LJType` .

##### Set_Unit_Transfer_Function

```python
Set_Unit_Transfer_Function(sometype)
```

This function is used to return a function to set a static method  `_unit_transfer` for a class. It is recommended used as a **decorator**.  The origin `_unit_transfer`  method will be replaced.

###### Input

- `sometype`: the class

- `somefunction`: the method

### Xponge.GlobalSetting

an **instance** of the class `Xponge._GlobalSetting`

# XpongeLib

XpongeLib is the C/C++ compiled library for Xponge.

## \_\_init\_\_

### XpongeLib.\_parmchk2

``` python
XpongeLib._parmchk2(i, iformat, o, datapath, print_all, print_dihedral_contain_X,  gaffORgaff2)
```
This **function** is used to check the force field parameters for gaff and gaff2
#### Input
- `i`: input file name
- `iformat`: input file format (prepi, prepc, ac, mol2, frcmod, leaplog)
- `o`: output file name
- `datapath`: the path to store the parmchk data
- `print_all`: print out all force field parameters including those in the parmfile or not. 0 for no and 1 for yes.
- `print_dihedral_contain_X`: print out parameters that matching improper dihedral parameters that contain 'X' in the force field parameter file. 0 for no and 1 for yes.
- `gaffORgaff2`: 1 for gaff and 2 for gaff2

#### Output:
- `None`

#### Example:
the implement of `Xponge.forcefield.AMBER.gaff.parmchk2_gaff`
``` python
def parmchk2_gaff(ifname, ofname, direct_load = True, keep = True):
    import XpongeLib as xlib
    import os
    datapath = os.path.split(xlib.__file__)[0]
    xlib._parmchk2(ifname, "mol2", ofname, datapath, 0, 1, 1)
    if direct_load:
        atoms, bonds, angles, propers, impropers, LJs, cmap = LOAD.frcmod(ofname)
        BOND.BondType.New_From_String(bonds)
        ANGLE.AngleType.New_From_String(angles)
        DIHEDRAL.ProperType.New_From_String(propers)
        DIHEDRAL.ImproperType.New_From_String(impropers)
    if not keep:
        os.remove(ofname)
```