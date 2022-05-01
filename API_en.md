[TOC]

# Xponge

The  **package** for molecular modelling

## \_\_Init\_\_

The basic **module** of Xponge

### Xponge._GlobalSetting

This **class** is used to set the global settings.

#### class variables 
- `farthest_bonded_force`: an **int**, the atom index distance of farthest bonded force. 4 for 'dihedral' and 5 for 'cmap'.
- `BondedForces`: a **list** to store all of the bonded force type.
- `BondedForcesMap`: a **dict** to store all of the name - bonded force type pairs. 
- `VirtualAtomTypes`: a **dict** to store all of the virtual atom type - dependent atom number pairs.
- `UnitMapping`: a **dict** to store all of the units.
- `PDBResidueNameMap`:  a **dict** to store all of the residue name mapping when loading pdb files.
- `HISMap`: a **dict** to store HIS residue information in different protonated state.
- `nocenter`: a **bool**ean value to ensure whether change the coordinate to the center when building.
- `boxspace`: a **float** to set the periodic box space between the inner molecules when building.
#### instance methods

##### Add_PDB_Residue_Name_Mapping

```python
Add_PDB_Residue_Name_Mapping(self, place, pdb_name, real_name)
```
This **function** is used to add the residue name mapping to the property `PDBResidueNameMap`.

###### Input

- `palce`: the **str** of the residue place, "head" or "tail".

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

- `types`: a **list** of **str**, the type names of the bonded forces.


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

- `types`: a **list** of **str**s, the type names of the bonded forces.

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

### Xponge.Type

This **class** is the abstract class of the types (atom types, bonded force types and so on).

#### Instance Variables

- `contents`: a **dict** to store the parameter - value pairs

#### class variables

- `name`: a **str**, the name of the class

- `parameters`: a **dict** to store the parameter - data type pairs for this class

- `types`: a **dict** to store all the name - instance pairs of the types for this class, one name for one instance.

- `types_different_name`: a **dict** to store all the name - instance pairs of the types for this class, some different names may be mapped to the same instance.

#### class methods

##### Add_Property
```python
Add_Property(cls, parm_fmt, parm_default = {})
```

This **function** is used to add a property to the class

###### Input

- `parm_fmt`: a **dict** to store the parameter-data type pair(s)
- `parm_default`: a **dict** to store the parameter-default value pair(s)

###### Output

- `None`

###### Example

In `Xponge.forcefield.BASE.LJ`

```python
AtomType.Add_Property({"LJtype":str})
```

##### Set_Property_Unit

```python
Set_Property_Unit(cls, prop, unit_type, base_unit)
```

This **function** is used to set the unit of the property of the class

###### Input

- `prop`: a **str**, the name of the property
- `unit_type`: a **str**, the name of the unit type
- `base_unit`: a **str**, the base unit of the property

###### Output

- `None`

###### Example

In `Xponge.forcefield.BASE.LJ`

```python
LJType.Set_Property_Unit("rmin", "distance", "A")
LJType.Set_Property_Unit("sigma", "distance", "A")
LJType.Set_Property_Unit("epsilon", "energy", "kcal/mol")
LJType.Set_Property_Unit("A", "energy·distance^6", "kcal/mol·A^6")
LJType.Set_Property_Unit("B", "energy·distance^12", "kcal/mol·A^12")
```

##### New_From_String

```python
New_From_String(cls, string, skip_lines=0)
```

This **function** is used to update the types of the class 

###### Input

- `string`: the **str** to update
- `skip_lines`: an **int**, the number of lines to skip

###### Output

- `None`

###### Example

In `Xponge.forcefield.AMBER.tip3p`

```python
AtomType.New_From_String(
"""
name mass    charge[e]  LJtype
HW   1.008    0.417       HW
OW   16      -0.834       OW
""")
```

##### New_From_File

```python
New_From_File(cls, filename, skip_lines=0)
```

This **function** is used to update the types of the class 

###### Input

- `filename`: the **str** of the file to update
- `skip_lines`: an **int**, the number of lines to skip

###### Output

- `None`

##### New_From_Dict

```python
New_From_String(cls, Dict)
```

This **function** is used to update the types of the class 

###### Input

- `Dict`: the **dict** to update

###### Output

- `None`

#### static methods

##### _unit_transfer

```python
_unit_transfer(self)
```

This **function** is used to convert units of the parameters.

#### instance methods

##### \_\_init\_\_

##### Update

```python
Update(self, **kwargs)
```

This **function** is used to update the properties of the instance

###### Input

- `kwargs`: the property - value pairs

###### Output

- `None`

### Xponge.AtomType

a sub**class** of Xponge.Type, for atom types

#### instance methods

##### \_\_init\_\_

 ### Xponge.ResidueType

a sub**class** of Xponge.Type, for residue types

#### instance variables

- `connectivity`: a **map** to store the atom - connected atom list pairs
- `builded`: a **bool**ean value to show whether the residue type is built or not
- `bonded_forces`: a **dict** to store the bonded force name - bonded force list pairs

- `atoms`: a **list** to store the atoms of the residue type
- `_name2atom`: a **dict** to store the name - Xponge.Atom pairs
- `_atom2name`: a **dict** to store the Xponge.Atom - name  pairs
- `_name2index`: a **dict** to store the name - index pairs
- `_atom2index`: a **dict** to store the Xponge.Atom - index pairs
- `connect_atoms`: a **dict** to store the distance - Xponge.Atom pairs
- `link`: a **dict** to store the information needed by linking

#### instance methods

##### \_\_init\_\_

##### Add_Atom

```python
Add_Atom(self, name, atom_type, x, y, z)
```

This **function** is used to add an atom to the residue type.

###### Input

- `name`: a **str**, the name of the atom to add
- `atom_type`: a **Xponge.AtomType**, the type of the atom to add
- `x`: a **float**, the position x of the atom to add
- `y`: a **float**, the position y of the atom to add
- `z`: a **float**, the position z of the atom to add

###### Output

- `None`

##### Add_Connectivity

```python
Add_Connectivity(self, atom0, atom1)
```

This **function** is used to add the connectivity between two atoms to the residue type.

###### Input

- `atom0`: the name of the atom to add
- `atom1`: the name of the atom to add

###### Output

- `None`

##### Add_Bonded_Force

```python
Add_Bonded_Force(self, bonded_force_entity, typename=None)
```

This **function** is used to add the bonded force to the residue type.

###### Input

- `bonded_force_entity`: a **Xponge.BondedForceEntity**, the bonded force entity to add
- `typename`: a **str**, the name of the bonded force entity. `None` to get the  type(bonded_force_entity).name

###### Output

- `None`

##### deepcopy

```python
deepcopy(self, name, forcopy=None)
```

This **function** is used to deep copy the instance

###### Input

- `name`: a **str**, the name for the new instance
- `forcopy`: a **str**, the string for temporary variables.  If forcopy == `None`, a hash string of the time will be used and the temporary variables will be deleted after copy.

###### Output

- `Xponge.ResidueType`: new copied residue type instance

### Xponge.Entity

This **class** is the abstract class of the entities (atoms, bonded forces, residues and so on).

#### class variables

- `count`: an **int** to record how many entities of this class there is
- `name`: a **str**, the name of the entity class

#### instance variables

- `contents`: a **dict**, which is the same as the one in Xponge.Type
- `count`: an **int** to record the entity instance index of the class
- `name`: a **str**, the name of the entity instance
- `type`: an Xponge.Type or a subclass instance of Xponge.Type, the type of the entity

#### instance methods

##### \_\_init\_\_

##### Update

```python
Update(self, **kwargs)
```

This **function** is used to update the properties of the instance

###### Input

- `kwargs`: the property - value pairs

###### Output

- `None`

### Xponge.Atom

a sub**class** of Xponge.Entity, for atoms

#### instance variables

- `residue`: an Xponge.Residue or an Xponge.ResidueType, the residue(type) which the atom belongs to
- `extra_excluded_atoms`: a **set** to store the extra excluded atoms for the atom
- `linked_atoms`: a **dict** to store distance - atom set pairs
- `copied`: a **dict** to store key - copied atom pairs

#### instance methods

##### \_\_init\_\_

##### deepcopy

```python
deepcopy(self, forcopy=None)
```

This **function** is used to deep copy the instance

###### Input

- `forcopy`: a **str**, the string for temporary variables.  If `forcopy` != `None`, the `forcopy` - new atom pair will be stored in the instance variable `copied`.

###### Output

- `Xponge.Atom`: new copied atom instance

##### Link_Atom

```python
Link_Atom(self, link_type, atom)
```

This **function** is used to link atoms for building

###### Input

- `link_type`:  an **int**, the distance between the `atom` atom and the instance atom
- `atom`: an Xponge.Atom, the atom to link

###### Output

- `None`

##### Extra_Exclude_Atom

```python
Extra_Exclude_Atom(self, atom)
```

This **function** is used to extra exclude one atom

###### Input

- `atom`: an Xponge.Atom, the atom to extra exclude

###### Output

- `None`

##### Extra_Exclude_Atoms

```python
Extra_Exclude_Atoms(self, lists)
```

This **function** is used to extra exclude a list of atoms

###### Input

- `lists`: a **list**, the atoms to extra exclude

###### Output

- `None`

### Xponge.Residue

a sub**class** of Xponge.Entity, for residues

#### instance variables

- `connectivity`: a **map** to store the atom - connected atom list pairs
- `builded`: a **bool**ean value to show whether the residue type is built or not
- `bonded_forces`: a **dict** to store the bonded force name - bonded force list pairs
- `atoms`: a **list** to store the atoms of the residue type
- `_name2atom`: a **dict** to store the name - Xponge.Atom pairs
- `_atom2name`: a **dict** to store the Xponge.Atom - name  pairs
- `_name2index`: a **dict** to store the name - index pairs
- `_atom2index`: a **dict** to store the Xponge.Atom - index pairs
- `connect_atoms`: a **dict** to store the distance - Xponge.Atom pairs

#### instance methods

##### \_\_init\_\_

##### Add_Atom

##### Add_Connectivity

##### Add_Bonded_Force

##### Add_Missing_Atoms

##### deepcopy

### Xponge.ResidueLink

a **class** for the link between residues

#### instance methods

##### \_\_init\_\_

##### Add_Bonded_Force

##### deepcopy

### Xponge.Molecule

a **class** for molecules

#### class methods

##### Set_Save_SPONGE_Input

##### Del_Save_SPONGE_Input

#### instance methods

##### \_\_init\_\_

##### Add_Residue

##### Add_Bonded_Force

##### Add_Residue_Link

##### Add_Missing_Atoms

##### deepcopy

### Xponge.Generate_New_Bonded_Force_Type

a **function** to generate the subclasses of the Xponge.Type and the Xponge.Entity for the bonded force

###  Xponge.Generate_New_Pairwise_Force_Type

a **function** to generate the subclasses of the Xponge.Type and the Xponge.Entity for the pairwise force

# Xponge.assign

the **package** to assign the properties for atoms, residues and molecules.

### \_\_init\_\_

the basic **module** for Xponge.assign

### Xponge.assign.Judge_Rule

the **class** of the rules to determinate the atom type for one atom

#### class variables

- `all`: a **dict** to store the rule name - rule pairs

####  instance methods

###### \_\_init\_\_

##### Add_Judge_Rule

### Xponge.assign._RING

the **class** of the rings

#### instance variables

- `atoms`
- `tohash`

#### instance methods

##### \_\_init\_\_

##### \_\_eq\_\_

##### get_3_neigbors

### Xponge.assign.Assign

the class to assign

#### class variables

#### instance variables

#### instance methods

### Xponge.assign.Guess_Element_From_Mass

### Xponge.assign.Get_Assignment_From_PubChem

this **function** is also stored in the main dict with the key `Get_Assignment_From_PubChem`

### Xponge.assign.Get_Assignment_From_PDB

this **function** is also stored in the main dict with the key `Get_Assignment_From_PDB`

### Xponge.assign.Get_Assignment_From_Mol2

this **function** is also stored in the main dict with the key `Get_Assignment_From_Mol2`

### Xponge.assign.Get_Assignment_From_ResidueType

this **function** is also stored in the main dict with the key `Get_Assignment_From_ResidueType`

## Xponge.assign.RESP

the **module** to calculate the RESP charge

## Xponge.assign.RDKit_tools

the **module** to use the RDKit tools

### Xponge.assign.RDKit_tools.Assign2RDKitMol

### Xponge.assign.RDKit_tools.Insert_Atom_Type_To_RDKitMol

### Xponge.assign.RDKit_tools.Find_Equal_Atoms

### Xponge.assign.RDKit_tools.Get_Conformer_Coordinate

### Xponge.assign.RDKit_tools.Get_Conformer_Coordinate_To_Residue

### Xponge.assign.RDKit_tools.Set_Conformer_Coordinate_From_Residue

### Xponge.assign.RDKit_tools.Set_Conformer_Coordinate

### Xponge.assign.RDKit_tools.Apply_Transform

### Xponge.assign.RDKit_tools.Get_Part_Align

# Xponge.BUILD

# Xponge.LOAD

# Xponge.IMPOSE

# Xponge.PROCESS

# Xponge.tools

# \_\_main\_\_

# Xponge.forcefield.BASE

# Xponge.forcefield.AMBER

# Xponge.forcefield.CHARMM27

# Xponge.forcefield.SPECIAL

# Xponge.Analysis

# XpongeLib

a **package** containing the C/C++ compiled library for Xponge.

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