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

```python
__init__(self, **kwargs)
```

###### Input

- `kwargs`: the property - value pairs

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

```python
__init__(self, **kwargs)
```

###### Input

- `kwargs`: the property - value pairs

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

```python
__init__(self, **kwargs)
```

###### Input

- `kwargs`: the property - value pairs

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

- `atom0`: a **str** or an Xponge.Atom instance, the name of the atom to add
- `atom1`: a **str** or an Xponge.Atom instance, the name of the atom to add

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

```python
__init__(self, entity_type, name=None)
```

###### Input

- `entity_type`: an instance of the (sub)class Xponge.Type, as the type of the entity instance
- `name`: a **str**, the name of the instance

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

```python
__init__(self, entity_type, name=None)
```

###### Input

- `entity_type`: an instance of the (sub)class Xponge.Type, as the type of the entity instance
- `name`: a **str**, the name of the instance

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

```python
__init__(self, entity_type, name=None, directly_copy=False)
```

###### Input

- `entity_type`: an instance of the (sub)class Xponge.Type, as the type of the entity instance
- `name`: a **str**, the name of the instance
- `directly_copy`: a **bool**ean value. If True, the instance will be initialized by directly copying from the `entity_type`

##### Add_Atom

```python
Add_Atom(self, name, atom_type=None, x=None, y=None, z=None)
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

This **function** is used to add the connectivity between two atoms to the residue entity.

###### Input

- `atom0`: a **str** or an Xponge.Atom instance, the name of the atom to add
- `atom1`: a **str** or an Xponge.Atom instance, the name of the atom to add

###### Output

- `None`

##### Add_Bonded_Force

```python
Add_Bonded_Force(self, bonded_force_entity)
```

This **function** is used to add the bonded force to the residue entity.

###### Input

- `bonded_force_entity`: a **Xponge.BondedForceEntity**, the bonded force entity to add

###### Output

- `None`

##### Add_Missing_Atoms

```python
Add_Missing_Atoms(self)
```

This **function** is used to add the missing atoms from the Xponge.ResidueType to the residue entity. 

##### deepcopy

```python
deepcopy(self, forcopy=None)
```

This **function** is used to deep copy the instance

###### Input

- `forcopy`: a **str**, the string for temporary variables.  If forcopy != `None`, the new atoms will be stored in the old atoms' attribute `copied`.

###### Output

- `Xponge.Residue`: new copied residue instance

### Xponge.ResidueLink

a **class** for the link between residues

#### instance methods

##### \_\_init\_\_

```python
__init__(self, atom1, atom2)
```

###### Input

- `atom1`: an Xponge.Atom instance, the name of the atom to add
- `atom2`: an Xponge.Atom instance, the name of the atom to add

##### Add_Bonded_Force

```python
Add_Bonded_Force(self, bonded_force_entity)
```

This **function** is used to add the bonded force to the residue link

###### Input

- `bonded_force_entity`: a **Xponge.BondedForceEntity**, the bonded force entity to add

###### Output

- `None`

##### deepcopy

```python
deepcopy(self, forcopy=None)
```

This **function** is used to deep copy the instance

###### Input

- `forcopy`: a **str**, the string for temporary variables.  If forcopy == `None`, a hash string of the time will be used.

###### Output

- `Xponge.ResidueLink`: new copied residue link instance

### Xponge.Molecule

a **class** for molecules

#### class methods

##### Set_Save_SPONGE_Input

```python
Set_Save_SPONGE_Input(cls, keyname)
```

This **function** is used to set the function when `Save_SPONGE_Input`. It is recommended used as a **decorator**. 

###### Input

- `keyname`: a **str**, which will be saved as "xxx\_`keyname`.txt" when saving SPONGE inputs.
- `somefunction`: the **function** that receives a Xponge.Molecule as the input and gives a str as the output.

###### Example

In `Xponge.__init__`

```python
@Molecule.Set_Save_SPONGE_Input("residue")
def write_residue(self):
    towrite = "%d %d\n" % (len(self.atoms), len(self.residues))
    towrite += "\n".join([str(len(res.atoms)) for res in self.residues])
    return towrite
```

##### Del_Save_SPONGE_Input

```python
Del_Save_SPONGE_Input(cls, keyname)
```

This **function** is used to delete the function when `Save_SPONGE_Input`. 

###### Input

- `keyname`: a **str**, which will be saved as "xxx\_`keyname`.txt" when saving SPONGE inputs.

#### instance methods

##### \_\_init\_\_

```python
__init__(self, name)
```

###### Input

- `name`: the name of the molecule

##### Add_Residue

```python
Add_Residue(self, residue)
```

This **function** is used to add a residue to the molecule

###### Input

- `residue`: a Xponge.Residue or Xponge.ResidueType instance, the residue to add

##### Add_Bonded_Force

```python
Add_Bonded_Force(self, bonded_force_entity)
```

This **function** is used to add the bonded force to the residue link

###### Input

- `bonded_force_entity`: a **Xponge.BondedForceEntity**, the bonded force entity to add

###### Output

- `None`

##### Add_Residue_Link

```python
Add_Residue_Link(self, atom1, atom2)
```

This **function** is used to add the connectivity between two atoms of two residues in the molecule.

###### Input

- `atom1`: a **str** or an Xponge.Atom instance, the name of the atom to add
- `atom2`: a **str** or an Xponge.Atom instance, the name of the atom to add

###### Output

- `None`

##### Add_Missing_Atoms

```python
Add_Missing_Atoms(self)
```

This **function** is used to add the missing atoms from the Xponge.ResidueType instances to the residues in the molecule. 

##### deepcopy

```python
deepcopy(self)
```

This **function** is used to deep copy the instance

### Xponge.Generate_New_Bonded_Force_Type

```python
Generate_New_Bonded_Force_Type(Type_Name, atoms, properties, Compulsory, Multiple=None) 
```

a **function** to generate the subclasses of the Xponge.Type and the Xponge.Entity for the bonded force

###### Input

- Type_Name: a **str**, the name of the bonded force type
- atoms: a **str**, the string to define the atom pattern of the bonded force type
- properties: a **dict**, the property name - type pairs
- Compulsory: a **bool**ean value to show whether every atom list which meet the atom pattern should be parametered.
- Multiple: a **list** to show the properties that could be multiple for one type.

###### Example

In `Xponge.forcefield.BASE.DIHEDRAL`

```python
ProperType = Generate_New_Bonded_Force_Type("dihedral", "1-2-3-4", {"k":float, "phi0": float, "periodicity":int}, True, ["k", "phi0", "periodicity"])
ImproperType = Generate_New_Bonded_Force_Type("improper", "1-3-2-3", {"k":float, "phi0": float, "periodicity":int}, False)
```

###  Xponge.Generate_New_Pairwise_Force_Type

```python
Generate_New_Pairwise_Force_Type(Type_Name, properties)
```

a **function** to generate the subclasses of the Xponge.Type and the Xponge.Entity for the pairwise force

###### Input

- Type_Name: a **str**, the name of the bonded force type
- properties: a **dict**, the property name - type pairs

###### Example

In `Xponge.forcefield.BASE.LJ`

```python
LJType = Generate_New_Pairwise_Force_Type("LJ", {"epsilon": float, "rmin": float, "sigma":float, "A":float, "B":float})
```

# Xponge.assign

the **package** to assign the properties for atoms, residues and molecules.

## \_\_init\_\_

the basic **module** for Xponge.assign

### Xponge.assign.Judge_Rule

the **class** of the rules to determine the atom type for one atom

#### class variables

- `all`: a **dict** to store the rule name - rule pairs

#### instance variables

- `rules`: an **OrderedDict** to store the atom type - judge function pairs 

####  instance methods

##### \_\_init\_\_

```python
__init__(self, name)
```

###### Input

- `name`: the name of the instance

###### Example

In `Xponge.forcefield.AMBER.gaff`

```python
GAFF = assign.Judge_Rule("GAFF")
```



##### Add_Judge_Rule

```python
Add_Judge_Rule(self, atomtype)
```

This **function** is used as an **decorator** to add the atom type - judge function

###### Input

- `atomtype`: a **str** or a Xponge.AtomType, the atom type to judge
- `somefunction`: the judge function that receives an index of the atom to judge and an Xponge.assign.Assign as the input and gives a boolean value as the output

###### Example

In `Xponge.forcefield.AMBER.gaff`

```python
@GAFF.Add_Judge_Rule("cx")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "C4") and "RG3" in Assign.atom_marker[i].keys()
```

### Xponge.assign._RING

the **class** of the rings

### Xponge.assign.Assign

the class to assign properties for atoms, which is called an "assignment"

#### class variables

- `XX`: set("CNOPS")
- `XA`: set("OS")
- `XB`: set("NP")
- `XC`: set(["F", "Cl", "Br", "I"])
- `XD`: set("SP")
- `XE`: set(["N", "O", "F", "Cl", "Br", "S", "I"])

#### instance variables

- `name`: a **str**, the name of the assignment
- `atom_numbers`: an **int**, the number of the atoms
- `atoms`: a **list**, the element names of the atoms
- `element_details`: a **list**, some detailed information on element of the atoms
- `coordinate`: a **numpy.array** with a shape of (`atom_numbers`,3), the coordinates of the atoms
- `charge`: a **numpy.array** with a shape of (`atom_numbers`, ), the partial charges of the atoms
- `atom_types`: a **dict**, the index of the atom - atom type pairs
- `atom_marker`: a **dict**, the index of the atom - atom marker **dict** (marker name - amount pairs) pairs
- `bonds`: a **dict**, the index of the first atom - **dict** (the index of the second atom - bond order pairs) pairs
- `ar_bonds`: a **dict**, the index of the first atom - the index of the second atoms **list** pairs for the aromatic bonds
- `am_bonds`: a **dict**, the index of the first atom - the index of the second atoms **list** pairs for the amide bonds

#### instance methods

##### \_\_init\_\_

```python
__init__(self, name = "ASN")
```

###### Input

- `name`: the name of the assignment

##### Atom_Judge

```python
Atom_Judge(self, atom, string)
```

this **function** is used to judge atom whether meet the pattern

###### Input

- `atom`: the index of the atom
- `string`: the string pattern

###### Output

- `True` or `False`

###### Example

In `Xponge.forcefield.AMBER.gaff`

```python
@GAFF.Add_Judge_Rule("c3")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "C4")
```

##### Add_Atom

```python
Add_Atom(self, element, x, y, z, name="", charge = 0.0)
```

this **function** is used to add an atom to the assignment

###### Input

- `element`: a **str**, the element name of the atom
- `x`: the x coordinate of the atom
- `y`: the y coordinate of the atom
- `z`: the z coordinate of the atom
- `name`: the name of the atom
- `charge`: the charge of the atom, in unit of e

##### Add_Atom_Marker

```python
Add_Atom_Marker(self, atom, marker)
```

This **function** is used to add a marker to the atom

###### Input

- `atom`: an **int**, the index of the atom
- `marker`: a `str`, the marker

##### Add_Bond

```python
Add_Bond(self, atom1, atom2, order=-1)
```

This **function** is used to add a bond to the assignment

###### Input

- `atom1`: an **int**, the index of the first atom
- `atom2`: an **int**, the index of the second atom
- `order`: an **int**, the bond order. `-1` for the unknown bond order

##### Add_Bond_Marker

```python
Add_Bond_Marker(self, atom1, atom2, marker, only1=False)
```

This **function** is used to add a marker to the bond

###### Input

- `atom1`: an **int**, the index of the first atom
- `atom2`: an **int**, the index of the second atom
- `marker`: a `str`, the marker
- `only1`: a **bool**ean value, if `True`, only the atom1 - atom2 bond will get the marker, while the atom2 - atom1 not

##### Determine_Equal_Atoms

```python
Determine_Equal_Atoms(self)
```

This **function** is used to get the equal atoms in the assignment

###### Output

- `list`: a list of sub lists, every sub list consist of equal atoms

###### Example

```python
import Xponge
t = Get_Assignment_From_PubChem("methane", "name")
print(t.Determine_Equal_Atoms())
#[[0],[1,2,3,4]]
```

##### Determine_Ring_And_Bond_Type

```python
Determine_Ring_And_Bond_Type(self)
```

This **function** is used to determine the rings, atom markers, bond types and bond markers in the assignment

##### Determine_Atom_Type

```python
Determine_Atom_Type(self, rule)
```

This **function** is used to determine the atom type from the `rule`

###### Input

- `rule`: a **str**, the rule to judge

##### To_ResidueType

```python
To_ResidueType(self, name, charge=None)
```

This **function** is used to get a Xponge.ResidueType from the assignment

###### Input

- `name`: a **str**, the name of the residue type
- `charge`: a `numpy.array` or `None`, the partial charges of the atoms. If `None`, `self.charge` will be used.

##### Calculate_Charge

```python
Calculate_Charge(self, method, **parameters)
```

This **function** is used to calculate the charges of the atoms

###### Input

- `method`: a **str**, the method to calculate the partial charges. "RESP" only now.
- `parameters`: the parameters to calculate the partial charges. See Xponge.assign.RESP.RESP_Fit for RESP parameters

##### Save_As_PDB

```python
Save_As_PDB(self, filename)
```

This **function** is used to save the assignment to a pdb file

###### Input

- `filename`: a **str**, the filename

##### Save_As_Mol2

```python
Save_As_Mol2(self, filename)
```

This **function** is used to save the assignment to a mol2 file

###### Input

- `filename`: a **str**, the filename

### Xponge.assign.Guess_Element_From_Mass

```python
Guess_Element_From_Mass(mass)
```

This **function** is used to guess element from mass

###### Input

- `mass`: a **float**, the mass

###### Output

- `str`: the name of the element

### Xponge.assign.Get_Assignment_From_PubChem

This **function** is also stored in the main dict with the key `Get_Assignment_From_PubChem`

```python
Get_Assignment_From_PubChem(parameter, keyword)
```

###### Input

- `parameter`: a **str**, the parameter
- `keyword`: a **str**, the type of the parameter

###### Output

- `Xponge.assign.Assign`: the assignment

###### Example

```python
import Xponge
methane = Get_Assignment_From_PubChem("methane", "name")
ethane = Get_Assignment_From_PubChem("CC", "smiles")
```

### Xponge.assign.Get_Assignment_From_PDB

this **function** is also stored in the main dict with the key `Get_Assignment_From_PDB`

```python
Get_Assignment_From_PDB(filename, determine_bond_order = True)
```

###### Input 

- `filename`: a **str**, the name of the file
- `determine_bond_order`: a **bool**ean value, whether determine the bond order after reading the pdb file

###### Output

- `Xponge.assign.Assign`: the assignment

### Xponge.assign.Get_Assignment_From_Mol2

this **function** is also stored in the main dict with the key `Get_Assignment_From_Mol2`

```python
Get_Assignment_From_Mol2(filename)
```

### Xponge.assign.Get_Assignment_From_ResidueType

this **function** is also stored in the main dict with the key `Get_Assignment_From_ResidueType`

```python
Get_Assignment_From_ResidueType(restype)
```

## Xponge.assign.RESP

the **module** to calculate the RESP charge

### Xponge.assign.RESP.RESP_Fit

```python
RESP_Fit(Assign, basis = "6-31g*", opt = False, opt_params = None, charge = 0, spin = 0, extra_equivalence = [], grid_density = 6, grid_cell_layer = 4, 
    radius = None, a1 = 0.0005, a2 = 0.001, two_stage = True, only_ESP  = False)
```

this **function** is uesd to calculate the RESP charge

###### Input

- `Assign`: an Xponge.assign.Assign, the assignment

- `basis`: a **str**, the basis for quantum quantum mechanics calculation

- `opt`: a **bool**eam value, whether do optimization before the charge fitting

- `charge`: an `int`, the total charge of the assignment

- `spin`: an `int`, the number of unpaired electrons *2S*, i.e. the difference between the number of alpha and beta electrons.

- `extra_equivalence`: a **list** of **list**s, the extra equal atoms besides the hydrogens in -CH3 and =CH2 groups, which is called the "extra" equivalence because the origin reference of the RESP charge only mentioned the equivalence of the hydrogens in -CH3 and =CH2 groups

- `grid_density`: a **float**, the grid density in unit of grid/$\AA^2$

- `grid_cell_layer`: an **int**, the number of grid cell layers

- `radius`: a **dict**, element - vdw radius pairs. 

  ```python
  default_radius = {"H": 1.2, "C":1.5, "N":1.5, 
                    "O": 1.4, "P":1.8, "S":1.75,
                    "F": 1.35,"Cl":1.7,"Br":2.3}
  ```
  
- `a1`: a **float**, the first stage restraint weight

- `a2`: a **float**, the second stage restraint weight

- `two_stage` : a **bool**ean value, whether to do the second stage of fitting 

- `only_ESP`: a **bool**ean value, whether to only do the ESP fitting instead of the RESP fitting
###### Output

`numpy.array`: the partial charges of the atoms

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