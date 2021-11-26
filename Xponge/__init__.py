##########################################################################
# Copyright 2021 Gao's lab, Peking University, CCME. All rights reserved.
#
# NOTICE TO LICENSEE:
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##########################################################################

__version__ = "alpha-test"
__author__ = "Yijie Xia"

from collections import OrderedDict
import os
import numpy as np
import types
from itertools import product, permutations
import time

##########################################################################
#Basic Classes
##########################################################################
class _GlobalSetting():
    #最远的成键距离，用于拓扑分析时最远分析多远
    farthest_bonded_force = 0
    #所有的成键类型力的Type
    BondedForces = []
    #所有虚拟原子的Type和对应的依赖的其他原子的数量
    VirtualAtomTypes = {}
    #单位换算
    UnitMapping = {"distance": {"nm": 1e-9, "A": 1e-10},
                   "energy"  : {"kcal/mol":4.184, "eV":96.4853, "kJ/mol": 1},
                   "charge"  : {"e":1,   "SPONGE":1.0/18.2223},
                   "angle"   : {"degree":3.141592654,  "rad": 180}
                   
                    }
    PDBResidueNameMap = {"head" : {}, "tail": {}}
    HISMap = {"DeltaH": "", "EpsilonH": "",   "HIS": {}}
    @staticmethod
    def Set_Unit_Transfer_Function(sometype):
        def wrapper(func):
            sometype._unit_transfer = func
        return wrapper
    @staticmethod
    def Add_Unit_Transfer_Function(sometype):
        func0 = sometype._unit_transfer
        def wrapper(func):
            def temp(self):
                func0(self)
                func(self)
            sometype._unit_transfer = temp
        return wrapper

GlobalSetting = _GlobalSetting()

class Type():
    name = None
    parameters = {"name":str }
    types = {}
    types_different_name = {}
    
    @classmethod
    def Add_Property(cls, parm_fmt):
        cls.parameters.update(parm_fmt)
    
    @classmethod
    def Set_Property_Unit(cls, prop, unit_type, base_unit):
        assert prop in cls.parameters.keys(), "Unknown property '%s' for type '%s'"%(prop, cls.name)
        temp_unit_lists = unit_type.split("·")
        temp_unit_power = []
        
        for i in range(len(temp_unit_lists)):
            unit = temp_unit_lists[i].split("^")
            if len(unit) == 2:
                temp_unit_lists[i] = unit[0]
                temp_unit_power.append(int(unit[1]))
            else:
                temp_unit_power.append(1)
            assert unit[0] in GlobalSetting.UnitMapping.keys(), "Unknown unit type '%s'"%unit_type
        
        temp_unit_lists_units = [list(GlobalSetting.UnitMapping[unit].keys()) for unit in temp_unit_lists]
        alls = {}
        for unit_combination in product(*temp_unit_lists_units):
            unit = []
            value = 1
            for i in range(len(unit_combination)):
                power = temp_unit_power[i]
                value *= GlobalSetting.UnitMapping[temp_unit_lists[i]][unit_combination[i]] ** power
                if power == 1:
                    unit.append(unit_combination[i])
                else:
                    unit.append(unit_combination[i] + "^" + str(power))
            alls["·".join(unit)] = value
        base_unit_rate = alls[base_unit]
        
        
        def temp_func(current_rate, base_unit_rate):
            return lambda x: float(x) * current_rate / base_unit_rate
        
        prop_alls = {}
        for unit, current_rate in alls.items():
            temp_prop = prop + '[' + unit + ']'
            prop_alls[temp_prop] = temp_func(current_rate, base_unit_rate)
        cls.Add_Property(prop_alls)
            
    @staticmethod
    def _unit_transfer(self):
        for prop in self.contents.keys():
            if "[" in prop and "]" in prop and self.contents[prop] != None:
                real_prop = prop.split('[')[0]
                self.contents[real_prop] = self.contents[prop]
                self.contents[prop] = None
                
                
    @classmethod
    def New_From_String(cls, string, skip_lines = 0):
        count = -1
        kwargs = OrderedDict()
        for line in string.split("\n"):
            if not line.strip() or line.strip()[0] in ("!", "#", ";", "/") :
                continue

            count += 1
            if count < skip_lines:
                continue
            elif count == skip_lines:
                kwargs = kwargs.fromkeys(line.split())
            else:
                words = line.split()
                i = 0
                tempkw = {}.fromkeys(kwargs.keys())
                for key in tempkw.keys():
                    tempkw[key] = words[i]
                    i += 1
                type_already_have = False
                
                if tempkw["name"] in cls.types_different_name.keys():
                    tempkw["name"] = cls.types_different_name[tempkw["name"]].name
                    type_already_have = True
                if not type_already_have:
                    if "reset" in tempkw.keys():
                        tempkw.pop("reset")
                    cls(**tempkw)
                else:
                    temp = cls.types[tempkw.pop("name")]
                    temp.update(**tempkw)

    @classmethod
    def New_From_File(cls, filename, skip_lines = 0):
        with open(filename) as f:
            New_From_String(f.read(), skip_lines)

    def __repr__(self):
        return "Type of " + type(self).name + ": " + self.name
        
    def __hash__(self):
        return hash(repr(self))
    
    def __getattribute__(self, attr):
        if attr != "contents" and attr in self.contents.keys():
            return self.contents[attr]
        else:
            return super().__getattribute__(attr)
    
    def __setattr__(self, attr, value):
        if attr != "contents" and attr in self.contents.keys():
            self.contents[attr] = value
        else:
            return super().__setattr__(attr, value)
    
    def __init__(self, **kwargs):
        prop_fmt = type(self).parameters
        
        self.contents = {}.fromkeys(prop_fmt.keys())
        self.name = kwargs.pop("name")
        assert self.name not in type(self).types.keys(), "The name '%s' has already existed in '%sType'"%(self.name, type(self).name)
        type(self).types[self.name] = self 
        type(self).types_different_name[self.name] = self 
        for key in kwargs.keys():
            assert key in self.contents.keys(), "The parameter '%s' is not one of the parameters of '%sType'"%(key, type(self).name)
            self.contents[key] = prop_fmt[key](kwargs[key])
        type(self)._unit_transfer(self)
        
    def update(self, **kwargs):
        for key in kwargs.keys():
            assert(key in self.contents.keys())
            self.contents[key] = type(self).parameters[key](kwargs[key])
        type(self)._unit_transfer(self)
    

            


class AtomType(Type):
    name =  "Atom"
    parameters = {"name":str,  "x":float, "y":float, "z":float}
    types = {}
    types_different_name = {}


class ResidueType(Type):
    name = "Residue"
    parameters = {"name":str,}
    types = {}
    types_different_name = {}
    
    @property
    def head(self):
        return self.connect_atoms["head"]
    
    @head.setter
    def head(self, atom):
        self.connect_atoms["head"] = atom

    @property
    def tail(self):
        return self.connect_atoms["tail"]
    
    @tail.setter
    def tail(self, atom):
        self.connect_atoms["tail"] = atom
        
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.atoms = []
        self._name2atom = {}
        self._atom2name = {}
        self.connectivity = {}
        self.bonded_forces = {}
        self.builded = False
        self.link = {}
        self.bonded_forces = {frc.name:[] for frc in GlobalSetting.BondedForces}
        self.connect_atoms = {"head": None, "tail":None}
        
    def Add_Atom(self, name, atom_type, x, y, z):
        new_atom = Atom(atom_type, name)
        self.atoms.append(new_atom)
        new_atom.x = float(x)
        new_atom.y = float(y)
        new_atom.z = float(z)
        self._name2atom[name] = new_atom
        self._atom2name[new_atom] = name
        self.connectivity[new_atom] = set([])

    def Add_Connectivity(self, atom0, atom1):
        if type(atom0) == str:
            atom0 = self._name2atom[atom0]
        if type(atom1) == str:
            atom1 = self._name2atom[atom1]
        self.connectivity[atom0].add(atom1)
        self.connectivity[atom1].add(atom0)
    
    def Add_Bonded_Force(self, bonded_force_entity):
        if type(bonded_force_entity).name not in self.bonded_forces.keys():
            self.bonded_forces[type(bonded_force_entity).name] = []
        self.bonded_forces[type(bonded_force_entity).name].append(bonded_force_entity)
    
    def Add_Connect_Atom(self, name, atom):
        self.connect_atoms[name] = atom


class Entity():
    count = 0
    name = None
    def __repr__(self):
        return "Entity of " + type(self).name + ": " + self.name  + "(" +str(self.count) + ")"
        
    def __hash__(self):
        return hash(repr(self))
        
    def __getattribute__(self, attr):
        if attr != "contents" and attr in self.contents.keys():
            return self.contents[attr]
        else:
            return super().__getattribute__(attr)
    
    def __setattr__(self, attr, value):
        if attr != "contents" and attr in self.contents.keys():
            self.contents[attr] = value
        else:
            return super().__setattr__(attr, value)
            
    def __init__(self, entity_type, name = None):
        self.contents = {**entity_type.contents}
        self.count = type(self).count
        if not name:
            name = entity_type.name
        type(self).count += 1
        self.name = name
        self.type = entity_type
    
    def update(self, **kwargs):
        for key in kwargs.keys():
            assert(key in self.contents.keys())
            self.contents[key] = type(self.type).parameters[key](kwargs[key])
        type(self.type)._unit_transfer(self)
   
        
        
        
class Atom(Entity):
    name = "Atom"
    count = 0
    def __init__(self, entity_type, name = None):
        super().__init__(entity_type, name)
        self.linked_atoms = { i+1 : [] for i in range(1,GlobalSetting.farthest_bonded_force)}
        self.residue = None
        self.copied = {}
    def deepcopy(self, forlink = None):
        new_atom = Atom(self.entity)
        if forlink:
            self.copied[forlink] = new_atom
        return new_atom


class Residue(Entity):
    name = "Residue"
    count = 0
    def __init__(self, entity_type, name = None):
        super().__init__(entity_type, name)
        self.atoms = []
        self._name2atom = {}
        self._atom2name = {}
        self.connectivity = {}
        self.bonded_forces = {frc.name:[] for frc in GlobalSetting.BondedForces}
        self.builded = False
    
    def Add_Atom(self, name, atom_type = None, x = None, y = None, z = None):
        if not atom_type: 
            atom_type = self.type._name2atom[name].type
            new_atom = Atom(atom_type, name)
            new_atom.contents = self.type._name2atom[name].contents
        else:
            new_atom = Atom(atom_type, name)
            new_atom.contents = self.type._name2atom[name].contents
        new_atom.residue = self
        self.atoms.append(new_atom)
        if x:
            new_atom.x = float(x)
        if y:
            new_atom.y = float(y)
        if z:
            new_atom.z = float(z)
        self._name2atom[name] = new_atom
        self._atom2name[new_atom] = name
        self.connectivity[new_atom] = set([])
        
    def Add_Connectivity(self, atom0, atom1):
        if type(atom0) == str:
            atom0 = self._name2atom[atom0]
        if type(atom1) == str:
            atom1 = self._name2atom[atom1]
        if atom0 in self.connectivity.keys():
            self.connectivity[atom0].add(atom1)
        if atom1 in self.connectivity.keys():
            self.connectivity[atom1].add(atom0)

    def Add_Bonded_Force(self, bonded_force_entity):
        if type(bonded_force_entity).name not in self.bonded_forces.keys():
            self.bonded_forces[type(bonded_force_entity).name] = []
        self.bonded_forces[type(bonded_force_entity).name].append(bonded_force_entity)
    
    def deepcopy(self, forlink = None):
        new_residue = Residue(self.type)
        for atom in self.atoms:
            new_residue.Add_Atom(atom.name)
            if forlink:
                atom.copied[forlink] = new_residue.atoms[-1]
        return new_residue

class ResidueLink():
    def __repr__(self):
        return "Entity of ResidueLink: " + repr(self.atom1) + "-" + repr(self.atom2)
    def __hash__(self):
        return hash(repr(self))
    def __init__(self, atom1, atom2):
        self.atom1 = atom1
        self.atom2 = atom2
        self.builded = False
        self.bonded_forces = {frc.name:[] for frc in GlobalSetting.BondedForces}
    def Add_Bonded_Force(self, bonded_force_entity):
        if type(bonded_force_entity).name not in self.bonded_forces.keys():
            self.bonded_forces[type(bonded_force_entity).name] = []
        self.bonded_forces[type(bonded_force_entity).name].append(bonded_force_entity)
    def deepcopy(self, forlink):
        if self.atom1.copied[forlink] and self.atom2.copied[forlink]:
            return ResidueLink(self.atom1.copied[forlink], self.atom2.copied[forlink])
        
class Molecule():
    all = {}
    save_functions = []
    
    @classmethod
    def Set_Save_SPONGE_Input(cls, func):
        cls.save_functions.append(func)
        
    def __repr__(self):
        return "Entity of Molecule: " + self.name
    def __init__(self, name):
        self.name = name
        Molecule.all[name] = self
        self.residues = []
        self.atoms = []
        self.residue_links = []
        self.bonded_forces = []
        self.builded = False
        
        
    def Add_Residue(self, residue):
        self.residues.append(residue)
        
    def Add_Bonded_Force(self, bonded_force_entity):
        if type(bonded_force_entity).name not in self.bonded_forces.keys():
            self.bonded_forces[type(bonded_force_entity).name] = []
        self.bonded_forces[type(bonded_force_entity).name].append(bonded_force_entity)
    
    def Add_Residue_Link(self, atom1, atom2):
        self.residue_links.append(ResidueLink(atom1, atom2))

    def deepcopy(self):
        new_molecule = Molecule(self.name)
        forlink = hash(str(time.time()))
        for res in self.residues:
            new_molecule.Add_Residue(res.deepcopy(forlink))
        for link in self.residue_links:
            new_molecule.residue_links.append(link.deepcopy(forlink))
        for res in self.residues:
            for atom in res.atoms:
                atom.copied.pop(forlink)
        return new_molecule
        
@Molecule.Set_Save_SPONGE_Input     
def write_residue(self, prefix, dirname):
    towrite = "%d %d\n"%(len(self.atoms), len(self.residues))
    towrite += "\n".join([str(len(res.atoms)) for res in self.residues])
    f = open(os.path.join(dirname, prefix + "_residue.txt"),"w")
    f.write(towrite)
    f.close()

@Molecule.Set_Save_SPONGE_Input     
def write_coordinate(self, prefix, dirname):
    towrite = "%d\n"%(len(self.atoms))
    boxlength = [0, 0, 0, 90, 90, 90]    
    for atom in self.atoms:
        towrite += "%f %f %f\n"%(atom.x, atom.y, atom.z)
        if atom.x > boxlength[0]:
            boxlength[0] = atom.x
        if atom.y > boxlength[1]:
            boxlength[1] = atom.y
        if atom.z > boxlength[2]:
            boxlength[2] = atom.z
    towrite += "\n".join(["%f %f %f"%(atom.x, atom.y, atom.z) for atom in self.atoms])
    boxlength[0] += 3
    boxlength[1] += 3
    boxlength[2] += 3
    towrite += "\n%f %f %f %f %f %f"%(boxlength[0], boxlength[1], boxlength[2], boxlength[3], boxlength[4], boxlength[5])
    f = open(os.path.join(dirname, prefix + "_coordinate.txt"),"w")
    f.write(towrite)
    f.close()

def Generate_New_Bonded_Force_Type(Type_Name, atoms, properties, Compulsory, Multiple = None):
    
    class BondedForceEntity(Entity):
        name = Type_Name
        count = 0
        def __init__(self, atoms, entity_type, name = None):
            super().__init__(entity_type, name)
            self.atoms = atoms
            
    temp = [int(i) for i in atoms.split("-")]
    class BondedForceType(Type):
        name = Type_Name
        topology_like = temp
        compulsory = Compulsory
        multiple = Multiple
        atom_numbers = len(atoms.split("-"))
        topology_matrix = [[temp[i]-j if i > j else 1 for i in range(len(atoms.split("-")))] for j in range(len(atoms.split("-")))]
        parameters = {
            "name":str,
        }
        entity = BondedForceEntity
        types = {}
        types_different_name = {}
            
        @classmethod
        def Same_Force(cls, atom_list):
            if type(atom_list) == str:
                atom_list_temp = [ atom.strip() for atom in atom_list.split("-")]
                temp = [atom_list, "-".join(atom_list_temp[::-1])]
            else:
                temp = [atom_list, atom_list[::-1]]
            return temp
            
        @classmethod
        def Set_Same_Force_Function(cls, func):
            cls.Same_Force = classmethod(func)
        
        def __init__(self, **kwargs):
            super().__init__(**kwargs)
            
            for name in type(self).Same_Force(self.name):
                type(self).types_different_name[name] = self
                    
            if type(self).multiple:
                for key in self.multiple:
                    self.contents[key + "s"] = [self.contents[key]]
                    self.contents[key] = None
                self.contents["multiple_numbers"] = 1
        
        def update(self, **kwargs):
            reset = 1
            if "reset" in kwargs.keys():
                reset = int(kwargs.pop("reset"))
            if reset:
                for key in self.contents.keys():
                    if key != "name":
                        self.contents[key] = None
                if type(self).multiple:
                    for key in self.multiple:
                        self.contents[key + "s"] = []
                    self.contents["multiple_numbers"] = 0 
            super().update(**kwargs)
            if type(self).multiple:
                for key in type(self).multiple:
                    self.contents[key + "s"].append(self.contents[key])
                    self.contents[key] = None
                self.multiple_numbers += 1
        
    BondedForceType.Add_Property(properties)
    
    global GlobalSetting
    GlobalSetting.BondedForces.append(BondedForceType)
    for i in atoms.split("-"):
        if int(i) > GlobalSetting.farthest_bonded_force:
            GlobalSetting.farthest_bonded_force = int(i)
    return BondedForceType

def Generate_New_Pairwise_Force_Type(Type_Name, properties):
    class PairwiseForceType(Type):
        name = Type_Name
        parameters = {
            "name":str,
        }
        types = {}
        
    PairwiseForceType.Add_Property(properties)
    
    return PairwiseForceType

from . import LOAD, BUILD

def ResidueType_Add(self, other):
    if type(other) == ResidueType:
        new_molecule = Molecule(self.name)
        resA = Residue(self)
        resB = Residue(other)
        for atom in self.atoms:
            resA.Add_Atom(atom.name)
        for atom in other.atoms:
            resB.Add_Atom(atom.name)
        new_molecule.Add_Residue(resA)
        new_molecule.Add_Residue(resB)
        if self.tail and other.head:
            new_molecule.Add_Residue_Link(resA._name2atom[self.tail], resB._name2atom[other.head])
        return new_molecule
    elif type(other) == Molecule:
        new_molecule = other.deepcopy()
        resA = Residue(self)
        for atom in self.atoms:
            resA.Add_Atom(atom.name)
        if self.tail and other.residues[0].type.head:
            new_molecule.Add_Residue_Link(resA._name2atom[self.tail], other.residues[0]._name2atom[other.residues[0].type.head])
        new_molecule.residues.insert(0, resA)
        return new_molecule
    elif type(other) == type(None):
        return self
    else:
        raise TypeError("unsupported operand type(s) for +: '%s' and '%s'"%(type(self), type(other)))

def Molecule_Add(self, other):
    if type(other) == ResidueType:
        new_molecule = self.deepcopy()
        resB = Residue(other)
        for atom in other.atoms:
            resB.Add_Atom(atom.name)
        new_molecule.Add_Residue(resB)
        if new_molecule.residues[-1].type.tail and other.head:
            new_molecule.Add_Residue_Link(new_molecule.residues[-1]._name2atom[new_molecule.residues[-1].type.tail], resB._name2atom[other.head])
        return new_molecule
    elif type(other) == Molecule:
        new_molecule = self.deepcopy()
        new_molecule2 = other.deepcopy()
        if new_molecule.residues[-1].type.tail and new_molecule2.residues[0].type.head:
            new_molecule.Add_Residue_Link(new_molecule.residues[-1]._name2atom[new_molecule.residues[-1].type.tail],
            new_molecule2.residues[0]._name2atom[new_molecule2.residues[0].type.head])
        for res in new_molecule2.residues:
            new_molecule.Add_Residue(res)
        return new_molecule
    elif type(other) == type(None):
        return self
    else:
        raise TypeError("unsupported operand type(s) for +: '%s' and '%s'"%(type(self), type(other)))

def iMolecule_Add(self, other):
    if type(other) == ResidueType:
        resB = Residue(other)
        for atom in other.atoms:
            resB.Add_Atom(atom.name)
        self.Add_Residue(resB)
        if self.residues[-1].type.tail and other.head:
            self.Add_Residue_Link(self.residues[-1]._name2atom[self.residues[-1].type.tail], resB._name2atom[other.head])
        return self
    elif type(other) == Molecule:
        new_molecule2 = other.deepcopy()
        if self.residues[-1].type.tail and new_molecule2.residues[0].type.head:
            self.Add_Residue_Link(self.residues[-1]._name2atom[self.residues[-1].type.tail],
            new_molecule2.residues[0]._name2atom[new_molecule2.residues[0].type.head])
        for res in new_molecule2.residues:
            self.Add_Residue(res)
        return self
    elif type(other) == type(None):
        return self
    else:
        raise TypeError("unsupported operand type(s) for +: '%s' and '%s'"%(type(self), type(other)))

def Muls(self, other):
    if type(other) == int:
        assert other >= 0
        t = self
        for i in range(other - 1):
            t = t + self
        return t
    else:
        raise TypeError("unsupported operand type(s) for +: '%s' and '%s'"%(type(self), type(other)))

def iMuls(self, other):
    if type(other) == int:
        assert other >= 0
        for i in range(other - 1):
            self += self
        return self
    else:
        raise TypeError("unsupported operand type(s) for +: '%s' and '%s'"%(type(self), type(other)))

ResidueType.__add__ = ResidueType_Add
ResidueType.__radd__ = ResidueType_Add
ResidueType.__mul__ = Muls
ResidueType.__rmul__ = Muls
Molecule.__add__ = Molecule_Add
Molecule.__radd__ = Molecule_Add
Molecule.__iadd__ = iMolecule_Add
Molecule.__mul__ = Muls
Molecule.__rmul__ = Muls
Molecule.__imul__ = iMuls

del ResidueType_Add
del Molecule_Add
del Muls
del iMuls
del iMolecule_Add