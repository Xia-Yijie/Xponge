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

__version__ = "beta-test"
__author__ = "Yijie Xia"

from collections import OrderedDict, deque
import os
import numpy as np
import types
from itertools import product, permutations
import time
import sys

##########################################################################
#Basic Classes
##########################################################################
class _GlobalSetting():
    #最远的成键距离，用于拓扑分析时最远分析多远
    farthest_bonded_force = 0
    #所有的成键类型力的Type
    BondedForces = []
    BondedForcesMap = {}
    #所有虚拟原子的Type和对应的依赖的其他原子的数量
    VirtualAtomTypes = {}
    #单位换算
    UnitMapping = {"distance": {"nm": 1e-9, "A": 1e-10},
                   "energy"  : {"kcal/mol":4.184, "eV":96.4853, "kJ/mol": 1},
                   "charge"  : {"e":1,   "SPONGE":1.0/18.2223},
                   "angle"   : {"degree":np.pi,  "rad": 180}
                  }
    PDBResidueNameMap = {"head" : {}, "tail": {}, "save":{}}
    HISMap = {"DeltaH": "", "EpsilonH": "",   "HIS": {}}
    
    def Add_PDB_Residue_Name_Mapping(self, place, pdb_name, real_name):
        assert place in ("head", "tail")
        self.PDBResidueNameMap[place][pdb_name] = real_name
        self.PDBResidueNameMap["save"][real_name] = pdb_name

    def Set_Invisible_Bonded_Forces(self, types):
        for typename in types:
            self.BondedForces.remove(self.BondedForcesMap[typename])
            
    def Set_Visible_Bonded_Forces(self, types):    
        self.BondedForces = []
        for typename in types:
            self.BondedForces.append(self.BondedForcesMap[typename])
            
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
        return self.link["head"]
    
    @head.setter
    def head(self, atom):
        self.link["head"] = atom
        
    @property
    def tail(self):
        return self.link["tail"]
    
    @tail.setter
    def tail(self, atom):
        self.link["tail"] = atom
        
    @property
    def head_next(self):
        return self.link["head_next"]
    
    @head_next.setter
    def head_next(self, atom):
        self.link["head_next"] = atom
        
    @property
    def tail_next(self):
        return self.link["tail_next"]
    
    @tail_next.setter
    def tail_next(self, atom):
        self.link["tail_next"] = atom 
    
    @property
    def head_length(self):
        return self.link["head_length"]
    
    @head_length.setter
    def head_length(self, atom):
        self.link["head_length"] = atom
        
    @property
    def tail_length(self):
        return self.link["tail_length"]
    
    @tail_length.setter
    def tail_length(self, atom):
        self.link["tail_length"] = atom   
    
    @property
    def head_link_conditions(self):
        return self.link["head_link_conditions"]

    @property
    def tail_link_conditions(self):
        return self.link["tail_link_conditions"]
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.atoms = []
        self._name2atom = {}
        self._atom2name = {}
        self.connectivity = {}
        self.bonded_forces = {}
        self.builded = False
        self.bonded_forces = {frc.name:[] for frc in GlobalSetting.BondedForces}
        self.link = {"head": None, "tail": None, "head_next":None, "tail_next":None, 
        "head_length":1.5, "tail_length":1.5, "head_link_conditions":[], "tail_link_conditions":[]}
        self.connect_atoms = {}
        
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
    def deepcopy(self, forcopy = None):
        new_atom = Atom(self.entity)
        new_atom.contents = {**self.contents}
        if forcopy:
            self.copied[forcopy] = new_atom
        return new_atom


class Residue(Entity):
    name = "Residue"
    count = 0
    
    def __getattribute__(self, attr):
        if attr not in  ("_name2atom", "contents") and attr in self._name2atom.keys():
            return self._name2atom[attr]
        else:
            return super().__getattribute__(attr)
    def __init__(self, entity_type, name = None):
        super().__init__(entity_type, name)
        self.atoms = []
        self._name2atom = {}
        self._atom2name = {}
        self.connectivity = {}
        self.bonded_forces = {frc.name:[] for frc in GlobalSetting.BondedForces}
        self.builded = False
    
    def Add_Atom(self, name, atom_type = None, x = None, y = None, z = None):
        if type(name) == Atom:
            assert atom_type == None
            new_atom = Atom(name.type, name.name)
            new_atom.contents = {**name.contents}
            name = name.name
        else:
            if not atom_type: 
                atom_type = self.type._name2atom[name].type
                new_atom = Atom(atom_type, name)
                new_atom.contents = {**self.type._name2atom[name].contents}
            else:
                new_atom = Atom(atom_type, name)
                new_atom.contents = {**self.type._name2atom[name].contents}
        
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
    
    def deepcopy(self, forcopy = None):
        new_residue = Residue(self.type)
        for atom in self.atoms:
            new_residue.Add_Atom(atom)
            if forcopy:
                atom.copied[forcopy] = new_residue.atoms[-1]

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
    def deepcopy(self, forcopy):
        if self.atom1.copied[forcopy] and self.atom2.copied[forcopy]:
            return ResidueLink(self.atom1.copied[forcopy], self.atom2.copied[forcopy])
        
class Molecule():
    all = {}
    save_functions = []
    
    @classmethod
    def Set_Save_SPONGE_Input(cls, func):
        cls.save_functions.append(func)
        
    def __repr__(self):
        return "Entity of Molecule: " + self.name
    def __init__(self, name):
        if type(name) == ResidueType:
            self.name = name.name
        else:
            self.name = name
        Molecule.all[self.name] = self
        self.residues = []
        self.atoms = []
        self.residue_links = []
        self.bonded_forces = []
        self.builded = False
        if type(name) == ResidueType:
            new_residue = Residue(name)
            for i in name.atoms:
                new_residue.Add_Atom(i)
            self.Add_Residue(new_residue)
            
            
        
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
        forcopy = hash(str(time.time()))
        for res in self.residues:
            new_molecule.Add_Residue(res.deepcopy(forcopy))
        for link in self.residue_links:
            new_molecule.residue_links.append(link.deepcopy(forcopy))
        for res in self.residues:
            for atom in res.atoms:
                atom.copied.pop(forcopy)
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
    maxi = [-float("inf"), -float("inf"), -float("inf")]
    mini = [ float("inf"),  float("inf"),  float("inf")]
    for atom in self.atoms:
        if atom.x > maxi[0]:
            maxi[0] = atom.x
        if atom.y > maxi[1]:
            maxi[1] = atom.y
        if atom.z > maxi[2]:
            maxi[2] = atom.z
        if atom.x < mini[0]:
            mini[0] = atom.x
        if atom.y < mini[1]:
            mini[1] = atom.y
        if atom.z < mini[2]:
            mini[2] = atom.z

    towrite += "\n".join(["%f %f %f"%(atom.x - mini[0] + 3, atom.y - mini[1] + 3, atom.z - mini[2] + 3) for atom in self.atoms])
    
    boxlength[0] = maxi[0] - mini[0] + 6
    boxlength[1] = maxi[1] - mini[1] + 6
    boxlength[2] = maxi[2] - mini[2] + 6
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
    GlobalSetting.BondedForcesMap[BondedForceType.name] = BondedForceType
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



from . import assign, LOAD, BUILD, IMPOSE, PROCESS
