import sys
import time
import numpy as np
from types import MethodType
from functools import partial
from collections import OrderedDict
from itertools import product


class _GlobalSetting():
    # 最远的成键距离，用于拓扑分析时最远分析多远
    farthest_bonded_force = 0
    # 所有的成键类型力的Type
    BondedForces = []
    BondedForcesMap = {}
    # 所有虚拟原子的Type和对应的依赖的其他原子的数量
    VirtualAtomTypes = {}
    # 单位换算
    UnitMapping = {"distance": {"nm": 1e-9, "A": 1e-10},
                   "energy": {"kcal/mol": 4.184, "eV": 96.4853, "kJ/mol": 1},
                   "charge": {"e": 1, "SPONGE": 1.0 / 18.2223},
                   "angle": {"degree": np.pi, "rad": 180}
                   }
    PDBResidueNameMap = {"head": {}, "tail": {}, "save": {}}
    HISMap = {"DeltaH": "", "EpsilonH": "", "HIS": {}}
    nocenter = False
    boxspace = 3  # space between the molecule and the box border

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


class Type:
    name = None
    parameters = {"name": str}
    types = {}
    types_different_name = {}

    @classmethod
    def Add_Property(cls, parm_fmt, parm_default=None):
        cls.parameters.update(parm_fmt)
        if parm_default is None:
            parm_default = {}
        for _type in cls.types.values():
            _type.contents.update({key: parm_default.get(key, None) for key in parm_fmt.keys()})

    @classmethod
    def Set_Property_Unit(cls, prop, unit_type, base_unit):
        assert prop in cls.parameters.keys(), "Unknown property '%s' for type '%s'" % (prop, cls.name)
        temp_unit_lists = unit_type.split("·")
        temp_unit_power = []

        for i in range(len(temp_unit_lists)):
            unit = temp_unit_lists[i].split("^")
            if len(unit) == 2:
                temp_unit_lists[i] = unit[0]
                temp_unit_power.append(int(unit[1]))
            else:
                temp_unit_power.append(1)
            assert unit[0] in GlobalSetting.UnitMapping.keys(), "Unknown unit type '%s'" % unit_type

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
    def New_From_String(cls, string, skip_lines=0):
        count = -1
        kwargs = OrderedDict()
        for line in string.split("\n"):
            if not line.strip() or line.strip()[0] in ("!", "#", ";", "/"):
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
                    temp.Update(**tempkw)

    @classmethod
    def New_From_File(cls, filename, skip_lines=0):
        with open(filename, encoding='utf-8') as f:
            cls.New_From_String(f.read(), skip_lines)

    @classmethod
    def New_From_Dict(cls, Dict):
        for name, values in Dict.items():
            type_already_have = False
            if name in cls.types_different_name.keys():
                name = cls.types_different_name[name].name
                type_already_have = True
            if not type_already_have:
                if "reset" in values.keys():
                    values.pop("reset")
                values["name"] = name
                cls(**values)
            else:
                temp = cls.types[name]
                temp.Update(**values)

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
        assert self.name not in type(self).types.keys(), "The name '%s' has already existed in '%sType'" % (
            self.name, type(self).name)
        type(self).types[self.name] = self
        type(self).types_different_name[self.name] = self
        for key in kwargs.keys():
            assert key in self.contents.keys(), "The parameter '%s' is not one of the parameters of '%sType'" % (
                key, type(self).name)
            self.contents[key] = prop_fmt[key](kwargs[key])
        type(self)._unit_transfer(self)

    def Update(self, **kwargs):
        for key in kwargs.keys():
            assert (key in self.contents.keys())
            self.contents[key] = type(self).parameters[key](kwargs[key])
        type(self)._unit_transfer(self)


class AtomType(Type):
    name = "Atom"
    parameters = {"name": str, "x": float, "y": float, "z": float}
    types = {}
    types_different_name = {}


AtomType.New_From_String("name\nUNKNOWN")


class ResidueType(Type):
    name = "Residue"
    parameters = {"name": str, }
    types = {}
    types_different_name = {}

    def __getattribute__(self, attr):
        if attr not in ("_name2atom", "contents") and attr in self._name2atom.keys():
            return self._name2atom[attr]
        elif attr in AtomType.parameters.keys() and AtomType.parameters[attr] == float:
            return np.sum([getattr(atom, attr) for atom in self.atoms])
        else:
            return super().__getattribute__(attr)

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
        # 力场构建相关
        self.contents = {}
        self.connectivity = {}
        self.built = False
        self.bonded_forces = {frc.name: [] for frc in GlobalSetting.BondedForces}

        # 索引相关
        self._name2atom = {}
        self.atoms = []
        self._atom2name = {}
        self._atom2index = {}
        self._name2index = {}

        super().__init__(**kwargs)
        # 连接功能
        self.link = {"head": None, "tail": None, "head_next": None, "tail_next": None,
                     "head_length": 1.5, "tail_length": 1.5, "head_link_conditions": [], "tail_link_conditions": []}
        self.connect_atoms = {}

    def name2atom(self, name):
        return self._name2atom[name]

    def Add_Atom(self, name, atom_type, x, y, z):
        new_atom = Atom(atom_type, name)
        self.atoms.append(new_atom)
        new_atom.residue = self
        new_atom.x = float(x)
        new_atom.y = float(y)
        new_atom.z = float(z)
        self._name2atom[name] = new_atom
        self._atom2name[new_atom] = name
        self._atom2index[new_atom] = len(self.atoms) - 1
        self._name2index[name] = len(self.atoms) - 1
        self.connectivity[new_atom] = set([])

    def Add_Connectivity(self, atom0, atom1):
        if type(atom0) == str:
            atom0 = self._name2atom[atom0]
        if type(atom1) == str:
            atom1 = self._name2atom[atom1]
        self.connectivity[atom0].add(atom1)
        self.connectivity[atom1].add(atom0)

    def Add_Bonded_Force(self, bonded_force_entity, typename=None):
        if typename is None:
            typename = type(bonded_force_entity).name
        if typename not in self.bonded_forces.keys():
            self.bonded_forces[typename] = []
        self.bonded_forces[typename].append(bonded_force_entity)

    def deepcopy(self, name, forcopy=None):
        new_restype = ResidueType(name=name)
        donot_delete = True
        if forcopy is None:
            donot_delete = False
            forcopy = hash(str(time.time()))

        for atom in self.atoms:
            new_restype.Add_Atom(atom.name, atom.type, atom.x, atom.y, atom.z)
            atom.copied[forcopy] = new_restype.atoms[-1]
            atom.copied[forcopy].contents = {key: value for key, value in atom.contents.items()}

        for atom, connect_set in self.connectivity.items():
            for aton in connect_set:
                new_restype.Add_Connectivity(atom.copied[forcopy], aton.copied[forcopy])

        for atom in self.atoms:
            atom.copied[forcopy].Extra_Exclude_Atoms(map(lambda aton: aton.copied[forcopy], atom.extra_excluded_atoms))

        if self.built:
            for bond_entities in self.bonded_forces.values():
                for bond_entity in bond_entities:
                    new_restype.Add_Bonded_Force(bond_entity.deepcopy(forcopy))
            new_restype.built = True
            for atom in self.atoms:
                atom.copied[forcopy].linked_atoms = {key: set(map(lambda _atom: _atom.copied[forcopy], value)) for
                                                     key, value in atom.linked_atoms.items()}

        if not donot_delete:
            for atom in self.atoms:
                atom.copied.pop(forcopy)

        return new_restype


class Entity:
    count = 0
    name = None

    def __repr__(self):
        return "Entity of " + type(self).name + ": " + self.name + "(" + str(self.count) + ")"

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

    def __init__(self, entity_type, name=None):
        self.contents = {**entity_type.contents}
        self.count = type(self).count
        if not name:
            name = entity_type.name
        type(self).count += 1
        self.name = name
        self.type = entity_type

    def Update(self, **kwargs):
        for key in kwargs.keys():
            assert (key in self.contents.keys())
            self.contents[key] = type(self.type).parameters[key](kwargs[key])
        type(self.type)._unit_transfer(self)


class Atom(Entity):
    name = "Atom"
    count = 0

    @property
    def extra_excluded_atoms(self):
        return self.linked_atoms["extra_excluded_atoms"]

    def __init__(self, entity_type, name=None):
        # 力场基本信息
        super().__init__(entity_type, name)
        # 父信息
        self.residue = None

        # 成键信息
        self.linked_atoms = {i + 1: set() for i in range(1, GlobalSetting.farthest_bonded_force)}
        self.linked_atoms["extra_excluded_atoms"] = set()

        # 复制信息
        self.copied = {}

    def deepcopy(self, forcopy=None):
        new_atom = Atom(self.type, self.name)
        new_atom.contents = {**self.contents}
        if forcopy:
            self.copied[forcopy] = new_atom
        return new_atom

    def Link_Atom(self, link_type, atom):
        if link_type not in self.linked_atoms.keys():
            self.linked_atoms[link_type] = set()
        self.linked_atoms[link_type].add(atom)

    def Extra_Exclude_Atom(self, atom):
        self.extra_excluded_atoms.add(atom)
        atom.extra_excluded_atoms.add(self)

    def Extra_Exclude_Atoms(self, lists):
        for atom in lists:
            self.Extra_Exclude_Atom(atom)


class Residue(Entity):
    name = "Residue"
    count = 0

    def __getattribute__(self, attr):
        if attr not in ("_name2atom", "contents") and attr in self._name2atom.keys():
            return self._name2atom[attr]
        elif attr in AtomType.parameters.keys() and AtomType.parameters[attr] == float:
            return np.sum([getattr(atom, attr) for atom in self.atoms])
        else:
            return super().__getattribute__(attr)

    def __init__(self, entity_type, name=None, directly_copy=False):
        super().__init__(entity_type, name)
        self.atoms = []
        self._name2atom = {}
        self._atom2name = {}
        self._atom2index = {}
        self._name2index = {}
        self.connectivity = {}
        self.bonded_forces = {frc.name: [] for frc in GlobalSetting.BondedForces}
        self.built = False
        if directly_copy:
            forcopy = hash(int(time.time()))
            for atom in self.type.atoms:
                self.Add_Atom(atom.name, atom.type, atom.x, atom.y, atom.z)
                atom.copied[forcopy] = self.atoms[-1]
            for atom in self.type.atoms:
                for aton in atom.extra_excluded_atoms:
                    atom.copied[forcopy].Extra_Exclude_Atom(aton.copied[forcopy])

    def name2atom(self, name):
        return self._name2atom[name]

    def Add_Atom(self, name, atom_type=None, x=None, y=None, z=None):
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
        self._atom2index[new_atom] = len(self.atoms) - 1
        self._name2index[name] = len(self.atoms) - 1
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

    def Add_Missing_Atoms(self):
        t = set([atom.name for atom in self.atoms])
        uncertified = set([atom.name for atom in self.type.atoms])
        for atom in self.type.atoms:
            if atom.name in t:
                uncertified.remove(atom.name)
        while uncertified:
            movedlist = []
            for atom_name in uncertified:
                temp_atom = getattr(self.type, atom_name)
                for connected_atom in self.type.connectivity[temp_atom]:
                    if connected_atom.name in t:
                        fact_connected_atom = self._name2atom[connected_atom.name]
                        _x = temp_atom.x - connected_atom.x + fact_connected_atom.x
                        _y = temp_atom.y - connected_atom.y + fact_connected_atom.y
                        _z = temp_atom.z - connected_atom.z + fact_connected_atom.z
                        t.add(atom_name)
                        movedlist.append(atom_name)
                        self.Add_Atom(atom_name, x=_x, y=_y, z=_z)
                        break
            for atom_name in movedlist:
                uncertified.remove(atom_name)

    def deepcopy(self, forcopy=None):
        new_residue = Residue(self.type)
        for atom in self.atoms:
            new_residue.Add_Atom(atom)
            new_residue.atoms[-1].contents = {key: value for key, value in atom.contents.items()}
            if forcopy:
                atom.copied[forcopy] = new_residue.atoms[-1]

        return new_residue


class ResidueLink:
    def __repr__(self):
        return "Entity of ResidueLink: " + repr(self.atom1) + "-" + repr(self.atom2)

    def __hash__(self):
        return hash(repr(self))

    def __init__(self, atom1, atom2):
        self.atom1 = atom1
        self.atom2 = atom2
        self.built = False
        self.bonded_forces = {frc.name: [] for frc in GlobalSetting.BondedForces}

    def Add_Bonded_Force(self, bonded_force_entity):
        if type(bonded_force_entity).name not in self.bonded_forces.keys():
            self.bonded_forces[type(bonded_force_entity).name] = []
        self.bonded_forces[type(bonded_force_entity).name].append(bonded_force_entity)

    def deepcopy(self, forcopy):
        if self.atom1.copied[forcopy] and self.atom2.copied[forcopy]:
            return ResidueLink(self.atom1.copied[forcopy], self.atom2.copied[forcopy])


class Molecule:
    all = {}
    save_functions = {}

    @staticmethod
    def _set_friends_in_different_residue(molecule, atom1, atom2):
        if molecule.atom_index[atom1.residue.atoms[0]] < molecule.atom_index[atom1.residue.atoms[1]]:
            res_index = molecule.atom_index[atom1.residue.atoms[-1]]
            atom1_friends = list(range(res_index + 1))
            atom2_friends = list(range(res_index + 1, len(molecule.atoms)))
        else:
            res_index = molecule.atom_index[atom2.residue.atoms[-1]]
            atom2_friends = list(range(res_index + 1))
            atom1_friends = list(range(res_index + 1, len(molecule.atoms)))
        return atom1_friends, atom2_friends

    @staticmethod
    def _get_head_and_tail(head, tail, molecule, typeatom1, typeatom2, restype, toset, atom1_friends, _resatom):
        assert typeatom2 in restype.connectivity[typeatom1]
        index_dict = {}.fromkeys(restype.connectivity[typeatom1], typeatom1)
        if typeatom2 in index_dict.keys():
            index_dict.pop(typeatom2)

        while index_dict:
            index_next = {}
            for atom0, from_atom in index_dict.items():
                if atom0.name == restype.head:
                    head = toset
                elif atom0.name == restype.tail:
                    tail = toset
                atom1_friends.append(molecule.atom_index[_resatom(atom0)])
                index_temp = {}.fromkeys(restype.connectivity[atom0], atom0)
                index_temp.pop(from_atom)
                if typeatom2 in index_temp.keys():
                    index_temp.pop(typeatom2)
                index_next.update(index_temp)
            index_dict = index_next
        return head, tail

    @classmethod
    def Set_Save_SPONGE_Input(cls, keyname):
        def wrapper(func):
            cls.save_functions[keyname] = func
            return func

        return wrapper

    @classmethod
    def Del_Save_SPONGE_Input(cls, keyname):
        cls.save_functions.pop(keyname)

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
        self.bonded_forces = {}
        self.built = False
        self.box_length = None
        self.box_angle = [90.0, 90.0, 90.0]
        if type(name) == ResidueType:
            new_residue = Residue(name)
            for i in name.atoms:
                new_residue.Add_Atom(i)
            self.Add_Residue(new_residue)

    def __getattribute__(self, attr):
        if attr in AtomType.parameters.keys() and AtomType.parameters[attr] == float:
            return np.sum([getattr(atom, attr) for res in self.residues for atom in res.atoms])
        else:
            return super().__getattribute__(attr)

    def Add_Residue(self, residue):
        if type(residue) == Residue:
            self.built = False
            self.residues.append(residue)
        elif type(residue) == ResidueType:
            self.built = False
            self.residues.append(Residue(residue, directly_copy=True))

    def Add_Bonded_Force(self, bonded_force_entity):
        if type(bonded_force_entity).name not in self.bonded_forces.keys():
            self.bonded_forces[type(bonded_force_entity).name] = []
        self.bonded_forces[type(bonded_force_entity).name].append(bonded_force_entity)

    def Add_Residue_Link(self, atom1, atom2):
        self.built = False
        self.residue_links.append(ResidueLink(atom1, atom2))

    def Add_Missing_Atoms(self):
        for residue in self.residues:
            residue.Add_Missing_Atoms()

    def deepcopy(self):
        new_molecule = Molecule(self.name)
        forcopy = hash(str(time.time()))
        for res in self.residues:
            new_molecule.Add_Residue(res.deepcopy(forcopy))

        for link in self.residue_links:
            new_molecule.residue_links.append(link.deepcopy(forcopy))

        for res in self.residues:
            for atom in res.atoms:
                atom.copied[forcopy].Extra_Exclude_Atoms(
                    map(lambda aton: aton.copied[forcopy], atom.extra_excluded_atoms))

        if self.built:
            for bond_entities in self.bonded_forces.values():
                for bond_entity in bond_entities:
                    new_molecule.Add_Bonded_Force(bond_entity.deepcopy(forcopy))
            new_molecule.built = True
            new_molecule.atoms = [atom for residue in new_molecule.residues for atom in residue.atoms]
            new_molecule.atom_index = {new_molecule.atoms[i]: i for i in range(len(new_molecule.atoms))}
            for atom in self.atoms:
                atom.copied[forcopy].linked_atoms = {key: set(map(lambda aton: aton.copied[forcopy], value)) for
                                                     key, value in atom.linked_atoms.items()}

        for res in self.residues:
            for atom in res.atoms:
                atom.copied.pop(forcopy)
        return new_molecule

    def get_atom_coordinates(self):
        self.atoms = []
        for res in self.residues:
            self.atoms.extend(res.atoms)

        self.atom_index = {self.atoms[i]: i for i in range(len(self.atoms))}
        return np.array([[atom.x, atom.y, atom.z] for atom in self.atoms])

    def divide_into_two_parts(self, atom1, atom2):
        if atom1.residue != atom2.residue:
            atom1_friends, atom2_friends = self._set_friends_in_different_residue(self, atom1, atom2)
        else:
            link_front = 0
            atom1_friends = []
            atom2_friends = []
            head = 0
            tail = 0

            restype = atom1.residue.type

            _restype_atom = lambda atom: restype._name2atom[atom.name]
            _resatom = lambda atom: atom1.residue._name2atom[atom.name]

            typeatom1 = _restype_atom(atom1)
            typeatom2 = _restype_atom(atom2)

            head, tail = self._get_head_and_tail(head, tail, self, typeatom1, typeatom2, restype, 1, atom1_friends,
                                                 _resatom)
            head, tail = self._get_head_and_tail(head, tail, self, typeatom2, typeatom1, restype, 2, atom2_friends,
                                                 _resatom)

            if atom1.name == restype.head:
                head = 1
            elif atom1.name == restype.tail:
                tail = 1
            if atom2.name == restype.head:
                head = 2
            elif atom2.name == restype.tail:
                tail = 2

            resindex_head = min(self.atom_index[atom1.residue.atoms[0]],
                                self.atom_index[atom2.residue.atoms[0]])
            resindex_tail = max(self.atom_index[atom1.residue.atoms[-1]],
                                self.atom_index[atom2.residue.atoms[-1]])

            if head == 1:
                atom1_friends.extend(list(range(resindex_head)))
            else:
                atom2_friends.extend(list(range(resindex_head)))
            if tail == 1:
                atom1_friends.extend(list(range(resindex_tail + 1, len(self.atoms))))
            else:
                atom2_friends.extend(list(range(resindex_tail + 1, len(self.atoms))))

            atom1_friends = set(atom1_friends)
            atom1_friends.add(self.atom_index[atom1])
            atom1_friends = np.array(list(atom1_friends))
            atom2_friends = set(atom2_friends)
            atom2_friends.add(self.atom_index[atom2])
            atom2_friends = np.array(list(atom2_friends))
        return atom1_friends, atom2_friends


def _link_residue_process_coordinate(molecule, atom1, atom2):
    resA = atom1.residue
    resB = atom2.residue
    crd = molecule.get_atom_coordinates()
    atom1_friends, atom2_friends = molecule.divide_into_two_parts(atom1, atom2)
    crd[atom2_friends] += 2000

    bond_length = (resA.type.tail_length + resB.type.head_length) / 2
    r0 = crd[molecule.atom_index[atom2]] - crd[molecule.atom_index[atom1]]
    L0 = np.linalg.norm(r0)
    dr = (bond_length / L0 - 1) * r0
    crd[atom2_friends] += dr

    res = resA
    atomA = atom1
    atomB = atom2
    atomB_friends = atom2_friends
    for link_conditions in res.type.tail_link_conditions:
        atoms = [res._name2atom[atom] for atom in link_conditions["atoms"]]
        parameter = link_conditions["parameter"]
        if len(atoms) == 1:
            r0 = crd[molecule.atom_index[atomB]] - crd[molecule.atom_index[atoms[0]]]
            L0 = np.linalg.norm(r0)
            dr = (parameter / L0 - 1) * r0
            crd[atomB_friends] += dr
        elif len(atoms) == 2:
            rAO = crd[molecule.atom_index[atoms[0]]] - crd[molecule.atom_index[atoms[1]]]
            rOB = crd[molecule.atom_index[atomB]] - crd[molecule.atom_index[atoms[1]]]
            angle0 = np.arccos(np.dot(rAO, rOB) / np.linalg.norm(rAO) / np.linalg.norm(rOB))
            deltaAngle = parameter - angle0
            crd[atomB_friends] = np.dot(crd[atomB_friends] - crd[molecule.atom_index[atoms[1]]],
                                        get_rotate_matrix(np.cross(rAO, rOB), deltaAngle)) + crd[
                                     molecule.atom_index[atoms[1]]]
        elif len(atoms) == 3:
            rOO = crd[molecule.atom_index[atoms[0]]] - crd[molecule.atom_index[atoms[1]]]
            rOA = crd[molecule.atom_index[atoms[1]]] - crd[molecule.atom_index[atoms[2]]]
            rAB = crd[molecule.atom_index[atoms[1]]] - crd[molecule.atom_index[atomB]]
            r12xr23 = np.cross(rOO, rOA)
            r23xr34 = np.cross(rAB, rOA)
            cos = np.dot(r12xr23, r23xr34) / np.linalg.norm(r12xr23) / np.linalg.norm(r23xr34)
            cos = max(-0.999999, min(cos, 0.999999))
            dihedral0 = np.arccos(cos)
            dihedral0 = np.pi - np.copysign(dihedral0, np.cross(r23xr34, r12xr23).dot(rOA))
            deltaAngle = parameter - dihedral0
            crd[atomB_friends] = np.dot(crd[atomB_friends] - crd[molecule.atom_index[atoms[2]]],
                                        get_rotate_matrix(rOA, deltaAngle)) + crd[molecule.atom_index[atoms[2]]]

    res = resB
    atomA = atom2
    atomB = atom1
    atomB_friends = atom1_friends
    for link_conditions in res.type.head_link_conditions:
        atoms = [res._name2atom[atom] for atom in link_conditions["atoms"]]
        parameter = link_conditions["parameter"]
        if len(atoms) == 1:
            r0 = crd[molecule.atom_index[atomB]] - crd[molecule.atom_index[atoms[0]]]
            L0 = np.linalg.norm(r0)
            dr = (parameter / L0 - 1) * r0
            crd[atomB_friends] += dr
        elif len(atoms) == 2:
            rAO = crd[molecule.atom_index[atoms[0]]] - crd[molecule.atom_index[atoms[1]]]
            rOB = crd[molecule.atom_index[atomB]] - crd[molecule.atom_index[atoms[1]]]
            angle0 = np.arccos(np.dot(rAO, rOB) / np.linalg.norm(rAO) / np.linalg.norm(rOB))
            deltaAngle = parameter - angle0
            crd[atomB_friends] = np.dot(crd[atomB_friends] - crd[molecule.atom_index[atoms[1]]],
                                        get_rotate_matrix(np.cross(rAO, rOB), deltaAngle)) + crd[
                                     molecule.atom_index[atoms[1]]]
        elif len(atoms) == 3:
            rOO = crd[molecule.atom_index[atoms[0]]] - crd[molecule.atom_index[atoms[1]]]
            rOA = crd[molecule.atom_index[atoms[1]]] - crd[molecule.atom_index[atoms[2]]]
            rAB = crd[molecule.atom_index[atoms[2]]] - crd[molecule.atom_index[atomB]]
            r12xr23 = np.cross(rOO, rOA)
            r23xr34 = np.cross(rAB, rOA)
            cos = np.dot(r12xr23, r23xr34) / np.linalg.norm(r12xr23) / np.linalg.norm(r23xr34)
            cos = max(-0.999999, min(cos, 0.999999))
            dihedral0 = np.arccos(cos)
            dihedral0 = np.pi - np.copysign(dihedral0, np.cross(r23xr34, r12xr23).dot(rOA))
            deltaAngle = parameter - dihedral0
            crd[atomB_friends] = np.dot(crd[atomB_friends] - crd[molecule.atom_index[atoms[2]]],
                                        get_rotate_matrix(rOA, deltaAngle)) + crd[molecule.atom_index[atoms[2]]]

    if resA.type.tail_next and resB.type.head_next:
        atomA = resA._name2atom[resA.type.tail_next]
        atomB = resB._name2atom[resB.type.head_next]
        rOO = crd[molecule.atom_index[atomA]] - crd[molecule.atom_index[atom1]]
        rOA = crd[molecule.atom_index[atom1]] - crd[molecule.atom_index[atom2]]
        rAB = crd[molecule.atom_index[atom2]] - crd[molecule.atom_index[atomB]]
        r12xr23 = np.cross(rOO, rOA)
        r23xr34 = np.cross(rAB, rOA)
        cos = np.dot(r12xr23, r23xr34) / np.linalg.norm(r12xr23) / np.linalg.norm(r23xr34)
        cos = max(-0.999999, min(cos, 0.999999))
        dihedral0 = np.arccos(cos)
        dihedral0 = np.pi - np.copysign(dihedral0, np.cross(r23xr34, r12xr23).dot(rOA))
        deltaAngle = np.pi - dihedral0
        crd[atom2_friends] = np.dot(crd[atom2_friends] - crd[molecule.atom_index[atom2]],
                                    get_rotate_matrix(rOA, deltaAngle)) + crd[molecule.atom_index[atom2]]

    for atom in molecule.atoms:
        i = molecule.atom_index[atom]
        atom.x = crd[i][0]
        atom.y = crd[i][1]
        atom.z = crd[i][2]


def ResidueType_Add(self, other):
    if isinstance(other, ResidueType):
        new_molecule = Molecule(self.name)
        resA = Residue(self)
        resB = Residue(other)
        for atom in self.atoms:
            resA.Add_Atom(atom)
        for atom in other.atoms:
            resB.Add_Atom(atom)
        new_molecule.Add_Residue(resA)
        new_molecule.Add_Residue(resB)
        if resA.type.tail and resB.type.head:
            atom1 = resA._name2atom[self.tail]
            atom2 = resB._name2atom[other.head]
            new_molecule.Add_Residue_Link(atom1, atom2)
            _link_residue_process_coordinate(new_molecule, atom1, atom2)
        return new_molecule
    elif isinstance(other, Molecule):
        new_molecule = other.deepcopy()
        resA = Residue(self)
        resB = new_molecule.residues[0]
        for atom in self.atoms:
            resA.Add_Atom(atom)
        new_molecule.residues.insert(0, resA)
        if resA.type.tail and resB.type.head:
            atom1 = resA._name2atom[resA.type.tail]
            atom2 = resB._name2atom[resB.type.head]
            new_molecule.Add_Residue_Link(atom1, atom2)
            _link_residue_process_coordinate(new_molecule, atom1, atom2)
        return new_molecule
    elif other is None:
        return self
    else:
        raise TypeError("unsupported operand type(s) for +: '%s' and '%s'" % (type(self), type(other)))


def Molecule_Add(self, other):
    if isinstance(other, ResidueType):
        new_molecule = self.deepcopy()
        resA = new_molecule.residues[-1]
        resB = Residue(other)
        for atom in other.atoms:
            resB.Add_Atom(atom)
        new_molecule.Add_Residue(resB)
        if resA.type.tail and resB.type.head:
            atom1 = resA._name2atom[resA.type.tail]
            atom2 = resB._name2atom[resB.type.head]
            new_molecule.Add_Residue_Link(atom1, atom2)
            _link_residue_process_coordinate(new_molecule, atom1, atom2)
        return new_molecule
    elif isinstance(other, Molecule):
        new_molecule = self.deepcopy()
        new_molecule2 = other.deepcopy()
        resA = new_molecule.residues[-1]
        resB = new_molecule2.residues[0]
        for res in new_molecule2.residues:
            new_molecule.Add_Residue(res)
        for reslink in new_molecule2.residue_links:
            new_molecule.Add_Residue_Link(reslink.atom1, reslink.atom2)
        if resA.type.tail and resB.type.head:
            atom1 = resA._name2atom[resA.type.tail]
            atom2 = resB._name2atom[resB.type.head]
            new_molecule.Add_Residue_Link(atom1, atom2)
            _link_residue_process_coordinate(new_molecule, atom1, atom2)
        return new_molecule
    elif other is None:
        return self
    else:
        raise TypeError("unsupported operand type(s) for +: '%s' and '%s'" % (type(self), type(other)))


def iMolecule_Add(self, other):
    if isinstance(other, ResidueType):
        resA = self.residues[-1]
        resB = Residue(other)
        for atom in other.atoms:
            resB.Add_Atom(atom)
        self.Add_Residue(resB)
        if resA.type.tail and resB.type.head:
            atom1 = resA._name2atom[resA.type.tail]
            atom2 = resB._name2atom[other.head]
            self.Add_Residue_Link(atom1, atom2)
            _link_residue_process_coordinate(self, atom1, atom2)
        return self
    elif isinstance(other, Molecule):
        new_molecule2 = other.deepcopy()
        resA = self.residues[-1]
        resB = new_molecule2.residues[0]
        for res in new_molecule2.residues:
            self.Add_Residue(res)
        for reslink in new_molecule2.residue_links:
            self.Add_Residue_Link(reslink.atom1, reslink.atom2)
        if resA.type.tail and resB.type.head:
            atom1 = resA._name2atom[resA.type.tail]
            atom2 = resB._name2atom[resB.type.head]
            self.Add_Residue_Link(atom1, atom2)
            _link_residue_process_coordinate(self, atom1, atom2)
        return self
    elif other is None:
        return self
    else:
        raise TypeError("unsupported operand type(s) for +: '%s' and '%s'" % (type(self), type(other)))


def Muls(self, other):
    if type(other) == int:
        assert other >= 1
        if type(self) == ResidueType:
            t = self
        else:
            t = self.deepcopy()
        for i in range(other - 1):
            t += self
        return t
    else:
        raise TypeError("unsupported operand type(s) for *: '%s' and '%s'" % (type(self), type(other)))


def iMuls(self, other):
    if isinstance(other, int):
        assert other >= 1
        for i in range(other - 1):
            self += self
        return self
    else:
        raise TypeError("unsupported operand type(s) for *: '%s' and '%s'" % (type(self), type(other)))


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


def ResidueType_Or(self, other):
    if isinstance(other, ResidueType):
        new_molecule = Molecule(self.name)
        resA = Residue(self)
        resB = Residue(other)
        for atom in self.atoms:
            resA.Add_Atom(atom)
        for atom in other.atoms:
            resB.Add_Atom(atom)
        new_molecule.Add_Residue(resA)
        new_molecule.Add_Residue(resB)
        return new_molecule
    elif isinstance(other, Molecule):
        new_molecule = other.deepcopy()
        resA = Residue(self)
        resB = new_molecule.residues[0]
        for atom in self.atoms:
            resA.Add_Atom(atom)
        new_molecule.residues.insert(0, resA)
        return new_molecule
    elif other is None:
        return self
    else:
        raise TypeError("unsupported operand type(s) for |: '%s' and '%s'" % (type(self), type(other)))


def Molecule_Or(self, other):
    if isinstance(other, ResidueType):
        new_molecule = self.deepcopy()
        resA = new_molecule.residues[-1]
        resB = Residue(other)
        for atom in other.atoms:
            resB.Add_Atom(atom)
        new_molecule.Add_Residue(resB)
        return new_molecule
    elif isinstance(other, Molecule):
        new_molecule = self.deepcopy()
        new_molecule2 = other.deepcopy()
        resA = new_molecule.residues[-1]
        resB = new_molecule2.residues[0]
        for res in new_molecule2.residues:
            new_molecule.Add_Residue(res)
        new_molecule.residue_links.extend(new_molecule2.residue_links)
        return new_molecule
    elif other is None:
        return self
    else:
        raise TypeError("unsupported operand type(s) for +: '%s' and '%s'" % (type(self), type(other)))


def iMolecule_Or(self, other):
    if isinstance(other, ResidueType):
        resA = self.residues[-1]
        resB = Residue(other)
        for atom in other.atoms:
            resB.Add_Atom(atom)
        self.Add_Residue(resB)
        return self
    elif isinstance(other, Molecule):
        new_molecule2 = other.deepcopy()
        resA = self.residues[-1]
        resB = new_molecule2.residues[0]
        for res in new_molecule2.residues:
            self.Add_Residue(res)
        self.residue_links.extend(new_molecule2.residue_links)
        return self
    elif other is None:
        return self
    else:
        raise TypeError("unsupported operand type(s) for +: '%s' and '%s'" % (type(self), type(other)))


ResidueType.__or__ = ResidueType_Or
ResidueType.__ror__ = ResidueType_Or

Molecule.__or__ = Molecule_Or
Molecule.__ror__ = Molecule_Or
Molecule.__ior__ = iMolecule_Or

del ResidueType_Or
del Molecule_Or
del iMolecule_Or


def get_rotate_matrix(r0, angle):
    cost = np.cos(angle)
    cost_one = 1 - cost
    sint = np.sin(angle)
    r0 /= np.linalg.norm(r0)
    return np.array([[r0[0] * r0[0] * cost_one + cost, r0[0] * r0[1] * cost_one - r0[2] * sint,
                      r0[0] * r0[2] * cost_one + r0[1] * sint],
                     [r0[0] * r0[1] * cost_one + r0[2] * sint, r0[1] * r0[1] * cost_one + cost,
                      r0[1] * r0[2] * cost_one - r0[0] * sint],
                     [r0[0] * r0[2] * cost_one - r0[1] * sint, r0[1] * r0[2] * cost_one + r0[0] * sint,
                      r0[2] * r0[2] * cost_one + cost]]).transpose()

def guess_element_from_mass(mass):
    """

    :param mass:
    :return:
    """
    elements = ["X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
                "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc",
                "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
                "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc",
                "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
                "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
                "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
                "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
                "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
                "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
                "Ds", "Rg"]
    masses = [0.00000, 1.00794, 4.00260, 6.941, 9.012182, 10.811,
              12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,
              22.989770, 24.3050, 26.981538, 28.0855, 30.973761,
              32.065, 35.453, 39.948, 39.0983, 40.078, 44.955910,
              47.867, 50.9415, 51.9961, 54.938049, 55.845, 58.9332,
              58.6934, 63.546, 65.409, 69.723, 72.64, 74.92160,
              78.96, 79.904, 83.798, 85.4678, 87.62, 88.90585,
              91.224, 92.90638, 95.94, 98.0, 101.07, 102.90550,
              106.42, 107.8682, 112.411, 114.818, 118.710, 121.760,
              127.60, 126.90447, 131.293, 132.90545, 137.327,
              138.9055, 140.116, 140.90765, 144.24, 145.0, 150.36,
              151.964, 157.25, 158.92534, 162.500, 164.93032,
              167.259, 168.93421, 173.04, 174.967, 178.49, 180.9479,
              183.84, 186.207, 190.23, 192.217, 195.078, 196.96655,
              200.59, 204.3833, 207.2, 208.98038, 209.0, 210.0, 222.0,
              223.0, 226.0, 227.0, 232.0381, 231.03588, 238.02891,
              237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0,
              258.0, 259.0, 262.0, 261.0, 262.0, 266.0, 264.0, 269.0,
              268.0, 271.0, 272.0]
    if mass > 0.0 and mass < 3.8:
        index = 1
    elif mass > 207.85 and mass < 208.99:
        index = 83
    elif mass > 56.50 and mass < 58.8133:
        index = 27
    else:
        index = 0
        for j in range(0, 111):
            if abs(mass - masses[j]) < 0.65:
                index = j
                break
    return elements[index]

SPECIAL_STRINGS = {"Pdb": "PDB", "Sponge": "SPONGE", "Nb14": "NB14", "Nb14extra": "NB14EXTRA", "Lj": "LJ",
                   "Residuetype": "ResidueType", "Pubchem": "PubChem"}


def set_real_global_variable(name, value):
    sys.modules["__main__"].__dict__[name] = value


def set_alternative_name(object, func, set_method):
    name = func.__name__
    new_name = "_".join([i.capitalize() for i in name.split("_")])
    second_new_name = "".join([i.capitalize() for i in name.split("_")])
    third_new_name = second_new_name[0].lower() + second_new_name[1:]

    set_method(object, new_name, func)
    new_new_name = new_name
    for t, newt in SPECIAL_STRINGS.items():
        new_new_name = new_new_name.replace(t, newt)
        second_new_name = second_new_name.replace(t, newt)
        third_new_name = third_new_name.replace(t, newt)
    if new_new_name != new_name:
        set_method(object, new_new_name, func)
        set_method(object, second_new_name.replace(t, newt), func)
        set_method(object, third_new_name.replace(t, newt), func)


set_attribute_alternative_name = partial(set_alternative_name, set_method=setattr)


def set_attribute_alternative_names(instance):
    for i in type(instance).__dict__:
        self_func = getattr(instance, i, None)
        if isinstance(self_func, MethodType) and not self_func.__name__.startswith("_"):
            set_attribute_alternative_name(instance, self_func)


def _dict_set_method(obj, name, func):
    obj[name] = func


set_dict_value_alternative_name = partial(set_alternative_name, set_method=_dict_set_method)


def set_global_alternative_names(dict, real_global=False):
    from types import FunctionType
    new_dict = {}
    for key, value in dict.items():
        if not isinstance(value, FunctionType) or value.__name__.startswith("_"):
            continue

        if real_global:
            def _global_set_method(obj, name, func):
                obj[name] = func
                set_real_global_variable(name, func)

            set_alternative_name(new_dict, value, _global_set_method)
        else:
            set_dict_value_alternative_name(new_dict, value)
    dict.update(new_dict)
