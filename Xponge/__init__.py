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

__version__ = "stable-1.2.5"
__author__ = "Yijie Xia"

from collections import OrderedDict, deque
import os
import numpy as np
from itertools import product, permutations
import time
import sys

from . import assign

from .helper import GlobalSetting, Type, AtomType, ResidueType, Entity, Atom, Residue, ResidueLink, Molecule, AtomType, \
    set_global_alternative_names
from .load import load_ffitp, load_mol2, load_rst7, load_frcmod, load_pdb, load_parmdat
from .build import save_mol2, save_pdb, save_sponge_input, save_gro
from .process import impose_bond, impose_angle, impose_dihedral, add_solvent_box, h_mass_repartition, solvent_replace, \
    main_axis_rotate


def initialize():
    set_global_alternative_names(globals(), True)
    AtomType.New_From_String("name\nUNKNOWN")

    @Molecule.Set_Save_SPONGE_Input("residue")
    def write_residue(self):
        towrite = "%d %d\n" % (len(self.atoms), len(self.residues))
        towrite += "\n".join([str(len(res.atoms)) for res in self.residues])
        return towrite

    @Molecule.Set_Save_SPONGE_Input("coordinate")
    def write_coordinate(self):
        towrite = "%d\n" % (len(self.atoms))
        boxlength = [0, 0, 0, self.box_angle[0], self.box_angle[1], self.box_angle[2]]
        maxi = [-float("inf"), -float("inf"), -float("inf")]
        mini = [float("inf"), float("inf"), float("inf")]
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
        if not GlobalSetting.nocenter and self.box_length is None:
            towrite += "\n".join(
                ["%f %f %f" % (atom.x - mini[0] + 3, atom.y - mini[1] + 3, atom.z - mini[2] + 3) for atom in
                 self.atoms])
        else:
            towrite += "\n".join(["%f %f %f" % (atom.x, atom.y, atom.z) for atom in self.atoms])
        if self.box_length is None:
            boxlength[0] = maxi[0] - mini[0] + 6
            boxlength[1] = maxi[1] - mini[1] + 6
            boxlength[2] = maxi[2] - mini[2] + 6
            self.box_length = [boxlength[0], boxlength[1], boxlength[2]]
        else:
            boxlength[0] = self.box_length[0]
            boxlength[1] = self.box_length[1]
            boxlength[2] = self.box_length[2]
        towrite += "\n%f %f %f %f %f %f" % (
            boxlength[0], boxlength[1], boxlength[2], boxlength[3], boxlength[4], boxlength[5])
        return towrite


initialize()


def Generate_New_Bonded_Force_Type(Type_Name, atoms, properties, Compulsory, Multiple=None):
    class BondedForceEntity(Entity):
        name = Type_Name
        count = 0

        def __init__(self, atoms, entity_type, name=None):
            super().__init__(entity_type, name)
            self.atoms = atoms

        def deepcopy(self, forcopy):
            _atoms = [atom.copied[forcopy] for atom in self.atoms]
            newone = type(self)(_atoms, self.type, self.name)
            newone.contents = {**self.contents}
            return newone

    temp = [int(i) for i in atoms.split("-")]

    class BondedForceType(Type):
        name = Type_Name
        topology_like = temp
        compulsory = Compulsory
        multiple = Multiple
        atom_numbers = len(atoms.split("-"))
        topology_matrix = [[temp[i] - j if i > j else 1 for i in range(len(atoms.split("-")))] for j in
                           range(len(atoms.split("-")))]
        parameters = {
            "name": str,
        }
        entity = BondedForceEntity
        types = {}
        types_different_name = {}

        @classmethod
        def Same_Force(cls, atom_list):
            if type(atom_list) == str:
                atom_list_temp = [atom.strip() for atom in atom_list.split("-")]
                _temp = [atom_list, "-".join(atom_list_temp[::-1])]
            else:
                _temp = [atom_list, atom_list[::-1]]
            return _temp

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

        def Update(self, **kwargs):
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
            super().Update(**kwargs)
            if type(self).multiple:
                for key in type(self).multiple:
                    self.contents[key + "s"].append(self.contents[key])
                    self.contents[key] = None
                self.multiple_numbers += 1

    BondedForceType.Add_Property(properties)
    BondedForceType.New_From_String("""name
    UNKNOWNS""")
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
            "name": str,
        }
        types = {}

    PairwiseForceType.Add_Property(properties)

    return PairwiseForceType
