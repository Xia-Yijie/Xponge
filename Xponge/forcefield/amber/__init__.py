"""
This **package** sets the basic configuration of amber force field
"""
import os
from ...helper import set_global_alternative_names, Generate_New_Bonded_Force_Type
from ... import AtomType, load_parmdat, load_frcmod

from ..base import charge_base, mass_base, lj_base, bond_base, angle_base, dihedral_base, nb14_base,\
    virtual_atom_base, exclude_base, cmap_base

AMBER_DATA_DIR = os.path.dirname(__file__)

lj_base.LJType.combining_method_A = lj_base.Lorentz_Berthelot_For_A
lj_base.LJType.combining_method_B = lj_base.Lorentz_Berthelot_For_B

nb14_base.NB14Type.New_From_String(r"""
name    kLJ     kee
X-X     0.5     0.833333
""")

exclude_base.Exclude(4)

# pylint: disable=invalid-name
AmberCMapType = Generate_New_Bonded_Force_Type("cmap", "1-2-3-4-5", {"resolution": int, "parameters": list}, False)

@AmberCMapType.Type_Name_Getter
def _(atoms):
    atom_names = [atom.type.name for atom in atoms]
    res_name = atoms[2].residue.name
    atom_names[2] = res_name + "@" + atom_names[2]
    return "-".join(atom_names)


@AmberCMapType.Set_Same_Force_Function
def cmap_same_force(_, atom_list):
    """
    This **function** is used to return the same force type for an atom list
    :param _:
    :param atom_list:
    :return:
    """
    return [atom_list]


def load_parameters_from_parmdat(filename, prefix=True):
    """
    This **function** is used to get amber force field parameters from parmdat files

    :param filename: the name of the input file
    :param prefix: whether add the AMBER_DATA_DIR to the filename
    :return:
    """
    if prefix:
        filename = os.path.join(AMBER_DATA_DIR, filename)
    atoms, bonds, angles, propers, impropers, ljs, nb14s = load_parmdat(filename)
    AtomType.New_From_String(atoms)
    bond_base.BondType.New_From_String(bonds)
    angle_base.AngleType.New_From_String(angles)
    dihedral_base.ProperType.New_From_String(propers)
    dihedral_base.ImproperType.New_From_String(impropers)
    lj_base.LJType.New_From_String(ljs)
    nb14_base.NB14Type.New_From_String(nb14s)


def load_parameters_from_frcmod(filename, include_cmap=False, prefix=True):
    """
    This **function** is used to get amber force field parameters from frcmod files

    :param filename: the name of the input file
    :param include_cmap: whether include cmap
    :param prefix: whether add the AMBER_DATA_DIR to the filename
    :return: None
    """
    if prefix:
        filename = os.path.join(AMBER_DATA_DIR, filename)
    atoms, bonds, angles, propers, impropers, ljs, cmap = load_frcmod(filename)

    AtomType.New_From_String(atoms)
    bond_base.BondType.New_From_String(bonds)
    angle_base.AngleType.New_From_String(angles)
    dihedral_base.ProperType.New_From_String(propers)
    dihedral_base.ImproperType.New_From_String(impropers)
    lj_base.LJType.New_From_String(ljs)

    if include_cmap:
        AmberCMapType.New_From_Dict(cmap)


set_global_alternative_names()
