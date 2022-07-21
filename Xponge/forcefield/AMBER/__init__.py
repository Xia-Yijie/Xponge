"""
This **package** sets the basic configuration of amber force field
"""
import os
from ...helper import set_global_alternative_names
from ... import AtomType, load_parmdat, load_frcmod

from ..BASE import CHARGE, MASS, LJ, BOND, ANGLE, DIHEDRAL, NB14, VATOM, EXCLUDE

AMBER_DATA_DIR = os.path.dirname(__file__)

LJ.LJType.combining_method_A = LJ.Lorentz_Berthelot_For_A
LJ.LJType.combining_method_B = LJ.Lorentz_Berthelot_For_B

NB14.NB14Type.New_From_String(r"""
name    kLJ     kee
X-X     0.5     0.833333
""")

EXCLUDE.Exclude(4)


def load_parameters_from_parmdat(filename):
    """
This **function** is used to get amber force field parameters from parmdat files
    :param filename:
    :return:
    """
    atoms, bonds, angles, propers, impropers, ljs = load_parmdat(os.path.join(AMBER_DATA_DIR, filename))
    AtomType.New_From_String(atoms)
    BOND.BondType.New_From_String(bonds)
    ANGLE.AngleType.New_From_String(angles)
    DIHEDRAL.ProperType.New_From_String(propers)
    DIHEDRAL.ImproperType.New_From_String(impropers)
    LJ.LJType.New_From_String(ljs)


def load_parameters_from_frcmod(filename, include_cmap=False):
    """
This **function** is used to get amber force field parameters from frcmod files
    :param filename:
    :param include_cmap:
    :return:
    """
    atoms, bonds, angles, propers, impropers, ljs, cmap = load_frcmod(os.path.join(AMBER_DATA_DIR, filename))

    AtomType.New_From_String(atoms)
    BOND.BondType.New_From_String(bonds)
    ANGLE.AngleType.New_From_String(angles)
    DIHEDRAL.ProperType.New_From_String(propers)
    DIHEDRAL.ImproperType.New_From_String(impropers)
    LJ.LJType.New_From_String(ljs)

    if include_cmap:
        from ..BASE import RCMAP
        RCMAP.CMapType.Residue_Map.update(cmap)


set_global_alternative_names(globals())
