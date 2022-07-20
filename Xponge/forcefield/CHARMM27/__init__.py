from ... import *
from ..BASE import CHARGE, MASS, LJ, BOND, UBANGLE, DIHEDRAL, NB14, NB14EXTRA, UBANGLE, IMPROPER, \
    VATOM, ACMAP, EXCLUDE
import os

CHARMM27_DATA_DIR = os.path.dirname(__file__)

LJ.LJType.combining_method_A = LJ.Lorentz_Berthelot_For_A
LJ.LJType.combining_method_B = LJ.Lorentz_Berthelot_For_B

GlobalSetting.Set_Invisible_Bonded_Forces(["improper"])

DIHEDRAL.ProperType.New_From_String(r"""
name        k reset  phi0 periodicity
X-X-X-X     0 0      0    0
""")
EXCLUDE.Exclude(4)

output = load_ffitp(os.path.join(CHARMM27_DATA_DIR, "forcefield.itp"))

AtomType.New_From_String(output["atomtypes"])
BOND.BondType.New_From_String(output["bonds"])
DIHEDRAL.ProperType.New_From_String(output["dihedrals"])
LJ.LJType.New_From_String(output["LJ"])
UBANGLE.UreyBradleyType.New_From_String(output["Urey-Bradley"])
IMPROPER.ImproperType.New_From_String(output["impropers"])
NB14EXTRA.NB14Type.New_From_String(output["nb14_extra"])
NB14.NB14Type.New_From_String(output["nb14"])
ACMAP.CMapType.New_From_Dict(output["cmaps"])
