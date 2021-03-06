from ... import *
from ..BASE import CHARGE, MASS, LJ, BOND, UREY_BRADLEY, DIHEDRAL, NB14, NB14_EXTRA, UREY_BRADLEY, IMPROPER, VIRTUAL_ATOM, ACMAP, EXCLUDE
import os

CHARMM27_DATA_DIR = os.path.dirname(__file__)

LJ.LJType.combining_method_A = LJ.Lorentz_Berthelot_For_A
LJ.LJType.combining_method_B = LJ.Lorents_Berthelot_For_B

GlobalSetting.Set_Invisible_Bonded_Forces(["improper"])



DIHEDRAL.ProperType.New_From_String(r"""
name        k reset  phi0 periodicity
X-X-X-X     0 0      0    0
""")
EXCLUDE.Exclude(4)


output = LOAD.ffitp(os.path.join(CHARMM27_DATA_DIR, "forcefield.itp"))

AtomType.New_From_String(output["atomtypes"])
BOND.BondType.New_From_String(output["bonds"])
DIHEDRAL.ProperType.New_From_String(output["dihedrals"])
LJ.LJType.New_From_String(output["LJ"])
UREY_BRADLEY.UreyBradleyType.New_From_String(output["Urey-Bradley"])
IMPROPER.ImproperType.New_From_String(output["impropers"])
NB14_EXTRA.NB14Type.New_From_String(output["nb14_extra"])
NB14.NB14Type.New_From_String(output["nb14"])
ACMAP.CMAP.New_From_Dict(output["cmaps"])
