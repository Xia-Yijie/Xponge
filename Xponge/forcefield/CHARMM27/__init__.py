from ... import *
from ..BASE import CHARGE, MASS, LJ, BOND, DIHEDRAL, NB14, NB14_EXTRA, UREY_BRADLEY, IMPROPER, VIRTUAL_ATOM, ACMAP, EXCLUDE
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

