from ... import *
from ..BASE import CHARGE, MASS, LJ, BOND, ANGLE, DIHEDRAL, NB14, VATOM, EXCLUDE
import os

AMBER_DATA_DIR = os.path.dirname(__file__)

LJ.LJType.combining_method_A = LJ.Lorentz_Berthelot_For_A
LJ.LJType.combining_method_B = LJ.Lorentz_Berthelot_For_B

NB14.NB14Type.New_From_String(r"""
name    kLJ     kee
X-X     0.5     0.833333
""")

EXCLUDE.Exclude(4)
