from . import *
import os
import sys

AtomType.New_From_String(
"""
name mass    charge[e]  LJtype
HW   1.008    0.417       HW
OW   16      -0.834       OW
""")

BOND.BondType.New_From_String(r"""
name k b
OW-HW  553  0.9572
HW-HW  553  1.5136
""")

ANGLE.AngleType.New_From_String(r"""
name   k    b
HW-OW-HW  0   0
OW-HW-HW  0   0
""")

LJ.LJType.New_From_String(r"""
name    epsilon    rmin
OW-OW   0.152      1.7683
HW-HW   0          0
""")

TIP3P = LOAD.mol2(os.path.join(AMBER_DATA_DIR, "tip3p.mol2"))


BUILD.Build_Bonded_Force(TIP3P)
sys.modules['__main__'].__dict__["WAT"] = TIP3P