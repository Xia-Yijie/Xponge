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
name   k[kcal/mol·A^-2]   b[A]
OW-HW  553                0.9572
HW-HW  553                1.5136
""")

UREY_BRADLEY.UreyBradleyType.New_From_String(r"""
name      k   kUB  
HW-OW-HW  0   0
OW-HW-HW  0   0
""")

LJ.LJType.New_From_String(r"""
name    sigma[nm]   epsilon[kJ/mol]
OW-OW   0.315057422683    0.6363864
HW-HW   0.0400013524445 0.192464
""")

TIPS3P = load_mol2(os.path.join(CHARMM27_DATA_DIR, "tip3p.mol2"))

load_mol2(os.path.join(CHARMM27_DATA_DIR, "atomic_ions.mol2"))

print("""Reference for CHARMM modified tip3p:
  to do
""")

sys.modules['__main__'].__dict__["WAT"] = TIPS3P
