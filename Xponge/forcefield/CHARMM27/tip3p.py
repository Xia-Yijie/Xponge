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

UBANGLE.UreyBradleyType.New_From_String(r"""
name      k   kUB  
HW-OW-HW  0   0
OW-HW-HW  0   0
""")

LJ.LJType.New_From_String(r"""
name    epsilon[kcal/mol]   rmin[A]
OW-OW   0.152               1.7683
HW-HW   0                   0
""")

TIP3P = load_mol2(os.path.join(CHARMM27_DATA_DIR, "tip3p.mol2"))

load_mol2(os.path.join(CHARMM27_DATA_DIR, "atomic_ions.mol2"))

print("""Reference for tip3p:
  William L. Jorgensen, Jayaraman Chandrasekhar, and Jeffry D. Madura
    Comparison of simple potential functions for simulating liquid water
    The Journal of Chemical Physics 1983 79, 926-935, 
    DOI: 10.1063/1.445869
""")

sys.modules['__main__'].__dict__["WAT"] = TIP3P
