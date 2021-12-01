from . import *
import os
import sys

AtomType.New_From_String(
"""
name mass    charge[e]  LJtype
HW   1.008    0.4238      HW
OW   16      -0.8476      OW
""")

BOND.BondType.New_From_String(r"""
name   k[kcal/mol·A^-2]   b[A]
OW-HW  553                1.0000
HW-HW  553                1.6330
""")

ANGLE.AngleType.New_From_String(r"""
name      k   b
HW-OW-HW  0   0
OW-HW-HW  0   0
""")

LJ.LJType.New_From_String(r"""
name    epsilon[kcal/mol]   rmin[A]
OW-OW   0.1553              1.7767
HW-HW   0                   0
""")

SPCE = LOAD.mol2(os.path.join(AMBER_DATA_DIR, "tip3p.mol2"))

atoms, bonds, angles, propers, impropers, LJs, cmap = LOAD.frcmod(os.path.join(AMBER_DATA_DIR, "ions1lm_126_spce.frcmod"))
AtomType.New_From_String(atoms)
BOND.BondType.New_From_String(bonds)
ANGLE.AngleType.New_From_String(angles)
DIHEDRAL.ProperType.New_From_String(propers)
DIHEDRAL.ImproperType.New_From_String(impropers)
LJ.LJType.New_From_String(LJs)

atoms, bonds, angles, propers, impropers, LJs, cmap = LOAD.frcmod(os.path.join(AMBER_DATA_DIR, "ionsjc_spce.frcmod"))
AtomType.New_From_String(atoms)
BOND.BondType.New_From_String(bonds)
ANGLE.AngleType.New_From_String(angles)
DIHEDRAL.ProperType.New_From_String(propers)
DIHEDRAL.ImproperType.New_From_String(impropers)
LJ.LJType.New_From_String(LJs)

atoms, bonds, angles, propers, impropers, LJs, cmap = LOAD.frcmod(os.path.join(AMBER_DATA_DIR, "ions234lm_126_spce.frcmod"))
AtomType.New_From_String(atoms)
BOND.BondType.New_From_String(bonds)
ANGLE.AngleType.New_From_String(angles)
DIHEDRAL.ProperType.New_From_String(propers)
DIHEDRAL.ImproperType.New_From_String(impropers)
LJ.LJType.New_From_String(LJs)

LOAD.mol2(os.path.join(AMBER_DATA_DIR, "atomic_ions.mol2"))



sys.modules['__main__'].__dict__["WAT"] = SPCE