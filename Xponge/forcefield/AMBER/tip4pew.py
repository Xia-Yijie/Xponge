from . import *
import os
import sys

AtomType.New_From_String(
"""
name mass    charge[e]      LJtype
HW   1.008   0.52422        HW
OW   16      0              OW
EPW  0      -1.04844        EPW
""")


BOND.BondType.New_From_String(r"""
name    k[kcal/mol·A^-2]   b[A]
OW-HW   553                0.9572
HW-HW   553                1.5136
""")


ANGLE.AngleType.New_From_String(r"""
name        k       b
HW-OW-HW    0       0
OW-HW-HW    0       0
""")


LJ.LJType.New_From_String(r"""
name    epsilon[kcal/mol]   sigma[A]
OW-OW   0.162750   	        3.16435
HW-HW   0                   0
EPW-EPW 0                   0
""")

VIRTUAL_ATOM.VirtualType2.New_From_String(r"""
name   atom0    atom1   atom2   k1         k2
EPW    -3       -2      -1      0.1066413  0.1066413
""")



TIP3P = LOAD.mol2(os.path.join(AMBER_DATA_DIR, "tip4pew.mol2"))

atoms, bonds, angles, propers, impropers, LJs, cmap = LOAD.frcmod(os.path.join(AMBER_DATA_DIR, "ions1lm_126_tip4pew.frcmod"))
AtomType.New_From_String(atoms)
BOND.BondType.New_From_String(bonds)
ANGLE.AngleType.New_From_String(angles)
DIHEDRAL.ProperType.New_From_String(propers)
DIHEDRAL.ImproperType.New_From_String(impropers)
LJ.LJType.New_From_String(LJs)

atoms, bonds, angles, propers, impropers, LJs, cmap = LOAD.frcmod(os.path.join(AMBER_DATA_DIR, "ionsjc_tip4pew.frcmod"))
AtomType.New_From_String(atoms)
BOND.BondType.New_From_String(bonds)
ANGLE.AngleType.New_From_String(angles)
DIHEDRAL.ProperType.New_From_String(propers)
DIHEDRAL.ImproperType.New_From_String(impropers)
LJ.LJType.New_From_String(LJs)

atoms, bonds, angles, propers, impropers, LJs, cmap = LOAD.frcmod(os.path.join(AMBER_DATA_DIR, "ions234lm_126_tip4pew.frcmod"))
AtomType.New_From_String(atoms)
BOND.BondType.New_From_String(bonds)
ANGLE.AngleType.New_From_String(angles)
DIHEDRAL.ProperType.New_From_String(propers)
DIHEDRAL.ImproperType.New_From_String(impropers)
LJ.LJType.New_From_String(LJs)

LOAD.mol2(os.path.join(AMBER_DATA_DIR, "atomic_ions.mol2"))

sys.modules['__main__'].__dict__["WAT"] = TIP3P