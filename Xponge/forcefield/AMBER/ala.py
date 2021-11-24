from . import *
import os
import sys

atoms, bonds, angles, propers, impropers, LJs = LOAD.parmdat(os.path.join(os.path.dirname(__file__), "parm10.dat"))

AtomType.New_From_String(atoms)
BOND.BondType.New_From_String(bonds)
ANGLE.AngleType.New_From_String(angles)
DIHEDRAL.ProperType.New_From_String(propers)
DIHEDRAL.ImproperType.New_From_String(impropers)
LJ.LJType.New_From_String(LJs)


atoms, bonds, angles, propers, impropers, LJs = LOAD.frcmod(os.path.join(os.path.dirname(__file__), "frcmod.ff14SB"))

AtomType.New_From_String(atoms)
BOND.BondType.New_From_String(bonds)
ANGLE.AngleType.New_From_String(angles)
DIHEDRAL.ProperType.New_From_String(propers)
DIHEDRAL.ImproperType.New_From_String(impropers)
LJ.LJType.New_From_String(LJs)

ALA = LOAD.mol2(os.path.join(os.path.dirname(__file__), "ALA.mol2"))


BUILD.Build_Bonded_Force(ALA)
sys.modules['__main__'].__dict__["ALA"] = ALA