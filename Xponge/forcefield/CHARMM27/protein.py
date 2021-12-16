from . import *

import sys

output = LOAD.itp(os.path.join(CHARMM27_DATA_DIR, "forcefield.itp"))

AtomType.New_From_String(output["atomtypes"])
BOND.BondType.New_From_String(output["bonds"])
DIHEDRAL.ProperType.New_From_String(output["dihedrals"])
LJ.LJType.New_From_String(output["LJ"])
UREY_BRADLEY.UreyBradleyType.New_From_String(output["Urey-Bradley"])
IMPROPER.ImproperType.New_From_String(output["impropers"])
NB14_EXTRA.NB14Type.New_From_String(output["nb14_extra"])
NB14.NB14Type.New_From_String(output["nb14"])


test = LOAD.mol2(os.path.join(CHARMM27_DATA_DIR, "test.mol2"))

