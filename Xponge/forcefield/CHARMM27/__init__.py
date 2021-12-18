from ... import *
from ..BASE import CHARGE, MASS, LJ, BOND, DIHEDRAL, NB14, NB14_EXTRA, UREY_BRADLEY, IMPROPER, VIRTUAL_ATOM, ACMAP
import os

CHARMM27_DATA_DIR = os.path.dirname(__file__)

LJ.LJType.combining_method_A = LJ.Lorentz_Berthelot_For_A
LJ.LJType.combining_method_B = LJ.Lorents_Berthelot_For_B

GlobalSetting.Set_Invisible_Bonded_Forces(["improper"])

@Molecule.Set_Save_SPONGE_Input
def write_exclude(self, prefix, dirname):
    exclude_numbers = 0
    excludes = []
    
    for atom in self.atoms:
        temp = set()
        atom_self_index = self.atom_index[atom]
        excludes.append([])
        for i in range(2, 5):
            for aton in atom.linked_atoms[i]:
                if self.atom_index[aton] > atom_self_index and aton not in temp:
                    temp.add(aton)
                    exclude_numbers += 1
                    excludes[-1].append(self.atom_index[aton])
                    
        if "v" in atom.linked_atoms.keys():
            for aton in atom.linked_atoms["v"]:
                if self.atom_index[aton] > atom_self_index and aton not in temp:
                    temp.add(aton)
                    exclude_numbers += 1
                    excludes[-1].append(self.atom_index[aton])
    
    towrite = "%d %d\n"%(len(self.atoms), exclude_numbers)
    for exclude in excludes:
        exclude.sort()
        towrite += "%d %s\n"%(len(exclude), " ".join([str(atom_index) for atom_index in exclude]))
    
    f = open(os.path.join(dirname, prefix + "_exclude.txt"),"w")
    f.write(towrite)
    f.close()


