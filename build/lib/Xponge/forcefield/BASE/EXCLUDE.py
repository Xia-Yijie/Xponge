from ... import Molecule
import os

class Exclude():
    def __init__(self, *args, **kwargs):
        if len(args) == 1:
            n = args[0]
        if "n" in kwargs.keys():
            n = kwargs["n"]
        @Molecule.Set_Save_SPONGE_Input
        def write_exclude(self, prefix, dirname):
            exclude_numbers = 0
            excludes = []
            
            for atom in self.atoms:
                temp = set()
                atom_self_index = self.atom_index[atom]
                excludes.append([])
                for i in range(2, n + 1):
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


