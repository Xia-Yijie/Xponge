from ... import Molecule
import os

class Exclude():
    current = None
    def __init__(self, *args, **kwargs):
        if len(args) == 1:
            n = args[0]
        if "n" in kwargs.keys():
            n = kwargs["n"]
        self.n = n
        Exclude.current = self
        @Molecule.Set_Save_SPONGE_Input("exclude")
        def write_exclude(self):
            exclude_numbers = 0
            excludes = []
            
            for atom in self.atoms:
                temp = atom.extra_excluded_atoms.copy()
                atom_self_index = self.atom_index[atom]
                excludes.append(list(map(lambda x: self.atom_index[x], filter(lambda x: self.atom_index[x] > atom_self_index, temp))))
                exclude_numbers += len(excludes[-1])
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
                excludes[-1].sort()
            towrite = "%d %d\n"%(len(self.atoms), exclude_numbers)
            for exclude in excludes:
                exclude.sort()
                towrite += "%d %s\n"%(len(exclude), " ".join([str(atom_index) for atom_index in exclude]))
            
            return towrite
    def Get_Excluded_Atoms(self, molecule):
        temp_dict = {}
        for atom in molecule.atoms:
            temp_dict[atom] = atom.extra_excluded_atoms.copy()
            for i in range(2, self.n + 1):
                for aton in atom.linked_atoms[i]:
                    temp_dict[atom].add(aton)
                        
            if "v" in atom.linked_atoms.keys():
                for aton in atom.linked_atoms["v"]:
                    temp_dict[atom].add(aton)
        return temp_dict

