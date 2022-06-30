from ... import Molecule


class Exclude:
    current = None

    def __init__(self, *args, **kwargs):
        n = 4
        if len(args) == 1:
            n = args[0]
        if "n" in kwargs.keys():
            n = kwargs["n"]
        self.n = n
        Exclude.current = self

        @Molecule.Set_Save_SPONGE_Input("exclude")
        def write_exclude(mol):
            exclude_numbers = 0
            excludes = []
            for atom in mol.atoms:
                temp = atom.extra_excluded_atoms.copy()
                atom_self_index = mol.atom_index[atom]
                excludes.append(
                    list(map(lambda x: mol.atom_index[x], filter(lambda x: mol.atom_index[x] > atom_self_index, temp))))
                exclude_numbers += len(excludes[-1])
                for i in range(2, n + 1):
                    for aton in atom.linked_atoms.get(i, []):
                        if mol.atom_index[aton] > atom_self_index and aton not in temp:
                            temp.add(aton)
                            exclude_numbers += 1
                            excludes[-1].append(mol.atom_index[aton])
                for aton in atom.linked_atoms.get("v", []):
                    if mol.atom_index[aton] > atom_self_index and aton not in temp:
                        temp.add(aton)
                        exclude_numbers += 1
                        excludes[-1].append(mol.atom_index[aton])
                excludes[-1].sort()
            towrite = "%d %d\n" % (len(mol.atoms), exclude_numbers)
            for exclude in excludes:
                exclude.sort()
                towrite += "%d %s\n" % (len(exclude), " ".join([str(atom_index) for atom_index in exclude]))

            return towrite

    def Get_Excluded_Atoms(self, molecule):
        temp_dict = {}
        for atom in molecule.atoms:
            temp_dict[atom] = atom.extra_excluded_atoms.copy()
            for i in range(2, self.n + 1):
                for aton in atom.linked_atoms.get(i, []):
                    temp_dict[atom].add(aton)
            for aton in atom.linked_atoms.get("v", []):
                temp_dict[atom].add(aton)
        return temp_dict
