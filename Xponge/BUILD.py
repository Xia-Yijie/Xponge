from . import *
import sys

def _analyze_connectivity(cls):
    for atom0, c in cls.connectivity.items():
        index_dict = {}.fromkeys(c, atom0)
        for i in range(2, GlobalSetting.farthest_bonded_force + 1):
            index_next = {}
            for atom1, from_atom in index_dict.items():
                atom0.Link_Atom(i, atom1)
                index_temp = {}.fromkeys(cls.connectivity[atom1], atom1)
                index_temp.pop(from_atom)
                index_next.update(index_temp)
            index_dict = index_next


def _build_bfrc(cls):
    _analyze_connectivity(cls)
    for frc in GlobalSetting.BondedForces:
        top = frc.topology_like
        top_matrix = frc.topology_matrix
        frc_all = []
        for atom0 in cls.atoms:
            backups = {i: [] for i in range(len(top))}
            backups[0].append([atom0])
            for i, d in enumerate(top):
                if i == 0:
                    continue
                for atom1 in atom0.linked_atoms[d]:
                    for backup in backups[i - 1]:
                        good_backup = True
                        for j, atomj in enumerate(backup):
                            if atomj == atom1 or abs(top_matrix[j][i]) <= 1 \
                                    or atom1 not in atomj.linked_atoms[abs(top_matrix[j][i])]:
                                good_backup = False
                                break
                            if top_matrix[j][i] <= -1:
                                for d2 in range(2, d):
                                    if atom1 in atomj.linked_atoms[d2]:
                                        good_backup = False
                                        break
                        if good_backup:
                            backups[i].append([*backup, atom1])
            frc_all.extend(backups[len(top) - 1])

        frc_all_final = []
        frc_keys = {}
        for frc_one in frc_all:
            frc_one_name = "".join([str(hash(atom)) for atom in frc_one])
            if frc_one_name in frc_keys.keys():
                frc_keys[frc_one_name].append(frc_one)
            else:
                temp_list = [frc_one]
                frc_all_final.append(temp_list)
                for atom_permutation in frc.Same_Force(frc_one):
                    frc_one_name = "".join([str(hash(atom)) for atom in atom_permutation])
                    frc_keys[frc_one_name] = temp_list

        for frc_ones in frc_all_final:
            finded = {}
            # 先直接找
            for frc_one in frc_ones:
                tofindname = "-".join([atom.type.name for atom in frc_one])
                if tofindname in frc.types.keys():
                    finded[tofindname] = [frc.types[tofindname], frc_one]
                    break
            # 没找到再找通用的
            if not finded:
                leastfindedX = 999
                for frc_one in frc_ones:
                    tofind = [[atom.type.name, "X"] for atom in frc_one]
                    for p in product(*tofind):
                        pcountx = p.count("X")
                        if pcountx > leastfindedX:
                            continue
                        tofindname = "-".join(p)
                        if tofindname in frc.types.keys():
                            finded = {tofindname: [frc.types[tofindname], frc_one]}
                            leastfindedX = pcountx
                            break

            assert (not frc.compulsory or len(finded) == 1), "None of %s type found for %s" % (
                frc.name, "-".join([atom.type.name for atom in frc_one]))

            if finded:
                for finded_type, finded_atoms in finded.values():
                    cls.Add_Bonded_Force(frc.entity(finded_atoms, finded_type))


def _build_bfrc_from_type(cls):
    if not cls.type.built:
        _build_bfrc(cls.type)

    res_type_atom_map = {}
    res_type_atom_map_inverse = {}
    clsatoms = {atom: None for atom in cls.atoms}
    for atom0 in cls.type.atoms:
        for atom in clsatoms.keys():
            if atom0.name == atom.name:
                res_type_atom_map[atom0] = atom
                res_type_atom_map_inverse[atom] = atom0
                clsatoms.pop(atom)
                break

    for atom in cls.atoms:
        atom0 = res_type_atom_map_inverse[atom]
        for key in atom0.linked_atoms.keys():
            for atomi in atom0.linked_atoms[key]:
                atom.Link_Atom(key, res_type_atom_map[atomi])

    for frc in GlobalSetting.BondedForces:
        frc_entities = cls.type.bonded_forces.get(frc.name, [])
        for frc_entity in frc_entities:
            finded_atoms = [res_type_atom_map[atom] for atom in frc_entity.atoms]
            finded_type = frc_entity.type
            cls.Add_Bonded_Force(frc.entity(finded_atoms, finded_type))
            cls.bonded_forces[frc.name][-1].contents = frc_entity.contents


def _build_bfrc_link(cls):
    atom1 = cls.atom1
    atom2 = cls.atom2

    atom1_friends = set([atom1])
    atom2_friends = set([atom2])

    far = GlobalSetting.farthest_bonded_force
    temp_atom1_linked = {i: set() for i in range(far, 2, -1)}
    temp_atom2_linked = {i: set() for i in range(far, 2, -1)}

    for i in range(far - 1, 1, -1):
        for atom in atom1.linked_atoms[i]:
            atom.Link_Atom(i + 1, atom2)
            temp_atom2_linked[i + 1].add(atom)
            atom1_friends.add(atom)
        for atom in atom2.linked_atoms[i]:
            atom.Link_Atom(i + 1, atom1)
            temp_atom1_linked[i + 1].add(atom)
            atom2_friends.add(atom)
    for i in range(far - 1, 1, -1):
        atom1.linked_atoms[i + 1] |= temp_atom1_linked[i + 1]
        atom2.linked_atoms[i + 1] |= temp_atom2_linked[i + 1]

    atom1.Link_Atom(2, atom2)
    atom2.Link_Atom(2, atom1)

    for i in range(2, far):
        for j in range(2, far + 1 - i):
            for atom1_linked_atom in atom1.linked_atoms[i]:
                for atom2_linked_atom in atom2.linked_atoms[j]:
                    if atom1_linked_atom not in atom2_friends and atom2_linked_atom not in atom1_friends:
                        atom1_linked_atom.Link_Atom(i + j, atom2_linked_atom)
                        atom2_linked_atom.Link_Atom(i + j, atom1_linked_atom)

    atom12_friends = atom1_friends | atom2_friends
    for frc in GlobalSetting.BondedForces:
        top = frc.topology_like
        top_matrix = frc.topology_matrix
        frc_all = []

        for atom0 in atom12_friends:
            backups = {i: [] for i in range(len(top))}
            backups[0].append([atom0])
            for i, d in enumerate(top):
                if i == 0:
                    continue
                for _atom1 in atom0.linked_atoms[d]:
                    for backup in backups[i - 1]:
                        good_backup = True
                        for j, atomj in enumerate(backup):
                            if atomj == _atom1 or abs(top_matrix[j][i]) <= 1 or \
                                    _atom1 not in atomj.linked_atoms[abs(top_matrix[j][i])]:
                                good_backup = False
                                break
                            if top_matrix[j][i] <= -1:
                                for d2 in range(2, d):
                                    if _atom1 in atomj.linked_atoms[d2]:
                                        good_backup = False
                                        break
                        if good_backup:
                            backups[i].append([*backup, _atom1])
            for backup in backups[len(top) - 1]:
                backupset = set(backup)
                if atom1_friends & backupset and backupset & atom2_friends:
                    frc_all.append(backup)

        frc_all_final = []
        frc_keys = {}
        for frc_one in frc_all:
            frc_one_name = "".join([str(hash(atom)) for atom in frc_one])
            if frc_one_name in frc_keys.keys():
                frc_keys[frc_one_name].append(frc_one)
            else:
                temp_list = [frc_one]
                frc_all_final.append(temp_list)
                for atom_permutation in frc.Same_Force(frc_one):
                    frc_one_name = "".join([str(hash(atom)) for atom in atom_permutation])
                    frc_keys[frc_one_name] = temp_list

        for frc_ones in frc_all_final:
           
            finded = {}
            # 先直接找
            for frc_one in frc_ones:
                tofindname = "-".join([atom.type.name for atom in frc_one])
                if tofindname in frc.types.keys():
                    finded[tofindname] = [frc.types[tofindname], frc_one]
                    break
            # 没找到再找通用的
            if not finded:
                leastfindedX = 999
                for frc_one in frc_ones:
                    tofind = [[atom.type.name, "X"] for atom in frc_one]
                    for p in product(*tofind):
                        pcountx = p.count("X")
                        if pcountx > leastfindedX:
                            continue
                        tofindname = "-".join(p)
                        if tofindname in frc.types.keys():
                            finded = {tofindname: [frc.types[tofindname], frc_one]}
                            leastfindedX = pcountx
                            break
            
            assert (not frc.compulsory or len(finded) > 0), "None of %s type found for %s" % (
                frc.name, "-".join([atom.type.name for atom in frc_one]))

            if finded:
                for finded_type, finded_atoms in finded.values():
                    cls.Add_Bonded_Force(frc.entity(finded_atoms, finded_type))


def Build_Bonded_Force(cls):
        if cls.built:
            return

        if type(cls) == ResidueType:
            _build_bfrc(cls)
            cls.built = True

        elif type(cls) == Residue:
            _build_bfrc_from_type(cls)
            cls.built = True

        elif type(cls) == ResidueLink:
            _build_bfrc_link(cls)
            cls.built = True

        elif type(cls) == Molecule:
            for res in cls.residues:
                if not res.type.built:
                    Build_Bonded_Force(res.type)
                Build_Bonded_Force(res)
            for link in cls.residue_links:
                Build_Bonded_Force(link)

            cls.atoms = []
            cls.bonded_forces = {frc.name: [] for frc in GlobalSetting.BondedForces}
            for res in cls.residues:
                cls.atoms.extend(res.atoms)
                for frc in GlobalSetting.BondedForces:
                    cls.bonded_forces[frc.name].extend(res.bonded_forces.get(frc.name, []))
            for link in cls.residue_links:
                for frc in GlobalSetting.BondedForces:
                    cls.bonded_forces[frc.name].extend(link.bonded_forces.get(frc.name, []))
            cls.atom_index = {cls.atoms[i]: i for i in range(len(cls.atoms))}

            for vatom_type_name, vatom_type_atom_numbers in GlobalSetting.VirtualAtomTypes.items():
                for vatom in cls.bonded_forces.get(vatom_type_name, []):
                    this_vatoms = [vatom.atoms[0]]
                    for i in range(vatom_type_atom_numbers):
                        this_vatoms.append(cls.atoms[cls.atom_index[vatom.atoms[0]] + getattr(vatom, "atom%d" % i)])
                    this_vatoms.sort(key=lambda x: cls.atom_index[x])
                    while this_vatoms:
                        tolink = this_vatoms.pop(0)
                        for i in this_vatoms:
                            tolink.Link_Atom("v", i)

            cls.built = True
        else:
            raise NotImplementedError


def Save_SPONGE_Input(molecule, prefix=None, dirname="."):
    if type(molecule) == Molecule:
        Build_Bonded_Force(molecule)

        if not prefix:
            prefix = molecule.name

        for key, func in Molecule.save_functions.items():
            towrite = func(molecule)
            if towrite:
                f = open(os.path.join(dirname, prefix + "_" + key + ".txt"), "w")
                f.write(towrite)
                f.close()

    elif type(molecule) == Residue:
        mol = Molecule(name=molecule.name)
        mol.Add_Residue(molecule)
        Save_SPONGE_Input(mol, prefix, dirname)
    elif type(molecule) == ResidueType:
        residue = Residue(molecule, name=molecule.name)
        for atom in molecule.atoms:
            residue.Add_Atom(atom)
        Save_SPONGE_Input(residue, prefix, dirname)


sys.modules['__main__'].__dict__["Save_SPONGE_Input"] = Save_SPONGE_Input


def Save_PDB(molecule, filename=None):
    if type(molecule) == Molecule:

        molecule.atoms = []
        for res in molecule.residues:
            molecule.atoms.extend(res.atoms)

        molecule.atom_index = {molecule.atoms[i]: i for i in range(len(molecule.atoms))}
        molecule.residue_index = {molecule.residues[i]: i for i in range(len(molecule.residues))}
        molecule.link_to_next = [False for res in molecule.residues]
        for link in molecule.residue_links:
            if molecule.residue_index[link.atom1.residue] - molecule.residue_index[link.atom2.residue] == 1:
                molecule.link_to_next[molecule.residue_index[link.atom2.residue]] = True
            elif molecule.residue_index[link.atom2.residue] - molecule.residue_index[link.atom1.residue] == 1:
                molecule.link_to_next[molecule.residue_index[link.atom1.residue]] = True

        towrite = "REMARK   Generated By Xponge (Molecule)\n"

        chain_atom0 = -1
        chain_residue0 = -1
        real_chain_residue0 = -1
        chain_counter = 0
        for atom in molecule.atoms:
            resname = atom.residue.name
            if resname in GlobalSetting.PDBResidueNameMap["save"].keys():
                resname = GlobalSetting.PDBResidueNameMap["save"][resname]
            towrite += "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%17s%2s\n" % (
                molecule.atom_index[atom] - chain_atom0, atom.name,
                resname, " ", (molecule.residue_index[atom.residue] - chain_residue0)%10000,
                atom.x, atom.y, atom.z, " ", " ")
            if atom == atom.residue.atoms[-1] and molecule.link_to_next[molecule.residue_index[atom.residue]] == False:
                towrite += "TER\n"
                chain_atom0 = molecule.atom_index[atom]
                if molecule.residue_index[atom.residue] - real_chain_residue0 != 1:
                    chain_residue0 = molecule.residue_index[atom.residue]
                    real_chain_residue0 = chain_residue0
                else:
                    real_chain_residue0 = molecule.residue_index[atom.residue]
        if not filename:
            filename = molecule.name + ".pdb"

        f = open(filename, "w")
        f.write(towrite)
        f.close()
    elif type(molecule) == Residue:
        mol = Molecule(name=molecule.name)
        mol.Add_Residue(molecule)
        Save_PDB(mol, filename)
    elif type(molecule) == ResidueType:
        residue = Residue(molecule, name=molecule.name)
        for atom in molecule.atoms:
            residue.Add_Atom(atom)
        Save_PDB(residue, filename)
    elif type(molecule) == assign.Assign:
        molecule.Save_As_PDB(filename)
    else:
        raise NotImplementedError


sys.modules['__main__'].__dict__["Save_PDB"] = Save_PDB


def Save_Mol2(molecule, filename=None):
    if type(molecule) == Molecule:
        molecule.atoms = []
        for res in molecule.residues:
            molecule.atoms.extend(res.atoms)
        molecule.atom_index = {molecule.atoms[i]: i for i in range(len(molecule.atoms))}
        bonds = []
        for res in molecule.residues:
            for atom1, atom1_con in res.type.connectivity.items():
                atom1_index = molecule.atom_index[res._name2atom[atom1.name]] + 1
                for atom2 in atom1_con:
                    atom2_index = molecule.atom_index[res._name2atom[atom2.name]] + 1
                    if atom1_index < atom2_index:
                        bonds.append("%6d %6d" % (atom1_index, atom2_index))

        for link in molecule.residue_links:
            atom1_index = molecule.atom_index[link.atom1] + 1
            atom2_index = molecule.atom_index[link.atom2] + 1
            if atom1_index < atom2_index:
                bonds.append("%6d %6d" % (atom1_index, atom2_index))
            else:
                bonds.append("%6d %6d" % (atom2_index, atom1_index))
        bonds.sort(key=lambda x: list(map(int, x.split())))
        towrite = "@<TRIPOS>MOLECULE\n"
        towrite += "%s\n" % (molecule.name)
        towrite += " %d %d %d 0 1\n" % (len(molecule.atoms), len(bonds), len(molecule.residues))
        towrite += "SMALL\n"
        towrite += "USER_CHARGES\n"

        towrite += "@<TRIPOS>ATOM\n"
        count = 0
        res_count = 0
        residue_start = []
        for atom in molecule.atoms:
            count += 1
            if atom == atom.residue.atoms[0]:
                res_count += 1
                residue_start.append(count)
            resname = atom.residue.name
            towrite += "%6d %4s %8.3f %8.3f %8.3f %4s %5d %8s %10.6f\n" % (
                count, atom.name, atom.x, atom.y, atom.z, atom.type.name, res_count, resname, atom.charge)

        towrite += "@<TRIPOS>BOND\n"
        for i, bond in enumerate(bonds):
            towrite += "%6d %s 1\n" % (i + 1, bond)
        towrite += "@<TRIPOS>SUBSTRUCTURE\n"
        for i, residue in enumerate(molecule.residues):
            towrite += "%5d %8s %6d ****               0 ****  **** \n" % (i + 1, residue.name, residue_start[i])

        if not filename:
            filename = molecule.name + ".mol2"

        f = open(filename, "w")
        f.write(towrite)
        f.close()
    elif type(molecule) == Residue:
        mol = Molecule(name=molecule.name)
        mol.Add_Residue(molecule)
        Save_Mol2(mol, filename)
    elif type(molecule) == ResidueType:
        residue = Residue(molecule, name=molecule.name)
        for atom in molecule.atoms:
            residue.Add_Atom(atom)
        Save_Mol2(residue, filename)
    elif type(molecule) == assign.Assign:
        molecule.Save_As_Mol2(filename)
    else:
        raise NotImplementedError


sys.modules['__main__'].__dict__["Save_Mol2"] = Save_Mol2

def Save_Gro(molecule, filename):
    towrite = "Generated By Xponge\n"
    molecule.atoms = []
    for res in molecule.residues:
        molecule.atoms.extend(res.atoms)
    molecule.residue_index = {molecule.residues[i]: i for i in range(len(molecule.residues))}
    
    boxlength = [0, 0, 0]
    maxi = [-float("inf"), -float("inf"), -float("inf")]
    mini = [float("inf"), float("inf"), float("inf")]
    for atom in molecule.atoms:
        if atom.x > maxi[0]:
            maxi[0] = atom.x
        if atom.y > maxi[1]:
            maxi[1] = atom.y
        if atom.z > maxi[2]:
            maxi[2] = atom.z
        if atom.x < mini[0]:
            mini[0] = atom.x
        if atom.y < mini[1]:
            mini[1] = atom.y
        if atom.z < mini[2]:
            mini[2] = atom.z

    towrite += "%d\n"%len(molecule.atoms)
    for i, atom in enumerate(molecule.atoms):
        residue = atom.residue
        if not GlobalSetting.nocenter:
            x = atom.x - mini[0] + GlobalSetting.boxspace
            y = atom.y - mini[1] + GlobalSetting.boxspace
            z = atom.z - mini[2] + GlobalSetting.boxspace
        else:
            x, y, z = atom.x. atom.y, atom.z
         
        towrite += "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"%(molecule.residue_index[residue] + 1, residue.name, atom.name, i+1, x/10, y/10, z/10)
    if molecule.box_length is None:
        boxlength[0] = maxi[0] - mini[0] + GlobalSetting.boxspace * 2
        boxlength[1] = maxi[1] - mini[1] + GlobalSetting.boxspace * 2
        boxlength[2] = maxi[2] - mini[2] + GlobalSetting.boxspace * 2
        molecule.box_length = [boxlength[0], boxlength[1], boxlength[2]]
    else:
        boxlength[0] = molecule.box_length[0]
        boxlength[1] = molecule.box_length[1]
        boxlength[2] = molecule.box_length[2]
    towrite += "%8.3f %8.3f %8.3f"%(molecule.box_length[0] /10, molecule.box_length[1]/10, molecule.box_length[2]/10)
    with open(filename, "w") as f:
        f.write(towrite)
    
sys.modules['__main__'].__dict__["Save_Gro"] = Save_Gro

def Save_NPZ(molecule, filename=None):
    import numpy as np
    bondtypeindex = {}

    def update_bonddict(dic_to, dic_from):
        for key, value in dic_from.items():
            if value is not None:
                if key in dic_to.keys():
                    dic_to[key].append(value)
                else:
                    dic_to[key] = [value]

    for bondedforce in GlobalSetting.BondedForces:
        bonddict = {}
        count = -1
        for bondtype in bondedforce.types.values():
            update_bonddict(bonddict, bondtype.contents)
            count += 1
            bondtypeindex[bondtype] = count
        max_len = {}
        for key, value in bonddict.items():
            max_len[key] = 1
            for vi in value:
                if isinstance(vi, list) and len(vi) > max_len[key]:
                    max_len[key] = len(vi)
        for key, value in bonddict.items():
            for i, vi in enumerate(value):
                if isinstance(vi, list) and len(vi) < max_len[key]:
                    bonddict[key][i] = np.array(vi + [0] * (max_len[key] - len(vi)))
        np.savez(bondedforce.name, **bonddict)

    if type(molecule) == Molecule:
        Build_Bonded_Force(molecule)

        if not filename:
            filename = molecule.name

        molecule.atoms = []
        molecule.bonded_forces = {frc.name: [] for frc in GlobalSetting.BondedForces}
        for res in molecule.residues:
            molecule.atoms.extend(res.atoms)
            for frc in GlobalSetting.BondedForces:
                molecule.bonded_forces[frc.name].extend(res.bonded_forces[frc.name])

        for link in molecule.residue_links:
            for frc in GlobalSetting.BondedForces:
                molecule.bonded_forces[frc.name].extend(link.bonded_forces[frc.name])

        molecule.atom_index = {molecule.atoms[i]: i for i in range(len(molecule.atoms))}

        bonddict = {}
        for bondname, bondtypes in molecule.bonded_forces.items():
            bonddict[bondname] = []
            for b in bondtypes:
                bonddict[bondname].append([molecule.atom_index[a] for a in b.atoms])
                bonddict[bondname][-1].append(bondtypeindex[b.type])
        np.savez(filename, **bonddict)
    elif type(molecule) == Residue:
        mol = Molecule(name=molecule.name)
        mol.Add_Residue(molecule)
        Save_NPZ(mol, filename)
    elif type(molecule) == ResidueType:
        residue = Residue(molecule, name=molecule.name)
        for atom in molecule.atoms:
            residue.Add_Atom(atom)
        Save_NPZ(residue, filename)


sys.modules['__main__'].__dict__["Save_NPZ"] = Save_NPZ
