from . import *
from time import time
import sys

def _build_bfrc(cls):
    #t = time()
    for atom0, c in cls.connectivity.items():
        index_dict = {}.fromkeys(c, atom0)
        for i in range(2,GlobalSetting.farthest_bonded_force+1):
            index_next = {}
            for atom1, from_atom in index_dict.items():
                atom0.linked_atoms[i].append(atom1)
                index_temp = {}.fromkeys(cls.connectivity[atom1], atom1)
                index_temp.pop(from_atom)
                index_next.update(index_temp)
            index_dict = index_next
    #print("analysis of connectivity: %f"%(time()-t))
    for frc in GlobalSetting.BondedForces:
        top = frc.topology_like
        top_matrix = frc.topology_matrix
        frc_all = []
        #t = time()
        for atom0 in cls.atoms:
            backups = {i:[] for i in range(len(top))}
            backups[0].append([atom0])
            for i,d in enumerate(top):
                if i == 0:
                    continue
                for atom1 in atom0.linked_atoms[d]:
                    for backup in backups[i-1]:
                        good_backup = True
                        for j, atomj in enumerate(backup):
                            if atomj == atom1 or top_matrix[j][i] <= 1 or atom1 not in atomj.linked_atoms[top_matrix[j][i]]:
                                good_backup = False
                                break
                        if  good_backup:
                            backups[i].append([*backup,atom1])
            frc_all.extend(backups[len(top)-1])
        #print(frc.name, "analysis of links:%f"%(time()-t))  
        #t = time()
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
                    
        #print(frc.name, "analysis of same force: %f"%(time()-t))  
        #t = time()
        for frc_ones in frc_all_final:
            finded = {}
            #先直接找
            for frc_one in frc_ones:
                tofindname = "-".join([atom.type.name for atom in frc_one])
                if tofindname in frc.types.keys():
                    finded[tofindname] = [frc.types[tofindname], frc_one]
                    break
            #没找到再找通用的
            if not finded:
                for frc_one in frc_ones:
                    tofind = [[atom.type.name, "X"] for atom in frc_one]
                    for p in product(*tofind):
                        tofindname = "-".join(p) 
                        if tofindname in frc.types.keys():
                            finded[tofindname] = [frc.types[tofindname], frc_one]
                            break                         
            
            assert (not frc.compulsory or len(finded) == 1), "None of %s type found for %s"%(frc.name, "-".join([atom.type.name for atom in frc_one]))
            
            if finded:
                for finded_type, finded_atoms in finded.values():
                    cls.Add_Bonded_Force(frc.entity(finded_atoms, finded_type))
                
        #print(frc.name, "analysis of force type: %f"%(time()-t))

def _build_bfrc_from_type(cls):
    if not cls.type.builded:
        _build_bfrc(cls.type)
    
    res_type_atom_map = {}
    res_type_atom_map_inverse = {}
    clsatoms = {atom:None for atom in cls.atoms}
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
            atom.linked_atoms[key] = []
            for i in range(len(atom0.linked_atoms[key])):
                atom.linked_atoms[key].append(res_type_atom_map[atom0.linked_atoms[key][i]])
        
    for frc in GlobalSetting.BondedForces:
        frc_entities = cls.type.bonded_forces[frc.name]
        for frc_entity in frc_entities:
            finded_atoms = [res_type_atom_map[atom] for atom in frc_entity.atoms]
            finded_type = frc_entity.type
            cls.Add_Bonded_Force(frc.entity(finded_atoms, finded_type))
        
def _build_bfrc_link(cls):
    atom1 = cls.atom1
    atom2 = cls.atom2

    atom1_friends = set([atom1])
    atom2_friends = set([atom2])
    
    far = GlobalSetting.farthest_bonded_force
    #t = time()
    for i in range(far-1, 1, -1):
        for atom in atom1.linked_atoms[i]:
            atom.linked_atoms[i+1].append(atom2)
            atom2.linked_atoms[i+1].append(atom)
            atom1_friends.add(atom)
        for atom in atom2.linked_atoms[i]:
            atom.linked_atoms[i+1].append(atom1)
            atom1.linked_atoms[i+1].append(atom)
            atom2_friends.add(atom)
    atom1.linked_atoms[2].append(atom2)
    atom2.linked_atoms[2].append(atom1)
    for i in range(2, far):
        for j in range(2, far + 1 - i):
            for atom1_linked_atom in atom1.linked_atoms[i]:
                for atom2_linked_atom in atom2.linked_atoms[j]:
                    atom1_linked_atom.linked_atoms[i+j].append(atom2_linked_atom)
                    atom2_linked_atom.linked_atoms[j+i].append(atom1_linked_atom)
    #print("analysis of connectivity: %f"%(time()-t))
    atom12_friends = atom1_friends | atom2_friends
    for frc in GlobalSetting.BondedForces:
        top = frc.topology_like
        top_matrix = frc.topology_matrix
        frc_all = []
        #t = time()
        for atom0 in atom12_friends:
            backups = {i:[] for i in range(len(top))}
            backups[0].append([atom0])
            for i,d in enumerate(top):
                if i == 0:
                    continue
                for atom1 in atom0.linked_atoms[d]:
                    for backup in backups[i-1]:
                        good_backup = True
                        for j, atomj in enumerate(backup):
                            if atomj == atom1 or atomj not in atom12_friends or top_matrix[j][i] <= 1 or atom1 not in atomj.linked_atoms[top_matrix[j][i]]:
                                good_backup = False
                                break
                        if  good_backup:
                            backups[i].append([*backup,atom1])
            for backup in backups[len(top)-1]:
                backupset = set(backup)
                if atom1_friends & backupset and backupset & atom2_friends:
                    frc_all.append(backup)
        #print(frc.name, "analysis of links:%f"%(time()-t))  
        #t = time()
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
                    
        #print(frc.name, "analysis of same force: %f"%(time()-t))  
        #t = time()
        for frc_ones in frc_all_final:
            finded = {}
            #先直接找
            for frc_one in frc_ones:
                tofindname = "-".join([atom.type.name for atom in frc_one])
                if tofindname in frc.types.keys():
                    finded[tofindname] = [frc.types[tofindname], frc_one]
                    break
            #没找到再找通用的
            if not finded:
                for frc_one in frc_ones:
                    tofind = [[atom.type.name, "X"] for atom in frc_one]
                    for p in product(*tofind):
                        tofindname = "-".join(p) 
                        if tofindname in frc.types.keys():
                            finded[tofindname] = [frc.types[tofindname], frc_one]
                            break                    

            assert (not frc.compulsory or len(finded) == 1), "None of %s type found for %s"%(frc.name, "-".join([atom.type.name for atom in frc_one]))
            
            if finded:
                for finded_type, finded_atoms in finded.values():
                    cls.Add_Bonded_Force(frc.entity(finded_atoms, finded_type))
        #print(frc.name, "analysis of force type: %f"%(time()-t))
    
def Build_Bonded_Force(cls):
    if cls.builded:
        return

    if type(cls) == ResidueType:
        _build_bfrc(cls)
        cls.builded = True  
        
    elif type(cls) == Residue:
        _build_bfrc_from_type(cls)
        cls.builded = True
        
    elif type(cls) == ResidueLink:
        _build_bfrc_link(cls)
        cls.builded = True
    elif type(cls) == Molecule:
        for res in cls.residues:
            if not res.type.builded:
                Build_Bonded_Force(res.type)
            Build_Bonded_Force(res)
        for link in cls.residue_links:
            Build_Bonded_Force(link)
            cls.builded = True
    
        
        
 
    else:
        raise NotImplementedError




def Save_SPONGE_Input(molecule, prefix = None, dirname = "."):
    if type(molecule)== Molecule:
        Build_Bonded_Force(molecule)
        
        if not prefix:
            prefix = molecule.name
        
        molecule.atoms = []
        molecule.bonded_forces = { frc.name:[] for frc in GlobalSetting.BondedForces}
        for res in molecule.residues:
            molecule.atoms.extend(res.atoms)
            for frc in GlobalSetting.BondedForces:
                molecule.bonded_forces[frc.name].extend(res.bonded_forces[frc.name])
                
        for link in molecule.residue_links:
            for frc in GlobalSetting.BondedForces:
                molecule.bonded_forces[frc.name].extend(link.bonded_forces[frc.name])
                
        molecule.atom_index = { molecule.atoms[i]: i for i in range(len(molecule.atoms))}
        
        for vatom_type_name, vatom_type_atom_numbers in GlobalSetting.VirtualAtomTypes.items():
            for vatom in molecule.bonded_forces[vatom_type_name]:
                this_vatoms = [vatom.atoms[0]]
                for i in range(vatom_type_atom_numbers):
                    this_vatoms.append(molecule.atoms[molecule.atom_index[vatom.atoms[0]] + getattr(vatom, "atom%d"%i)]) 
                this_vatoms.sort(key = lambda x:molecule.atom_index[x])
                while this_vatoms:
                    tolink = this_vatoms.pop(0)
                    for i in this_vatoms:
                        if "v" not in tolink.linked_atoms.keys():
                            tolink.linked_atoms["v"] = []
                        tolink.linked_atoms["v"].append(i) 
                        
        for func in Molecule.save_functions:
            func(molecule, prefix, dirname)
    elif type(molecule) == Residue:
        mol = Molecule(name = molecule.name)
        mol.Add_Residue(molecule)
        Save_SPONGE_Input(mol)
    elif type(molecule) == ResidueType:
        residue = Residue(molecule, name = molecule.name)
        for atom in molecule.atoms:
            residue.Add_Atom(atom.name, x = atom.x, y = atom.y, z = atom.z)
        Save_SPONGE_Input(residue)

sys.modules['__main__'].__dict__["Save_SPONGE_Input"] = Save_SPONGE_Input 