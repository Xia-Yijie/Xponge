def Assign2RDKitMol(assign, ignore_bond_type = False):
    from rdkit import Chem
    molA = Chem.RWMol()
    for atom in assign.atoms:
        molA.AddAtom(Chem.Atom(atom))
    for atom, bonds in assign.bonds.items():
        for aton, n in bonds.items():
            if aton < atom:
                continue
            if  ignore_bond_type or n == -1:
                temp_bond = Chem.BondType.UNSPECIFIED
            elif n == 1:
                temp_bond = Chem.BondType.SINGLE
            elif n == 2:
                temp_bond = Chem.BondType.DOUBLE
            elif n == 3:
                temp_bond = Chem.BondType.TRIPLE
            else:
                raise NotImplementedError
            molA.AddBond(atom, aton, temp_bond)
    conf = Chem.Conformer(assign.atom_numbers)
    for i in range(assign.atom_numbers):
        conf.SetAtomPosition(i, assign.coordinate[i])
    mol = molA.GetMol()
    mol.AddConformer(conf)
    Chem.SanitizeMol(mol)
    return mol

def Insert_Atom_Type_To_RDKitMol(mol, res, assign, atom_type_dict = None):
    i = 0
    if atom_type_dict is None:
        atom_type_dict = {}
    for a in mol.GetAtoms():
        atom_type = res._name2atom[assign.names[i]].type.name
        if atom_type not in atom_type_dict.keys():
            atom_type_dict[atom_type] = len(atom_type_dict)
        a.SetIsotope(atom_type_dict[atom_type])
        i += 1

def Find_Equal_Atoms(assign):
    from rdkit import Chem
    mols = []
    CanonSmiles = []
    for i in range(len(assign.atoms)):
        mols.append(Assign2RDKitMol(assign))
        mols[-1].GetAtoms()[i].SetIsotope(1)
        CanonSmiles.append(Chem.MolToSmiles(mols[-1]))
    group = {i:i for i in range(len(assign.atoms))}
    for i in range(len(assign.atoms)):
        if group[i] == i:
            for j in range(i+1, len(assign.atoms)):
                if CanonSmiles[j] == CanonSmiles[i]:
                    group[j] = i
    ret = []
    realmap = {}
    for i in group:
        if group[i] == i:
            ret.append([i])
            realmap[i] = len(realmap)
        else:
            ret[realmap[group[i]]].append(i)
    return list(filter(lambda x:len(x) > 1, ret))

def Get_Conformer_Coordinate(mol):
    import numpy as np
    crd = np.zeros((mol.GetNumAtoms(), 3))
    for i in range(mol.GetNumAtoms()):
        crd[i] = mol.GetConformer(0).GetAtomPosition(i)
    return crd

def Get_Conformer_Coordinate_To_Residue(mol, res, assign):
    for i in range(mol.GetNumAtoms()):
        j = res._name2index[assign.names[i]]
        res.atoms[j].x, res.atoms[j].y, res.atoms[j].z = mol.GetConformer(0).GetAtomPosition(i)

def Set_Conformer_Coordinate_From_Residue(mol, res, assign):
    for i in range(mol.GetNumAtoms()):    
        j = res._name2index[assign.names[i]]
        mol.GetConformer(0).SetAtomPosition(i, [res.atoms[j].x, res.atoms[j].y, res.atoms[j].z])

def Set_Conformer_Coordinate(mol, crd):
    for i in range(mol.GetNumAtoms()):
        mol.GetConformer(0).SetAtomPosition(i, crd[i])

def Apply_Transform(crd, matrix):
    import numpy as np
    temp = np.dot(crd, matrix[:3, :3]) + matrix[:3,3]
    return temp[:, :3]

def Get_Part_Align(molA, molB, partA, partB):
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    import numpy as np
    RWmolA = Chem.RWMol(molA)
    RWmolB = Chem.RWMol(molB)
    for i in range(molA.GetNumAtoms())[::-1]:
        if i not in partA:
            RWmolA.RemoveAtom(i)
    for i in range(molB.GetNumAtoms())[::-1]:
        if i not in partB:
            RWmolB.RemoveAtom(i)

    molAn = RWmolA.GetMol()
    molBn = RWmolB.GetMol()
    rmsd, transformer = rdMolAlign.GetAlignmentTransform(molBn, molAn)
    crd = Apply_Transform(Get_Conformer_Coordinate(molB), transformer)
    Set_Conformer_Coordinate(molB, crd)



