from . import *

def _get_rotate_matrix(r0, angle):
    cost = np.cos(angle)
    cost_one = 1 - cost
    sint = np.sin(angle)
    r0 /= np.linalg.norm(r0)    
    return np.array([[r0[0] * r0[0] * cost_one + cost, r0[0] * r0[1] * cost_one - r0[2] * sint, r0[0] * r0[2] * cost_one + r0[1] * sint],
                     [r0[0] * r0[1] * cost_one + r0[2] * sint, r0[1] * r0[1] * cost_one + cost, r0[1] * r0[2] * cost_one - r0[0] * sint],
                     [r0[0] * r0[2] * cost_one - r0[1] * sint, r0[1] * r0[2] * cost_one + r0[0] * sint, r0[2] * r0[2] * cost_one + cost]]).transpose()

def _get_crd(molecule):
    molecule.atoms = []
    for res in molecule.residues:
        molecule.atoms.extend(res.atoms)
            
    molecule.atom_index = { molecule.atoms[i]: i for i in range(len(molecule.atoms))}
    return np.array([[atom.x, atom.y, atom.z] for atom in molecule.atoms])

def _get_friends(molecule, atom1, atom2):    
    if atom1.residue != atom2.residue:
        if molecule.atom_index[atom1.residue.atoms[0]] < molecule.atom_index[atom1.residue.atoms[1]]:
            res_index = molecule.atom_index[atom1.residue.atoms[-1]]
            atom1_friends = list(range(res_index+1))
            atom2_friends = list(range(res_index+1, len(molecule.atoms)))
        else:
            res_index = molecule.atom_index[atom2.residue.atoms[-1]]
            atom2_friends = list(range(res_index+1))
            atom1_friends = list(range(res_index+1, len(molecule.atoms)))
    else:
        link_front = 0
        atom1_friends = []
        atom2_friends = []
        head = 0
        tail = 0
        
        restype = atom1.residue.type
        def restype_atom(atom):
            return restype._name2atom[atom.name]
        
        def resatom(atom):
            return atom1.residue._name2atom[atom.name]
        
        typeatom1 = restype_atom(atom1)
        typeatom2 = restype_atom(atom2)
        assert typeatom2 in restype.connectivity[typeatom1] and typeatom1 in restype.connectivity[typeatom2]
        index_dict = {}.fromkeys(restype.connectivity[typeatom1], typeatom1)
        if typeatom2 in index_dict.keys():
            index_dict.pop(typeatom2)
        while index_dict:
            index_next = {}
            for atom0, from_atom in index_dict.items():
                if atom0.name == restype.head:
                    head = 1
                elif atom0.name == restype.tail:
                    tail = 1
                atom1_friends.append(molecule.atom_index[resatom(atom0)])
                index_temp = {}.fromkeys(restype.connectivity[atom0], atom0)
                index_temp.pop(from_atom)
                if typeatom2 in index_temp.keys():
                    index_temp.pop(typeatom2)
                index_next.update(index_temp)
            index_dict = index_next
            
        
        index_dict = {}.fromkeys(restype.connectivity[typeatom2], typeatom2)
        if typeatom1 in index_dict.keys():
            index_dict.pop(typeatom1)

        while index_dict:
            index_next = {}
            for atom0, from_atom in index_dict.items():
                if atom0.name == restype.head:
                    head = 2
                elif atom0.name == restype.tail:
                    tail = 2
                atom2_friends.append(molecule.atom_index[resatom(atom0)])
                index_temp = {}.fromkeys(restype.connectivity[atom0], atom0)
                index_temp.pop(from_atom)
                if typeatom1 in index_temp.keys():
                    index_temp.pop(typeatom1)
                index_next.update(index_temp)
            index_dict = index_next
        
        if atom1.name == restype.head:
            head = 1
        elif atom1.name == restype.tail:
            tail = 1
        if atom2.name == restype.head:
            head = 2
        elif atom2.name == restype.tail:
            tail = 2
        atom1_friends = set(atom1_friends)
        atom1_friends.add(molecule.atom_index[atom1])
        atom1_friends = np.array(list(atom1_friends))
        atom2_friends = set(atom2_friends)
        atom2_friends.add(molecule.atom_index[atom2])
        atom2_friends = np.array(list(atom2_friends))
    return  atom1_friends, atom2_friends
        
    

def Impose_Bond(molecule, atom1, atom2, length):
    crd = _get_crd(molecule)
    atom1_friends, atom2_friends = _get_friends(molecule, atom1, atom2)
    r0 = crd[molecule.atom_index[atom2]] - crd[molecule.atom_index[atom1]]
    L0 = np.linalg.norm(r0)
    if L0 == 0:
        crd[molecule.atom_index[atom2]] += (1/3)**(0.5)
        r0 = crd[molecule.atom_index[atom2]] - crd[molecule.atom_index[atom1]]
        L0 = np.linalg.norm(r0)
    dr = (L0 - length) * r0
    crd[atom2_friends] += dr
    for atom in molecule.atoms:
        i = molecule.atom_index[atom]
        atom.x = crd[i][0]
        atom.y = crd[i][1]
        atom.z = crd[i][2]

sys.modules['__main__'].__dict__["Impose_Bond"] = Impose_Bond 

def Impose_Angle(molecule, atom1, atom2, atom3, angle):
    crd = _get_crd(molecule)
    atom2_friends, atom3_friends = _get_friends(molecule, atom1, atom2)
    r12 = crd[molecule.atom_index[atom1]] - crd[molecule.atom_index[atom2]]
    r23 = crd[molecule.atom_index[atom3]] - crd[molecule.atom_index[atom2]]
    angle0 = np.arccos(np.dot(r12, r23)  / np.linalg.norm(r23) / np.linalg.norm(r12))
    deltaAngle = angle - angle0
    crd[atom3_friends] = np.dot(crd[atom3_friends] - crd[molecule.atom_index[atom2]], _get_rotate_matrix(np.cross(r12, r23), deltaAngle)) + crd[molecule.atom_index[atom2]]
        
sys.modules['__main__'].__dict__["Impose_Angle"] = Impose_Angle 

def Impose_Dihedral(molecule, atom1, atom2, atom3, atom4, dihedral):
    crd = _get_crd(molecule)
    atom3_friends, atom4_friends = _get_friends(molecule, atom2, atom3)
    r12 = crd[molecule.atom_index[atom1]] - crd[molecule.atom_index[atom2]]
    r23 = crd[molecule.atom_index[atom2]] - crd[molecule.atom_index[atom3]]
    r34 = crd[molecule.atom_index[atom3]] - crd[molecule.atom_index[atom4]]
    r12xr23 = np.cross(r12, r23)
    r23xr34 = np.cross(r23, r34)
    dihedral0 = np.arccos(np.dot(r12xr23, r23xr34)  / np.linalg.norm(r12xr23) / np.linalg.norm(r23xr34))
    deltaAngle = dihedral - dihedral0
    crd[atom4_friends] = np.dot(crd[atom4_friends] - crd[molecule.atom_index[atom3]], _get_rotate_matrix(r23, deltaAngle)) + crd[molecule.atom_index[atom3]]
    

sys.modules['__main__'].__dict__["Impose_Dihedral"] = Impose_Dihedral 

def _link_residue_process_coordinate(molecule, atom1, atom2, atom10, atom100, linkbond, linkangle, linkdihedral):
    crd = _get_crd(molecule)
    atom1_friends, atom2_friends = _get_friends(molecule, atom1, atom2)
    crd[atom2_friends] += 2000
    r0 = crd[molecule.atom_index[atom2]] - crd[molecule.atom_index[atom1]]
    L0 = np.linalg.norm(r0)
    dr = (linkbond/L0 - 1) * r0
    crd[atom2_friends] += dr
    
    if atom10:        
        r101 = crd[molecule.atom_index[atom10]] - crd[molecule.atom_index[atom1]]
        angle0 = np.arccos(np.dot(r0, r101)  / np.linalg.norm(r0) / np.linalg.norm(r101))
        deltaAngle =  linkangle - angle0 
        crd[atom2_friends] = np.dot(crd[atom2_friends] - crd[molecule.atom_index[atom1]], _get_rotate_matrix( np.cross(r101, r0), deltaAngle)) + crd[molecule.atom_index[atom1]]
        if atom100:
            atom1_friends, atom2_friends = _get_friends(molecule, atom10, atom1)
            r000 = crd[molecule.atom_index[atom100]] - crd[molecule.atom_index[atom10]]
            r12xr23 = np.cross(r000, r101)
            r23xr34 = np.cross(r101, r0)
            dihedral0 = np.arccos(np.dot(r12xr23, r23xr34)  / np.linalg.norm(r12xr23) / np.linalg.norm(r23xr34))
            deltaAngle = linkdihedral - dihedral0
            crd[atom2_friends] = np.dot(crd[atom2_friends] - crd[molecule.atom_index[atom1]], _get_rotate_matrix(r101, deltaAngle)) + crd[molecule.atom_index[atom1]]
    
    for atom in molecule.atoms:
        i = molecule.atom_index[atom]
        atom.x = crd[i][0]
        atom.y = crd[i][1]
        atom.z = crd[i][2]
        
def ResidueType_Add(self, other):
    if type(other) == ResidueType:
        new_molecule = Molecule(self.name)
        resA = Residue(self)
        resB = Residue(other)
        for atom in self.atoms:
            resA.Add_Atom(atom)
        for atom in other.atoms:
            resB.Add_Atom(atom)
        new_molecule.Add_Residue(resA)
        new_molecule.Add_Residue(resB)
        assert resA.type.tail and resB.type.head, "+ and * can only be used to link head and tail"
        atom1 = resA._name2atom[self.tail]
        atom2 = resB._name2atom[other.head] 
        new_molecule.Add_Residue_Link(atom1, atom2)
        atom10 = resA._name2atom[self.tail_second]
        atom100 = resA._name2atom[self.tail_third] 
        _link_residue_process_coordinate(new_molecule, atom1, atom2, atom10, atom100, 
            resA.type.tail_bond, resA.type.tail_angle, resA.type.tail_dihedral)            
        return new_molecule
    elif type(other) == Molecule:
        new_molecule = other.deepcopy()
        resA = Residue(self)
        resB = other.residues[0]
        for atom in self.atoms:
            resA.Add_Atom(atom)
        new_molecule.residues.insert(0, resA)
        assert resA.type.tail and resB.type.head, "+ and * can only be used to link head and tail"
        atom1 = resA._name2atom[self.tail]
        atom2 = resB._name2atom[other.residues[0].type.head] 
        new_molecule.Add_Residue_Link(atom1, atom2)
        atom10 = resA._name2atom[self.tail_second]
        atom100 = resA._name2atom[self.tail_third] 
        _link_residue_process_coordinate(new_molecule, atom1, atom2, atom10, atom100, 
            resA.type.tail_bond, resA.type.tail_angle, resA.type.tail_dihedral)   
        
        return new_molecule
    elif type(other) == type(None):
        return self
    else:
        raise TypeError("unsupported operand type(s) for +: '%s' and '%s'"%(type(self), type(other)))

def Molecule_Add(self, other):
    if type(other) == ResidueType:
        new_molecule = self.deepcopy()
        resA = new_molecule.residues[-1]
        resB = Residue(other)
        for atom in other.atoms:
            resB.Add_Atom(atom)
        new_molecule.Add_Residue(resB)
        assert resA.type.tail and resB.type.head, "+ and * can only be used to link head and tail"
        atom1 = resA._name2atom[resA.type.tail]
        atom2 = resB._name2atom[other.head]
        new_molecule.Add_Residue_Link(atom1, atom2)
        atom10 = resA._name2atom[resA.type.tail_second]
        atom100 = resA._name2atom[resA.type.tail_third] 
        _link_residue_process_coordinate(new_molecule, atom1, atom2, atom10, atom100, 
            resA.type.tail_bond, resA.type.tail_angle, resA.type.tail_dihedral)    
        
        return new_molecule
    elif type(other) == Molecule:
        new_molecule = self.deepcopy()
        new_molecule2 = other.deepcopy()
        resA = new_molecule.residues[-1]
        resB = new_molecule2.residues[0]
        for res in new_molecule2.residues:
            new_molecule.Add_Residue(res)
        assert resA.type.tail and resB.type.head, "+ and * can only be used to link head and tail"
        atom1 = resA._name2atom[resA.type.tail]
        atom2 = resB._name2atom[resB.type.head]
        new_molecule.Add_Residue_Link(atom1, atom2)
        atom10 = resA._name2atom[resA.type.tail_second]
        atom100 = resA._name2atom[resA.type.tail_third] 
        _link_residue_process_coordinate(new_molecule, atom1, atom2, atom10, atom100, 
            resA.type.tail_bond, resA.type.tail_angle, resA.type.tail_dihedral)   

        return new_molecule
    elif type(other) == type(None):
        return self
    else:
        raise TypeError("unsupported operand type(s) for +: '%s' and '%s'"%(type(self), type(other)))

def iMolecule_Add(self, other):
    if type(other) == ResidueType:
        resA = self.residues[-1]
        resB = Residue(other)
        for atom in other.atoms:
            resB.Add_Atom(atom)
        self.Add_Residue(resB)
        if resA.type.tail and other.head:
            atom1 = resA._name2atom[resA.type.tail]
            atom2 = resB._name2atom[other.head]
            self.Add_Residue_Link(atom1, atom2)
            atom10 = resA._name2atom[resA.type.tail_second]
            atom100 = resA._name2atom[resA.type.tail_third] 
            _link_residue_process_coordinate(self, atom1, atom2, atom10, atom100, 
                resA.type.tail_bond, resA.type.tail_angle, resA.type.tail_dihedral)  
        return self
    elif type(other) == Molecule:
        self = self.deepcopy()
        new_molecule2 = other.deepcopy()
        resA = self.residues[-1]
        resB = new_molecule2.residues[0]
        for res in new_molecule2.residues:
            self.Add_Residue(res)
        if resA.type.tail and resB.type.head:
            atom1 = resA._name2atom[resA.type.tail]
            atom2 = resB._name2atom[resB.type.head]
            self.Add_Residue_Link(atom1, atom2)
            atom10 = resA._name2atom[resA.type.tail_second]
            atom100 = resA._name2atom[resA.type.tail_third] 
            _link_residue_process_coordinate(self, atom1, atom2, atom10, atom100, 
                resA.type.tail_bond, resA.type.tail_angle, resA.type.tail_dihedral)   
        return self
    elif type(other) == type(None):
        return self
    else:
        raise TypeError("unsupported operand type(s) for +: '%s' and '%s'"%(type(self), type(other)))

def Muls(self, other):
    if type(other) == int:
        assert other >= 0
        t = self
        for i in range(other - 1):
            t = t + self
        return t
    else:
        raise TypeError("unsupported operand type(s) for +: '%s' and '%s'"%(type(self), type(other)))

def iMuls(self, other):
    if type(other) == int:
        assert other >= 0
        for i in range(other - 1):
            self += self
        return self
    else:
        raise TypeError("unsupported operand type(s) for +: '%s' and '%s'"%(type(self), type(other)))

ResidueType.__add__ = ResidueType_Add
ResidueType.__radd__ = ResidueType_Add
ResidueType.__mul__ = Muls
ResidueType.__rmul__ = Muls
Molecule.__add__ = Molecule_Add
Molecule.__radd__ = Molecule_Add
Molecule.__iadd__ = iMolecule_Add
Molecule.__mul__ = Muls
Molecule.__rmul__ = Muls
Molecule.__imul__ = iMuls

del ResidueType_Add
del Molecule_Add
del Muls
del iMuls
del iMolecule_Add