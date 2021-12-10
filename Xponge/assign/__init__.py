from .. import *
import pubchempy as pcp
from itertools import groupby

class Judge_Rule():
    all = {}
    def __init__(self, name):
        self.name = name
        Judge_Rule.all[name] = self
        self.rules = OrderedDict()
    def Add_Judge_Rule(self, atomtype):
        if type(atomtype) == type(""):
            atomtype = AtomType.types[atomtype]
        def wrapper(rule_function):
            self.rules[atomtype] = rule_function
        return wrapper


class _RING():
    def __repr__(self):
        return self.tohash
    def __hash__(self):
        return hash(self.tohash)
    def __eq__(self, other):
        return type(other) == _RING and self.tohash == other.tohash
    def __init__(self, atom_list):
        min_index = np.argmin(atom_list)
        self.atoms = atom_list[min_index:] + atom_list[:min_index]
        reverse_list = self.atoms[::-1]
        reverse_list = reverse_list[-1:] + reverse_list[:-1]
        if reverse_list[1] < self.atoms[1]:
            self.atoms = reverse_list
        self.tohash = "-".join(["%d"%atom for atom in self.atoms])
    def get_3_neighbors(self):
        for i in range(len(self.atoms)):
            yield self.atoms[i-2],self.atoms[i-1],self.atoms[i]
            
class ASSIGN():
    XX = set("CNOPS")
    XA = set("OS")
    XB = set("NP")
    XC = set(["F", "Cl", "Br", "I"])
    XD = set("SP")
    XE = set(["N", "O", "F", "Cl", "Br"])
    def __init__(self):
        self.atom_numbers = 0
        self.atoms = []
        self.names = []
        self.coordinate = None
        self.atom_types = {}
        self.atom_marker = {}
        self.bonds = {}
        self.bond_marker = {}
    
    def Atom_Judge(self, atom, string):
        element, links = [''.join(list(g)) for k, g in groupby(string, key=lambda x: x.isdigit())]
        return self.atoms[atom] == element and len(self.bonds[atom]) == int(links)
    
    def Add_Atom(self, element, x, y, z, name = ""):
        self.atoms.append(element)
        self.bonds[self.atom_numbers] = {}
        self.bond_marker[self.atom_numbers] = {}
        self.atom_marker[self.atom_numbers] = {}
        self.atom_types[self.atom_numbers] = None
        self.atom_numbers += 1
        self.names.append(name)
        if type(self.coordinate) == type(None):
            self.coordinate = np.array([[x,y,z]])
        else:
            self.coordinate = np.vstack((self.coordinate, np.array([x,y,z])))
    
    def Add_Atom_Marker(self, atom, marker):
        if marker in self.atom_marker[atom].keys():
            self.atom_marker[atom][marker] += 1
        else:
            self.atom_marker[atom][marker] = 1
    
    def Add_Bond(self, atom1, atom2, order = -1):
        self.bonds[atom1][atom2] = order
        self.bond_marker[atom1][atom2] = set([])
        self.bonds[atom2][atom1] = order   
        self.bond_marker[atom2][atom1] = set([])
        
    
    def Add_Bond_Marker(self, atom1, atom2, marker, only1 = False):
        self.bond_marker[atom1][atom2].add(marker)
        if marker in self.atom_marker[atom1]:
            self.atom_marker[atom1][marker] += 1
        else:
            self.atom_marker[atom1][marker] = 1
        if not only1:
            self.bond_marker[atom2][atom1].add(marker)
            if marker in self.atom_marker[atom2]:
                self.atom_marker[atom2][marker] += 1
            else:
                self.atom_marker[atom2][marker] = 1
    def Determine_Bond_Order(self):
        raise NotImplementedError

    def Determine_Ring_And_Bond_Type(self):
        current_path = []
        current_path_sons = {}
        current_work = []
        current_path_father = {}
        have_found_rings = set([])
        for atom0 in range(len(self.atoms)):
            current_path.append(atom0)
            current_work.extend([[atom, atom0] for atom in self.bonds[atom0].keys()])
            current_path_sons[atom0] = len(self.bonds[atom0])
            current_path_father = []
            while current_path:    
                work_atom, from_atom = current_work.pop()
                current_path.append(work_atom)
                current_path_father.append(from_atom)
                bond_atom = []
                for atom in self.bonds[work_atom].keys():
                    if atom != from_atom:
                        try:
                            index = current_path.index(atom)                           
                            have_found_rings.add(_RING(current_path[index:]))
                        except ValueError:
                            bond_atom.append([atom, work_atom])
                           
                if len(current_path) < 9:                    
                    current_path_sons[work_atom] = len(bond_atom)
                    current_work.extend(bond_atom)
                    
                else:
                    current_path_sons[work_atom] = 0
                    
                for atom in current_path[::-1]:
                    if current_path_sons[atom] == 0:
                        pop_atom = current_path.pop()
                        current_path_sons.pop(pop_atom)
                        if current_path_father:
                            father_atom = current_path_father.pop()
                            current_path_sons[father_atom] -= 1
        
        for ring in have_found_rings:
            if len(ring.atoms) == 6:
                ring.is_pure_aromatic_ring = True
                for atom in ring.atoms:
                    if not self.Atom_Judge(atom, "C3") and not self.Atom_Judge(atom, "N2") and not self.Atom_Judge(atom, "N3"):
                        ring.is_pure_aromatic_ring = False
                        break
            else:
                ring.is_pure_aromatic_ring = False      
            ring.is_pure_aliphatic_ring = True
            ring.is_planar_ring = True
            for atom in ring.atoms:
                self.Add_Atom_Marker(atom, "RG")
                self.Add_Atom_Marker(atom, "RG%d"%len(ring.atoms))
                if ring.is_pure_aromatic_ring:
                    self.Add_Atom_Marker(atom, "AR1")
                    for i in range(6):
                        self.Add_Bond_Marker(ring.atoms[i-1], ring.atoms[i], "AB")
                if not self.Atom_Judge(atom, "C4"):
                    ring.is_pure_aliphatic_ring = False
                if (not self.Atom_Judge(atom, "C3") and not self.Atom_Judge(atom, "N2") 
                and not self.Atom_Judge(atom, "N3") and not self.Atom_Judge(atom, "O2")
                and not self.Atom_Judge(atom, "S2") and not self.Atom_Judge(atom, "P2")
                and not self.Atom_Judge(atom, "P3")):
                    ring.is_planar_ring = False
            
            ring.out_plane_double_bond = False
            for atom in ring.atoms:
                for bonded_atom, order in self.bonds[atom].items():
                    if self.atoms[bonded_atom] != "C" and order == 2 and bonded_atom not in ring.atoms:
                        ring.out_plane_double_bond = True
                        
            if not ring.is_pure_aromatic_ring:
                for atom in ring.atoms:
                    if ring.is_pure_aliphatic_ring:
                        self.Add_Atom_Marker(atom, "AR5")
                    elif ring.is_planar_ring:
                        if ring.out_plane_double_bond:
                            self.Add_Atom_Marker(atom, "AR3")
                        else:
                            self.Add_Atom_Marker(atom, "AR2")       
                    else:
                        self.Add_Atom_Marker(atom, "AR4")
        
        for atom in range(len(self.atoms)):
            for atom2, order in self.bonds[atom].items():
                if order == 1:
                    self.Add_Bond_Marker(atom, atom2, "SB", True)
                    if "AB" not in self.bond_marker[atom][atom2]:
                        self.Add_Bond_Marker(atom, atom2, "sb", True)
                elif order == 2:
                    self.Add_Bond_Marker(atom, atom2, "DB", True)
                    if "AB" not in self.bond_marker[atom][atom2]:
                        self.Add_Bond_Marker(atom, atom2, "db", True)
                else:
                    self.Add_Bond_Marker(atom, atom2, "tb", True)
                
        #print(self.atom_marker, "\n\n", self.atom_marker, "\n\n", self.bond_marker)
            
            
    def Determine_Atom_Type(self, rule):
        if type(rule) == type(""):
            rule = Judge_Rule.all[rule]
        #print(self.atoms)
        for i in range(len(self.atoms)):
            find_type = False
            for atom_type, type_rule in rule.rules.items():
                if type_rule(i, self):
                    self.atom_types[i] = atom_type
                    find_type = True
                    break
            
            assert find_type
        #print(self.atom_types)
                
    def To_ResidueType(self, name, charge):
        temp = ResidueType(name = name)
        count = {}
        for i in range(self.atom_numbers):
            assert self.atom_types[i] != None
            if self.names[i]:
                atom_name = self.names[i]
            elif self.atoms[i] in count.keys():
                atom_name = self.atoms[i] + "%d"%count[self.atoms[i]]
                count[self.atoms[i]] += 1
            else:
                count[self.atoms[i]] = 1
                atom_name = self.atoms[i]
            temp.Add_Atom(atom_name, self.atom_types[i], x = self.coordinate[i][0],
                y = self.coordinate[i][1],  z = self.coordinate[i][2])
            temp.atoms[-1].charge = charge[i] * 18.2223
        for i, bondi in self.bonds.items():
            for j in bondi.keys():
                temp.Add_Connectivity(temp.atoms[i], temp.atoms[j])
        sys.modules["__main__"].__dict__[name] = temp
        return temp
        
def Get_Molecule_From_PubChem(parameter, keyword):
    cs = pcp.get_compounds(parameter, keyword, record_type='3d')
    if len(cs) == 0:
        raise pcp.NotFoundError
    elif len(cs) == 1:
        assign = ASSIGN()
        c = cs[0]
        for atom in c.atoms:
            assign.Add_Atom(atom.element, atom.x, atom.y, atom.z)
        for bond in c.bonds:
            assign.Add_Bond(bond.aid1 - 1, bond.aid2 - 1, bond.order)
        assign.Determine_Ring_And_Bond_Type()
        return assign
    else:
        raise NotImplementedError


from . import RESP

sys.modules["__main__"].__dict__["Get_Molecule_From_PubChem"] = Get_Molecule_From_PubChem