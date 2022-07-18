from .. import *
import pubchempy as pcp
from itertools import groupby
import sys


class Judge_Rule:
    all = {}

    def __init__(self, name):
        self.name = name
        Judge_Rule.all[name] = self
        self.rules = OrderedDict()

    def Add_Judge_Rule(self, atomtype):
        if isinstance(atomtype, str):
            atomtype = AtomType.types[atomtype]

        def wrapper(rule_function):
            self.rules[atomtype] = rule_function

        return wrapper


class _RING():
    @staticmethod
    def Add_Rings_Basic_Marker(assign, rings):
        for ring in rings:
            for atom in ring.atoms:
                assign.Add_Atom_Marker(atom, "RG")
                assign.Add_Atom_Marker(atom, "RG%d" % len(ring.atoms))

    @staticmethod
    def Check_Rings_Type(assign, rings):
        for ring in rings:
            ring.check_pure_aromatic(assign)
            ring.check_pure_aliphatic_and_planar(assign)
            ring.check_out_plane_double_bond(assign)

            if not ring.is_pure_aromatic_ring:
                for atom in ring.atoms:
                    if ring.is_pure_aliphatic_ring:
                        assign.Add_Atom_Marker(atom, "AR5")
                    elif ring.is_planar_ring:
                        if ring.out_plane_double_bond:
                            assign.Add_Atom_Marker(atom, "AR3")
                        else:
                            assign.Add_Atom_Marker(atom, "AR2")
                    else:
                        assign.Add_Atom_Marker(atom, "AR4")

    @staticmethod
    def Get_Rings(assign):
        current_path = []
        current_path_sons = {}
        current_work = []
        current_path_father = {}
        have_found_rings = set([])
        for atom0 in range(len(assign.atoms)):
            current_path.append(atom0)
            current_work.extend([[atom, atom0] for atom in assign.bonds[atom0].keys()])
            current_path_sons[atom0] = len(assign.bonds[atom0])
            current_path_father = []
            while current_path:
                work_atom, from_atom = current_work.pop()
                current_path.append(work_atom)
                current_path_father.append(from_atom)
                bond_atom = []
                for atom in assign.bonds[work_atom].keys():
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
        return have_found_rings

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
        self.tohash = "-".join(["%d" % atom for atom in self.atoms])

    def get_3_neighbors(self):
        for i in range(len(self.atoms)):
            yield self.atoms[i - 2], self.atoms[i - 1], self.atoms[i]

    def check_pure_aromatic(self, assign):
        if len(self.atoms) == 6:
            self.is_pure_aromatic_ring = True
            for atom in self.atoms:
                if not assign.Atom_Judge(atom, "C3") and not assign.Atom_Judge(atom, "N2") and not assign.Atom_Judge(
                        atom, "N3"):
                    self.is_pure_aromatic_ring = False
                    break
                if assign.Atom_Judge(atom, "N3"):
                    temp = 0
                    for bonded_atom, bond_order in assign.bonds[atom].items():
                        temp += bond_order
                    if temp == 3:
                        self.is_pure_aromatic_ring = False
                        break
                for bonded_atom, bond_order in assign.bonds[atom].items():
                    if bond_order == 2 and "RG" not in assign.atom_marker[bonded_atom].keys():
                        self.is_pure_aromatic_ring = False
                        break
                if self.is_pure_aromatic_ring == False:
                    break
        else:
            self.is_pure_aromatic_ring = False

    def check_pure_aliphatic_and_planar(self, assign):
        self.is_pure_aliphatic_ring = True
        self.is_planar_ring = True
        for atom in self.atoms:
            if self.is_pure_aromatic_ring:
                assign.Add_Atom_Marker(atom, "AR1")
                for i in range(6):
                    assign.Add_Bond_Marker(self.atoms[i - 1], self.atoms[i], "AB")
            if not assign.Atom_Judge(atom, "C4"):
                self.is_pure_aliphatic_ring = False
            if (not assign.Atom_Judge(atom, "C3") and not assign.Atom_Judge(atom, "N2")
                    and not assign.Atom_Judge(atom, "N3") and not assign.Atom_Judge(atom, "O2")
                    and not assign.Atom_Judge(atom, "S2") and not assign.Atom_Judge(atom, "P2")
                    and not assign.Atom_Judge(atom, "P3")):
                self.is_planar_ring = False

    def check_out_plane_double_bond(self, assign):
        self.out_plane_double_bond = False
        for atom in self.atoms:
            for bonded_atom, order in assign.bonds[atom].items():
                if assign.atoms[bonded_atom] != "C" and order == 2 and bonded_atom not in self.atoms:
                    self.out_plane_double_bond = True


class Assign():
    XX = set("CNOPS")
    XA = set("OS")
    XB = set("NP")
    XC = set(["F", "Cl", "Br", "I"])
    XD = set("SP")
    XE = set(["N", "O", "F", "Cl", "Br", "S", "I"])

    def __init__(self, name="ASN"):
        self.name = name
        self.atom_numbers = 0
        self.atoms = []
        self.names = []
        self.element_details = []
        self.coordinate = None
        self.charge = None
        self.atom_types = {}
        self.atom_marker = {}
        self.bonds = {}
        self.ar_bonds = {}
        self.am_bonds = {}
        self.bond_marker = {}

    def Add_Index_To_Name(self):
        for i in range(self.atom_numbers):
            self.names[i] += str(i)

    def Atom_Judge(self, atom, string):
        assert isinstance(string, list) or isinstance(string, str)
        if isinstance(string, str):
            todo = [string]
        else:
            todo = string
        judge = False
        for s in todo:
            element, links = [''.join(list(g)) for k, g in groupby(s, key=lambda x: x.isdigit())]
            if self.atoms[atom] == element and len(self.bonds[atom]) == int(links):
                judge = True
                break
        return judge

    def Add_Atom(self, element, x, y, z, name="", charge=0.0):
        if "." in element:
            element, element_detail = element.split(".")
            element_detail = "." + element_detail
        else:
            element_detail = ""
        self.element_details.append(element_detail)
        self.atoms.append(element)
        self.bonds[self.atom_numbers] = {}
        self.bond_marker[self.atom_numbers] = {}
        self.atom_marker[self.atom_numbers] = {}
        self.atom_types[self.atom_numbers] = None
        self.atom_numbers += 1
        self.names.append(name)
        if self.coordinate is None:
            self.coordinate = np.array([[x, y, z]])
        else:
            self.coordinate = np.vstack((self.coordinate, np.array([x, y, z])))
        if self.charge is None:
            self.charge = np.array([charge])
        else:
            self.charge = np.hstack((self.charge, np.array([charge])))

    def Add_Atom_Marker(self, atom, marker):
        if marker in self.atom_marker[atom].keys():
            self.atom_marker[atom][marker] += 1
        else:
            self.atom_marker[atom][marker] = 1

    def Add_Bond(self, atom1, atom2, order=-1):
        self.bonds[atom1][atom2] = order
        self.bond_marker[atom1][atom2] = set([])
        self.bonds[atom2][atom1] = order
        self.bond_marker[atom2][atom1] = set([])

    def Add_Bond_Marker(self, atom1, atom2, marker, only1=False):
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

    def Determine_Equal_Atoms(self):
        from .RDKit_tools import Find_Equal_Atoms
        return Find_Equal_Atoms(self)

    def Determine_Bond_Order(self):
        raise NotImplementedError

    def Determine_Ring_And_Bond_Type(self):
        have_found_rings = _RING.Get_Rings(self)
        _RING.Add_Rings_Basic_Marker(self, have_found_rings)
        _RING.Check_Rings_Type(self, have_found_rings)

        for atom in range(len(self.atoms)):
            DLO = 0
            NOTO = 0
            for atom2, order in self.bonds[atom].items():
                if self.Atom_Judge(atom2, "O1"):
                    DLO += 1
                else:
                    NOTO += 1
            if DLO >= 1 and NOTO <= 1:
                for atom2, order in self.bonds[atom].items():
                    if self.Atom_Judge(atom2, "O1"):
                        self.Add_Bond_Marker(atom, atom2, "DLB")
            for atom2, order in self.bonds[atom].items():
                if "DLB" in self.bond_marker[atom][atom2]:
                    self.Add_Bond_Marker(atom, atom2, "DL", True)
                    self.Add_Bond_Marker(atom, atom2, "sb", True)
                elif order == 1:
                    self.Add_Bond_Marker(atom, atom2, "sb", True)
                    if "AB" not in self.bond_marker[atom][atom2]:
                        self.Add_Bond_Marker(atom, atom2, "SB", True)
                elif order == 2:
                    self.Add_Bond_Marker(atom, atom2, "db", True)
                    if "AB" not in self.bond_marker[atom][atom2]:
                        self.Add_Bond_Marker(atom, atom2, "DB", True)
                else:
                    self.Add_Bond_Marker(atom, atom2, "tb", True)

    def Determine_Atom_Type(self, rule):
        if isinstance(rule, str):
            rule = Judge_Rule.all[rule]

        for i in range(len(self.atoms)):
            find_type = False
            for atom_type, type_rule in rule.rules.items():
                if type_rule(i, self):
                    self.atom_types[i] = atom_type
                    find_type = True
                    break

            assert find_type, "No atom type found for assignment %s of atom #%d" % (self.name, i)

    def To_ResidueType(self, name, charge=None):
        temp = ResidueType(name=name)
        if not charge:
            if self.charge is None:
                charge = np.zeros(self.atom_numbers)
            else:
                charge = self.charge
        count = {}
        for i in range(self.atom_numbers):
            assert self.atom_types[i] != None
            if self.names[i]:
                atom_name = self.names[i]
            elif self.atoms[i] in count.keys():
                atom_name = self.atoms[i] + "%d" % count[self.atoms[i]]
                self.names[i] = atom_name
                count[self.atoms[i]] += 1
            else:
                count[self.atoms[i]] = 1
                atom_name = self.atoms[i]
                self.names[i] = atom_name
            temp.Add_Atom(atom_name, self.atom_types[i], x=self.coordinate[i][0],
                          y=self.coordinate[i][1], z=self.coordinate[i][2])
            temp.atoms[-1].charge = charge[i]
        for i, bondi in self.bonds.items():
            for j in bondi.keys():
                temp.Add_Connectivity(temp.atoms[i], temp.atoms[j])
        sys.modules["__main__"].__dict__[name] = temp
        return temp

    def Calculate_Charge(self, method, **parameters):
        if method == "RESP":
            from . import RESP
            self.charge = RESP.RESP_Fit(self, basis=parameters.get("basis", "6-31g*"), opt=parameters.get("opt", False),
                                        charge=parameters.get("charge", 0), spin=parameters.get("spin", 0),
                                        extra_equivalence=parameters.get("extra_equivalence", []),
                                        grid_density=parameters.get("grid_density", 6),
                                        grid_cell_layer=parameters.get("grid_cell_layer", 4),
                                        a1=parameters.get("a1", 0.0005),
                                        a2=parameters.get("a2", 0.001), two_stage=parameters.get("two_stage", True),
                                        only_ESP=parameters.get("only_ESP", False),
                                        radius=parameters.get("radius", None))

    def Save_As_PDB(self, filename):
        towrite = towrite = "REMARK   Generated By Xponge (Assignment)\n"
        count = {}
        for i in range(self.atom_numbers):
            if self.names[i]:
                atom_name = self.names[i]
            elif self.atoms[i] in count.keys():
                atom_name = self.atoms[i] + "%d" % count[self.atoms[i]]
                self.names[i] = atom_name
                count[self.atoms[i]] += 1
            else:
                count[self.atoms[i]] = 1
                atom_name = self.atoms[i]
                self.names[i] = atom_name
        for i, atom in enumerate(self.atoms):
            towrite += "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%17s%2s\n" % (i + 1, self.names[i],
                                                                                 self.name, " ", 1,
                                                                                 self.coordinate[i][0],
                                                                                 self.coordinate[i][1],
                                                                                 self.coordinate[i][2], " ", atom)

        for i in range(self.atom_numbers):
            bonded_atoms = list(self.bonds[i].keys())
            bonded_atoms.sort()
            bonded_atoms = [bonded_atoms[i:i + 4] for i in range(0, len(bonded_atoms), 4)]
            if len(bonded_atoms) > 0:
                for atoms in bonded_atoms:
                    towrite += "CONECT %4d" % (i + 1)
                    for atom in atoms:
                        towrite += " %4d" % (atom + 1)
                    towrite += "\n"

        f = open(filename, "w")
        f.write(towrite)
        f.close()

    def Save_As_Mol2(self, filename):
        bonds = []
        for i in range(self.atom_numbers):
            for j, order in self.bonds[i].items():
                if i < j:
                    if i in self.ar_bonds.keys() and j in self.ar_bonds[i]:
                        bonds.append("%6d %6d ar\n" % (i + 1, j + 1))
                    elif i in self.am_bonds.keys() and j in self.am_bonds[i]:
                        bonds.append("%6d %6d am\n" % (i + 1, j + 1))
                    else:
                        bonds.append("%6d %6d %1d\n" % (i + 1, j + 1, order))
        bonds.sort(key=lambda x: list(map(int, x.split()[:2])))
        count = {}
        for i in range(self.atom_numbers):
            if self.names[i]:
                atom_name = self.names[i]
            elif self.atoms[i] in count.keys():
                atom_name = self.atoms[i] + "%d" % count[self.atoms[i]]
                self.names[i] = atom_name
                count[self.atoms[i]] += 1
            else:
                count[self.atoms[i]] = 1
                atom_name = self.atoms[i]
                self.names[i] = atom_name
        towrite = "@<TRIPOS>MOLECULE\n%s\n %d %d 1 0 1\nSMALL\nUSER_CHARGES\n" % (
        self.name, self.atom_numbers, len(bonds))
        towrite += "@<TRIPOS>ATOM\n"
        for i, atom in enumerate(self.atoms):
            towrite += "%6d %4s %8.4f %8.4f %8.4f   %-8s %5d %8s %10.6f\n" % (
                i + 1, self.names[i], self.coordinate[i][0], self.coordinate[i][1], self.coordinate[i][2],
                atom + self.element_details[i], 1, self.name, self.charge[i])

        towrite += "@<TRIPOS>BOND\n"
        for i, bond in enumerate(bonds):
            towrite += "%6d %s" % (i + 1, bond)
        towrite += "@<TRIPOS>SUBSTRUCTURE\n"
        towrite += "%5d %8s %6d ****               0 ****  **** \n" % (1, self.name, 1)

        f = open(filename, "w")
        f.write(towrite)
        f.close()


def Guess_Element_From_Mass(mass):
    elements = ["X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
                "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc",
                "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
                "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc",
                "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
                "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
                "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
                "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
                "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
                "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
                "Ds", "Rg"]
    masses = [0.00000, 1.00794, 4.00260, 6.941, 9.012182, 10.811,
              12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,
              22.989770, 24.3050, 26.981538, 28.0855, 30.973761,
              32.065, 35.453, 39.948, 39.0983, 40.078, 44.955910,
              47.867, 50.9415, 51.9961, 54.938049, 55.845, 58.9332,
              58.6934, 63.546, 65.409, 69.723, 72.64, 74.92160,
              78.96, 79.904, 83.798, 85.4678, 87.62, 88.90585,
              91.224, 92.90638, 95.94, 98.0, 101.07, 102.90550,
              106.42, 107.8682, 112.411, 114.818, 118.710, 121.760,
              127.60, 126.90447, 131.293, 132.90545, 137.327,
              138.9055, 140.116, 140.90765, 144.24, 145.0, 150.36,
              151.964, 157.25, 158.92534, 162.500, 164.93032,
              167.259, 168.93421, 173.04, 174.967, 178.49, 180.9479,
              183.84, 186.207, 190.23, 192.217, 195.078, 196.96655,
              200.59, 204.3833, 207.2, 208.98038, 209.0, 210.0, 222.0,
              223.0, 226.0, 227.0, 232.0381, 231.03588, 238.02891,
              237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0,
              258.0, 259.0, 262.0, 261.0, 262.0, 266.0, 264.0, 269.0,
              268.0, 271.0, 272.0]
    if mass > 0.0 and mass < 3.8:
        index = 1;
    elif mass > 207.85 and mass < 208.99:
        index = 83;
    elif mass > 56.50 and mass < 58.8133:
        index = 27;
    else:
        index = 0;
        for j in range(0, 111):
            if abs(mass - masses[j]) < 0.65:
                index = j
                break
    return elements[index]


def Get_Assignment_From_PubChem(parameter, keyword):
    cs = pcp.get_compounds(parameter, keyword, record_type='3d')
    if len(cs) == 0:
        raise pcp.NotFoundError
    elif len(cs) == 1:
        assign = Assign()
        c = cs[0]
        for atom in c.atoms:
            assign.Add_Atom(atom.element, atom.x, atom.y, atom.z)
        for bond in c.bonds:
            assign.Add_Bond(bond.aid1 - 1, bond.aid2 - 1, bond.order)
        assign.Determine_Ring_And_Bond_Type()
        return assign
    else:
        raise NotImplementedError


def Get_Assignment_From_PDB(filename, determine_bond_order=True, only_residue=""):
    assign = Assign()
    index_atom_map = {}
    with open(filename) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if only_residue:
                    resname = line[17:20].strip()
                    if resname != only_residue:
                        continue
                index = int(line[6:11])
                index_atom_map[index] = assign.atom_numbers
                atom_name = line[12:16].strip()
                element = line[76:78].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                assign.Add_Atom(element, x, y, z, atom_name)
            if line.startswith("CONECT"):
                atom = int(line[6:11])
                if atom not in index_atom_map.keys():
                    continue
                for bonded_atom_i in range(11, 31, 5):
                    try:
                        temp = line[bonded_atom_i:bonded_atom_i + 5]
                        bonded_atom = int(temp)
                    except:
                        break
                    if bonded_atom in index_atom_map.keys():
                        assign.Add_Bond(index_atom_map[atom], index_atom_map[int(bonded_atom)])
    if determine_bond_order:
        assign.Determine_Bond_Order()
        assign.Determine_Ring_And_Bond_Type()
    return assign


def Get_Assignment_From_ResidueType(restype):
    assign = Assign()
    for atom in restype.atoms:
        assign.Add_Atom(Guess_Element_From_Mass(atom.mass), atom.x, atom.y, atom.z, atom.name)
    for atom in restype.atoms:
        i = restype._atom2index[atom]
        for atomb in restype.connectivity[atom]:
            j = restype._atom2index[atomb]
            if i < j:
                assign.Add_Bond(i, j)
    return assign


def _deal_with_ar_bonds(assign):
    ar_bonds_atoms = list(assign.ar_bonds.keys())
    ar_bonds_atoms.sort(key=lambda x: (x, len(assign.ar_bonds[x])))
    doubled = {}
    checked = {}
    for ar_atom in ar_bonds_atoms:
        assign.ar_bonds[ar_atom].sort(key=lambda x: (x, len(assign.ar_bonds[x])))
        doubled[ar_atom] = False
        checked[ar_atom] = False

    working_space = []
    while ar_bonds_atoms:
        working_space.append(ar_bonds_atoms.pop())
        while working_space:
            work_atom = working_space.pop()
            if checked[work_atom]:
                continue
            checked[work_atom] = True
            for neighbor in assign.ar_bonds[work_atom]:
                if not checked[neighbor]:
                    working_space.append(neighbor)
            for neighbor in assign.ar_bonds[work_atom][::-1]:
                if doubled[work_atom]:
                    break
                if not doubled[neighbor] and not doubled[work_atom]:
                    assign.bonds[neighbor][work_atom] = 2
                    doubled[neighbor] = True
                    doubled[work_atom] = True


def Get_Assignment_From_Mol2(filename):
    with open(filename) as f:
        flag = None
        subflag = None
        for line in f:
            if not line.strip():
                continue
            if line.startswith("@<TRIPOS>"):
                flag = line[9:].strip()
            elif flag == "MOLECULE":
                if subflag is None:
                    assign = Assign(line.strip())
                    subflag = "0"
            # 处理原子信息
            elif flag == "ATOM":
                words = line.split()
                atom_name = words[1]
                element = words[5]
                x = float(words[2])
                y = float(words[3])
                z = float(words[4])
                charge = float(words[8])
                assign.Add_Atom(element, x, y, z, atom_name, charge)
            elif flag == "BOND":
                words = line.split()
                if words[3] in "1234567890":
                    assign.Add_Bond(int(words[1]) - 1, int(words[2]) - 1, int(words[3]))
                elif words[3] == "ar":
                    atom1 = int(words[1]) - 1
                    atom2 = int(words[2]) - 1
                    assign.Add_Bond(atom1, atom2, 1)
                    if atom1 not in assign.ar_bonds.keys():
                        assign.ar_bonds[atom1] = [atom2]
                    else:
                        assign.ar_bonds[atom1].append(atom2)
                    if atom2 not in assign.ar_bonds.keys():
                        assign.ar_bonds[atom2] = [atom1]
                    else:
                        assign.ar_bonds[atom2].append(atom1)
                elif words[3] == "am":
                    atom1 = int(words[1]) - 1
                    atom2 = int(words[2]) - 1
                    assign.Add_Bond(atom1, atom2, 1)
                    if atom1 not in assign.am_bonds.keys():
                        assign.am_bonds[atom1] = [atom2]
                    else:
                        assign.am_bonds[atom1].append(atom2)
                    if atom2 not in assign.am_bonds.keys():
                        assign.am_bonds[atom2] = [atom1]
                    else:
                        assign.am_bonds[atom2].append(atom1)
                else:
                    raise NotImplementedError(
                        "No implemented method to process bond #%s type %s" % (words[0], words[3]))

    _deal_with_ar_bonds(assign)
    assign.Determine_Ring_And_Bond_Type()
    return assign


def Get_Assignment_From_SDF(filename):
    with open(filename) as f:
        assign = Assign(f.readline().strip())
        f.readline()
        f.readline()
        line = f.readline()
        atom_numbers = int(line[:3].strip())
        bond_numbers = int(line[3:6].strip())
        for i in range(atom_numbers):
            line = f.readline()
            element = line[31:33].strip()
            x = float(line[:10])
            y = float(line[10:20])
            z = float(line[20:30])
            assign.Add_Atom(element, x, y, z)
        for i in range(bond_numbers):
            line = f.readline()
            atom1 = int(line[:3]) - 1
            atom2 = int(line[3:6]) - 1
            bond_order = int(line[6:9])
            assign.Add_Bond(atom1, atom2, bond_order)

    assign.Determine_Ring_And_Bond_Type()
    return assign


sys.modules["__main__"].__dict__["Get_Assignment_From_PubChem"] = Get_Assignment_From_PubChem
sys.modules["__main__"].__dict__["Get_Assignment_From_PDB"] = Get_Assignment_From_PDB
sys.modules["__main__"].__dict__["Get_Assignment_From_Mol2"] = Get_Assignment_From_Mol2
sys.modules["__main__"].__dict__["Get_Assignment_From_ResidueType"] = Get_Assignment_From_ResidueType
