from ... import *


NB14Type = Generate_New_Bonded_Force_Type("nb14_extra", "1-4", {"epsilon": float, "rmin": float, "sigma":float, "A":float, "B":float, "kee": float}, False)

NB14Type.Set_Property_Unit("rmin", "distance", "A")
NB14Type.Set_Property_Unit("sigma", "distance", "A")
NB14Type.Set_Property_Unit("epsilon", "energy", "kcal/mol")
NB14Type.Set_Property_Unit("A", "energy·distance^6", "kcal/mol·A^6")
NB14Type.Set_Property_Unit("B", "energy·distance^12", "kcal/mol·A^12")

NB14Type.topology_matrix = [[1, -4],
                            [1, 1]]


def Get_NB14EXTRA_AB(atom1, atom2):
    from . import LJ 
    LJType = LJ.LJType

    A = 0
    B = 0
        
    LJ_i = LJType.types[atom1.LJtype + "-" + atom1.LJtype]
    LJ_j = LJType.types[atom2.LJtype + "-" + atom2.LJtype]
    finded = False
    findnames = [atom1.LJtype + "-" + atom2.LJtype, 
        atom2.LJtype + "-" + atom1.LJtype]
    
    for findname in findnames:
        if findname in LJType.types.keys():
            finded = True
            LJ_ij = LJType.types[findname]
            A = LJType.combining_method_A(LJ_ij.epsilon, LJ_ij.rmin, LJ_ij.epsilon, LJ_ij.rmin)
            B = LJType.combining_method_B(LJ_ij.epsilon, LJ_ij.rmin, LJ_ij.epsilon, LJ_ij.rmin)
            break
    if not finded:
            A = LJType.combining_method_A(LJ_i.epsilon, LJ_i.rmin, LJ_j.epsilon, LJ_j.rmin)
            B = LJType.combining_method_B(LJ_i.epsilon, LJ_i.rmin, LJ_j.epsilon, LJ_j.rmin)
    return A, B

def Exclude_To_NB14EXTRA(molecule, atom1, atom2):
    new_force = NB14Type.entity([atom1, atom2], NB14Type.types["UNKNOWNS"])
    A, B = Get_NB14EXTRA_AB(new_force.atoms[0], new_force.atoms[1])
    new_force.A = nb14_bond.kLJ * A
    new_force.B = nb14_bond.kLJ * B
    new_force.kee = nb14_bond.kee
    atom1.Extra_Exclude_Atom(atom2)
    
    molecule.Add_Bonded_Force(new_force)

def NB14_To_NB14EXTRA(molecule):

    #处理nb14
    #A、B中的nb14全部变为nb14_extra
    while molecule.bonded_forces.get("nb14",[]):
        nb14_bond = molecule.bonded_forces["nb14"].pop()
        new_force = NB14Type.entity(nb14_bond.atoms, NB14Type.types["UNKNOWNS"], nb14_bond.name)
        
        A, B = Get_NB14EXTRA_AB(new_force.atoms[0], new_force.atoms[1])
                
        new_force.A = nb14_bond.kLJ * A
        new_force.B = nb14_bond.kLJ * B
        new_force.kee = nb14_bond.kee
        molecule.Add_Bonded_Force(new_force)
        

@GlobalSetting.Add_Unit_Transfer_Function(NB14Type)
def LJ_Unit_Transfer(self):
    if self.sigma != None:
        self.rmin = self.sigma * (4 ** (1/12) / 2)
        self.sigma = None
    if self.rmin != None and self.epsilon != None:
        self.A = self.epsilon * (2 * self.rmin) ** 12
        self.B = self.epsilon * 2 * ((2 * self.rmin) ** 6)
        self.rmin = None
        self.epsilon = None
        

@Molecule.Set_Save_SPONGE_Input("nb14_extra")
def write_nb14(self):
    bonds = []
    for bond in self.bonded_forces.get("nb14_extra",[]):
        order = list(range(2))
        if bond.A != 0 and bond.B != 0 and  bond.kee != 0:
            if self.atom_index[bond.atoms[order[0]]] > self.atom_index[bond.atoms[order[-1]]]:
                temp_order = order[::-1]
            else:
                temp_order = order

            bonds.append("%d %d %.6e %.6e %.6e"%(self.atom_index[bond.atoms[temp_order[0]]]
            , self.atom_index[bond.atoms[temp_order[1]]], bond.A, bond.B, bond.kee))
    
    if (bonds):
        towrite = "%d\n"%len(bonds)
        bonds.sort(key = lambda x: list(map(int, x.split()[:2])))
        towrite += "\n".join(bonds)
        
        return towrite