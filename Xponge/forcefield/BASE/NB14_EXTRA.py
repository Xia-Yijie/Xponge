from ... import *

NB14Type = Generate_New_Bonded_Force_Type("nb14_extra", "1-4", {"epsilon": float, "rmin": float, "sigma":float, "A":float, "B":float, "kee": float}, False)

NB14Type.Set_Property_Unit("rmin", "distance", "A")
NB14Type.Set_Property_Unit("sigma", "distance", "A")
NB14Type.Set_Property_Unit("epsilon", "energy", "kcal/mol")
NB14Type.Set_Property_Unit("A", "energy·distance^6", "kcal/mol·A^6")
NB14Type.Set_Property_Unit("B", "energy·distance^12", "kcal/mol·A^12")

NB14Type.topology_matrix = [[1, -4],
                            [1, 1]]

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
        

@Molecule.Set_Save_SPONGE_Input
def write_nb14(self, prefix, dirname):
    bonds = []
    for bond in self.bonded_forces["nb14_extra"]:
        order = list(range(2))
        if bond.A != 0 and bond.B != 0 and  bond.kee != 0:
            if self.atom_index[bond.atoms[order[0]]] > self.atom_index[bond.atoms[order[-1]]]:
                temp_order = order[::-1]
            else:
                temp_order = order
            
            bonds.append("%d %d %f %f %f"%(self.atom_index[bond.atoms[temp_order[0]]]
            , self.atom_index[bond.atoms[temp_order[1]]], bond.A, bond.B, bond.kee))
    
    if (bonds):
        towrite = "%d\n"%len(bonds)
        bonds.sort(key = lambda x: list(map(int, x.split()[:2])))
        towrite += "\n".join(bonds)
        
        f = open(os.path.join(dirname, prefix + "_nb14_extra.txt"),"w")
        f.write(towrite)
        f.close()