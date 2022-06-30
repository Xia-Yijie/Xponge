from ... import *

BondType = Generate_New_Bonded_Force_Type("bond", "1-2", {"k":float, "b":float}, True)

BondType.Set_Property_Unit("k", "energy·distance^-2", "kcal/mol·A^-2")
BondType.Set_Property_Unit("b", "distance", "A")

@Molecule.Set_Save_SPONGE_Input("bond")
def write_bond(self):
    bonds = []
    for bond in self.bonded_forces.get("bond",[]):
        order = list(range(2))
        if bond.k != 0:
            if self.atom_index[bond.atoms[order[0]]] > self.atom_index[bond.atoms[order[-1]]]:
                temp_order = order[::-1]
            else:
                temp_order = order
            bonds.append("%d %d %f %f"%(self.atom_index[bond.atoms[temp_order[0]]]
            , self.atom_index[bond.atoms[temp_order[1]]], bond.k, bond.b))
    
    if (bonds):
        towrite = "%d\n"%len(bonds)
        bonds.sort(key = lambda x: list(map(int, x.split()[:2])))
        towrite += "\n".join(bonds)
        
        return towrite