from ... import *

VirtualType2 = Generate_New_Bonded_Force_Type("vatom2", "1", {"atom0":int, "atom1":int, "atom2":int, "k1":float, "k2":float}, False)
GlobalSetting.VirtualAtomTypes["vatom2"] = 3

@Molecule.Set_Save_SPONGE_Input("virtual_atom")
def write_virtual_atoms(self):
    vatoms = []
    for vatom in self.bonded_forces.get("vatom2",[]):
        vatoms.append("2 %d %d %d %d %f %f"%(self.atom_index[vatom.atoms[0]],
            self.atom_index[vatom.atoms[0]] + vatom.atom0, 
            self.atom_index[vatom.atoms[0]] + vatom.atom1,
            self.atom_index[vatom.atoms[0]] + vatom.atom2,
            vatom.k1, vatom.k2))
    
    if (vatoms):
        towrite = ""
        towrite += "\n".join(vatoms)
        
        return towrite