from ... import *

AtomType.Add_Property({"mass":float})

@Molecule.Set_Save_SPONGE_Input("mass")  
def write_mass(self):
    towrite = "%d\n"%(len(self.atoms))
    towrite += "\n".join(["%.3f"%(atom.mass) for atom in self.atoms])
    return towrite


