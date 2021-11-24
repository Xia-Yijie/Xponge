from ... import *

AtomType.Add_Property({"mass":float})

@Molecule.Set_Save_SPONGE_Input      
def write_mass(self, prefix, dirname):
    towrite = "%d\n"%(len(self.atoms))
    towrite += "\n".join(["%.3f"%(atom.mass) for atom in self.atoms])
    f = open(os.path.join(dirname, prefix + "_mass.txt"),"w")
    f.write(towrite)
    f.close()
