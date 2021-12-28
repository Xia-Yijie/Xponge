from ... import *

AtomType.Add_Property({"charge":float})

AtomType.Set_Property_Unit("charge", "charge", "e")

@Molecule.Set_Save_SPONGE_Input      
def write_charge(self, prefix, dirname):
    towrite = "%d\n"%(len(self.atoms))
    towrite += "\n".join(["%.6f"%(atom.charge * 18.2223) for atom in self.atoms])
    f = open(os.path.join(dirname, prefix + "_charge.txt"),"w")
    f.write(towrite)
    f.close()
