from ... import *

AtomType.Add_Property({"charge": float})

AtomType.Set_Property_Unit("charge", "charge", "e")


@Molecule.Set_Save_SPONGE_Input("charge")
def write_charge(self):
    towrite = "%d\n" % (len(self.atoms))
    towrite += "\n".join(["%.6f" % (atom.charge * 18.2223) for atom in self.atoms])
    return towrite
