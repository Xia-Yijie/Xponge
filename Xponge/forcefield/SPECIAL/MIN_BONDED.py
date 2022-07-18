from ... import *


def write_mass(self):
    towrite = "%d\n" % (len(self.atoms))
    towrite += "\n".join(
        ["%.3f" % (1) if (atom.mass < 3.999 or atom.LJtype == "ZERO_LJ_ATOM") else "%.3f" % (0) for atom in self.atoms])
    return towrite


def write_LJ(self):
    towrite = "%d %d\n\n" % (len(self.atoms), 1)
    for i in range(1):
        for j in range(1):
            towrite += "%16.7e" % 0 + " "
        towrite += "\n"
    towrite += "\n"

    for i in range(1):
        for j in range(1):
            towrite += "%16.7e" % 0 + " "
        towrite += "\n"
    towrite += "\n"
    towrite += "\n".join(["%d" % (0) for atom in self.atoms])
    return towrite


def write_charge(self):
    towrite = "%d\n" % (len(self.atoms))
    towrite += "\n".join(["%.6f" % (0) for atom in self.atoms])
    return towrite


def Save_Min_Bonded_Parameters():
    Molecule.Set_Save_SPONGE_Input("fake_mass")(write_mass)
    Molecule.Set_Save_SPONGE_Input("fake_LJ")(write_LJ)
    Molecule.Set_Save_SPONGE_Input("fake_charge")(write_charge)


def Do_Not_Save_Min_Bonded_Parameters():
    Molecule.Del_Save_SPONGE_Input("fake_mass")
    Molecule.Del_Save_SPONGE_Input("fake_LJ")
    Molecule.Del_Save_SPONGE_Input("fake_charge")
