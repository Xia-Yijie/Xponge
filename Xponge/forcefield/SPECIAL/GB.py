from ... import *
from ...assign import Guess_Element_From_Mass

AtomType.Add_Property({"born_radii": float})
AtomType.Set_Property_Unit("born_radii", "distance", "A")

bonded_radii = {"H":1.2, "C":1.7, "N":1.55, "O":1.52, "F": 1.47, 
                         "P": 1.8, "S": 1.8, "Cl": 1.75, "Br":1.85, "I":1.98, 
                         "Reference": "Reference todo"}

modified_bonded_radii = {"H":1.3, "C":1.7, "N":1.55, "O":1.52, "F": 1.47, 
                         "P": 1.8, "S": 1.8, "Cl": 1.75, "Br":1.85, "I":1.98, 
                         "Reference": "Reference todo"}

available_radii = {"bonded_radii": bonded_radii, "modified_bonded_radii":modified_bonded_radii}

def write_born_radii(self):
    towrite = "%d\n"%(len(self.atoms))
    towrite += "\n".join(["%.4f"%(atom.born_radii) for atom in self.atoms])
    return towrite

def Set_Born_Radii(mol, radii = modified_bonded_radii):
    if radii.get("Reference", None) is not None:
        print(radii["Reference"])
    if isinstance(mol, Molecule):
        for residue in mol.residues:
            for atom in residue.atoms:
               element = Guess_Element_From_Mass(atom.mass)
               atom.born_radii = modified_bonded_radii.get(element)
    else:
        for atom in mol.atoms:
            element = Guess_Element_From_Mass(atom.mass)
            atom.born_radii = modified_bonded_radii.get(element)
    Molecule.Set_Save_SPONGE_Input("born_radii")(write_born_radii)    


