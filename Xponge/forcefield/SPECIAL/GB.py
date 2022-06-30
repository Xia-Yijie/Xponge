from ... import *
from ...assign import Guess_Element_From_Mass

AtomType.Add_Property({"GB_radii": float})
AtomType.Set_Property_Unit("GB_radii", "distance", "A")
AtomType.Add_Property({"GB_scaler": float})


def Bondi_radii(atom):
    if isinstance(atom,str) and atom == "ref":
        print("""Reference for Bondi radii:
    A. Bondi
    van der Waals Volumes and Radii
    Journal of Physical Chemistry 1964 68 (3) 441-451
    DOI: 10.1021/j100785a001
""")
    else:
        temp_dict = {"H":1.2, "C":1.7, "N":1.55, "O":1.52, "F": 1.47, 
                     "P": 1.8, "S": 1.8, "Cl": 1.75, "Br":1.85, "I":1.98}
        element = Guess_Element_From_Mass(atom.mass)
        atom.GB_radii = temp_dict.get(element, 1.5)
        temp_dict2 = {"H":0.85, "C":0.72, "N":0.79, "O":0.85, "F": 0.88, 
                     "P": 0.86, "S": 0.96}
        atom.GB_scaler = temp_dict2.get(element, 0.8)

def modified_Bondi_radii(atom):
    if isinstance(atom,str) and atom == "ref":
        print("""Reference for modified Bondi radii:
    Vickie Tsui, David A. Case
    Theory and Applications of the Generalized Born Solvation Model in Macromolecular Simulations
    Biopolymers 2001 56 (4) 275-291
    DOI: 10.1002/1097-0282(2000)56:4<275::AID-BIP10024>3.0.CO;2-E
""")
    else:
        temp_dict = {"H":1.2, "C":1.7, "N":1.55, "O":1.5, "F": 1.5, 
     "Si":2.1, "P": 1.85, "S": 1.8, "Cl": 1.7, "Br":1.85, "I":1.98}
        element = Guess_Element_From_Mass(atom.mass)
        atom.GB_radii = temp_dict.get(element, 1.5)
        temp_dict2 = {"H":0.85, "C":0.72, "N":0.79, "O":0.85, "F": 0.88, 
                     "P": 0.86, "S": 0.96}
        atom.GB_scaler = temp_dict2.get(element, 0.8)
        if element == "H":
            for atomA in atom.residue.connectivity[atom]:
                elementA = Guess_Element_From_Mass(atomA.mass)
                if elementA == "C" or elementA == "N":
                    atom.GB_radii = 1.3
                elif elementA == "S" or elementA == "O" or elementA == "H":
                    atom.GB_radii = 0.8
                break

available_radius_set = {"Bondi_radii": Bondi_radii, "modified_Bondi_radii": modified_Bondi_radii}

def write_GB_radii_and_scaler(self):
    towrite = "%d\n"%(len(self.atoms))
    towrite += "\n".join(["%.4f %.4f"%(atom.GB_radii, atom.GB_scaler) for atom in self.atoms])
    return towrite

def Set_GB_Radii(mol, radius_set = modified_Bondi_radii):
    radius_set("ref")
    if isinstance(mol, Molecule):
        for residue in mol.residues:
            for atom in residue.atoms:
               radius_set(atom)
    else:
        for atom in mol.atoms:
            radius_set(atom)
    mol.box_length = [999,999,999]
    Molecule.Set_Save_SPONGE_Input("gb")(write_GB_radii_and_scaler)    


