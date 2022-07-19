from . import *
import sys

atoms, bonds, angles, propers, impropers, LJs = load_parmdat(os.path.join(AMBER_DATA_DIR, "lipid14.dat"))

AtomType.New_From_String(atoms)
BOND.BondType.New_From_String(bonds)
ANGLE.AngleType.New_From_String(angles)
DIHEDRAL.ProperType.New_From_String(propers)
DIHEDRAL.ImproperType.New_From_String(impropers)
LJ.LJType.New_From_String(LJs)

lipid14 = load_mol2(os.path.join(AMBER_DATA_DIR, "lipid14.mol2"))

for res in "LA PA MY OL".split():
    ResidueType.types[res].head = "C12"
    ResidueType.types[res].tail = "C12"
    ResidueType.types[res].head_next = "C13"
    ResidueType.types[res].tail_next = "C13"
    ResidueType.types[res].head_length = 1.5
    ResidueType.types[res].tail_length = 1.5
    ResidueType.types[res].head_link_conditions.append({"atoms": ["H2R", "C12"], "parameter": 109.5 / 180 * np.pi})
    ResidueType.types[res].head_link_conditions.append(
        {"atoms": ["H2S", "H2R", "C12"], "parameter": -120 / 180 * np.pi})
    ResidueType.types[res].tail_link_conditions.append({"atoms": ["H2R", "C12"], "parameter": 109.5 / 180 * np.pi})
    ResidueType.types[res].tail_link_conditions.append(
        {"atoms": ["H2S", "H2R", "C12"], "parameter": -120 / 180 * np.pi})

for res in "PC PE".split():
    ResidueType.types[res].head = "C11"
    ResidueType.types[res].tail = "C21"
    ResidueType.types[res].head_next = "O11"
    ResidueType.types[res].tail_next = "O21"
    ResidueType.types[res].head_length = 1.5
    ResidueType.types[res].tail_length = 1.5
    ResidueType.types[res].head_link_conditions.append({"atoms": ["O11", "C11"], "parameter": 120 / 180 * np.pi})
    ResidueType.types[res].head_link_conditions.append({"atoms": ["O12", "O11", "C11"], "parameter": np.pi})
    ResidueType.types[res].tail_link_conditions.append({"atoms": ["O21", "C21"], "parameter": 120 / 180 * np.pi})
    ResidueType.types[res].tail_link_conditions.append({"atoms": ["O22", "O21", "C21"], "parameter": np.pi})

print("""Reference for lipid14:
  Dickson, C.J., Madej, B.D., Skjevik, A.A., Betz, R.M., Teigen, K., Gould, I.R., Walker, R.C. 
    Lipid14: The Amber Lipid Force Field.
    Journal of Chemical Theory and Computation 2014 10(2), 865-879,
    DOI: 10.1021/ct4010307
""")
