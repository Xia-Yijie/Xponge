from . import *

import sys

atoms, bonds, angles, propers, impropers, LJs = LOAD.parmdat(os.path.join(AMBER_DATA_DIR, "parm19.dat"))

AtomType.New_From_String(atoms)
BOND.BondType.New_From_String(bonds)
ANGLE.AngleType.New_From_String(angles)
DIHEDRAL.ProperType.New_From_String(propers)
DIHEDRAL.ImproperType.New_From_String(impropers)
LJ.LJType.New_From_String(LJs)


atoms, bonds, angles, propers, impropers, LJs, cmap = LOAD.frcmod(os.path.join(AMBER_DATA_DIR, "ff19SB.frcmod"))

AtomType.New_From_String(atoms)
BOND.BondType.New_From_String(bonds)
ANGLE.AngleType.New_From_String(angles)
DIHEDRAL.ProperType.New_From_String(propers)
DIHEDRAL.ImproperType.New_From_String(impropers)
LJ.LJType.New_From_String(LJs)

from ..BASE import RCMAP
RCMAP.CMAP.Residue_Map.update(cmap)

RCMAP.CMAP.New_From_String(r"""
name
C-N-XC-C-N
""")


ff19SB = LOAD.mol2(os.path.join(AMBER_DATA_DIR, "ff19SB.mol2"))
ResidueType.types["HIS"] = ResidueType.types["HIE"]
ResidueType.types["NHIS"] = ResidueType.types["NHIE"]
ResidueType.types["CHIS"] = ResidueType.types["CHIE"]
sys.modules['__main__'].__dict__["HIS"] = ResidueType.types["HIS"] 
sys.modules['__main__'].__dict__["NHIS"] = ResidueType.types["NHIS"] 
sys.modules['__main__'].__dict__["CHIS"] = ResidueType.types["CHIS"] 

residues = "ALA ARG ASN ASP CYS CYX GLN GLU GLY HID HIE HIP ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HIS".split()

for res in residues:
    ResidueType.types[res].head_next = "CA"
    ResidueType.types[res].head_length = 1.3
    ResidueType.types[res].tail_next = "CA"
    ResidueType.types[res].tail_length = 1.3
    ResidueType.types[res].head_link_conditions.append({"atoms":["CA", "N"], "parameter": 120/180 * np.pi})
    ResidueType.types[res].head_link_conditions.append({"atoms":["H", "CA", "N"], "parameter": -np.pi})
    ResidueType.types[res].tail_link_conditions.append({"atoms":["CA", "C"], "parameter": 120/180 * np.pi})
    ResidueType.types[res].tail_link_conditions.append({"atoms":["O", "CA", "C"], "parameter": -np.pi})    


    ResidueType.types["N" + res].tail_next = "CA"
    ResidueType.types["N" + res].tail_length = 1.3
    ResidueType.types["N" + res].tail_link_conditions.append({"atoms":["CA", "C"], "parameter": 120/180 * np.pi})
    ResidueType.types["N" + res].tail_link_conditions.append({"atoms":["O", "CA", "C"], "parameter": -np.pi})   
    
    ResidueType.types["C" + res].head_next = "CA"
    ResidueType.types["C" + res].head_length = 1.3
    ResidueType.types["C" + res].head_link_conditions.append({"atoms":["CA", "N"], "parameter": 120/180 * np.pi})
    ResidueType.types["C" + res].head_link_conditions.append({"atoms":["H", "CA", "N"], "parameter": -np.pi})
    
    GlobalSetting.Add_PDB_Residue_Name_Mapping("head", res, "N" + res)
    GlobalSetting.Add_PDB_Residue_Name_Mapping("tail", res, "C" + res)

ResidueType.types["ACE"].tail_next = "CH3"
ResidueType.types["ACE"].tail_length = 1.3
ResidueType.types["ACE"].tail_link_conditions.append({"atoms":["CH3", "C"], "parameter": 120/180 * np.pi})
ResidueType.types["ACE"].tail_link_conditions.append({"atoms":["O", "CH3", "C"], "parameter": -np.pi})   

ResidueType.types["NME"].head_next = "CH3"
ResidueType.types["NME"].head_length = 1.3
ResidueType.types["NME"].head_link_conditions.append({"atoms":["CH3", "N"], "parameter": 120/180 * np.pi})
ResidueType.types["NME"].head_link_conditions.append({"atoms":["H", "CH3", "N"], "parameter": -np.pi})   

GlobalSetting.HISMap["DeltaH"] = "HD1"
GlobalSetting.HISMap["EpsilonH"] = "HE2"
GlobalSetting.HISMap["HIS"].update({"HIS": {"HID":"HID", "HIE":"HIE", "HIP":"HIP"}, 
                                    "CHIS":{"HID":"CHID", "HIE":"CHIE", "HIP":"CHIP"},
                                    "NHIS":{"HID":"NHID", "HIE":"NHIE", "HIP":"NHIP"}})

ResidueType.types["CYX"].connect_atoms["ssbond"] = "SG"


print("""Reference for ff19SB.py:
  Chuan Tian, Koushik Kasavajhala, Kellon A. A. Belfon, Lauren Raguette, He Huang, Angela N. Migues, John Bickel, Yuzhang Wang, Jorge Pincay, Qin Wu, and Carlos Simmerling
    ff19SB: Amino-Acid-Specific Protein Backbone Parameters Trained against Quantum Mechanics Energy Surfaces in Solution
    The Journal of Chemical Physics 2020 16 (1), 528-552, 
    DOI: 10.1021/acs.jctc.9b00591
""")