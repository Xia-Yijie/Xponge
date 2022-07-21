"""
This **module** set the basic configuration for ff14sb
"""
from ...helper import source, set_real_global_variable, Xprint

source("....")
amber = source("...amber")

amber.load_parameters_from_parmdat("parm10.dat")
amber.load_parameters_from_frcmod("ff14SB.frcmod")

load_mol2(os.path.join(AMBER_DATA_DIR, "ff14SB.mol2"))

ResidueType.types["HIS"] = ResidueType.types["HIE"]
ResidueType.types["NHIS"] = ResidueType.types["NHIE"]
ResidueType.types["CHIS"] = ResidueType.types["CHIE"]

set_real_global_variable("HIS", ResidueType.types["HIS"])
set_real_global_variable("NHIS", ResidueType.types["NHIS"])
set_real_global_variable("CHIS", ResidueType.types["CHIS"])

residues = "ALA ARG ASN ASP CYS CYX GLN GLU GLY HID HIE HIP ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HIS".split()

for res in residues:
    ResidueType.types[res].head_next = "CA"
    ResidueType.types[res].head_length = 1.3
    ResidueType.types[res].tail_next = "CA"
    ResidueType.types[res].tail_length = 1.3
    ResidueType.types[res].head_link_conditions.append({"atoms": ["CA", "N"], "parameter": 120 / 180 * np.pi})
    ResidueType.types[res].head_link_conditions.append({"atoms": ["H", "CA", "N"], "parameter": -np.pi})
    ResidueType.types[res].tail_link_conditions.append({"atoms": ["CA", "C"], "parameter": 120 / 180 * np.pi})
    ResidueType.types[res].tail_link_conditions.append({"atoms": ["O", "CA", "C"], "parameter": -np.pi})

    ResidueType.types["N" + res].tail_next = "CA"
    ResidueType.types["N" + res].tail_length = 1.3
    ResidueType.types["N" + res].tail_link_conditions.append({"atoms": ["CA", "C"], "parameter": 120 / 180 * np.pi})
    ResidueType.types["N" + res].tail_link_conditions.append({"atoms": ["O", "CA", "C"], "parameter": -np.pi})

    ResidueType.types["C" + res].head_next = "CA"
    ResidueType.types["C" + res].head_length = 1.3
    ResidueType.types["C" + res].head_link_conditions.append({"atoms": ["CA", "N"], "parameter": 120 / 180 * np.pi})
    ResidueType.types["C" + res].head_link_conditions.append({"atoms": ["H", "CA", "N"], "parameter": -np.pi})

    GlobalSetting.Add_PDB_Residue_Name_Mapping("head", res, "N" + res)
    GlobalSetting.Add_PDB_Residue_Name_Mapping("tail", res, "C" + res)

ResidueType.types["ACE"].tail_next = "CH3"
ResidueType.types["ACE"].tail_length = 1.3
ResidueType.types["ACE"].tail_link_conditions.append({"atoms": ["CH3", "C"], "parameter": 120 / 180 * np.pi})
ResidueType.types["ACE"].tail_link_conditions.append({"atoms": ["O", "CH3", "C"], "parameter": -np.pi})

ResidueType.types["NME"].head_next = "CH3"
ResidueType.types["NME"].head_length = 1.3
ResidueType.types["NME"].head_link_conditions.append({"atoms": ["CH3", "N"], "parameter": 120 / 180 * np.pi})
ResidueType.types["NME"].head_link_conditions.append({"atoms": ["H", "CH3", "N"], "parameter": -np.pi})

GlobalSetting.HISMap["DeltaH"] = "HD1"
GlobalSetting.HISMap["EpsilonH"] = "HE2"
GlobalSetting.HISMap["HIS"].update({"HIS": {"HID": "HID", "HIE": "HIE", "HIP": "HIP"},
                                    "CHIS": {"HID": "CHID", "HIE": "CHIE", "HIP": "CHIP"},
                                    "NHIS": {"HID": "NHID", "HIE": "NHIE", "HIP": "NHIP"}})

ResidueType.types["CYX"].connect_atoms["ssbond"] = "SG"

Xprint("""Reference for ff14SB:
  James A. Maier, Carmenza Martinez, Koushik Kasavajhala, Lauren Wickstrom, Kevin E. Hauser, and Carlos Simmerling
    ff14SB: Improving the accuracy of protein side chain and backbone parameters from ff99SB
    Journal of Chemical Theory and Computation 2015 11 (8), 3696-3713
    DOI: 10.1021/acs.jctc.5b00255
""")
