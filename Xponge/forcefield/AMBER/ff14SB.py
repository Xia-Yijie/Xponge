from . import *

import sys

atoms, bonds, angles, propers, impropers, LJs = LOAD.parmdat(os.path.join(AMBER_DATA_DIR, "parm10.dat"))

AtomType.New_From_String(atoms)
BOND.BondType.New_From_String(bonds)
ANGLE.AngleType.New_From_String(angles)
DIHEDRAL.ProperType.New_From_String(propers)
DIHEDRAL.ImproperType.New_From_String(impropers)
LJ.LJType.New_From_String(LJs)


atoms, bonds, angles, propers, impropers, LJs, cmap = LOAD.frcmod(os.path.join(AMBER_DATA_DIR, "ff14SB.frcmod"))

AtomType.New_From_String(atoms)
BOND.BondType.New_From_String(bonds)
ANGLE.AngleType.New_From_String(angles)
DIHEDRAL.ProperType.New_From_String(propers)
DIHEDRAL.ImproperType.New_From_String(impropers)
LJ.LJType.New_From_String(LJs)

ff14SB = LOAD.mol2(os.path.join(AMBER_DATA_DIR, "ff14SB.mol2"))

residues = "ALA ARG ASN ASP CYS CYX GLN GLU GLY HID HIE HIP ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HIS".split()

GlobalSetting.PDBResidueNameMap["head"].update({resname:"N" + resname for resname in residues})
GlobalSetting.PDBResidueNameMap["tail"].update({resname:"C" + resname for resname in residues})

GlobalSetting.HISMap["DeltaH"] = "HD1"
GlobalSetting.HISMap["EpsilonH"] = "HE2"
GlobalSetting.HISMap["HIS"].update({"HIS": {"HID":"HID", "HIE":"HIE", "HIP":"HIP"}, 
                                    "CHIS":{"HID":"CHID", "HIE":"CHIE", "HIP":"CHIP"},
                                    "NHIS":{"HID":"NHID", "HIE":"NHIE", "HIP":"NHIP"}})

ResidueType.types["CYX"].connect_atoms["ssbond"] = "SG"


ResidueType.types["HIS"] = ResidueType.types["HIE"]