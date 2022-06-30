
"""
import Xponge
t = Get_Assignment_From_PubChem("ethylene glycol", "name")
#t = Get_Assignment_From_Mol2("PHE_ASN.mol2")
equivalence = [[0,1],[2,3],[4,5,6,7],[8,9]] #for butane [[0,1], [2,3],[5,6,7,4],[9,10,11,12,13,8]]
q = t.Calculate_Charge("RESP", extra_equivalence = equivalence, opt = False)

import Xponge.forcefield.AMBER.gaff
t.Determine_Atom_Type("GAFF")

PHE = t.To_ResidueType("PHE", q)

Save_Mol2(PHE)
Save_Mol2(t, "PHE_ASN.mol2")
Save_PDB(t, "PHE.pdb")
"""

"""
import Xponge.forcefield.CHARMM27.protein
t = NLYS + LYS + CLYS
Save_SPONGE_Input(t)


from Xponge.analysis.MDAnalysis import Universe
u = Universe("/home/dhx/Desktop/validation/sponge1.2/speed_test/system3/originace2.parm7",
"/home/dhx/Desktop/validation/sponge1.2/speed_test/system3/mdcrd.dat", "/home/dhx/Desktop/validation/sponge1.2/speed_test/system3/mdbox.txt")
O = u.select_atoms("resname WAT and name O")
from MDAnalysis.analysis import rdf as RDF
rdf = RDF.InterRDF(O,O)
rdf.run(verbose = True)
import matplotlib.pyplot as plt
plt.plot(rdf.results.bins[1:], rdf.results.rdf[1:])
plt.show()
"""

"""
import Xponge
import Xponge.forcefield.AMBER.ff14SB
import Xponge.forcefield.AMBER.tip3p

t = ACE + ALA + NME
t |= WAT

import Xponge.forcefield.SPECIAL.FEP as FEP
FEP.Intramolecule_NB_To_NB14(t, t.residues[:3])

free_t = FEP.Get_Free_Molecule(t, t.residues[:3])
tt = FEP_Merge(t, free_t, 0.5)

Save_SPONGE_Input(t, "FEP_test/A")
Save_SPONGE_Input(free_t, "FEP_test/B")

FEP.Save_Soft_Core_LJ()
Save_SPONGE_Input(tt, "FEP_test/AB")
"""
"""
import Xponge
import Xponge.forcefield.AMBER.gaff

from rdkit import  Chem
from rdkit.Chem import Draw, rdFMCS, rdMolAlign, AllChem

from Xponge.assign.RDKit_tools import * 
from Xponge.forcefield.SPECIAL import FEP


t = Get_Assignment_From_PubChem("ethane", "name")
t2 = Get_Assignment_From_PubChem("methane", "name")

import Xponge.forcefield.AMBER.gaff
import Xponge.forcefield.AMBER.tip3p

t.Determine_Atom_Type("GAFF")
t2.Determine_Atom_Type("GAFF")
restype = t.To_ResidueType("MOL")
restype2 = t2.To_ResidueType("MOL2")

rmol = Xponge.Molecule("test")

rmol.Add_Residue(restype)

rmol |= WAT

restypeR, restypeR2 = FEP.Temp_Function_Name(rmol, rmol.residues[0], restype2, t, t2)
"""

#力场
import Xponge
import Xponge.forcefield.AMBER.gaff
import Xponge.forcefield.SPECIAL.FEP as FEP
import Xponge.forcefield.AMBER.tip3p

#获得乙烷和甲烷单个residue的top
t = Get_Assignment_From_PubChem("Phenol", "name")
t2 = Get_Assignment_From_PubChem("Benzoic acid", "name")
t.Determine_Atom_Type("GAFF")
t2.Determine_Atom_Type("GAFF")
restype2 = t.To_ResidueType("MOL")
restype = t2.To_ResidueType("MOL2")

#建模，正常的话应该是loadpdb的，这儿简单举例
rmol = Xponge.Molecule("test")
rmol.Add_Residue(restype)

#溶剂化，简单举例，加了一个水
rmol |= WAT
#把水和乙烷坐标错开10埃
for atom in rmol.residues[1].atoms:
    atom.x += 10

#把分子、需要变的残基、需要变成的残疾类，上面的两个top信息（正常可以从Mol2得）
restypeR, restypeR2 = FEP.Merge_Dual_Topology(rmol, rmol.residues[0], restype2, t2, t)


#或者是获得”自由“的分子，乙烷直接消失的那种
#free_t = FEP.Get_Free_Molecule(rmol, rmol.residues[0])

#先保存原来的信息
Save_SPONGE_Input(rmol)
Save_SPONGE_Input(restypeR)
Save_SPONGE_Input(restypeR2)
Save_Mol2(restypeR)
Save_Mol2(restypeR2)

#获得lambda = 0.5的分子
FEP.Save_Soft_Core_LJ()
tt = FEP.Merge_Force_Field(restypeR, restypeR2, 0.5)
tt.name =  "testAB"

#保存lambda = 0.5的分子
Save_SPONGE_Input(tt)


