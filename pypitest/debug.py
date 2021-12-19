
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
"""

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
