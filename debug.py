import Xponge
t = Get_Assignment_From_PubChem("butane", "name")
#t = Get_Assignment_From_Mol2("PHE_ASN.mol2")
q = t.Calculate_Charge("RESP", extra_equivalence = [[0,1], [2,3],[5,6,7,4],[9,10,11,12,13,8]])

import Xponge.forcefield.AMBER.gaff
t.Determine_Atom_Type("GAFF")

PHE = t.To_ResidueType("PHE", q)

Save_Mol2(PHE)
Save_Mol2(t, "PHE_ASN.mol2")
Save_PDB(t, "PHE.pdb")