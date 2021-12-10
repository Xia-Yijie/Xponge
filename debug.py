import Xponge.forcefield.AMBER.gaff

t = Get_Molecule_From_PubChem("ethyl acetate", "name")
t.Determine_Atom_Type("GAFF")

#import Xponge.assign as assign

#qq = assign.ASSIGN()
#qq.Add_Atom("C", x=0, y=0, z=0)
q = RESP_Fit(t, "6-31g*", opt = True)

PHE = t.To_ResidueType("PHE", q)

Save_Mol2(PHE)