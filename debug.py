import Xponge.forcefield.AMBER.gaff

t = Get_Molecule_From_PubChem("Phenylalanine", "name")
t.Determine_Atom_Type("GAFF")

q = RESP_Fit(t, "6-31+g*", opt = False, charge = -1)
print(q)