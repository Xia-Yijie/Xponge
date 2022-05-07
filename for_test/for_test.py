import Xponge
import Xponge.forcefield.AMBER.ff14SB
import Xponge.forcefield.AMBER.tip3p

t = ACE + ALA * 10 + NME

for i in range(1,len(t.residues)-1):
    head = t.residues[i-1]
    res = t.residues[i]
    tail = t.residues[i+1]
    Impose_Dihedral(t, head.C, res.N, res.CA, res.C, -3.1415926/3)
    Impose_Dihedral(t, res.N, res.CA, res.C, tail.N, -3.1415926/3)

Save_Mol2(t, "test.mol2")
c = int(round(t.charge))
Process_Box(t, WAT, 10)
Ion_Replace(t, lambda res: res.type.name == "WAT", {CL:30 + c, K:30})
t.residues.sort(key = lambda residue: {"CL":2, "K":1, "WAT":3}.get(residue.type.name, 0))
Save_PDB(t, "test.pdb")
Save_SPONGE_Input(t, "test")