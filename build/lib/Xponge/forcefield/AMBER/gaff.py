from . import *
import sys

from ... import assign

atoms, bonds, angles, propers, impropers, LJs = LOAD.parmdat(os.path.join(AMBER_DATA_DIR, "gaff.dat"))

AtomType.New_From_String(atoms)
BOND.BondType.New_From_String(bonds)
ANGLE.AngleType.New_From_String(angles)
DIHEDRAL.ProperType.New_From_String(propers)
DIHEDRAL.ImproperType.New_From_String(impropers)
LJ.LJType.New_From_String(LJs)

#print(AtomType.types)
GAFF = assign.Judge_Rule("GAFF")

@GAFF.Add_Judge_Rule("cx")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "C4") and "RG3" in Assign.atom_marker[i].keys()

@GAFF.Add_Judge_Rule("cy")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "C4") and "RG4" in Assign.atom_marker[i].keys()

@GAFF.Add_Judge_Rule("c3")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "C4")

@GAFF.Add_Judge_Rule("c")
def temp(i, Assign):
    tofind = False
    for bonded_atom in Assign.bonds[i].keys():
        if Assign.Atom_Judge(bonded_atom, "O1") or Assign.Atom_Judge(bonded_atom, "S1"):
            tofind = True
            break
    return Assign.Atom_Judge(i, "C3") and tofind

@GAFF.Add_Judge_Rule("cz")
def temp(i, Assign):
    tofind = False
    for bonded_atom in Assign.bonds[i].keys():
        if Assign.Atom_Judge(bonded_atom, "N3"):
            tofind = True
            break
            
    return Assign.Atom_Judge(i, "C3") and tofind

@GAFF.Add_Judge_Rule("cp")
def temp(i, Assign):
    tofind = Assign.Atom_Judge(i, "C3") and "AR1" in Assign.atom_marker[i] and Assign.atom_marker[i].get("RG6") == 1
    for bonded_atom in Assign.bonds[i].keys():
        if not tofind:
            break
        if Assign.atoms[bonded_atom] not in Assign.XX or "AR1" not in Assign.atom_marker[bonded_atom].keys():
            tofind = False
            
    return tofind

@GAFF.Add_Judge_Rule("ca")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "C3") and "AR1" in Assign.atom_marker[i]

@GAFF.Add_Judge_Rule("cc")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "C3") and "sb" in Assign.atom_marker[i] and "db" in Assign.atom_marker[i] and ("AR2" in Assign.atom_marker[i] or "AR3" in Assign.atom_marker[i])

@GAFF.Add_Judge_Rule("ce")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "C3") and "sb" in Assign.atom_marker[i] and "db" in Assign.atom_marker[i]

@GAFF.Add_Judge_Rule("cu")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "C3") and "RG3" in Assign.atom_marker[i].keys()

@GAFF.Add_Judge_Rule("cv")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "C3") and "RG4" in Assign.atom_marker[i].keys()

@GAFF.Add_Judge_Rule("c2")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "C3")

@GAFF.Add_Judge_Rule("cg")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "C3") and "sb" in Assign.atom_marker[i] and "tb" in Assign.atom_marker[i]

@GAFF.Add_Judge_Rule("c1")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "C2") or Assign.Atom_Judge(i, "C1")

@GAFF.Add_Judge_Rule("hn")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "H1") and Assign.atoms[list(Assign.bonds[i].keys())[0]] == "N"

@GAFF.Add_Judge_Rule("ho")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "H1") and Assign.atoms[list(Assign.bonds[i].keys())[0]] == "O"

@GAFF.Add_Judge_Rule("hs")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "H1") and Assign.atoms[list(Assign.bonds[i].keys())[0]] == "S"

@GAFF.Add_Judge_Rule("hp")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "H1") and Assign.atoms[list(Assign.bonds[i].keys())[0]] == "P"

@GAFF.Add_Judge_Rule("hx")
def temp(i, Assign):
    tofind = False
    for bonded_atom in Assign.bonds[i].keys():
        if Assign.atoms[bonded_atom] == "C":
            for bonded_atom_bonded in Assign.bonds[bonded_atom].keys():
                if Assign.Atom_Judge(bonded_atom_bonded, "N4"):
                    tofind = True
                    break
    return Assign.Atom_Judge(i, "H1") and tofind

@GAFF.Add_Judge_Rule("hw")
def temp(i, Assign):
    tofind = False
    for bonded_atom in Assign.bonds[i].keys():
        if Assign.atoms[bonded_atom] == "O":
            for bonded_atom_bonded in Assign.bonds[bonded_atom].keys():
                if Assign.Atom_Judge(bonded_atom_bonded, "H1"):
                    tofind = True
                    break
    return Assign.Atom_Judge(i, "H1") and tofind

@GAFF.Add_Judge_Rule("h3")
def temp(i, Assign):
    tofind = 0
    for bonded_atom in Assign.bonds[i].keys():
        if Assign.Atom_Judge(bonded_atom, "C4"):
            for bonded_atom_bonded in Assign.bonds[bonded_atom].keys():
                if Assign.atoms[bonded_atom_bonded] in Assign.XE:
                    tofind += 1
                    
    return Assign.Atom_Judge(i, "H1") and tofind == 3

@GAFF.Add_Judge_Rule("h2")
def temp(i, Assign):
    tofind = 0
    for bonded_atom in Assign.bonds[i].keys():
        if Assign.Atom_Judge(bonded_atom, "C4"):
            for bonded_atom_bonded in Assign.bonds[bonded_atom].keys():
                if Assign.atoms[bonded_atom_bonded] in Assign.XE:
                    tofind += 1                    
    return Assign.Atom_Judge(i, "H1") and tofind == 2

@GAFF.Add_Judge_Rule("h1")
def temp(i, Assign):
    tofind = 0
    for bonded_atom in Assign.bonds[i].keys():
        if Assign.Atom_Judge(bonded_atom, "C4"):
            for bonded_atom_bonded in Assign.bonds[bonded_atom].keys():
                if Assign.atoms[bonded_atom_bonded] in Assign.XE:
                    tofind += 1                    
    return Assign.Atom_Judge(i, "H1") and tofind == 1

@GAFF.Add_Judge_Rule("hc")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "H1") and Assign.Atom_Judge(list(Assign.bonds[i].keys())[0], "C4")
 
@GAFF.Add_Judge_Rule("h5")
def temp(i, Assign):
    tofind = 0
    for bonded_atom in Assign.bonds[i].keys():
        if Assign.Atom_Judge(bonded_atom, "C3"):
            for bonded_atom_bonded in Assign.bonds[bonded_atom].keys():
                if Assign.atoms[bonded_atom_bonded] in Assign.XE:
                    tofind += 1                    
    return Assign.Atom_Judge(i, "H1") and tofind == 2

@GAFF.Add_Judge_Rule("h4")
def temp(i, Assign):
    tofind = 0
    for bonded_atom in Assign.bonds[i].keys():
        if Assign.Atom_Judge(bonded_atom, "C3"):
            for bonded_atom_bonded in Assign.bonds[bonded_atom].keys():
                if Assign.atoms[bonded_atom_bonded] in Assign.XE:
                    tofind += 1                    
    return Assign.Atom_Judge(i, "H1") and tofind == 1

@GAFF.Add_Judge_Rule("ha")
def temp(i, Assign):                  
    return Assign.Atom_Judge(i, "H1")

@GAFF.Add_Judge_Rule("f")
def temp(i, Assign):                  
    return Assign.atoms[i] == "F"

@GAFF.Add_Judge_Rule("cl")
def temp(i, Assign):                  
    return Assign.atoms[i] == "Cl"

@GAFF.Add_Judge_Rule("br")
def temp(i, Assign):                  
    return Assign.atoms[i] == "Br"

@GAFF.Add_Judge_Rule("i")
def temp(i, Assign):                  
    return Assign.atoms[i] == "I"
    
@GAFF.Add_Judge_Rule("o")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "O1")

@GAFF.Add_Judge_Rule("op")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "O2") and "RG3" in Assign.atom_marker[i].keys()

@GAFF.Add_Judge_Rule("oq")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "O2") and "RG4" in Assign.atom_marker[i].keys()

@GAFF.Add_Judge_Rule("oh")
def temp(i, Assign):
    tofind = False
    if Assign.Atom_Judge(i, "O2") or Assign.Atom_Judge(i, "O3"):
        for bonded_atom in Assign.bonds[i].keys():
            if Assign.Atom_Judge(bonded_atom, "H1"):
                tofind = True
                break    
    return tofind

@GAFF.Add_Judge_Rule("os")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "O2")

@GAFF.Add_Judge_Rule("ni")
def temp(i, Assign):
    tofind = False
    if Assign.Atom_Judge(i, "N3") and "RG3" in Assign.atom_marker[i].keys():
        for bonded_atom in Assign.bonds[i].keys():
            if Assign.Atom_Judge(bonded_atom, "C3"):
                for bonded_atom_bonded in Assign.bonds[bonded_atom].keys():
                    if Assign.atoms[bonded_atom_bonded] in Assign.XA:
                        tofind = True
                        break
            if tofind:
                break
    return tofind

@GAFF.Add_Judge_Rule("nj")
def temp(i, Assign):
    tofind = False
    if Assign.Atom_Judge(i, "N3") and "RG4" in Assign.atom_marker[i].keys():
        for bonded_atom in Assign.bonds[i].keys():
            if Assign.Atom_Judge(bonded_atom, "C3"):
                for bonded_atom_bonded in Assign.bonds[bonded_atom].keys():
                    if Assign.atoms[bonded_atom_bonded] in Assign.XA:
                        tofind = True
                        break
            if tofind:
                break
    return tofind

@GAFF.Add_Judge_Rule("n")
def temp(i, Assign):
    tofind = False
    if Assign.Atom_Judge(i, "N3"):
        for bonded_atom in Assign.bonds[i].keys():
            if Assign.Atom_Judge(bonded_atom, "C3"):
                for bonded_atom_bonded in Assign.bonds[bonded_atom].keys():
                    if Assign.atoms[bonded_atom_bonded] in Assign.XA:
                        tofind = True
                        break
            if tofind:
                break
    return tofind

@GAFF.Add_Judge_Rule("nk")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "N4") and "RG3" in Assign.atom_marker[i].keys()

@GAFF.Add_Judge_Rule("nl")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "N4") and "RG4" in Assign.atom_marker[i].keys()

@GAFF.Add_Judge_Rule("n4")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "N4")

@GAFF.Add_Judge_Rule("no")
def temp(i, Assign):
    tofind = 0
    if Assign.Atom_Judge(i, "N3"):
         for bonded_atom in Assign.bonds[i].keys():
            if Assign.Atom_Judge(bonded_atom, "O1"):
                tofind += 1
    return tofind == 2

@GAFF.Add_Judge_Rule("na")
def temp(i, Assign):
    tofind = 0
    if Assign.Atom_Judge(i, "N3"):
         for bonded_atom in Assign.bonds[i].keys():
            if Assign.Atom_Judge(bonded_atom, "O1"):
                tofind += 1
    return Assign.Atom_Judge(i, "N3") and ("AR1" in Assign.atom_marker[i].keys() or "AR2" in Assign.atom_marker[i].keys() or "AR3" in Assign.atom_marker[i].keys())

@GAFF.Add_Judge_Rule("nm")
def temp(i, Assign):
    pass

@GAFF.Add_Judge_Rule("nn")
def temp(i, Assign):
    pass

@GAFF.Add_Judge_Rule("nh")
def temp(i, Assign):
    pass

@GAFF.Add_Judge_Rule("np")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "N3") and "RG3" in Assign.atom_marker[i].keys()

@GAFF.Add_Judge_Rule("nq")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "N3") and "RG4" in Assign.atom_marker[i].keys()

@GAFF.Add_Judge_Rule("n3")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "N3")

@GAFF.Add_Judge_Rule("nb")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "N2") and "AR1" in Assign.atom_marker[i].keys()

@GAFF.Add_Judge_Rule("nc")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "N2") and "sb" in Assign.atom_marker[i] and "db" in Assign.atom_marker[i] and ("AR2" in Assign.atom_marker[i] or "AR3" in Assign.atom_marker[i])

@GAFF.Add_Judge_Rule("ne")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "N2") and "sb" in Assign.atom_marker[i] and "db" in Assign.atom_marker[i]

@GAFF.Add_Judge_Rule("n1")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "N1") or (Assign.Atom_Judge(i, "N2") and (("sb" in Assign.atom_marker[i] and "tb" in Assign.atom_marker[i]) or (Assign.atom_marker[i].get("db") == 2) ))

@GAFF.Add_Judge_Rule("n2")
def temp(i, Assign):
    return Assign.Atom_Judge(i, "N2")

@GAFF.Add_Judge_Rule("s")
def temp(i, Assign):
    pass

@GAFF.Add_Judge_Rule("pc")
def temp(i, Assign):
    pass


print("""Reference for gaff.py:
  Wang, J., Wolf, R.M., Caldwell, J.W., Kollman, P.A. and Case, D.A.
    Development and testing of a general amber force field.
    Journal of Computational Chemistry 2004 25, 1157-1174
    DOI: 10.1002/jcc.20035
""")
