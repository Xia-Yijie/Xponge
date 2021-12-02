from . import *

def Box(molecule, solutions, distance):
    try:
        distance[0]
    except:
        distance = [distance, distance, distance]
    
    if type(molecule) == ResidueType:
        new_molecule = Molecule(molecule.name)
        resA = Residue(molecule)
        for atom in molecule.atoms:
            resA.Add_Atom(atom)
        new_molecule.Add_Residue(resA)
        for key, value in sys.modules['__main__'].__dict__.items():
            if value == molecule:
                sys.modules['__main__'].__dict__[key] = new_molecule 
        molecule = new_molecule
    molcrd = IMPOSE._get_crd(molecule)
    molmin = np.min(molcrd, axis = 0)
    molmax = np.max(molcrd, axis = 0)
    if type(solutions) == ResidueType:
        new_molecule = Molecule(solutions.name)
        resA = Residue(solutions)
        for atom in solutions.atoms:
            resA.Add_Atom(atom)
        new_molecule.Add_Residue(resA)
    else:
        new_molecule = solutions.deepcopy()
    solcrd = IMPOSE._get_crd(new_molecule)
    solmin = np.min(solcrd, axis = 0)
    minindex = np.argmin(solcrd, axis = 0)
    solmax = np.max(solcrd, axis = 0)
    maxindex = np.argmin(solcrd, axis = 0)
    solshape = solmax - solmin + 3
    
    x0 = molmin[0] - solshape[0] - distance[0]
    while x0 < molmax[0] + distance[0] + solshape[0]:
        y0 = molmin[1] - solshape[1] - distance[1]
        while y0 < molmax[1] + distance[1] + solshape[1]:
            z0 = molmin[2] - solshape[2] - distance[2]
            while z0 < molmax[2] + distance[2] + solshape[2]:
                if (x0 > molmin[0] - 3 - solshape[0] and x0 < molmax[0] + 3 + solshape[0] and 
                   y0 > molmin[1] - 3 - solshape[1] and y0 < molmax[1] + 3 + solshape[1] and
                   z0 > molmin[2] - 3 - solshape[2] and z0 < molmax[2] + 3 + solshape[2]):
                   z0 += solshape[2]
                   continue
                for atom in new_molecule.atoms:
                    i = new_molecule.atom_index[atom]
                    atom.x = solcrd[i][0] + x0
                    atom.y = solcrd[i][1] + y0
                    atom.z = solcrd[i][2] + z0
                molecule |= new_molecule
                z0 += solshape[2]
            y0 += solshape[1]
        x0 += solshape[0]

sys.modules['__main__'].__dict__["Process_Box"] = Box

def Pack(molecules, box):
    pass
