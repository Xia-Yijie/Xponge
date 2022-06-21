from . import *
import sys


def Box(molecule, solvent, distance, tolerance = 3):
    if isinstance(distance, float) or isinstance(distance, int):
        distance = [distance] * 6
    elif not isinstance(distance, list):
        raise Exception("parameter distance should be a list, an int or a float")
        

    if len(distance) == 3:
        distance = distance + distance
    elif len(distance) != 6:
        raise Exception("the length of parameter distance should be 3 or 6")

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
    molmin = np.min(molcrd, axis=0)
    molmax = np.max(molcrd, axis=0)
    if type(solvent) == ResidueType:
        new_molecule = Molecule(solvent.name)
        resA = Residue(solvent)
        for atom in solvent.atoms:
            resA.Add_Atom(atom)
        new_molecule.Add_Residue(resA)
    else:
        new_molecule = solvent.deepcopy()
    solcrd = IMPOSE._get_crd(new_molecule)
    solmin = np.min(solcrd, axis=0)
    minindex = np.argmin(solcrd, axis=0)
    solmax = np.max(solcrd, axis=0)
    maxindex = np.argmin(solcrd, axis=0)
    solshape = solmax - solmin + tolerance

    x0 = molmin[0] - solshape[0] - distance[0]
    while x0 < molmax[0] + distance[3] + solshape[0]:
        y0 = molmin[1] - solshape[1] - distance[1]
        while y0 < molmax[1] + distance[4] + solshape[1]:
            z0 = molmin[2] - solshape[2] - distance[2]
            while z0 < molmax[2] + distance[5] + solshape[2]:
                if (x0 > molmin[0] - tolerance - solshape[0] and x0 < molmax[0] + tolerance + solshape[0] and
                        y0 > molmin[1] - tolerance - solshape[1] and y0 < molmax[1] + tolerance + solshape[1] and
                        z0 > molmin[2] - tolerance - solshape[2] and z0 < molmax[2] + tolerance + solshape[2]):
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


def HMass_Repartition(molecules, repartition_mass=1.1, repartition_rate=3, exclude_residue_name="WAT"):
    for res in molecules.residues:
        if res.name == exclude_residue_name:
            continue
        for atom in res.atoms:
            if atom.mass <= repartition_mass:
                connect_atoms = res.type.connectivity[res.type._name2atom[atom.name]]
                assert len(connect_atoms) == 1
                origin_mass = atom.mass
                atom.mass *= repartition_rate
                delta_mass = atom.mass - origin_mass
                for heavy_atom in connect_atoms:
                    res._name2atom[heavy_atom.name].mass -= delta_mass


sys.modules['__main__'].__dict__["HMass_Repartition"] = HMass_Repartition


def Replace(molecule, select, toreplace):
    solutions = []
    for i in range(len(molecule.residues)):
        if select(molecule.residues[i]):
            solutions.append(i)

    np.random.shuffle(solutions)
    count = 0
    for key, value in toreplace.items():
        assert type(key) == ResidueType or (type(key) == Molecule and len(key.residues) == 1)
        if type(key) == Molecule:
            key = key.residues[0].type

        tempi = solutions[count:count + value]
        count += value
        for i in tempi:
            new_residue = Residue(key)
            crd_o = [molecule.residues[i].atoms[0].x, molecule.residues[i].atoms[0].y, molecule.residues[i].atoms[0].z]
            crd0 = [key.atoms[0].x, key.atoms[0].y, key.atoms[0].z]
            for atom in key.atoms:
                new_residue.Add_Atom(atom, x=atom.x + crd_o[0] - crd0[0],
                                     y=atom.y + crd_o[1] - crd0[1], z=atom.z + crd_o[2] - crd0[2])
            molecule.residues[i] = new_residue


sys.modules['__main__'].__dict__["Ion_Replace"] = Replace


def Rotate(molecule, direction_long = None, direction_middle  = None, direction_short = None):
    if direction_long is None:
        direction_long = [0,0,1]
    if direction_middle is None:
        direction_middle = [0, 1, 0]     
    if direction_short is None:
        direction_short = [1, 0, 0]
    molcrd = IMPOSE._get_crd(molecule)
    I = np.zeros((3, 3))
    mass_of_center = np.zeros(3)
    total_mass = 0
    for i, atom in enumerate(molecule.atoms):
        xi, yi, zi = molcrd[i]
        total_mass += atom.mass
        mass_of_center += atom.mass * np.array([xi, yi, zi])
    mass_of_center /= total_mass

    for i, atom in enumerate(molecule.atoms):
        xi, yi, zi = molcrd[i] - mass_of_center
        I += atom.mass * np.array([[yi * yi + zi * zi, -xi * yi, -xi * zi],
                                   [-xi * yi, xi * xi + zi * zi, -yi * zi],
                                   [-xi * zi, -yi * zi, xi * xi + yi * yi]])

    eigval, eigvec = np.linalg.eig(I)
    t = np.argsort(eigval)
    matrix0 = np.vstack([direction_short, direction_middle, direction_long])
    rotation_matrix = np.dot(matrix0, np.linalg.inv(np.vstack((eigvec[:, t[2]], eigvec[:, t[1]], eigvec[:, t[0]]))))
    molcrd = np.dot(molcrd - mass_of_center, rotation_matrix) + mass_of_center
    for i, atom in enumerate(molecule.atoms):
        atom.x = molcrd[i][0]
        atom.y = molcrd[i][1]
        atom.z = molcrd[i][2]


sys.modules['__main__'].__dict__["Molecule_Rotate"] = Rotate


def Mutation(molecule, residue, tomutate):
    pass


def Repeat(molecule, repeat_unit, repeat_place, repeat_range):
    pass


def Pack(molecules, box, conditions):
    pass
