try:
    import pyscf
except:
    raise Exception("To calculate RESP charge, 'pyscf' package needed. Maybe you need 'pip install pyscf'")


import numpy as np
from . import *
from .. import *

def Get_Fibonacci_Grid(N, origin, radius):
    n = np.arange(1,N+1)
    factorn = (np.sqrt(5)-1) * np.pi * n
    out = np.zeros((N,3))
    out[:,2] = (2*n - 1) / N - 1
    sqrtz = np.sqrt(1 - out[:,2] * out[:,2])
    out[:,0] = sqrtz * np.cos(factorn)
    out[:,1] = sqrtz * np.sin(factorn)
    out *= radius
    out += origin
    return out

#Pay Attention To !!!UNIT!!!
default_radius = {"H": 1.2, "C":1.5, "N":1.5, 
                  "O": 1.4, "P":1.8, "S":1.75}

#Pay Attention To !!!UNIT!!!
def Get_MK_Grid(molecule, crd, area_density = 1.0, radius = default_radius, ):
    grids = []
    for atom in molecule.atoms:
        if getattr(atom, "element", -1) == -1:
            if getattr(atom, "mass", -1) == -1:
                raise Exception("unknown atom element")
            else:
                raise NotImplementedError("From Mass To Element")
        elif atom.element not in radius.keys():
            raise Exception("unknown atom vdw radius")
    
    factor = area_density * 4 * np.pi
    lists0 = np.array([1.4, 1.6, 1.8, 2.0])
    for i, atom in enumerate(molecule.atoms):
        r0 = radius[atom.element] / 0.5291772094723385
        for r in r0 * lists0:
            grids.extend([*Get_Fibonacci_Grid(int(factor * r * r), crd[i], r)])   
    grids = np.array(grids).reshape(-1, 3)
    for i, atom in enumerate(molecule.atoms):
        r0 = 1.4 * radius[atom.element] / 0.5291772094723385
        t = np.linalg.norm(grids - crd[i], axis = 1)
        grids = grids[t > r0, :]
    return grids



#Pay Attention To !!!UNIT!!!
def RESP_Fit(molecule, basis = "6-31g*", opt = True, charge = 0, spin = 0, method = None):
    from pyscf import gto, scf
    mols = ""
    for atom in molecule.atoms:
        mols += "%s %f %f %f\n"%(atom.element, atom.x, atom.y, atom.z)
    mol = gto.M(atom = mols, verbose = 0, basis = basis, charge = charge, spin = spin, symmetry = True)
    if spin == 0:
        fun = scf.RHF(mol)
    else:
        fun = scf.UHF(mol)
    
    if opt:
        from pyscf.geomopt.geometric_solver import optimize as geometric_opt
        mol = geometric_opt(fun)
        if spin == 0:
            fun = scf.RHF(mol)
        else:
            fun = scf.UHF(mol)
    
    fun.run()
    
    grids = Get_MK_Grid(molecule, mol.atom_coords(), 6)

    #step1
    a = 0.0005
    b = 0.1
    step = 0
    q_last_step = np.ones(mol.natm)
    q_out = np.zeros(mol.natm)
    while np.sum(np.abs(q_out - q_last_step)) > 1e-3:
        step += 1
        if step == 1:
            q_last_step = np.zeros(mol.natm)
        else:
            q_last_step = q_out
        Vnuc = 0
        A = np.zeros((mol.natm, mol.natm))
        for i in range(mol.natm):
            r = mol.atom_coord(i)
            Z = mol.atom_charge(i)
            rp = r - grids
            for j in range(mol.natm):
                rpj = mol.atom_coord(j) - grids
                A[i][j] = np.sum(1.0 / np.linalg.norm(rp, axis = 1) / np.linalg.norm(rpj, axis = 1))
            if molecule.atoms[i].element != "H":
                A[i][i] += a * (np.sqrt(q_last_step[i] ** 2 + b * b) - b )
            Vnuc += Z / np.einsum('xi,xi->x', rp, rp)**.5
        
        A = np.hstack((A,np.ones(mol.natm).reshape(-1,1)))
        temp = np.ones(mol.natm+1)
        temp[-1] = 0
        A = np.vstack((A,temp.reshape(1,-1)))
        Ainv = np.linalg.inv(A)
        
        from pyscf import df
        fakemol = gto.fakemol_for_charges(grids)
        Vele = np.einsum('ijp,ij->p', df.incore.aux_e2(mol, fakemol), fun.make_rdm1())
        
        MEP = Vnuc - Vele
        
        B = np.zeros((mol.natm + 1))
        B[-1] = charge
        for i in range(mol.natm):
            r = mol.atom_coord(i)
            rp = np.linalg.norm(r - grids, axis = 1)
            B[i] = np.sum(MEP / rp)

        q = np.dot(Ainv, B.reshape(-1,1))
        
        q_out = q[:-1]

    #step2
    a = 0.001
    b = 0.1
    step = 0
    q_last_step = np.ones(mol.natm)
    while np.sum(np.abs(q_out - q_last_step)) > 1e-3:
        step += 1
        q_last_step = q_out

        Vnuc = 0
        A = np.zeros((mol.natm, mol.natm))
        for i in range(mol.natm):
            r = mol.atom_coord(i)
            Z = mol.atom_charge(i)
            rp = r - grids
            for j in range(mol.natm):
                rpj = mol.atom_coord(j) - grids
                A[i][j] = np.sum(1.0 / np.linalg.norm(rp, axis = 1) / np.linalg.norm(rpj, axis = 1))
            if molecule.atoms[i].element != "H":
                A[i][i] += a * (np.sqrt(q_last_step[i] ** 2 + b * b) -b)
            Vnuc += Z / np.einsum('xi,xi->x', rp, rp)**.5
        
        A = np.hstack((A,np.ones(mol.natm).reshape(-1,1)))
        temp = np.ones(mol.natm+1)
        temp[-1] = 0
        A = np.vstack((A,temp.reshape(1,-1)))
        Ainv = np.linalg.inv(A)
        
        from pyscf import df
        fakemol = gto.fakemol_for_charges(grids)
        Vele = np.einsum('ijp,ij->p', df.incore.aux_e2(mol, fakemol), fun.make_rdm1())
        
        MEP = Vnuc - Vele
        
        B = np.zeros((mol.natm + 1))
        B[-1] = charge
        for i in range(mol.natm):
            r = mol.atom_coord(i)
            rp = np.linalg.norm(r - grids, axis = 1)
            B[i] = np.sum(MEP / rp)

        q = np.dot(Ainv, B.reshape(-1,1))
        
        q_out = q[:-1]
        
    return q_out

t = Get_Molecule_From_PubChem(175,"cid")
q = RESP_Fit(t, "6-31+g*", opt = False, charge = -1)
print(q)
