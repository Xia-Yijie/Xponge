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
def Get_MK_Grid(Assign, crd, area_density = 1.0, layer = 4, radius = default_radius):
    grids = []
    factor = area_density * 0.52918 * 0.52918 * 4 * np.pi
    lists0 = np.array([1.4 + 0.2 * i for i in range(layer)])
    for i, atom in enumerate(Assign.atoms):
        r0 = radius[atom] / 0.52918
        for r in r0 * lists0:
            grids.extend([*Get_Fibonacci_Grid(int(factor * r * r), crd[i], r)])   
    grids = np.array(grids).reshape(-1, 3)
    for i, atom in enumerate(Assign.atoms):
        r0 = 1.39 * radius[atom] / 0.52918
        t = np.linalg.norm(grids - crd[i], axis = 1)
        grids = grids[t >= r0, :]
    return grids



#Pay Attention To !!!UNIT!!!
def RESP_Fit(Assign, basis = "6-31g*", opt = False, charge = 0, spin = 0, method = None):
    from pyscf import gto, scf
    mols = ""
    for i, atom in enumerate(Assign.atoms):
        mols += "%s %f %f %f\n"%(atom, Assign.coordinate[i][0], Assign.coordinate[i][1], Assign.coordinate[i][2])
    mol = gto.M(atom = mols, verbose = 0, basis = basis, charge = charge, spin = spin, symmetry = True)
    from pyscf import symm

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
        Assign.coordinate = mol.atom_coords() * 0.52918
    
    fun.run()
    grids = Get_MK_Grid(Assign, mol.atom_coords(), 1)

    #step1
    #fit all atoms    
    Vnuc = 0
    A0 = np.zeros((mol.natm, mol.natm))
    for i in range(mol.natm):
        r = mol.atom_coord(i)
        Z = mol.atom_charge(i)
        rp = r - grids
        for j in range(mol.natm):
            rpj = mol.atom_coord(j) - grids
            A0[i][j] = np.sum(1.0 / np.linalg.norm(rp, axis = 1) / np.linalg.norm(rpj, axis = 1))

        Vnuc += Z / np.einsum('xi,xi->x', rp, rp)**.5
    
    A0 = np.hstack((A0,np.ones(mol.natm).reshape(-1,1)))
    temp = np.ones(mol.natm+1)
    temp[-1] = 0
    A0 = np.vstack((A0,temp.reshape(1,-1)))
    A = np.zeros_like(A0)
    A[:] = A0
    Ainv = np.linalg.inv(A)

    
    from pyscf import df
    fakemol = gto.fakemol_for_charges(grids)
    Vele = np.einsum('ijp,ij->p', df.incore.aux_e2(mol, fakemol), fun.make_rdm1())
    
    MEP = Vnuc - Vele
    
    B = np.zeros((mol.natm + 1))
    for i in range(mol.natm):
        r = mol.atom_coord(i)
        rp = np.linalg.norm(r - grids, axis = 1)
        B[i] = np.sum(MEP / rp) 
    B[-1] = charge
    B = B.reshape(-1,1)
    
    q = np.dot(Ainv, B.reshape(-1,1))[:-1]
    step = 0
    a = 0.0005
    b = 0.1
    while step == 0 or np.max(np.abs(q - q_last_step)) > 1e-4:
        step += 1
        q_last_step = q
        for i in range(mol.natm):
            if Assign.atoms[i] != "H":
                A[i][i] = A0[i][i] + a / np.sqrt(q_last_step[i] *q_last_step[i] + b * b)
        Ainv = np.linalg.inv(A)
        q = np.dot(Ainv, B)
        q = q[:-1]
    
    #step2
    #fit the sp3 C and the hydrogen connected to it (pay attension to the symmetry!)
    tofit_second = []
    notH_group = []
    fit_group = {i : -1 for i in range(mol.natm)}
    sublength = 0
    for i in range(mol.natm):
        if Assign.Atom_Judge(i, "C4"):
            fit_group[i] = len(tofit_second)
            notH_group.append(len(tofit_second))
            tofit_second.append([i])
            temp = []
            for j in Assign.bonds[i].keys():
                if Assign.atoms[j] == "H":
                    temp.append(j)
            if temp:
                for j in temp:
                    fit_group[j] = len(tofit_second)
                tofit_second.append(temp)
                sublength += len(temp) - 1
                    
        if Assign.Atom_Judge(i, "C3"):
            temp = []
            for j in Assign.bonds[i].keys():
                if Assign.atoms[j] == "H":
                    temp.append(j)
            if len(temp) == 2:
                fit_group[i] = len(tofit_second)
                notH_group.append(len(tofit_second))
                tofit_second.append([i])
                for j in temp:
                    fit_group[j] = len(tofit_second)
                tofit_second.append(temp)
                sublength += 1

    if tofit_second:
        total_length = mol.natm - sublength + 1 + mol.natm - sublength - len(tofit_second)
        A20 = np.zeros((total_length, total_length))
        count = len(tofit_second)
        for i in range(mol.natm):
            if fit_group[i] == -1:
                fit_group[i] = count
                count += 1
            A20[mol.natm - sublength][fit_group[i]] += 1
            A20[fit_group[i]][mol.natm - sublength] += 1
            
        B20 = np.zeros(total_length)
        for i in range(mol.natm):
            B20[fit_group[i]] += B[i]
            for j in range(mol.natm):
                A20[fit_group[i]][fit_group[j]] += A0[i][j]
        
        
        B20[mol.natm - sublength] = charge
        count = 0
        for i in range(mol.natm):
            if fit_group[i] >= len(tofit_second):
                B20[mol.natm - sublength + count + 1] = q[i]
                A20[mol.natm - sublength + count + 1][len(tofit_second) + count] = 1
                A20[len(tofit_second) + count][mol.natm - sublength + count + 1] = 1
                count += 1
        
        A = np.zeros_like(A20)
        A[:] = A20[:]
        B = B20.reshape(-1,1)
        Ainv = np.linalg.inv(A)
        q_temp = np.dot(Ainv, B)[:-1]
   
        a = 0.001
        b = 0.1
        step = 0
        while step == 0 or np.max(np.abs(q_temp - q_last_step)) > 1e-4:
            step += 1
            q_last_step = q_temp
            for i in notH_group:
                A[i][i] = A20[i][i] + a / np.sqrt(q_last_step[i] * q_last_step[i] + b * b)
            for i in range(mol.natm - sublength):
                if Assign.atoms[i] != "H" and fit_group[i] >= len(tofit_second):
                    A[i][i] = A20[i][i] + a / np.sqrt(q_last_step[i] * q_last_step[i] + b * b)
            Ainv = np.linalg.inv(A)
            q_temp = np.dot(Ainv, B)[:-1]

        for i, group in enumerate(tofit_second):
            for j in group:
                q[j] = q_temp[i]
            
    return q

sys.modules["__main__"].__dict__["RESP_Fit"] = RESP_Fit