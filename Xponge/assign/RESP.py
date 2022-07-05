try:
    import pyscf
except ModuleNotFoundError:
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
                  "O": 1.4, "P":1.8, "S":1.75,
                  "F": 1.35,"Cl":1.7,"Br":2.3}

#Pay Attention To !!!UNIT!!!
def Get_MK_Grid(Assign, crd, area_density = 1.0, layer = 4, radius = None):
    grids = []
    factor = area_density * 0.52918 * 0.52918 * 4 * np.pi
    real_radius = {}
    real_radius.update(default_radius)
    if radius:
        real_radius.update(radius)
    
    lists0 = np.array([1.4 + 0.2 * i for i in range(layer)])
    for i, atom in enumerate(Assign.atoms):
        r0 = real_radius[atom] / 0.52918
        for r in r0 * lists0:
            grids.extend([*Get_Fibonacci_Grid(int(factor * r * r), crd[i], r)])   
    grids = np.array(grids).reshape(-1, 3)
    for i, atom in enumerate(Assign.atoms):
        r0 = 1.39 * real_radius[atom] / 0.52918
        t = np.linalg.norm(grids - crd[i], axis = 1)
        grids = grids[t >= r0, :]
    return grids

def force_equivalence_q(q, extra_equivalence):
    for eq_group in extra_equivalence:
        q_mean = np.mean(q[eq_group])
        q[eq_group] = q_mean
    return q

def _get_pyscf_mol(Assign, basis, charge, spin, opt):
    from pyscf import gto, scf
    mols = ""
    for i, atom in enumerate(Assign.atoms):
        mols += "%s %f %f %f\n"%(atom, Assign.coordinate[i][0], Assign.coordinate[i][1], Assign.coordinate[i][2])
    mol = gto.M(atom = mols, verbose = 0, basis = basis, charge = charge, spin = spin)

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
    return mol, fun

def _resp_scf_kernel(mol, Assign, a, b, A, A0, B, q):
    step = 0
    while step == 0 or np.max(np.abs(q - q_last_step)) > 1e-4:
        step += 1
        q_last_step = q
        for i in range(mol.natm):
            if Assign.atoms[i] != "H":
                A[i][i] = A0[i][i] + a / np.sqrt(q_last_step[i] *q_last_step[i] + b * b)

        Ainv = np.linalg.inv(A)
        q = np.dot(Ainv, B)
        q = q[:-1]
        
    return q

def _find_tofit_second(mol, Assign):
    tofit_second = []
    fit_group = {i : -1 for i in range(mol.natm)}
    sublength = 0
    for i in range(mol.natm):
        if Assign.Atom_Judge(i, "C4"):
            fit_group[i] = len(tofit_second)
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
                tofit_second.append([i])
                for j in temp:
                    fit_group[j] = len(tofit_second)
                tofit_second.append(temp)
                sublength += 1
    return tofit_second, fit_group, sublength

def _correct_extra_equivalence(tofit_second, fit_group, sublength, extra_equivalence):
    if extra_equivalence:
        equi_group = [set() for i in extra_equivalence]    
        for i, eq in enumerate(extra_equivalence):
            for eq_atom in eq:
                if fit_group[eq_atom] != -1:
                    equi_group[i].add(fit_group[eq_atom])
            equi_group[i] = list(equi_group[i])
            equi_group[i].sort()
            
        all_groups = set()
        new_all_groups = set()
        for atom in range(len(Assign.atoms)):
            all_groups.add(fit_group[atom])
        all_groups = list(all_groups)
        all_groups.sort()
        
        
        group_map = {i:i for i in all_groups}
        for eq in equi_group:
            for group in eq:
                group_map[group] = eq[0]
        
        temp_max = 0
        for group in all_groups:
            if group == -1:
                continue
            elif group_map[group] == group:
                group_map[group] = temp_max
                temp_max += 1
            else:
                group_map[group] = group_map[group_map[group]]

        
        temp = tofit_second
        tofit_second = [[] for i in range(temp_max)]
        for i, group in enumerate(temp):
            tofit_second[group_map[i]].extend(group)
            sublength -= len(group) - 1
        
        for group in tofit_second:
            sublength += len(group) - 1
        
        for atom in range(len(Assign.atoms)):
            fit_group[atom] = group_map[fit_group[atom]]
            
    return tofit_second, fit_group, sublength

def _get_A20_and_B20(total_length, tofit_second, fit_group, sublength, mol, A0, B, charge, q):
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
    return A20, B20
        

#Pay Attention To !!!UNIT!!!
def RESP_Fit(Assign, basis = "6-31g*", opt = False, opt_params = None, charge = None, spin = 0, extra_equivalence = None, grid_density = 6, grid_cell_layer = 4, 
    radius = None, a1 = 0.0005, a2 = 0.001, two_stage = True, only_ESP  = False):
    if extra_equivalence is None:
        extra_equivalence = []
    if charge is None:
        charge = int(round(np.sum(Assign.charge)))

    mol, fun = _get_pyscf_mol(Assign, basis, charge, spin, opt)
    grids = Get_MK_Grid(Assign, mol.atom_coords(), grid_density, grid_cell_layer, radius)
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

    try:
        from pyscf import df
        fakemol = gto.fakemol_for_charges(grids)
        Vele = np.einsum('ijp,ij->p', df.incore.aux_e2(mol, fakemol), fun.make_rdm1())
    except:
        dm = fun.make_rdm1()
        Vele = []
        for p in grids:
            mol.set_rinv_orig_(p)
            Vele.append(np.einsum('ij,ij', mol.intor('int1e_rinv'),dm))
        Vele = np.array(Vele)

    MEP = Vnuc - Vele
    
    B = np.zeros((mol.natm + 1))
    for i in range(mol.natm):
        r = mol.atom_coord(i)
        rp = np.linalg.norm(r - grids, axis = 1)
        B[i] = np.sum(MEP / rp) 
    B[-1] = charge
    B = B.reshape(-1,1)
    
    q = np.dot(Ainv, B.reshape(-1,1))[:-1]
    
    if only_ESP:
        return force_equivalence_q(q, extra_equivalence)

    q = _resp_scf_kernel(mol, Assign, a1, 0.1, A, A0, B, q)
        
    if not two_stage:
        return force_equivalence_q(q, extra_equivalence)
    
    #step2
    #fit the sp3 C and the hydrogen connected to it (pay attention to the symmetry!)
    tofit_second, fit_group, sublength = _find_tofit_second(mol, Assign)
    tofit_second, fit_group, sublength = _correct_extra_equivalence(tofit_second, fit_group, sublength, extra_equivalence)

    if tofit_second:
        total_length = mol.natm - sublength + 1 + mol.natm - sublength - len(tofit_second)

        A20, B20 = _get_A20_and_B20(total_length, tofit_second, fit_group, sublength, mol, A0, B, charge, q)

        A = np.zeros_like(A20)
        A[:] = A20[:]
        B = B20.reshape(-1,1)
        Ainv = np.linalg.inv(A)
        q_temp = np.dot(Ainv, B)[:-1]
   
        a = a2
        b = 0.1
        step = 0
        while step == 0 or np.max(np.abs(q_temp - q_last_step)) > 1e-4:
            step += 1
            q_last_step = q_temp
            for i in range(mol.natm - sublength):
                if Assign.atoms[i] != "H":
                    A[i][i] = A20[i][i] + a / np.sqrt(q_last_step[i] * q_last_step[i] + b * b)
            Ainv = np.linalg.inv(A)
            q_temp = np.dot(Ainv, B)[:-1]

        for i, group in enumerate(tofit_second):
            for j in group:
                q[j] = q_temp[i]
            
    return force_equivalence_q(q, extra_equivalence)

print("""Reference for RESP.py:
1. pyscf
  Q. Sun, T. C. Berkelbach, N. S. Blunt, G. H. Booth, S. Guo, Z. Li, J. Liu, J. McClain, S. Sharma, S. Wouters, and G. K.-L. Chan
    PySCF: the Python-based simulations of chemistry framework
    WIREs Computational Molecular Science 2018 8(e1340) 
    DOI: 10.1002/wcms.1340
    
2. ESP MK grid generation
  Brent H. Besler, Kenneth M. Merz Jr., Peter A. Kollman
    Atomic charges derived from semiempirical methods
    Journal of Computational Chemistry 1990 11 431-439
    DOI: 10.1002/jcc.540110404
    
3. RESP
  Christopher I. Bayly, Piotr Cieplak, Wendy Cornell, and Peter A. Kollman
    A well-behaved electrostatic potential-based method using charge restraints for deriving atomic char
    Journal of Physical Chemistry 1993 97(40) 10269-10280
    DOI: 10.1021/j100142a004
    
""")
