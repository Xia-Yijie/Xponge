from . import *
import os
import sys

##########################################################################
#General Format
##########################################################################
def _mol2_atom(line, current_residue_index, current_residue, ignore_atom_type, temp, current_molecule, atom_residue_map):
    words = line.split()
    if current_residue_index is None or int(words[6]) != current_residue_index:
        current_residue_index = int(words[6])
        if words[7] not in ResidueType.types.keys():
            sys.modules['__main__'].__dict__[words[7]] = ResidueType(name = words[7])
            temp = True
        else:
            temp = False
        if current_residue:
            current_molecule.Add_Residue(current_residue)
        current_residue = Residue(ResidueType.types[words[7]])
    if ignore_atom_type:
        temp_atom_type = AtomType.types["UNKNOWN"]
    else:
        temp_atom_type = AtomType.types[words[5]]
    
    if temp:
        current_residue.type.Add_Atom(words[1], temp_atom_type, *words[2:5])
        current_residue.type.atoms[-1].Update(**{"charge[e]": float(words[8])})
    current_residue.Add_Atom(words[1], temp_atom_type, *words[2:5])
    current_residue.atoms[-1].Update(**{"charge[e]": float(words[8])})
    atom_residue_map[words[0]]=[words[1], current_residue, current_residue_index, temp, current_residue.atoms[-1], current_residue.type.atoms[-1]]
    return current_residue_index, current_residue, temp

def _mol2_bond(line, current_molecule, atom_residue_map):
    words = line.split()
    if atom_residue_map[words[1]][1] == atom_residue_map[words[2]][1]:
        if atom_residue_map[words[1]][3]:
            atom_residue_map[words[1]][1].type.Add_Connectivity(atom_residue_map[words[1]][0], atom_residue_map[words[2]][0])
        atom_residue_map[words[1]][1].Add_Connectivity(atom_residue_map[words[1]][0], atom_residue_map[words[2]][0])
    else:
        current_molecule.Add_Residue_Link(atom_residue_map[words[1]][4], atom_residue_map[words[2]][4])
        index_diff = atom_residue_map[words[1]][2] - atom_residue_map[words[2]][2]
        if  abs(index_diff) == 1:
            if atom_residue_map[words[1]][3]:
                if index_diff < 0:
                    atom_residue_map[words[1]][1].type.tail = atom_residue_map[words[1]][0]
                else:
                    atom_residue_map[words[1]][1].type.head = atom_residue_map[words[1]][0]
            if atom_residue_map[words[2]][3]:
                if index_diff < 0:
                    atom_residue_map[words[2]][1].type.head = atom_residue_map[words[2]][0]
                else:
                    atom_residue_map[words[2]][1].type.tail = atom_residue_map[words[2]][0]
                    
def mol2(filename, ignore_atom_type = False):
    with open(filename) as f:
        #存储读的时候的临时信息，key是编号
        #value是list：原子名(0)、residue(1)、residue编号(2)、是否是新的residue type(3)、该原子(4)、residue type的最新原子(5)
        atom_residue_map = {}  
        flag = None
        nline = 0
        current_molecule = None
        current_residue = None
        current_residue_index = None
        temp = None
        for line in f:
            if line.strip():
                nline += 1
            else:
                continue
            
            if line.startswith("@<TRIPOS>"):
                if flag == "ATOM":
                    current_molecule.Add_Residue(current_residue)
                flag = line[9:].strip()
                nline = 0
            elif flag == "MOLECULE":
                if nline == 1:
                    current_molecule = Molecule(line.strip())
            elif flag == "ATOM":
                current_residue_index, current_residue, temp = _mol2_atom(line, current_residue_index, current_residue, ignore_atom_type, temp, current_molecule, atom_residue_map)
            elif flag == "BOND":
                _mol2_bond(line, current_molecule, atom_residue_map)
    return current_molecule
    
sys.modules['__main__'].__dict__["loadmol2"] = mol2

def _pdb_SSBONDS(chain, residue_type_map, SSBOND, molecule):
    for ssbond in SSBOND:
        resA = chain[ssbond[15]][int(ssbond[17:21])]
        residue_type_map[resA] = "CYX"
        resB = chain[ssbond[29]][int(ssbond[31:35])]
        residue_type_map[resB] = "CYX"
        if resA > resB:
            resA, resB = (resB, resA)
        ResA = molecule.residues[resA]
        ResB = molecule.residues[resB]
        molecule.Add_Residue_Link(ResA._name2atom[ResA.type.connect_atoms["ssbond"]], ResB._name2atom[ResB.type.connect_atoms["ssbond"]])

    
def _pdb_add_residue(filename, molecule, position_need, residue_type_map, ignoreH):
    current_residue_count = -1
    current_residue_index = None
    current_residue = None
    current_resname = None
    links = []
    with open(filename) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                extra = line[16]
                resindex = int(line[22:26])
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                if current_residue_index is None or current_residue_index != resindex or \
                   current_resname != resname:
                    current_residue_count += 1
                    if current_residue:
                        molecule.Add_Residue(current_residue)
                        if current_residue_index is not None and current_residue.type.tail and ResidueType.types[residue_type_map[current_residue_count]].head:
                            links.append(len(molecule.residues))
                    current_residue = Residue(ResidueType.types[residue_type_map[current_residue_count]])
                    current_residue_index = resindex
                    current_resname = resname
                if extra not in (" ", position_need):
                    continue
                if ignoreH and atomname.startswith("H"):
                    continue
                current_residue.Add_Atom(atomname, x = x, y = y, z = z)
            elif line.startswith("TER"):
                current_residue_index = None
                current_resname = None
                
    if current_residue:
        molecule.Add_Residue(current_residue)
    for count in links:
        ResA = molecule.residues[count]
        ResB = molecule.residues[count-1]
        molecule.Add_Residue_Link(ResA._name2atom[ResA.type.head], ResB._name2atom[ResB.type.tail])
        
def _pdb_judge_histone(judge_HIS, residue_type_map, current_HIS):
    if judge_HIS and residue_type_map and residue_type_map[-1] in GlobalSetting.HISMap["HIS"].keys():
        if current_HIS["DeltaH"]: 
            if current_HIS["EpsilonH"]:
                residue_type_map[-1] = GlobalSetting.HISMap["HIS"][residue_type_map[-1]]["HIP"]
            else:
                residue_type_map[-1] = GlobalSetting.HISMap["HIS"][residue_type_map[-1]]["HID"]
        else:
            residue_type_map[-1] = GlobalSetting.HISMap["HIS"][residue_type_map[-1]]["HIE"]
        current_HIS = {"DeltaH":False, "EpsilonH":False}

def pdb(filename, judge_HIS = True, position_need = "A", ignoreH = False):
    molecule = Molecule(os.path.splitext(os.path.basename(filename))[0])
    chain = {}
    SSBOND = []
    residue_type_map = []
    current_residue_count = -1
    current_residue_index = None
    current_resname = None
    current_HIS = {"DeltaH":False, "EpsilonH":False}
    with open(filename) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resindex = int(line[22:26])
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                if current_residue_index == None:
                    _pdb_judge_histone(judge_HIS, residue_type_map, current_HIS)
                    current_residue_count += 1
                    current_resname = resname
                    residue_type_map.append(resname)
                    current_residue_index = resindex
                    chain[chr(ord("A") + len(chain.keys()))] = {resindex:current_residue_count}
                    if resname in GlobalSetting.PDBResidueNameMap["head"].keys():
                        resname = GlobalSetting.PDBResidueNameMap["head"][resname]
                elif current_residue_index != resindex or current_resname != resname:
                    _pdb_judge_histone(judge_HIS, residue_type_map, current_HIS)
                    current_residue_count += 1
                    current_resname = resname
                    current_residue_index = resindex
                    chain[chr(ord("A") + len(chain.keys()) - 1)][resindex] = current_residue_count
                    
                    residue_type_map.append(resname)

                if judge_HIS and resname in GlobalSetting.HISMap["HIS"].keys():
                    if atomname == GlobalSetting.HISMap["DeltaH"]:
                        current_HIS["DeltaH"] = True
                    elif atomname == GlobalSetting.HISMap["EpsilonH"]:
                        current_HIS["EpsilonH"] = True
            elif line.startswith("TER"):
                current_residue_index = None
                current_resname = None
                if residue_type_map[-1] in GlobalSetting.PDBResidueNameMap["tail"].keys():
                   residue_type_map[-1] = GlobalSetting.PDBResidueNameMap["tail"][residue_type_map[-1]]
                _pdb_judge_histone(judge_HIS, residue_type_map, current_HIS)

            elif line.startswith("SSBOND"):
                SSBOND.append(line)
    
    current_residue_index = None
    if residue_type_map[-1] in GlobalSetting.PDBResidueNameMap["tail"].keys():
       residue_type_map[-1] = GlobalSetting.PDBResidueNameMap["tail"][residue_type_map[-1]]
    
    _pdb_add_residue(filename, molecule, position_need, residue_type_map, ignoreH)
    _pdb_SSBONDS(chain, residue_type_map, SSBOND, molecule)
    
    return molecule
    
sys.modules['__main__'].__dict__["loadpdb"] = pdb     

##########################################################################
#AMBER Format
##########################################################################
def frcmod(filename, nbtype = "RE"):
    with open(filename) as f:
        f.readline()
        flag = None
        atom_types = {}  #元素符号和质量
        bonds = "name  k[kcal/mol·A^-2]    b[A]\n"
        angles = "name  k[kcal/mol·rad^-2]    b[degree]\n"
        propers = "name  k[kcal/mol]    phi0[degree]    periodicity    reset\n"
        reset = 1
        impropers = "name  k[kcal/mol]    phi0[degree]    periodicity\n"
        cmap = {}
        cmap_flag = None
        temp_cmp = {}
        if nbtype == "SK":
            raise NotImplementedError
        elif nbtype == "AC":
            LJs = "name A[kcal/mol·A^-12]   B[kcal/mol·A^-6]\n"
        elif nbtype == "RE":
            LJs = "name rmin[A]   epsilon[kcal/mol]\n"
        
        for line in f:
            if not line.strip():
                continue
            words = line.split()
            if flag != "CMAP" and len(words) == 1:
                flag = line.strip()
            elif flag[:4] == "MASS":
                atom_types[words[0]] = words[1]
            elif flag[:4] == "BOND":
                atoms = [ word.strip() for word in line[:5].split("-") ]
                words = line[5:].split()
                bonds += "-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\n"
            elif flag[:4] == "ANGL":
                atoms = [ word.strip() for word in line[:8].split("-") ]
                words = line[8:].split()
                angles += "-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\n"
            elif flag[:4] == "DIHE":
                atoms = [ word.strip() for word in line[:11].split("-") ]
                words = line[11:].split()
                propers += "-".join(atoms) + "\t" + str(float(words[1]) / int(words[0])) + "\t" + words[2] + "\t" + str(abs(int(float(words[3])))) + "\t" + str(reset) + "\n"
                if int(float(words[3])) < 0:
                    reset = 0
                else:
                    reset = 1
            elif flag[:4] == "IMPR":
                atoms = [ word.strip() for word in line[:11].split("-") ]
                words = line[11:].split()
                impropers += "-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\t" + str(int(float(words[2]))) + "\n"
            elif flag[:4] == "NONB":
                words = line.split()
                LJs += words[0] + "-" + words[0] + "\t" + words[1] + "\t" + words[2] + "\n"
            elif flag[:4] == "CMAP":
                if line.startswith("%FLAG"):
                    if "CMAP_COUNT" in line:
                        if temp_cmp:
                            for res in temp_cmp["residues"]:
                                cmap[res] = temp_cmp["info"]
                        temp_cmp = {"residues":[], "info": {"resolution":24, "count":int(line.split()[-1]), "parameters": []}}
                        cmap_flag = "CMAP_COUNT"
                    elif "CMAP_RESOLUTION" in line:
                        temp_cmp["info"]["resolution"] = int(line.split()[-1])
                        cmap_flag = "CMAP_RESOLUTION"
                    elif  "CMAP_RESLIST" in line:
                        cmap_flag = "CMAP_RESLIST"
                    elif "CMAP_TITLE" in line:
                        cmap_flag = "CMAP_TITLE"
                    elif "CMAP_PARAMETER" in line:
                        cmap_flag = "CMAP_PARAMETER"
                elif cmap_flag == "CMAP_RESLIST":
                    temp_cmp["residues"].extend(line.split())
                elif cmap_flag == "CMAP_PARAMETER":
                    temp_cmp["info"]["parameters"].extend([float(x) for x in line.split()])
                    
                    
    if temp_cmp:
        for res in temp_cmp["residues"]:
            cmap[res] = temp_cmp["info"]
    atoms = "name  mass  LJtype\n"
    for atom, mass in atom_types.items():
        atoms += atom + "\t" + mass + "\t" + atom + "\n"
        
    return atoms, bonds, angles, propers, impropers, LJs, cmap

sys.modules['__main__'].__dict__["loadfrcmod"] = frcmod    

def parmdat(filename):
    with open(filename) as f:
        f.readline()
        #读原子
        atom_types = {}  #元素符号和质量
        LJ_types = {}    #元素符号和LJ类型
        for line in f:
            if not line.strip():
                break
            else:
                words = line.split()
                atom_types[words[0]] = words[1]
                LJ_types[words[0]] = words[0]
        f.readline()
        #读键长
        bonds = "name  k[kcal/mol·A^-2]    b[A]\n"
        for line in f:
            if not line.strip():
                break
            else:
                atoms = [ word.strip() for word in line[:5].split("-") ]
                words = line[5:].split()

                bonds += "-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\n"
                
        #读键角
        angles = "name  k[kcal/mol·rad^-2]    b[degree]\n"
        for line in f:
            if not line.strip():
                break
            else:
                atoms = [ word.strip() for word in line[:8].split("-") ]
                words = line[8:].split()
                angles += "-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\n"

        #读恰当二面角
        reset = 1
        propers = "name  k[kcal/mol]    phi0[degree]    periodicity    reset\n"
        for line in f:
            if not line.strip():
                break
            else:
                atoms = [ word.strip() for word in line[:11].split("-") ]
                words = line[11:].split()
                propers += "-".join(atoms) + "\t" + str(float(words[1]) / int(words[0])) + "\t" + words[2] + "\t" + str(abs(int(float(words[3])))) + "\t" + str(reset) + "\n"
                if int(float(words[3])) < 0:
                    reset = 0
                else:
                    reset = 1
        #读非恰当二面角
        impropers = "name  k[kcal/mol]    phi0[degree]    periodicity\n"
        for line in f:
            if not line.strip():
                break
            else:
                atoms = [ word.strip() for word in line[:11].split("-") ]
                words = line[11:].split()
                impropers += "-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\t" + str(int(float(words[2]))) + "\n"
        
        #跳过水的信息
        f.readline()
        f.readline()
        
        #读LJ种类
        for line in f:
            if not line.strip():
                break
            else:
                atoms = line.split()
                atom0 = atoms.pop(0)
                for atom in atoms:
                    LJ_types[atom] = atom0
        
        #读LJ信息
        word = f.readline().split()[1]
        if word == "SK":
            raise NotImplementedError
        elif word == "AC":
            LJs = "name A[kcal/mol·A^-12]   B[kcal/mol·A^-6]\n"
        elif word == "RE":
            LJs = "name rmin[A]   epsilon[kcal/mol]\n"
        
        for line in f:
            if not line.strip():
                break    
            else:
                words = line.split()
                LJs += words[0] + "-" + words[0] + "\t" + words[1] + "\t" + words[2] + "\n"
        
        

    atoms = "name  mass  LJtype\n"
    for atom, mass in atom_types.items():
        atoms += atom + "\t" + mass + "\t" + LJ_types[atom] + "\n"
        
    return atoms, bonds, angles, propers, impropers, LJs
        
sys.modules['__main__'].__dict__["loadparmdat"] = parmdat 

def _read_parm7_someflag(length, word_process, offset, f):
    tempvar = []
    line = f.readline()
    while line.startswith("%"):
        line = f.readline()
    while 1:
        i = 0
        while i + offset < len(line):
            word = line[i:i+offset].strip()
            if word:
                tempvar.append(word_process(word))
            else:
                break
            i += offset
        if len(tempvar) == length:
            break
        line = f.readline()
    return tempvar

def _read_parm7_build(mol, BONDTYPE, bonds, parameter_head, parameter_lists, add_connectivity = False):
    bond_names = {}
    newBONDs = ""
    newBONDs += parameter_head
    for i in bonds.keys():
        bond_names[i] = "-".join([mol.atoms[atype].type.name for atype in bonds[i][0]])
        newBONDs += bond_names[i]
        for fmt, parm_list in parameter_lists:
            newBONDs += " " + fmt%parm_list[i]
        newBONDs += "\n"
    BONDTYPE.New_From_String(newBONDs)
    
    for typei, bondi in bonds.items():
        for atom_indexes in bondi:
            atoms = [mol.atoms[i] for i in atom_indexes ]
            residues = set([atom.residue for atom in atoms])
            if len(residues) == 1:
                residue = atoms[0].residue
                resatoms = [residue.type._name2atom[atom.name] for atom in atoms]
                residue.type.Add_Bonded_Force(BONDTYPE.entity(resatoms, BONDTYPE.types_different_name[bond_names[typei]]))
                residue.Add_Bonded_Force(BONDTYPE.entity(atoms, BONDTYPE.types_different_name[bond_names[typei]]))
                if add_connectivity:
                    residue.Add_Connectivity(atoms[0], atoms[-1])
            else:
                raise NotImplementedError
                
def _read_parm7_read(filename):   
    bonds = {}
    angles = {}
    propers = {}
    nb14s = {}
    impropers = {}
    with open(filename) as f:
        for line in f:
            if line.startswith("%FLAG TITLE"):
                line = f.readline()
                while line.startswith("%"):
                    line = f.readline()
                try:
                    name = line.split()[0].strip()
                except IndexError:
                    from random import randint
                    name = "Molecule%d"%randint(0, 10086)
                mol =  Molecule(name = name)
            elif line.startswith("%FLAG POINTERS"):
                line = f.readline()
                while line.startswith("%"):
                    line = f.readline()
                atom_numbers = int(line[0:8])
                atom_type_numbers = int(line[8:16])
                bond_numbers_H = int(line[16:24])
                bond_numbers_N = int(line[24:32])
                bond_numbers = bond_numbers_H + bond_numbers_N
                angle_numbers_H = int(line[32:40])
                angle_numbers_N = int(line[40:48])
                angle_numbers = angle_numbers_H + angle_numbers_N
                dihedral_numbers_H = int(line[48:56])
                dihedral_numbers_N = int(line[56:64])            
                dihedral_numbers = dihedral_numbers_H + dihedral_numbers_N
                line = f.readline()
                total_excluded_numbers = int(line[0:8])
                residue_numbers = int(line[8:16])
                bond_type_numbers = int(line[40:48])
                angle_type_numbers = int(line[48:56])
                dihedral_type_numbers = int(line[56:64])
            elif line.startswith("%FLAG DIHEDRAL_FORCE_CONSTANT"):
                dihedral_k = _read_parm7_someflag(dihedral_type_numbers, lambda x:float(x), 16, f)
            elif line.startswith("%FLAG DIHEDRAL_PERIODICITY"):
                dihedral_n = _read_parm7_someflag(dihedral_type_numbers, lambda x:float(x), 16, f)
            elif line.startswith("%FLAG DIHEDRAL_PHASE"):
                dihedral_phase = _read_parm7_someflag(dihedral_type_numbers, lambda x:float(x), 16, f)
            elif line.startswith("%FLAG DIHEDRALS_INC_HYDROGEN"):
                tempvar = _read_parm7_someflag(5 * dihedral_numbers_H, lambda x:int(x), 8, f)
                for i in range(0,len(tempvar),5):
                    if tempvar[i+3] < 0:
                        if tempvar[i+4] - 1 in impropers.keys():
                            impropers[tempvar[i+4] - 1].append([tempvar[i] // 3, tempvar[i+1] // 3, abs(tempvar[i+2] // 3), abs(tempvar[i+3]) // 3])
                        else:
                            impropers[tempvar[i+4] - 1] = [[tempvar[i] // 3, tempvar[i+1] // 3, abs(tempvar[i+2]) // 3, abs(tempvar[i+3]) // 3]]
                    else:
                        if tempvar[i+4] - 1 in propers.keys():
                            propers[tempvar[i+4] - 1].append([tempvar[i] // 3, tempvar[i+1] // 3, abs(tempvar[i+2] // 3), tempvar[i+3] // 3])
                        else:
                            propers[tempvar[i+4] - 1] = [[tempvar[i] // 3, tempvar[i+1] // 3, abs(tempvar[i+2] // 3), tempvar[i+3] // 3]]
                        if tempvar[i+2] > 0:
                            if tempvar[i+4] - 1 in nb14s.keys():
                                nb14s[tempvar[i+4] - 1].append([tempvar[i] // 3, tempvar[i+3] // 3])
                            else:
                                nb14s[tempvar[i+4] - 1] = [[tempvar[i] // 3, tempvar[i+3] // 3]]
            elif line.startswith("%FLAG DIHEDRALS_WITHOUT_HYDROGEN"):
                tempvar = _read_parm7_someflag(5 * dihedral_numbers_N, lambda x:int(x), 8, f)
                for i in range(0,len(tempvar),5):
                    if tempvar[i+3] < 0:
                        if tempvar[i+4] - 1 in impropers.keys():
                            impropers[tempvar[i+4] - 1].append([tempvar[i] // 3, tempvar[i+1] // 3, abs(tempvar[i+2] // 3), abs(tempvar[i+3]) // 3])
                        else:
                            impropers[tempvar[i+4] - 1] = [[tempvar[i] // 3, tempvar[i+1] // 3, abs(tempvar[i+2]) // 3, abs(tempvar[i+3]) // 3]]
                    else:
                        if tempvar[i+4] - 1 in propers.keys():
                            propers[tempvar[i+4] - 1].append([tempvar[i] // 3, tempvar[i+1] // 3, abs(tempvar[i+2] // 3), tempvar[i+3] // 3])
                        else:
                            propers[tempvar[i+4] - 1] = [[tempvar[i] // 3, tempvar[i+1] // 3, abs(tempvar[i+2] // 3), tempvar[i+3] // 3]]
                        if tempvar[i+2] > 0:
                            if tempvar[i+4] - 1 in nb14s.keys():
                                nb14s[tempvar[i+4] - 1].append([tempvar[i] // 3, tempvar[i+3] // 3])
                            else:
                                nb14s[tempvar[i+4] - 1] = [[tempvar[i] // 3, tempvar[i+3] // 3]]

            elif line.startswith("%FLAG ANGLE_FORCE_CONSTANT"):
                angle_k = _read_parm7_someflag(angle_type_numbers, lambda x:float(x), 16, f)
            elif line.startswith("%FLAG ANGLE_EQUIL_VALUE"):
                angle_b = _read_parm7_someflag(angle_type_numbers, lambda x:float(x), 16, f)
            elif line.startswith("%FLAG ANGLES_INC_HYDROGEN"):
                tempvar = _read_parm7_someflag(4 * angle_numbers_H, lambda x:int(x), 8, f)
                for i in range(0,len(tempvar),4):
                    if tempvar[i+3] - 1 in angles.keys():
                        angles[tempvar[i+3] - 1].append([tempvar[i] // 3, tempvar[i+1] // 3, tempvar[i+2] // 3])
                    else:
                        angles[tempvar[i+3] - 1] = [[tempvar[i] // 3, tempvar[i+1] // 3, tempvar[i+2] // 3]]
            elif line.startswith("%FLAG ANGLES_WITHOUT_HYDROGEN"):
                tempvar = _read_parm7_someflag(4 * angle_numbers_N, lambda x:int(x), 8, f)
                for i in range(0,len(tempvar),4):
                    if tempvar[i+3] - 1 in angles.keys():
                        angles[tempvar[i+3] - 1].append([tempvar[i] // 3, tempvar[i+1] // 3, tempvar[i+2] // 3])
                    else:
                        angles[tempvar[i+3] - 1] = [[tempvar[i] // 3, tempvar[i+1] // 3, tempvar[i+2] // 3]]
            elif line.startswith("%FLAG BOND_FORCE_CONSTANT"):
                bond_k = _read_parm7_someflag(bond_type_numbers, lambda x:float(x), 16, f)
            elif line.startswith("%FLAG BOND_EQUIL_VALUE"):
                bond_b = _read_parm7_someflag(bond_type_numbers, lambda x:float(x), 16, f)
            elif line.startswith("%FLAG BONDS_INC_HYDROGEN"):
                tempvar = _read_parm7_someflag(3 * bond_numbers_H, lambda x:int(x), 8, f)
                for i in range(0,len(tempvar),3):
                    if tempvar[i+2] - 1 in bonds.keys():
                        bonds[tempvar[i+2] - 1].append([tempvar[i] // 3, tempvar[i+1] // 3])
                    else:
                        bonds[tempvar[i+2] - 1] = [[tempvar[i] // 3, tempvar[i+1] // 3]]
            elif line.startswith("%FLAG BONDS_WITHOUT_HYDROGEN"):
                tempvar = _read_parm7_someflag(3 * bond_numbers_N, lambda x:int(x), 8, f)
                for i in range(0,len(tempvar),3):
                    if tempvar[i+2] - 1 in bonds.keys():
                        bonds[tempvar[i+2] - 1].append([tempvar[i] // 3, tempvar[i+1] // 3])
                    else:
                        bonds[tempvar[i+2] - 1] = [[tempvar[i] // 3, tempvar[i+1] // 3]]
            elif line.startswith("%FLAG NUMBER_EXCLUDED_ATOMS"):
                exclude_numbers = _read_parm7_someflag(atom_numbers, lambda x:int(x), 8, f)
            elif line.startswith("%FLAG EXCLUDED_ATOMS_LIST"):
                excluded_list = _read_parm7_someflag(total_excluded_numbers, lambda x:int(x) - 1, 8, f)
            elif line.startswith("%FLAG ATOM_NAME"):
                atom_names = _read_parm7_someflag(atom_numbers, lambda x:x, 4, f)
            elif line.startswith("%FLAG AMBER_ATOM_TYPE"):
                atom_types = _read_parm7_someflag(atom_numbers, lambda x:x, 4, f)
            elif line.startswith("%FLAG LENNARD_JONES_ACOEF"):
                LJ_A = _read_parm7_someflag(atom_type_numbers * (atom_type_numbers + 1) // 2, lambda x:float(x) , 16, f)
            elif line.startswith("%FLAG LENNARD_JONES_BCOEF"):
                LJ_B = _read_parm7_someflag(atom_type_numbers * (atom_type_numbers + 1) // 2, lambda x:float(x) , 16, f)
            elif line.startswith("%FLAG ATOM_TYPE_INDEX"):
                atom_type_index = _read_parm7_someflag(atom_numbers, lambda x:int(x) - 1 , 8, f)
            elif line.startswith("%FLAG CHARGE"):
                charges = _read_parm7_someflag(atom_numbers, lambda x:float(x)/18.2223, 16, f)
            elif line.startswith("%FLAG MASS"):
                masses = _read_parm7_someflag(atom_numbers, lambda x:float(x), 16, f)
            elif line.startswith("%FLAG RESIDUE_LABEL"):
                residues = _read_parm7_someflag(residue_numbers, lambda x:x, 4, f)
            elif line.startswith("%FLAG RESIDUE_POINTER"):
                residue_starts = _read_parm7_someflag(residue_numbers, lambda x: int(x) - 1, 8, f)
                residue_starts.append(atom_numbers)
    return locals()

def _read_parm7_construct(mol):
    mol.atoms = []
    mol.bonded_forces = {frc.name: [] for frc in GlobalSetting.BondedForces}
    for res in mol.residues:
        mol.atoms.extend(res.atoms)
        for frc in GlobalSetting.BondedForces:
            mol.bonded_forces[frc.name].extend(res.bonded_forces.get(frc.name, []))
        for link in mol.residue_links:
            for frc in GlobalSetting.BondedForces:
                mol.bonded_forces[frc.name].extend(link.bonded_forces.get(frc.name, []))
    mol.atom_index = {mol.atoms[i]: i for i in range(len(mol.atoms))}


def parm7(filename):
    import Xponge.forcefield.AMBER
    local = _read_parm7_read(filename)
    
    atom_types = local["atom_types"]
    atom_type_index = local["atom_type_index"]
    masses = local["masses"]
    residues = local["residues"]
    residue_starts = local["residue_starts"]
    atom_names = local["atom_names"]
    charges = local["charges"]
    mol = local["mol"]
    atom_type_numbers = local["atom_type_numbers"]
    LJ_A = local["LJ_A"]
    LJ_B = local["LJ_B"]
    atom_numbers = local["atom_numbers"]
    exclude_numbers = local["exclude_numbers"]
    excluded_list = local["excluded_list"]
    bonds = local["bonds"]
    angles = local["angles"]
    bond_k = local["bond_k"]
    bond_b = local["bond_b"]
    angle_k = local["angle_k"]
    angle_b = local["angle_b"]
    propers = local["propers"]
    dihedral_k = local["dihedral_k"]
    dihedral_phase = local["dihedral_phase"]
    dihedral_n = local["dihedral_n"]
    impropers = local["impropers"]
    nb14s = local["nb14s"]
    
    atom_type_map = {}
    for i, atype in enumerate(atom_types):
        if atom_type_index[i] not in atom_type_map.keys():
            atom_type_map[atom_type_index[i]] = atype
        if atype not in AtomType.types.keys():
            AtomType(name = atype, mass = masses[i], LJtype = atom_type_map[atom_type_index[i]])
    
    for ri, residue in enumerate(residues):
        if residue not in ResidueType.types.keys():
            New_ResidueType = ResidueType(name = residue)
            New_ResidueType.built = True
            for i in range(residue_starts[ri],residue_starts[ri+1]):
                New_ResidueType.Add_Atom(atom_names[i], AtomType.types[atom_types[i]], x = 0, y = 0, z = 0)
                New_ResidueType.atoms[-1].mass = masses[i]
                New_ResidueType.atoms[-1].charge = charges[i]
            sys.modules['__main__'].__dict__[residue] = New_ResidueType
        New_Residue = Residue(ResidueType.types[residue])
        for i in range(residue_starts[ri],residue_starts[ri+1]):
            New_Residue.Add_Atom(atom_names[i], AtomType.types[atom_types[i]])
            New_Residue.atoms[-1].mass = masses[i]
            New_Residue.atoms[-1].charge = charges[i]
        New_Residue.built = True
        mol.Add_Residue(New_Residue)
    
    _read_parm7_construct(mol)

    for vatom_type_name, vatom_type_atom_numbers in GlobalSetting.VirtualAtomTypes.items():
        for vatom in mol.bonded_forces.get(vatom_type_name, []):
            this_vatoms = [vatom.atoms[0]]
            for i in range(vatom_type_atom_numbers):
                this_vatoms.append(cls.atoms[cls.atom_index[vatom.atoms[0]] + getattr(vatom, "atom%d" % i)])
                this_vatoms.sort(key=lambda x: cls.atom_index[x])
                while this_vatoms:
                    tolink = this_vatoms.pop(0)
                    for i in this_vatoms:
                        tolink.Link_Atom("v", i)
    newLJs = "name A B\n"
    for i in range(atom_type_numbers):
        for j in range(0, i + 1):
            atype = atom_type_map[i]
            btype = atom_type_map[j]
            index = i * (i + 1) // 2 + j
            newLJs += "%s-%s %f %f\n"%(atype, btype, LJ_A[index], LJ_B[index])
    Xponge.forcefield.BASE.LJ.LJType.New_From_String(newLJs)
    mol.built = True

    count = -1
    for i in range(atom_numbers):
        atom_i = mol.atoms[i]
        for j in range(exclude_numbers[i]):
            count += 1
            if excluded_list[count] == -1:
                break
            atom_i.Extra_Exclude_Atom(mol.atoms[excluded_list[count]])
        atom_i.residue.type._name2atom[atom_i.name].Extra_Exclude_Atoms([ atom_i.residue.type._name2atom[aton.name] for aton in atom_i.extra_excluded_atoms])

    _read_parm7_build(mol, Xponge.forcefield.BASE.BOND.BondType, bonds, "name k b\n", [["%f", bond_k], ["%f", bond_b]], True)           
    _read_parm7_build(mol, Xponge.forcefield.BASE.ANGLE.AngleType, angles, "name k b\n", [["%f", angle_k], ["%f", angle_b]])
    _read_parm7_build(mol, Xponge.forcefield.BASE.DIHEDRAL.ProperType, propers, "name k phi0 periodicity reset\n", [["%f", dihedral_k], ["%f", dihedral_phase], ["%d", dihedral_n], ["%d", [0 for i in dihedral_n]]]) 
    _read_parm7_build(mol, Xponge.forcefield.BASE.DIHEDRAL.ImproperType, impropers, "name k phi0 periodicity\n", [["%f", dihedral_k], ["%f", dihedral_phase], ["%d", dihedral_n]]) 
    
    NB14TYPE = Xponge.forcefield.BASE.NB14.NB14Type
    for typei, bondi in nb14s.items():
        for ai, bi in bondi:
            aAtom = mol.atoms[ai]
            bAtom = mol.atoms[bi]
            if aAtom.residue == bAtom.residue:
                atoms = [aAtom, bAtom]
                resatoms = [aAtom.residue.type._name2atom[atom.name] for atom in atoms]
                aAtom.residue.type.Add_Bonded_Force(NB14TYPE.entity(resatoms, NB14TYPE.types["X-X"]))
                aAtom.residue.Add_Bonded_Force(NB14TYPE.entity(atoms, NB14TYPE.types["X-X"]))
            else:
                raise NotImplementedError

    _read_parm7_construct(mol)

    return mol

sys.modules['__main__'].__dict__["loadparm7"] = parm7 

def rst7(filename, mol = None):
    crds = []
    with open(filename) as f:
        f.readline()
        words = f.readline().split()
        atom_numbers = int(words[0])
        for line in f:
            words = line.split()
            while words and len(crds) < atom_numbers * 3:
                crds.append(float(words.pop(0)))
        box = [float(i) for i in line.split()]
    crds = np.array(crds).reshape((-1,3))
    mol.box_length = box[:3]
    count = -1
    for ri, residue in enumerate(mol.residues):
       for ai, atom in enumerate(residue.atoms):
           count += 1  
           atom.x = crds[count][0]
           atom.y = crds[count][1]
           atom.z = crds[count][2]

    return crds, box
     
sys.modules['__main__'].__dict__["loadrst7"] = rst7 

##########################################################################
#GROMACS Format
##########################################################################

class GROMACS_TOPOLOGY_ITERATOR():
    def __init__(self, filename = None, macros = None):
        self.files = []
        self.filenames = []
        if macros:
            self.defined_macros = macros 
        else:
            self.defined_macros = {}
        if filename:
            self.Add_Iterator_File(filename)
    def Add_Iterator_File(self, filename):
        if self.files:
            filename = os.path.abspath(os.path.join(os.path.dirname(self.filenames[-1]), filename.replace('"', '')))
        else:
            filename = os.path.abspath(filename.replace('"', ''))
            
        f = open(filename)
        self.files.append(f)
        self.filenames.append(filename)
    def __iter__(self):
        self.flag = ""
        self.macro_define_stat = []
        return self
    def __next__(self):
        while self.files:
            f = self.files[-1]
            line = f.readline()
            
            if line:
                line = line.strip()
                comment = line.find(";")
                if comment >= 0:
                    line = line[:comment]
                while line and line[-1] == "\\":
                    line = line[:-1] + " " + next(self).strip()        
                if not line:
                    line = next(self)
                if line[0] == "#":         
                    words = line.split()
                    if words[0] == "#ifdef":
                        macro = words[1]
                        if self.macro_define_stat and self.macro_define_stat[-1] == False:
                            self.macro_define_stat.append(False)
                        elif macro in self.defined_macros.keys():
                            self.macro_define_stat.append(True)
                        else:
                            self.macro_define_stat.append(False)
                    elif words[0] == "#else":
                        if len(self.macro_define_stat) <= 1 or self.macro_define_stat[-2] != False:
                            self.macro_define_stat[-1] = not self.macro_define_stat[-1]
                    elif words[0] == "#endif":
                        self.macro_define_stat.pop()    
                    elif self.macro_define_stat and self.macro_define_stat[-1] == False:
                        next(self)
                    elif words[0] == "#define":
                        if len(words) > 2:
                            self.defined_macros[words[1]] = line[line.find(words[1]) + len(words[1]):].strip()
                        else:    
                            self.defined_macros[words[1]] = ""
                    elif words[0] == "#include":
                        self.Add_Iterator_File(words[1])
                    elif words[0] == "#undef":
                        self.defined_macros.pop(words[1])
                    elif words[0] == "#error":
                        raise Exception(line)
                    line = next(self)
                elif self.macro_define_stat and self.macro_define_stat[-1] == False:
                    line = next(self)
                elif "[" in line and "]" in line:
                    self.flag = line[1:-1].strip()
                    line = next(self)
                for macro, tobecome in self.defined_macros.items():
                    line = line.replace(macro, tobecome)
                return line
            else:
                f.close()
                self.files.pop()
                self.filenames.pop()
                continue
        raise StopIteration
            
def ffitp(filename, macros = None):
    iterator = GROMACS_TOPOLOGY_ITERATOR(filename, macros)
    output = {}
    output["nb14"] = "name  kLJ  kee\n"
    output["atomtypes"] = "name mass charge[e] LJtype\n"
    output["bonds"] = "name b[nm] k[kJ/mol·nm^-2]\n"
    output["angles"] = "name b[degree] k[kJ/mol·rad^-2]\n"
    output["Urey-Bradley"] = "name b[degree] k[kJ/mol·rad^-2] r13[nm] kUB[kJ/mol·nm^-2]\n"
    output["dihedrals"] = "name phi0[degree] k[kJ/mol] periodicity  reset\n"
    output["impropers"] = "name phi0[degree] k[kJ/mol·rad^-2]\n"
    output["cmaps"] = {}
    for line in iterator:
        if iterator.flag == "":
            continue
        elif iterator.flag == "defaults":
            words = line.split()
            assert int(words[0]) == 1, "SPONGE Only supports Lennard-Jones now"
            if int(words[1]) == 1:
                output["LJ"] = "name A[kJ/mol·nm^6] B[kJ/mol·nm^12]\n"
                output["nb14_extra"] = "name A[kJ/mol·nm^6] B[kJ/mol·nm^12] kee\n"
            else:
                output["LJ"] = "name sigma[nm] epsilon[kJ/mol] \n"
                output["nb14_extra"] = "name sigma[nm] epsilon[kJ/mol] kee\n"
            fudgeLJ = float(words[3])
            fudgeQQ = float(words[4])
            if words[2] == "yes":
                output["nb14"] += "X-X {fudgeLJ} {fudgeQQ}\n".format(fudgeLJ=fudgeLJ, fudgeQQ=fudgeQQ)
                
        elif iterator.flag == "atomtypes":
            words = line.split()
            output["atomtypes"] += "{type} {mass} {charge} {type}\n".format(type=words[0], mass = float(words[2]), charge = float(words[3]))
            output["LJ"] += "{type}-{type} {V} {W}\n".format(type=words[0], V=float(words[5]), W=float(words[6]))
        elif iterator.flag == "pairtypes":
            words = line.split()
            if len(words) <= 3:
                output["nb14"] +=  "{atom1}-{atom2} {kLJ} {kee}\n".format(atom1 = words[0], atom2 = words[1], kLJ = fudgeLJ, kee = fudgeQQ)
            elif words[2] == "1":
                output["nb14_extra"] += "{atom1}-{atom2} {V} {W} {kee}\n".format(atom1 = words[0], atom2 = words[1], V = float(words[3]), W = float(words[4]), kee = fudgeQQ)
                output["nb14"] += "{atom1}-{atom2} 0 0\n".format(atom1 = words[0], atom2 = words[1])
            elif words[2] == "2":
                raise NotImplementedError
        elif iterator.flag == "bondtypes":
            words = line.split()
            func = words[2]
            if func == "1":
                output["bonds"] += "{atom1}-{atom2} {b} {k}\n".format(atom1 = words[0], atom2 = words[1], b = float(words[3]), k = float(words[4])/2)
            else:
                raise NotImplementedError
        elif iterator.flag == "angletypes":
            words = line.split()
            func = words[3]
            if func == "1":
                output["angles"] += "-".join(words[:3]) + " {b} {k}".format(b = float(words[4]), k = float(words[5])/2)  + "\n"
            elif func == "5":
                output["Urey-Bradley"] += "-".join(words[:3]) + " {b} {k} {b2} {k2}".format(b = float(words[4]), k = float(words[5])/2, b2 = float(words[6]), k2 = float(words[7])/2)  + "\n"
            else:
                raise NotImplementedError
        elif iterator.flag == "dihedraltypes":
            words = line.split()
            func = words[4]
            if func == "1":
                output["dihedrals"] += "-".join(words[:4]) + " " + " ".join(words[5:]) + " 0\n"
            elif func == "2":
                temp1 = [words[1], words[2], words[0], words[3]]
                temp2 = [words[1], words[2], words[3], words[0]]
                if words[0][0].upper() in ("C", "N", "S"):
                    output["impropers"] += "-".join(temp1) + " {b} {k}".format(b = float(words[5]), k = float(words[6])/2)  + "\n"
                else:
                    output["impropers"] += "-".join(temp2) + " {b} {k}".format(b = float(words[5]), k = float(words[6])/2)  + "\n"
            elif func == "9":
                for i in range(5,len(words),20):
                    output["dihedrals"] += "-".join(words[:4]) + " " + " ".join(words[i:i+3]) + " 0\n"
            else:
                raise NotImplementedError
        elif iterator.flag == "cmaptypes":
            words = line.split()
            output["cmaps"]["-".join(words[:5])] = {"resolution": int(words[7]), "parameters":list(map(lambda x:float(x) / 4.184, words[8:]))}
        elif iterator.flag == "nonbond_params":
            words = line.split()
            output["LJ"] += "{type1}-{type2} {V} {W}\n".format(type1=words[0], type2=words[1], V=float(words[3]), W=float(words[4]))
    return output
    
sys.modules['__main__'].__dict__["loadffitp"] = ffitp
