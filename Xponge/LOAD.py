from . import *
import os
import sys

##########################################################################
#General Format
##########################################################################
def mol2(filename):
    with open(filename) as f:
        #存储读的时候的临时信息，key是编号
        #value是list：原子名(0)、residue(1)、residue编号(2)、是否是新的residue type(3)、该原子(4)、residue type的最新原子(5)
        atom_residue_map = {}  
        flag = None
        nline = 0
        current_molecule = None
        current_residue = None
        current_residue_index = None
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
            #处理分子信息    
            elif flag == "MOLECULE":
                if nline == 1:
                    current_molecule = Molecule(line.strip())
            #处理原子信息
            elif flag == "ATOM":
                words = line.split()
                #新的残基
                if not current_residue_index or int(words[6]) != current_residue_index:
                    current_residue_index = int(words[6])
                    if words[7] in ResidueType.types.keys():
                        temp = False
                    else:
                        sys.modules['__main__'].__dict__[words[7]] = ResidueType(name = words[7])
                        temp = True
                    if current_residue:
                        current_molecule.Add_Residue(current_residue)
                    current_residue = Residue(ResidueType.types[words[7]])
                if temp:
                    current_residue.type.Add_Atom(words[1], AtomType.types[words[5]], *words[2:5])
                    current_residue.type.atoms[-1].update(**{"charge[e]": float(words[8])})
                current_residue.Add_Atom(words[1], AtomType.types[words[5]], *words[2:5])
                current_residue.atoms[-1].update(**{"charge[e]": float(words[8])})
                atom_residue_map[words[0]]=[words[1], current_residue, current_residue_index, temp, current_residue.atoms[-1], current_residue.type.atoms[-1]]
            elif flag == "BOND":
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
    return current_molecule
    
sys.modules['__main__'].__dict__["loadmol2"] = mol2    
                

def pdb(filename, judge_HIS = True):
    molecule = Molecule(os.path.splitext(os.path.basename(filename))[0])
    chain = {}
    SSBOND = []
    residue_type_map = []
    current_residue_count = -1
    current_residue_index = None
    current_HIS = {"DeltaH":False, "EpsilonH":False}
    with open(filename) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resindex = int(line[22:26])
                resname = line[17:20]
                atomname = line[12:16].strip()
                if current_residue_index == None:
                    current_residue_count += 1
                    if judge_HIS and residue_type_map and residue_type_map[-1] in GlobalSetting.HISMap["HIS"].keys():
                        if current_HIS["DeltaH"]: 
                            if current_HIS["EpsilonH"]:
                                residue_type_map[-1] = GlobalSetting.HISMap["HIS"][residue_type_map[-1]]["HIP"]
                            else:
                                residue_type_map[-1] = GlobalSetting.HISMap["HIS"][residue_type_map[-1]]["HID"]
                        else:
                            residue_type_map[-1] = GlobalSetting.HISMap["HIS"][residue_type_map[-1]]["HIE"]
                        current_HIS = {"DeltaH":False, "EpsilonH":False}
                    if resname in GlobalSetting.PDBResidueNameMap["head"].keys():
                        resname = GlobalSetting.PDBResidueNameMap["head"][resname]
                    residue_type_map.append(resname)
                    current_residue_index = resindex
                    chain[chr(ord("A") + len(chain.keys()))] = {resindex:current_residue_count}
                elif current_residue_index != resindex:
                    if judge_HIS and residue_type_map and residue_type_map[-1] in GlobalSetting.HISMap["HIS"].keys():
                        if current_HIS["DeltaH"]: 
                            if current_HIS["EpsilonH"]:
                                residue_type_map[-1] = GlobalSetting.HISMap["HIS"][residue_type_map[-1]]["HIP"]
                            else:
                                residue_type_map[-1] = GlobalSetting.HISMap["HIS"][residue_type_map[-1]]["HID"]
                        else:
                            residue_type_map[-1] = GlobalSetting.HISMap["HIS"][residue_type_map[-1]]["HIE"]
                        current_HIS = {"DeltaH":False, "EpsilonH":False}
                    current_residue_count += 1
                    residue_type_map.append(resname)
                    current_residue_index = resindex
                    chain[chr(ord("A") + len(chain.keys()) - 1)][resindex] = current_residue_count
                if judge_HIS and resname in GlobalSetting.HISMap["HIS"].keys():
                    if atomname == GlobalSetting.HISMap["DeltaH"]:
                        current_HIS["DeltaH"] = True
                    elif atomname == GlobalSetting.HISMap["EpsilonH"]:
                        current_HIS["EpsilonH"] = True
            elif line.startswith("TER"):
                current_residue_index = None
                if residue_type_map[-1] in GlobalSetting.PDBResidueNameMap["tail"].keys():
                   residue_type_map[-1] = GlobalSetting.PDBResidueNameMap["tail"][residue_type_map[-1]]
                if judge_HIS and residue_type_map and residue_type_map[-1] in GlobalSetting.HISMap["HIS"].keys():
                    if current_HIS["DeltaH"]: 
                        if current_HIS["EpsilonH"]:
                            residue_type_map[-1] = GlobalSetting.HISMap["HIS"][residue_type_map[-1]]["HIP"]
                        else:
                            residue_type_map[-1] = GlobalSetting.HISMap["HIS"][residue_type_map[-1]]["HID"]
                    else:
                        residue_type_map[-1] = GlobalSetting.HISMap["HIS"][residue_type_map[-1]]["HIE"]
                    current_HIS = {"DeltaH":False, "EpsilonH":False}
            elif line.startswith("SSBOND"):
                SSBOND.append(line)
    
    current_residue_index = None
    if residue_type_map[-1] in GlobalSetting.PDBResidueNameMap["tail"].keys():
       residue_type_map[-1] = GlobalSetting.PDBResidueNameMap["tail"][residue_type_map[-1]]
    
    SSBONDS = {}
    for ssbond in SSBOND:
        resA = chain[ssbond[15]][int(ssbond[17:21])]
        residue_type_map[resA] = "CYX"
        resB = chain[ssbond[29]][int(ssbond[31:35])]
        residue_type_map[resB] = "CYX"
        if resA > resB:
            resA, resB = (resB, resA)
        SSBONDS[resB] = resA

    current_residue_count = -1
    current_residue = None
    links = []
    with open(filename) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resindex = int(line[22:26])
                atomname = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                if not current_residue_index or current_residue_index != resindex:
                    current_residue_count += 1
                    if current_residue:
                        molecule.Add_Residue(current_residue)
                        if current_residue.type.tail and ResidueType.types[residue_type_map[current_residue_count]].head:
                            links.append(current_residue_count)
                    current_residue = Residue(ResidueType.types[residue_type_map[current_residue_count]])
                    current_residue_index = resindex
                current_residue.Add_Atom(atomname, x = x, y = y, z = z)
            elif line.startswith("TER"):
                current_residue_index = None
                
    if current_residue:
        molecule.Add_Residue(current_residue)
    for resA, resB in SSBONDS.items():
        ResA = molecule.residues[resA]
        ResB = molecule.residues[resB]
        molecule.Add_Residue_Link(ResA._name2atom[ResA.type.connect_atoms["ssbond"]], ResB._name2atom[ResB.type.connect_atoms["ssbond"]])
    for count in links:
        ResA = molecule.residues[count]
        ResB = molecule.residues[count-1]
        molecule.Add_Residue_Link(ResA._name2atom[ResA.type.head], ResB._name2atom[ResB.type.tail])
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
            elif flag == "MASS":
                atom_types[words[0]] = words[1]
            elif flag == "BOND":
                atoms = [ word.strip() for word in line[:5].split("-") ]
                words = line[5:].split()
                bonds += "-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\n"
            elif flag == "ANGL":
                atoms = [ word.strip() for word in line[:8].split("-") ]
                words = line[8:].split()
                angles += "-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\n"
            elif flag == "DIHE":
                atoms = [ word.strip() for word in line[:11].split("-") ]
                words = line[11:].split()
                propers += "-".join(atoms) + "\t" + str(float(words[1]) / int(words[0])) + "\t" + words[2] + "\t" + str(abs(int(float(words[3])))) + "\t" + str(reset) + "\n"
                if int(float(words[3])) < 0:
                    reset = 0
                else:
                    reset = 1
            elif flag == "IMPR":
                atoms = [ word.strip() for word in line[:11].split("-") ]
                words = line[11:].split()
                impropers += "-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\t" + str(int(float(words[2]))) + "\n"
            elif flag == "NONB" or flag == "NONBON":
                words = line.split()
                LJs += words[0] + "-" + words[0] + "\t" + words[1] + "\t" + words[2] + "\n"
            elif flag == "CMAP":
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

##########################################################################
#GROMACS Format
##########################################################################

class GROMACS_TOPOLOGY_ITERATOR():
    def __init__(self, filename = None):
        self.files = []
        if filename:
            self.Add_Iterator_File(filename)
    def Add_Iterator_File(self, filename):
        f = open(filename)
        self.files.append(f)
    def __iter__(self):
        return self
    def __next__(self):
        while self.files:
            f = self.files[-1]
            line = f.readline()
            if line:
                return line
            else:
                f.close()
                self.files.pop()
                continue
        raise StopIteration
            
def itp(filename):
    defined_macros = {}
    macro_define_stat = []
    flag = ""
    iterator = GROMACS_TOPOLOGY_ITERATOR(filename)
    output = {}
    output["nb14"] = "name  kLJ  kee\n"
    output["atomtypes"] = "name mass charge[e] LJtype\n"
    output["bonds"] = "name b[nm] k[kJ/mol·nm^-2]\n"
    output["angles"] = "name b[degree] k[kJ/mol·rad^-2]\n"
    output["Urey-Bradley"] = "name b[degree] k[kJ/mol·rad^-2] r13[nm] kUB[kJ/mol·nm^-2]\n"
    for line in iterator:
        line = line.strip()
        if not line or line[0] == ";":
            continue
        elif line[0] == "#":
            words = line.split()
            if words[0] == "#ifdef":
                macro = words[1]
                if macro in defined_macros.keys():
                    macro_define_stat.append(True)
                else:
                    macro_define_stat.append(False)
            elif words[0] == "#endif":
                macro_define_stat.pop()                 
            elif macro_define_stat and macro_define_stat[-1] == False:
                continue
            elif words[0] == "#define":
                if len(words) > 2:
                    defined_macros[words[1]] = words[2]
                else:
                    defined_macros[words[1]] = None
            elif words[0] == "#include":
                iterator.Add_Iterator_File(words[1])
            elif words[0] == "#undef":
                defined_macros.pop(words[1])
            elif words[0] == "#error":
                raise Exception(line)
        elif macro_define_stat and macro_define_stat[-1] == False:
            continue
        elif "[" in line and "]" in line:
            flag = line[1:-1].strip()
        elif flag == "":
            raise Exception("No flag! Check your file.")
        elif flag == "defaults":
            words = line.split()
            assert int(words[0]) == 1, "SPONGE Only supports Lennard-Jones now"
            if int(words[1]) == 1:
                output["LJ"] = "name A[kJ/mol·nm^6] B[kJ/mol·nm^12]\n"
                output["nb14_extra"] = "name A[kJ/mol·nm^6] B[kJ/mol·nm^12] kee\n"
            else:
                output["LJ"] = "name epsilon[kJ/mol] sigma[nm]\n"
                output["nb14_extra"] = "name epsilon[kJ/mol] sigma[nm] kee\n"
            fudgeLJ = float(words[3])
            fudgeQQ = float(words[4])
            if words[2] == "yes":
                output["nb14"] += "X-X {fudgeLJ} {fudgeQQ}\n".format(fudgeLJ=fudgeLJ, fudgeQQ=fudgeQQ)
                
        elif flag == "atomtypes":
            words = line.split()
            output["atomtypes"] += "{type} {mass} {charge} {type}\n".format(type=words[0], mass = float(words[2]), charge = float(words[3]))
            output["LJ"] += "{type}-{type} {V} {W}\n".format(type=words[0], V=float(words[5]), W=float(words[6]))
        elif flag == "pairtypes":
            words = line.split()
            if len(words) <= 3:
                output["nb14"] +=  "{atom1}-{atom2} {kLJ} {kee}\n".format(atom1 = words[0], atom2 = words[1], kLJ = fudgeLJ, kee = fudgeQQ)
            elif words[2] == "1":
                output["nb14_extra"] += "{atom1}-{atom2} {V} {W} {kee}\n".format(atom1 = words[0], atom2 = words[1], V = float(words[3]), W = float(words[4]), kee = fudgeQQ)
            elif words[2] == "2":
                raise NotImplementedError
        elif flag == "bondtypes":
            words = line.split()
            func = words[2]
            if func == "1":
                output["bonds"] += "{atom1}-{atom2} {b} {k}\n".format(atom1 = words[0], atom2 = words[1], b = float(words[3]), k = float(words[4]))
            else:
                raise NotImplementedError
        elif flag == "angletypes":
            words = line.split()
            func = words[3]
            if func == "1":
                output["angles"] += "-".join(words[:3]) + " " + " ".join(words[4:]) + "\n"
            elif func == "5":
                output["Urey-Bradley"] += "-".join(words[:3]) + " " + " ".join(words[4:])  + "\n"
            else:
                raise NotImplementedError
    return output
       

sys.modules['__main__'].__dict__["loaditp"] = itp 