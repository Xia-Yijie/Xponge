from . import *

##########################################################################
#General Format
##########################################################################
def mol2(filename):
    with open(filename) as f:
        flag = None
        nline = 0
        current_molecule = None
        current_residue = None
        new_residue_type = {}
        atom_residue_map = {}   #存储读的时候的临时信息，key是编号，value是list包括：原子名0、residue1、residue编号2、是否是新的residue type3、该原子4、residue type的最新原子5
        current_residue_index = None
        temp_link_count = 0
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
                        ResidueType(name = words[7])
                        temp = True
                    if current_residue:
                        current_molecule.Add_Residue(current_residue)
                    current_residue = Residue(ResidueType.types[words[7]])
                    new_residue_type[current_residue] = temp
                if temp:
                    current_residue.type.Add_Atom(words[1], AtomType.types[words[5]], *words[2:5])
                    current_residue.type.atoms[-1].update(**{"charge[e]": float(words[8])})
                current_residue.Add_Atom(words[1], AtomType.types[words[5]], *words[2:5])
                current_residue.atoms[-1].update(**{"charge[e]": float(words[8])})
                atom_residue_map[words[0]]=[words[1], current_residue, current_residue_index, temp, current_residue.atoms[-1], current_residue.type.atoms[-1]]
            elif flag == "BOND":
                words = line.split()
                if atom_residue_map[words[1]][1] == atom_residue_map[words[2]][1]:
                    if new_residue_type[atom_residue_map[words[1]][1]]:
                        atom_residue_map[words[1]][1].type.Add_Connectivity(atom_residue_map[words[1]][0], atom_residue_map[words[2]][0])
                    atom_residue_map[words[1]][1].Add_Connectivity(atom_residue_map[words[1]][0], atom_residue_map[words[2]][0])
                else:
                    pass
    return current_molecule
                
                

def pdb(filename):
    with open(filename) as f:
        pass


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
            if len(words) == 1:
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
            elif flag == "NONB":
                words = line.split()
                LJs += words[0] + "-" + words[0] + "\t" + words[1] + "\t" + words[2] + "\n"

    atoms = "name  mass  LJtype\n"
    for atom, mass in atom_types.items():
        atoms += atom + "\t" + mass + "\t" + atom + "\n"
        
    return atoms, bonds, angles, propers, impropers, LJs
            
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
        
        
        