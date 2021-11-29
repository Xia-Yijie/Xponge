from ... import *

AtomType.Add_Property({"LJtype":str})

LJType = Generate_New_Pairwise_Force_Type("LJ", {"epsilon": float, "rmin": float, "sigma":float, "A":float, "B":float})

LJType.Set_Property_Unit("rmin", "distance", "A")
LJType.Set_Property_Unit("sigma", "distance", "A")
LJType.Set_Property_Unit("epsilon", "energy", "kcal/mol")

@GlobalSetting.Add_Unit_Transfer_Function(LJType)
def LJ_Unit_Transfer(self):
    if self.A != None and self.B != None:
        self.sigma = (self.A / self.B) ** (1/6)
        self.epsilon = 0.25 * B * self.sigma ** (-6)
        self.A = None
        self.B = None
    if self.sigma != None:
        self.rmin = self.sigma * (4 ** (1/12) / 2)
        self.sigma = None


def Lorentz_Berthelot_For_A(epsilon1, rmin1, epsilon2, rmin2):
    return np.sqrt(epsilon1 * epsilon2) * ((rmin1 + rmin2) ** 12)

def Lorents_Berthelot_For_B(epsilon1, rmin1, epsilon2, rmin2):
    return np.sqrt(epsilon1 * epsilon2) * 2 * ((rmin1 + rmin2) ** 6)
    


@Molecule.Set_Save_SPONGE_Input      
def write_LJ(self, prefix, dirname):
    LJtypes = []
    LJtypemap = {}
    for atom in self.atoms:
        if atom.LJtype not in LJtypemap.keys():
            LJtypemap[atom.LJtype] = len(LJtypes)
            LJtypes.append(atom.LJtype)
             
    As = []
    Bs = []
    for i in range(len(LJtypes)):
        LJ_i = LJType.types[LJtypes[i] + "-" + LJtypes[i]]
        for j in range(len(LJtypes)):
            LJ_j = LJType.types[LJtypes[j] + "-" + LJtypes[j]]
            finded = False
            findnames = [LJtypes[i] + "-" + LJtypes[j], LJtypes[j] + "-" + LJtypes[i]]
            for findname in findnames:
                if findname in LJType.types.keys():
                    finded = True
                    LJ_ij = LJType.types[findname]
                    As.append(LJType.combining_method_A(LJ_ij.epsilon, LJ_ij.rmin, LJ_ij.epsilon, LJ_ij.rmin))
                    Bs.append(LJType.combining_method_B(LJ_ij.epsilon, LJ_ij.rmin, LJ_ij.epsilon, LJ_ij.rmin))
                    break
            if not finded:
                    As.append(LJType.combining_method_A(LJ_i.epsilon, LJ_i.rmin, LJ_j.epsilon, LJ_j.rmin))
                    Bs.append(LJType.combining_method_B(LJ_i.epsilon, LJ_i.rmin, LJ_j.epsilon, LJ_j.rmin))
    
    checks = {}
    count = 0
    for i in range(len(LJtypes)):
        check_string_A = ""
        check_string_B = ""
        for j in range(len(LJtypes)):
            check_string_A += "%16.7e"%As[count] + " "
            check_string_B += "%16.7e"%Bs[count] + " "
            count += 1
            
        checks[i] = check_string_A+check_string_B
    
    same_type = { i: i for i in range(len(LJtypes))}
    for i in range(len(LJtypes)-1, -1, -1):
        for j in range(i+1, len(LJtypes)):
            if checks[i] == checks[j]:
                same_type[j] = i

    real_LJtypes = []
    real_As = []
    real_Bs = []
    tosub = 0
    for i in range(len(LJtypes)):
        
        if same_type[i] == i:
            real_LJtypes.append(LJtypes[i])
            same_type[i] -= tosub
        else:
            same_type[i] = same_type[same_type[i]]
            tosub += 1

    for i in range(len(real_LJtypes)):
        LJ_i = LJType.types[real_LJtypes[i] + "-" + real_LJtypes[i]]
        for j in range(i+1):
            LJ_j = LJType.types[real_LJtypes[j] + "-" + real_LJtypes[j]]
            finded = False
            findnames = [real_LJtypes[i] + "-" + real_LJtypes[j], real_LJtypes[j] + "-" + real_LJtypes[i]]
            for findname in findnames:
                if findname in LJType.types.keys():
                    finded = True
                    LJ_ij = LJType.types[findname]
                    real_As.append(LJType.combining_method_A(LJ_ij.epsilon, LJ_ij.rmin, LJ_ij.epsilon, LJ_ij.rmin))
                    real_Bs.append(LJType.combining_method_B(LJ_ij.epsilon, LJ_ij.rmin, LJ_ij.epsilon, LJ_ij.rmin))
                    break
            if not finded:
                    real_As.append(LJType.combining_method_A(LJ_i.epsilon, LJ_i.rmin, LJ_j.epsilon, LJ_j.rmin))
                    real_Bs.append(LJType.combining_method_B(LJ_i.epsilon, LJ_i.rmin, LJ_j.epsilon, LJ_j.rmin))
    
                
    
    towrite = "%d %d\n\n"%(len(self.atoms), len(real_LJtypes))
    count = 0
    for i in range(len(real_LJtypes)):
        for j in range(i+1):
            towrite += "%16.7e"%real_As[count] + " "
            count += 1
        towrite +="\n"
    towrite += "\n"
    
    count = 0
    for i in range(len(real_LJtypes)):
        for j in range(i+1):
            towrite += "%16.7e"%real_Bs[count] + " "
            count += 1
        towrite +="\n"
    towrite += "\n"
    towrite += "\n".join(["%d"%(same_type[LJtypemap[atom.LJtype]]) for atom in self.atoms])
    f = open(os.path.join(dirname, prefix + "_LJ.txt"),"w")
    f.write(towrite)
    f.close()
