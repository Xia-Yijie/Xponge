from ... import *

AtomType.Add_Property({"LJtype":str})

LJType = Generate_New_Pairwise_Force_Type("LJ", {"epsilon": float, "rmin": float, "sigma":float, "A":float, "B":float})

LJType.Set_Property_Unit("rmin", "distance", "A")
LJType.Set_Property_Unit("sigma", "distance", "A")
LJType.Set_Property_Unit("epsilon", "energy", "kcal/mol")
LJType.Set_Property_Unit("A", "energy·distance^6", "kcal/mol·A^6")
LJType.Set_Property_Unit("B", "energy·distance^12", "kcal/mol·A^12")

@GlobalSetting.Add_Unit_Transfer_Function(LJType)
def LJ_Unit_Transfer(self):
    if self.A != None and self.B != None:
        if self.B == 0 or self.A == 0:
            self.sigma = 0
            self.epsilon = 0
        else:
            self.sigma = (self.A / self.B) ** (1/6)
            self.epsilon = 0.25 * self.B * self.sigma ** (-6)
        self.A = None
        self.B = None
    if self.sigma != None:
        self.rmin = self.sigma * (4 ** (1/12) / 2)
        self.sigma = None


def Lorentz_Berthelot_For_A(epsilon1, rmin1, epsilon2, rmin2):
    return np.sqrt(epsilon1 * epsilon2) * ((rmin1 + rmin2) ** 12)

def Lorents_Berthelot_For_B(epsilon1, rmin1, epsilon2, rmin2):
    return np.sqrt(epsilon1 * epsilon2) * 2 * ((rmin1 + rmin2) ** 6)

def find_AB_LJ(LJtypes, stat = True):
    As = []
    Bs = []
    
    for i in range(len(LJtypes)):
        LJ_i = LJType.types[LJtypes[i] + "-" + LJtypes[i]]
        if stat:
            j_max = len(LJtypes)
        else:
            j_max = i + 1
        for j in range(j_max):
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
    return As, Bs

def get_checks(LJtypes, As, Bs):
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
    return checks

def judge_same_type(LJtypes, checks):
    same_type = { i: i for i in range(len(LJtypes))}
    for i in range(len(LJtypes)-1, -1, -1):
        for j in range(i+1, len(LJtypes)):
            if checks[i] == checks[j]:
                same_type[j] = i
    return same_type

def get_real_LJ(LJtypes, same_type):
    real_LJtypes = []
    tosub = 0
    for i in range(len(LJtypes)):
        if same_type[i] == i:
            real_LJtypes.append(LJtypes[i])
            same_type[i] -= tosub
        else:
            same_type[i] = same_type[same_type[i]]
            tosub += 1
    return real_LJtypes

def write_LJ(self):
    LJtypes = []
    LJtypemap = {}
    for atom in self.atoms:
        if atom.LJtype not in LJtypemap.keys():
            LJtypemap[atom.LJtype] = len(LJtypes)
            LJtypes.append(atom.LJtype)
             
    As, Bs = find_AB_LJ(LJtypes)
    
    checks = get_checks(LJtypes, As, Bs)
    
    same_type = judge_same_type(LJtypes, checks)
    
    real_LJtypes = get_real_LJ(LJtypes, same_type)

    real_As, real_Bs = find_AB_LJ(real_LJtypes, False)
    
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
    return towrite

Molecule.Set_Save_SPONGE_Input("LJ")(write_LJ)