#this is used to modify mol2 more conviniently


def rtp2mol2_part(rtp, mol2, oldresname, newresname):
    towrite = ""
    lines = rtp.split("\n")
    mol2lines = mol2.split("\n")
    for i, line in enumerate(lines):
        words = line.split()
        front = mol2lines[i][:45]
        middle = mol2lines[i][48:60].replace(oldresname, newresname)
        towrite += front + " %3s"%words[1] + middle + " {:>5s}".format(words[2])+"\n"
    return towrite

def reorder(order):
    def func(rtp):
        lines = rtp.split("\n")
        return "\n".join([lines[i] for i in order])
    return func
    

def replace_standart_COO(rtp):
    lines = rtp.split("\n")[:-2]
    lines = "\n".join(lines) + """\n	C	CC	0.34	8
	OT1	OC	-0.67	9
    OT2 OC  -0.67   10"""
    return lines
        

def replace_standart_NH3(rtp):
    lines = rtp.split("\n")[4:]
    lines = """    N   NH3 -0.3
    H1  HC  0.33
    H2  HC  0.33
    H3  HC  0.33
    CA  CT1 0.21
    HA  HB  0.10\n""" + "\n".join(lines)
    return lines


def rtp2mol2(rtp, mol2, old_new_resnamemap = {}, gmx_sponge_resnamemap = {}, gmx_special_process = {}, to = ""):
    lines_front = ""
    lines_back = ""
    with open(mol2) as fmol2:
        t = fmol2.read()
    lines_front, t = t.split("@<TRIPOS>ATOM")
    lines_front += "@<TRIPOS>ATOM\n"
    t, lines_back = t.split("@<TRIPOS>BOND")
    lines_back = "@<TRIPOS>BOND" + lines_back
    lastresid = -1
    mol2line = ""
    mol2lines = []
    mol2resname = []
    for line in t.split("\n"):
        if not line.strip():
            continue
        words = line.split()
        if int(words[6]) != lastresid:
            if mol2line:
                mol2lines.append(mol2line)
            mol2line = line + "\n"
            mol2resname.append(words[7].strip())
            lastresid = int(words[6])
        else:
            mol2line += line + "\n"
    mol2lines.append(mol2line)
    real_old_new_resnamemap = { i:i for i in mol2resname}
    real_old_new_resnamemap.update(old_new_resnamemap)
    real_gmx_sponge_resnamemap = {i:i for i in mol2resname}
    real_gmx_sponge_resnamemap.update(gmx_sponge_resnamemap)
    with open(rtp) as frtp:
        t = [ i for i in frtp.read().split("\n")]
    for i in range(len(t)-1, -1, -1):
        index = t[i].find(";")
        if index >= 0:
            t[i] = t[i][:index]
        t[i] = t[i].strip()
        if not t[i]:
            t.pop(i)
    newlines = []
    for i, old_resname in enumerate(mol2resname):
        start = t.index("[ %s ]"%real_gmx_sponge_resnamemap[old_resname]) + 2
        end = start
        while 1:
            end += 1
            line = t[end]
            if "[" in line:
                break
        templine = "\n".join(t[start:end])
        templine = gmx_special_process.get(old_resname, lambda x:x)(templine)
        newlines.append(rtp2mol2_part(templine, mol2lines[i], old_resname, real_old_new_resnamemap[old_resname]))
        #print(newlines[-1], end="")
    f = open(to,"w")
    f.write(lines_front + "\n".join(newlines) + lines_back)
    f.close()


gmx_sponge_resnamemap = {"ASH":"ASPP", "CYM":"CYS2", "GLH":"GLUP", "LYN":"LSN", "NME":"CT3", "CYX":"CYS2", "NCYX":"CYS2", "CCYX": "CYS2",
                         "NGLY":"GLY", "CGLY":"GLY", "NHE":"CT3", "HYP":"PRO", "CHYP":"PRO", "PRO":"PRO", "NPRO":"PRO", "CPRO":"PRO"}


def NGLY(rtp):
    lines = rtp.split("\n")[3:]
    lines = """    N   NH3 -0.3
    H1  HC  0.33
    H2  HC  0.33
    H3  HC  0.33
    CA  CT2 0.13\n""" + "\n".join(lines)
    return lines

def NPRO(rtp):
    lines = rtp.split("\n")[4:]
    lines = """    N   NP -0.07
    H1  HC  0.24
    H2  HC  0.24
    CD  CP3 0.16
    HD1 HA  0.09
    HD2 HA  0.09
    CG  CP2 -0.18
    HG1 HA  0.09
    HG2 HA  0.09
    CB  CP2 -0.18
    HB1 HA  0.09
    HB2 HA  0.09
    CA  CP1 0.16
    HA  HB  0.09
    C   C   0.51
    O   O   -0.51"""
    return lines

gmx_special_process = {"ACE":reorder([1,0,2,3,4,5]), 
                      "NCYX":replace_standart_NH3, 
                      "CCYX":replace_standart_COO,
                      "NGLY": NGLY,
                      "CGLY":replace_standart_COO,
                      "NPRO": NPRO,
                      "CPRO": replace_standart_COO,
                      "NHE": lambda x: "    N   NH2 -0.62\n    H1  H   0.31\n    H2  H   0.31"}

for res in ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "ILE", "LEU", "LYS", "MET", "PHE", "SER", "THR", "TRP", "TYR", "VAL"]:
    gmx_sponge_resnamemap["N" + res] = res
    gmx_special_process["N"+res] = replace_standart_NH3
    gmx_sponge_resnamemap["C" + res] = res
    gmx_special_process["C"+res] = replace_standart_COO  


for res in ["HID", "HIE", "HIP"]:
    res_gmx = res.replace("I","S")
    gmx_sponge_resnamemap[res] = res_gmx
    gmx_sponge_resnamemap["N" + res] = res_gmx
    gmx_special_process["N" + res] = replace_standart_NH3
    gmx_sponge_resnamemap["C" + res] = res_gmx
    gmx_special_process["C"+res] = replace_standart_COO  

#No HYP for CHARMM27
gmx_special_process["HYP"] = lambda rtp: "\n".join(["    N   NH2 -0.62" for i in range(15)])
gmx_special_process["CHYP"] = lambda rtp: "\n".join(["    N   NH2 -0.62" for i in range(16)])


rtp2mol2("./aminoacids.rtp", "./ff19SB.mol2", {"CYM":"CYX"}, gmx_sponge_resnamemap, gmx_special_process, "../Xponge/forcefield/CHARMM27/protein.mol2")