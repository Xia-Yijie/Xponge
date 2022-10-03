"""
This **package** sets the basic configuration of RSFF2C
"""
from ...helper import source, Xdict
source("..ff14sb")
source("...base.cmap_base")

def _init():
    """
    Initialize the module
    """
    rsff2c_dihedral = Generate_New_Bonded_Force_Type("dihedral", "1-2-3-4",
                                                     {"k": float, "phi0": float, "periodicity": int},
                                                     False, ["k", "phi0", "periodicity"])
    rsff2c_dihedral.Set_Property_Unit("k", "energy", "kcal/mol")
    rsff2c_dihedral.Set_Property_Unit("phi0", "angle", "rad")
    @rsff2c_dihedral.Type_Name_Getter
    def build_check(atoms):
        return "-".join([atom.residue.name + "@" + atom.name if len(atom.residue.name) <= 3
                         else atom.residue.name[-3:] + "@" + atom.name for atom in atoms])
    source_data = """
    dih_E_chi2 = '-0.457   0.453  -0.105   0.269'
    dih_E_chi3 = ' 0.028   0.201   0.052   0.025'
    dih_Q_chi2 = ' 0.603   0.672   0.341   0.311'
    dih_Q_chi3 = '-0.187   0.012  -0.440   0.032'
    dih_N_chi2 = '-0.343   0.484  -0.154  -0.108'
    """.strip().replace("\n", "=").split("=")
    coeffs = Xdict({source_data[i].strip() : [float(j) for j in source_data[i + 1].split("'")[1].split()]
                    for i in range(0, len(source_data), 2)})
    phase = Xdict({key: [ 3.1415926 if x > 0 else 0 for x in coeff] for key, coeff in coeffs.items()})
    to_update = "name k phi0 periodicity reset\n"
    names = Xdict({"GLU@CA-GLU@CB-GLU@CG-GLU@CD": "dih_E_chi2",
                   "GLU@CB-GLU@CG-GLU@CD-GLU@OE1": "dih_E_chi3",
                   "GLU@CB-GLU@CG-GLU@CD-GLU@OE2": "dih_E_chi3",
                   "GLN@CA-GLN@CB-GLN@CG-GLN@CD": "dih_Q_chi2",
                   "GLN@CB-GLN@CG-GLN@CD-GLN@OE1": "dih_Q_chi3",
                   "ASN@CA-ASN@CB-ASN@CG-ASN@OD1": "dih_N_chi2"})
    for name, pname in names.items():
        to_update += f"{name} {abs(coeffs[pname][0])} {phase[pname][0]} 1 1\n"
        to_update += f"{name} {abs(coeffs[pname][1])} {phase[pname][1]} 2 0\n"
        to_update += f"{name} {abs(coeffs[pname][2])} {phase[pname][2]} 3 0\n"
        to_update += f"{name} {abs(coeffs[pname][3])} {phase[pname][3]} 4 0\n"

    rsff2c_dihedral.New_From_String(to_update)

    rsff2c_cmap = Generate_New_Bonded_Force_Type("cmap",
                                                 "1-2-3-4-5",
                                                 {"resolution": int, "parameters": list}, False)
    @rsff2c_cmap.Type_Name_Getter
    def _(atoms):
        atom_names = ["G" if ((atom.residue.name == "VAL" and atom.name == "CG2") or
                              (atom.residue.name != "VAL" and atom.name in ('CG', 'OG', 'SG', 'CG1', 'OG1')))
                      else atom.name for atom in atoms]
        res_name = atoms[2].residue.name
        if "-".join(atom_names) == "C-N-CA-C-N" and res_name != "GLY":
            return "C-N-ALA@CA-C-N"
        atom_names[2] = res_name + "@" + atom_names[2]
        return "-".join(atom_names)

    with open(os.path.join(AMBER_DATA_DIR, "RSFF2C.dat")) as f:
        source_data = f.read().strip().replace("%FLAG", "%FORMAT(8(F9.5))").split("%FORMAT(8(F9.5))")[1:]
    names = [line.split("COMMENT")[1].split() for line in source_data[::2]]
    temp_map = Xdict({"phi/psi": "C-N-%s@CA-C-N",
                      "chi/phi": "G-CB-%s@CA-N-C",
                      "chi/psi": "G-CB-%s@CA-C-N"})
    names = [temp_map[name[1]]%name[0] for name in names]
    parameters = [[float(word) for word in line.split()] for line in source_data[1::2]]
    to_update = Xdict()
    for name, parameter in zip(names, parameters):
        to_update[name] = {"resolution":24, "parameters": parameter}
    rsff2c_cmap.New_From_Dict(to_update)


_init()

Xprint("""Reference for RSFF2C:
  Wei Kang, Fan Jiang, Yun-Dong Wu
    Universal Implementation of a Residue-Specific Force Field Based on CMAP Potentials and Free Energy Decomposition
    Journal of Chemical Theory and Computation 2018, 14(8), 4474–4486
    DOI: 10.1021/acs.jctc.8b00285
""")
