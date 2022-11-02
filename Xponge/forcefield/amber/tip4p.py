"""
This **module** set the basic configuration for ff14sb
"""
from ...helper import source, Xprint, set_real_global_variable

source("....")
amber = source("...amber")

AtomType.New_From_String(
    """
    name mass    charge[e]      LJtype
    HW   1.008   0.52000        HW
    OW   16      0              OW
    EP   0      -1.04000        EPW
    """)

bond_base.BondType.New_From_String(r"""
name    k[kcal/mol·A^-2]   b[A]
OW-HW   553                0.9572
HW-HW   553                1.5136
""")

angle_base.AngleType.New_From_String(r"""
name        k       b
HW-OW-HW    0       0
OW-HW-HW    0       0
""")

lj_base.LJType.New_From_String(r"""
name    epsilon[kcal/mol]   sigma[A]
OW-OW   0.155               3.154
HW-HW   0                   0
EPW-EPW 0                   0
""")

virtual_atom_base.VirtualType2.New_From_String(r"""
name   atom0    atom1   atom2   k1         k2
EP     -3       -2      -1      0.1279696  0.1279696
""")

TIP4P = load_mol2(os.path.join(AMBER_DATA_DIR, "tip4pew.mol2"), as_template=True)

ResidueType.set_type("H2O", ResidueType.get_type("WAT"))
ResidueType.set_type("HOH", ResidueType.get_type("WAT"))

set_real_global_variable("H2O", ResidueType.get_type("WAT"))
set_real_global_variable("HOH", ResidueType.get_type("WAT"))

Xprint("""Reference for tip4p:
  William L. Jorgensen, Jayaraman Chandrasekhar, Jeffry D. Madura, Roger W. Impey and Michael L. Klein
    Comparison of simple potential functions for simulating liquid water
    The Journal of Chemical Physics 1983, 79, 926-935
    DOI: 10.1063/1.445869

""")
# pylint:disable=undefined-variable
