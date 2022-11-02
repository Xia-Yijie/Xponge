"""
This **module** set the basic configuration for ff14sb
"""
from ...helper import source, Xprint, set_real_global_variable

source("....")
amber = source("...amber")

AtomType.New_From_String(
    """
    name mass    charge[e]      LJtype
    HW   1.008   0.52422        HW
    OW   16      0              OW
    EP   0      -1.04844        EPW
    """)

bond_base.BondType.New_From_String(r"""
name    k[kcal/mol·A^-2]   b[A]
OW-HW   553                0.87243313
HW-HW   553                1.37120510
""")

angle_base.AngleType.New_From_String(r"""
name        k       b
HW-OW-HW    0       0
OW-HW-HW    0       0
""")

lj_base.LJType.New_From_String(r"""
name    epsilon[kcal/mol]   rmin[A]
OW-OW   0.2128008130           1.777167268
HW-HW   0                   0
EPW-EPW 0                   0
""")

virtual_atom_base.VirtualType2.New_From_String(r"""
name   atom0    atom1   atom2   k1         k2
EP     -3       -2      -1      0.1477206  0.1477206
""")

load_mol2(os.path.join(AMBER_DATA_DIR, "opc.mol2"), as_template=True)

ResidueType.set_type("H2O", ResidueType.get_type("WAT"))
ResidueType.set_type("HOH", ResidueType.get_type("WAT"))

set_real_global_variable("H2O", ResidueType.get_type("WAT"))
set_real_global_variable("HOH", ResidueType.get_type("WAT"))

amber.load_parameters_from_frcmod("ions1lm_126_tip4pew.frcmod")
amber.load_parameters_from_frcmod("ionsjc_tip4pew.frcmod")
amber.load_parameters_from_frcmod("ions234lm_126_tip4pew.frcmod")

load_mol2(os.path.join(AMBER_DATA_DIR, "atomic_ions.mol2"), as_template=True)

Xprint("""Reference for opc:
1. Water:
  Saeed Izadi, Ramu Anandakrishnan, and Alexey V. Onufriev
    Building Water Models: A Different Approach
    The Journal of Physical Chemistry Letters 2014 5 (21), 3863-3871, 
    DOI: 10.1021/jz501780a

2. Li+, Na+, K+, Rb+, Cs+, F-, Cl-, Br-, I-:
  In Suk Joung and Thomas E. Cheatham
    Determination of Alkali and Halide Monovalent Ion Parameters for Use in Explicitly Solvated Biomolecular Simulations
    The Journal of Physical Chemistry B 2008 112 (30), 9020-9041
    DOI: 10.1021/jp8001614

3. Ag+, Tl+, Cu+:
  Pengfei Li, Lin Frank Song, and Kenneth M. Merz
    Systematic Parameterization of Monovalent Ions Employing the Nonbonded Model
    Journal of Chemical Theory and Computation 2015 11 (4), 1645-1657, 
    DOI: 10.1021/ct500918t
    
4. Divalent Ions(Ba2+, Mg2+...)
  Pengfei Li and Kenneth M. Merz
    Taking into Account the Ion-Induced Dipole Interaction in the Nonbonded Model of Ions
    Journal of Chemical Theory and Computation 2014 10 (1), 289-297
    DOI: 10.1021/ct400751u

5. Trivalent and Tetravalent Cations(Al3+, Fe3+, Hf4+...)
  Pengfei Li, Lin Frank Song, and Kenneth M. Merz
    Parameterization of Highly Charged Metal Ions Using the 12-6-4 LJ-Type Nonbonded Model in Explicit Water
    The Journal of Physical Chemistry B 2015 119 (3), 883-895
    DOI: 10.1021/jp505875v  
""")
# pylint:disable=undefined-variable
