forcefield development
----------------------
.. include:: namespace.rst

new parameters
====================

string/text file/python dict
################################

Take an example of nitro-methane in the form of an OPLSAA force field constructed in the combination of  AMBER force field. The data were obtained from `GMXTPP <https://jerkwin.github.io/prog/gmxtop.html>`_.

.. code-block::

    import Xponge
    import Xponge.forcefield.amber

    # New atom types from string
    Xponge.AtomType.New_From_String("""name  mass  charge[e] LJtype
    opls_760  14.01   0.54   opls_760
    opls_761  16.00  -0.37   opls_761
    opls_762  12.01   0.02   opls_762
    opls_763  1.008   0.06   opls_763""")

    # New LJ information from string
    Xponge.forcefield.base.lj_base.LJType.New_From_String("""name  sigma[nm]  epsilon[kJ/mol]
    opls_760-opls_760 0.3250000    0.5020800
    opls_761-opls_761 0.2960000    0.7112800
    opls_762-opls_762 0.3500000    0.2761440
    opls_763-opls_763 0.2500000    0.0627600  
    """)

    # New bond information from dict
    temp_dict =  {"opls_762-opls_763": {"b[nm]":0.109, "k[kJ/mol·nm^-2]": 284512},
                  "opls_762-opls_760": {"b[nm]":0.14900, "k[kJ/mol·nm^-2]": 313800},
                  "opls_760-opls_761": {"b[nm]":0.1225, "k[kJ/mol·nm^-2]": 460240}}
    Xponge.forcefield.base.bond_base.BondType.New_From_Dict(temp_dict)

    # New angle information from a text file
    Xponge.forcefield.base.angle_base.AngleType.New_From_File("test_angle.txt")

    # New proper dihedral information from string
    Xponge.forcefield.base.dihedral_base.ProperType.New_From_String("""name phi0[degree] k[kJ/mol] periodicity  reset
    opls_763-opls_762-opls_760-opls_761  0 0 0 1
    """)

    # New improper dihedral information from string
    Xponge.forcefield.base.dihedral_base.ImproperType.New_From_String("""name phi0[degree] k[kJ/mol] periodicity
    opls_761-opls_761-opls_760-X 180.0     43.93200   2
    """)

    # New non bonded 1-4 energy parameter from string
    Xponge.forcefield.base.nb14_base.NB14Type.New_From_String("""name kLJ kee
    opls_763-opls_761 0.5 0.5
    """)

    # New Residue Type
    NIM = Xponge.ResidueType(name = "NIM")
    # Add atoms to the new residue type
    NIM.Add_Atom(name = "CT", atom_type = Xponge.AtomType.types["opls_762"], x = 1.107, y = 1.13, z = 1.117)
    NIM.Add_Atom(name = "HC1", atom_type = Xponge.AtomType.types["opls_763"], x = 1.143, y = 1.029, z = 1.117)
    NIM.Add_Atom(name = "HC2", atom_type = Xponge.AtomType.types["opls_763"], x = 1.143, y = 1.18, z = 1.03)
    NIM.Add_Atom(name = "HC3", atom_type = Xponge.AtomType.types["opls_763"], x = 1., y = 1.13, z = 1.117)
    NIM.Add_Atom(name = "NO", atom_type = Xponge.AtomType.types["opls_760"], x = 1.156, y = 1.199, z = 1.237)
    NIM.Add_Atom(name = "ON1", atom_type = Xponge.AtomType.types["opls_761"], x = 1.177, y = 1.321, z = 1.234)
    NIM.Add_Atom(name = "ON2", atom_type = Xponge.AtomType.types["opls_761"], x = 1.177, y = 1.135, z = 1.341)

    # Add connectivity to the new residue type
    NIM.Add_Connectivity(NIM.CT, NIM.HC1)
    NIM.Add_Connectivity(NIM.CT, NIM.HC2)
    NIM.Add_Connectivity(NIM.CT, NIM.HC3)
    NIM.Add_Connectivity(NIM.CT, NIM.NO)
    NIM.Add_Connectivity(NIM.ON1, NIM.NO)
    NIM.Add_Connectivity(NIM.ON2, NIM.NO)

    Save_PDB(NIM)
    Save_SPONGE_Input(NIM)

new forcefield format
=======================

non-bonded interactions
##########################

The information of non bonded forces can be seen as the property of the atom. Take the electrostatic force as an example, which is implemented in `Xponge.forcefield.base.charge_base`::

    import Xponge
    # Add a new property to the AtomType class
    Xponge.AtomType.Add_Property({"charge":float})

    # Set the unit of charge
    # The meaning of the three parameters is that the unit of the proper "charge" (the first parameter) is the same as "charge" (the second parameter)，and the basic unit is "e" (the third parameter).
    Xponge.AtomType.Set_Property_Unit("charge", "charge", "e")

    #The decorator "@Xponge.Molecule.Set_Save_SPONGE_Input" is used to append the following function to the list which will be used when "Save_SPONGE_Input"
    #See the API document for detailed information
    @Xponge.Molecule.Set_Save_SPONGE_Input      
    def write_charge(self):
        towrite = "%d\n"%(len(self.atoms))
        towrite += "\n".join(["%.6f"%(atom.charge * 18.2223) for atom in self.atoms])
        return towrite

bonded interactions
##########################

The simplest example is the harmonic bond in `Xponge.forcefield.base.bond_base`::

    from ... import *

    # See the API document for detailed information
    BondType = Generate_New_Bonded_Force_Type("bond", "1-2", {"k":float, "b":float}, True)

    # Set units
    BondType.Set_Property_Unit("k", "energy·distance^-2", "kcal/mol·A^-2")
    BondType.Set_Property_Unit("b", "distance", "A")

    @Molecule.Set_Save_SPONGE_Input("bond")
    def write_bond(self):
        bonds = []
        for bond in self.bonded_forces.get("bond",[]):
            order = list(range(2))
            if bond.k != 0:
                if self.atom_index[bond.atoms[order[0]]] > self.atom_index[bond.atoms[order[-1]]]:
                    temp_order = order[::-1]
                else:
                    temp_order = order
                bonds.append("%d %d %f %f"%(self.atom_index[bond.atoms[temp_order[0]]]
                , self.atom_index[bond.atoms[temp_order[1]]], bond.k, bond.b))
        
        if (bonds):
            towrite = "%d\n"%len(bonds)
            bonds.sort(key = lambda x: list(map(int, x.split()[:2])))
            towrite += "\n".join(bonds)
            
            return towrite

new forcefield combination
============================

A combination of force field formats and parameters is necessary for a new force field.

Take the AMBER-style force field as an example::


    import os
    from ...helper import set_global_alternative_names
    from ... import AtomType, load_parmdat, load_frcmod

    from ..base import charge_base, mass_base, lj_base, bond_base, angle_base, dihedral_base, nb14_base,\
        virtual_atom_base, exclude_base

    AMBER_DATA_DIR = os.path.dirname(__file__)

    lj_base.LJType.combining_method_A = lj_base.Lorentz_Berthelot_For_A
    lj_base.LJType.combining_method_B = lj_base.Lorentz_Berthelot_For_B

    nb14_base.NB14Type.New_From_String(r"""
    name    kLJ     kee
    X-X     0.5     0.833333
    """)

    exclude_base.Exclude(4)


new type assignment
=======================

Take gaff in `Xponge.forcefield.amber.gaff` as an example::

    from Xponge import assign

    # Register the new assignment judging rules
    GAFF = assign.Judge_Rule("gaff")

    # Using "GAFF.Add_Judge_Rule" to add new judging functions
    # Xponge will use the judging functions one by one, and if the function returns "True", Xponge will stop and returns the atom type

    # The following codes are the complete judgement of oxygen atoms
    @GAFF.Add_Judge_Rule("o")
    def temp(i, Assign):
        return Assign.Atom_Judge(i, "O1")

    @GAFF.Add_Judge_Rule("op")
    def temp(i, Assign):
        return Assign.Atom_Judge(i, "O2") and "RG3" in Assign.atom_marker[i].keys()

    @GAFF.Add_Judge_Rule("oq")
    def temp(i, Assign):
        return Assign.Atom_Judge(i, "O2") and "RG4" in Assign.atom_marker[i].keys()

    @GAFF.Add_Judge_Rule("oh")
    def temp(i, Assign):
        tofind = False
        if Assign.Atom_Judge(i, "O2") or Assign.Atom_Judge(i, "O3"):
            for bonded_atom in Assign.bonds[i].keys():
                if Assign.Atom_Judge(bonded_atom, "H1"):
                    tofind = True
                    break    
        return tofind

    @GAFF.Add_Judge_Rule("os")
    def temp(i, Assign):
        return Assign.Atom_Judge(i, "O2")

After the definition of the judging functions, then you can use the Xponge.assign.Assign instance method ``Determine_Atom_Type("gaff")`` to determine the atom type.
