type assignment
--------------------
.. include:: namespace.rst

Type Assignment
    There are many types in a force field. For example, every atom has its own atom type. It is important to assign the correct atom type for every atom.


from PubChem
=============================

The basic information can be directly obtained from the PubChem database by name or smiles. Here, take the ethyl 2,6-dimethylbenzoate under gaff as an example::

    import Xponge
    import Xponge.forcefield.amber.gaff

    # Get the basic information from the PubChem database by name
    assign = Get_Assignment_From_PubChem("ethyl 2,6-dimethylbenzoate", "name")
    # by smiles
    #assign = Get_Assignment_From_PubChem("CCOC(=O)C1=C(C=CC=C1C)C", "smiles")

    # Determin the atom type for each atom
    assign.Determine_Atom_Type("GAFF")

    # Save the assignment to a mol2 file
    Save_Mol2(assign, "EDF_ASN.mol2")

    # Get the chemical equivalent atoms and calculate the RESP charge
    equivalence = assign.Determine_Equal_Atoms()
    # Here in fact, "q" is None, the real charges are stored in assign.charge.
    q = assign.Calculate_Charge("RESP", basis = "6-311g*", grid_density = 1, 
        extra_equivalence = equivalence, opt = True)
    # You can give your own charge, too
    #q = np.array([0 for i in range(assign.atom_numbers)])

    # Convert the assignment to a Xponge.ResidueType
    EDF = assign.To_ResidueType("EDF", q)

    # If you want to use GAFF and do not know the parameters in the default GAFF module, you need to run the following line. This is a standard step for GAFF.
    Xponge.forcefield.AMBER.gaff.parmchk2_gaff(EDF, "EDF.frcmod")

    Save_Mol2(EDF, "EDF.mol2")


from a mol2 file
=================

You can also get an assignment from a standard mol2 file.

.. ATTENTION::

    The mol2 file for the assignments is in the standard format, pay attention to the difference in the bond orders and the atom types.

Simply replace the  ``Get_Assignment_From_PubChem`` with ``Get_Assignment_From_Mol2``, and do what we do in `from PubChem`_::

    assign = Get_Assignment_From_Mol2("EDF_ASN.mol2")
    # ^   ^   ^   ^   ^   ^   ^   ^   ^   ^   ^   ^   ^
    # |   |   |   |   |   |   |   |   |   |   |   |   |
    #assign = Get_Assignment_From_PubChem("ethyl 2,6-dimethylbenzoate", "name")

from a python script
=====================

Here, take water under gaff as an example::

    import Xponge
    import Xponge.forcefield.amber.gaff as gaff
    assign = Xponge.Assign()
    assign.Add_Atom(element = "O", x = 0, y = 0, z = 0)
    assign.Add_Atom(element = "H", x = 0.9572, y = 0, z = 0)
    assign.Add_Atom(element = "H", x = -0.239988, y = 0.926627, z = 0)
    assign.Add_Bond(0,1,1)
    assign.Add_Bond(0,2,1)

    # ATTENTION, detailed information on rings and bond types should be determined before the atom type assigning, though no influence on water
    assign.Determine_Ring_And_Bond_Type()

    assign.Determine_Atom_Type("gaff")
    assign.Calculate_Charge("RESP", extra_equivalence = [[1,2]])

    WAT = assign.To_ResidueType("WAT")

    # If you want to use GAFF and do not know the parameters in the default GAFF module, you need to run the following line. This is a standard step for GAFF.
    gaff.parmchk2_gaff(WAT, "WAT.frcmod")

    Save_Mol2(WAT, "WAT.mol2")

