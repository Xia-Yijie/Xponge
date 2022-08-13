building systems
--------------------
.. include:: namespace.rst

Residue
    Literally, a residue is a single molecular unit within a polymer. Here, the word is more general, which represents a single molecular unit within a molecule, no matter the molecule is a small organic molecule or a polymer.  

**Residue** is the basic unit in ``Xponge``, and there is no significant difference to build systems with one residue or a large number of residues.

existing force fields
=======================

Here are the existing force fields in ``Xponge`` now:

.. list-table:: 
    :header-rows: 1
    :stub-columns: 1
    
    * - Force Field Style
      - Molecule Type
      - Force Field
      - Module
    * - AMBER
      - proteins
      - ff14SB
      - Xponge.forcefield.amber.ff14sb
    * - 
      - proteins
      - ff19SB
      - Xponge.forcefield.amber.ff19sb
    * - 
      - lipids
      - lipid14
      - Xponge.forcefield.amber.lipid14
    * - 
      - lipids
      - lipid17
      - Xponge.forcefield.amber.lipid17
    * - 
      - organic molecules
      - gaff
      - Xponge.forcefield.amber.gaff
    * - 
      - water/ions
      - tip3p
      - Xponge.forcefield.amber.tip3p
    * - 
      - water/ions
      - spce
      - Xponge.forcefield.amber.spce
    * - 
      - water/ions
      - tip4pew
      - Xponge.forcefield.amber.tip4pew
    * - 
      - water/ions
      - opc
      - Xponge.forcefield.amber.opc
    * - CHARMM27
      - proteins
      - protein
      - Xponge.forcefield.charmm27.protein
    * - 
      - water/ions
      - tip3p
      - Xponge.forcefield.charmm27.tip3p
    * - 
      - water/ions
      - tips3p
      - Xponge.forcefield.charmm27.tips3p

building known residue types
==============================================

from a python script
######################

Here, let's build an alanine tetrapeptide molecule in vacuum under the ff14SB force field.

First of all, you need to import the force field module,

.. code-block:: python

    import Xponge
    import Xponge.forcefield.amber.ff14sb

You can get some information on the references::

    Reference for ff14SB.py:
      James A. Maier, Carmenza Martinez, Koushik Kasavajhala, Lauren Wickstrom, Kevin E. Hauser, and Carlos Simmerling
        ff14SB: Improving the accuracy of protein side chain and backbone parameters from ff99SB
        Journal of Chemical Theory and Computation 2015 11 (8), 3696-3713
        DOI: 10.1021/acs.jctc.5b00255

Of course, you can set the verbose level less than 0 to disable the reference print by

.. code-block:: python

    Xponge.GlobalSetting.verbose = -1

You can print all of the residue types imported. (This is not compulsory to build the system, but helpful for you to understand Xponge)::

    print(Xponge.ResidueType.get_all_types())

Then you can get::

    {'ACE': Type of Residue: ACE, 'ASH': Type of Residue: ASH, 'CYM': Type of Residue: CYM, 'GLH': Type of Residue: GLH, 'LYN': Type of Residue: LYN, 'NME': Type of Residue: NME, 'NALA': Type of Residue: NALA, 'ALA': Type of Residue: ALA, 'CALA': Type of Residue: CALA, 'NARG': Type of Residue: NARG, 'ARG': Type of Residue: ARG, 'CARG': Type of Residue: CARG, 'NASN': Type of Residue: NASN, 'ASN': Type of Residue: ASN, 'CASN': Type of Residue: CASN, 'NASP': Type of Residue: NASP, 'ASP': Type of Residue: ASP, 'CASP': Type of Residue: CASP, 'NCYS': Type of Residue: NCYS, 'CYS': Type of Residue: CYS, 'CCYS': Type of Residue: CCYS, 'NCYX': Type of Residue: NCYX, 'CYX': Type of Residue: CYX, 'CCYX': Type of Residue: CCYX, 'NGLN': Type of Residue: NGLN, 'GLN': Type of Residue: GLN, 'CGLN': Type of Residue: CGLN, 'NGLU': Type of Residue: NGLU, 'GLU': Type of Residue: GLU, 'CGLU': Type of Residue: CGLU, 'NGLY': Type of Residue: NGLY, 'GLY': Type of Residue: GLY, 'CGLY': Type of Residue: CGLY, 'NHID': Type of Residue: NHID, 'HID': Type of Residue: HID, 'CHID': Type of Residue: CHID, 'NHIE': Type of Residue: NHIE, 'HIE': Type of Residue: HIE, 'CHIE': Type of Residue: CHIE, 'NHIP': Type of Residue: NHIP, 'HIP': Type of Residue: HIP, 'CHIP': Type of Residue: CHIP, 'NHE': Type of Residue: NHE, 'HYP': Type of Residue: HYP, 'CHYP': Type of Residue: CHYP, 'NILE': Type of Residue: NILE, 'ILE': Type of Residue: ILE, 'CILE': Type of Residue: CILE, 'NLEU': Type of Residue: NLEU, 'LEU': Type of Residue: LEU, 'CLEU': Type of Residue: CLEU, 'NLYS': Type of Residue: NLYS, 'LYS': Type of Residue: LYS, 'CLYS': Type of Residue: CLYS, 'NMET': Type of Residue: NMET, 'MET': Type of Residue: MET, 'CMET': Type of Residue: CMET, 'NPHE': Type of Residue: NPHE, 'PHE': Type of Residue: PHE, 'CPHE': Type of Residue: CPHE, 'NPRO': Type of Residue: NPRO, 'PRO': Type of Residue: PRO, 'CPRO': Type of Residue: CPRO, 'NSER': Type of Residue: NSER, 'SER': Type of Residue: SER, 'CSER': Type of Residue: CSER, 'NTHR': Type of Residue: NTHR, 'THR': Type of Residue: THR, 'CTHR': Type of Residue: CTHR, 'NTRP': Type of Residue: NTRP, 'TRP': Type of Residue: TRP, 'CTRP': Type of Residue: CTRP, 'NTYR': Type of Residue: NTYR, 'TYR': Type of Residue: TYR, 'CTYR': Type of Residue: CTYR, 'NVAL': Type of Residue: NVAL, 'VAL': Type of Residue: VAL, 'CVAL': Type of Residue: CVAL, 'HIS': Type of Residue: HIE, 'NHIS': Type of Residue: NHIE, 'CHIS': Type of Residue: CHIE}

All of the residue type will be stored in the main dict too (which means you can directly use it as a global variable in Python), and the residues can be linked by adding and multiplying, so you can use the following simple codes to get an alanine tetrapeptide molecule::

    protein = ACE + ALA * 3 + NME
    # this sets the periodic conditions 
    protein.box_length = [50.0, 50.0, 50.0]

New in 1.2.16, we can get a peptide molecule from its sequence now::

    protein = get_peptide_from_sequence("AAAAA")
    # protein = NALA + ALA * 3 + CALA

Now what you need to do is just to save your ``protein`` in a format you like. The functions to save are  ``Save_PDB``, ``Save_SPONGE_Input`` and ``Save_Mol2`` now::

    Save_Mol2(protein, "protein.mol2")
    Save_PDB(protein, "protein.pdb")
    Save_SPONGE_Input(protein, "protein")

You can use `VMD <https://www.ks.uiuc.edu/Research/vmd>`_ to visualize the pdb file or the mol2 file::

    vmd protein.pdb
    # or 
    # vmd protein.mol2

You can also install the `SPONGE VMD plugin <https://spongemm.cn/download/files/tools/sponge_vmd_v1.2.1.zip>`_ to VMD to visualize the SPONGE input file::

    vmd -sponge_mass ./protein_mass.txt -sponge_crd ./protein_coordinate.txt

Then you can see,

.. image:: https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/1.png

Also, you can do the molecular dynamics simulations by `SPONGE <https://spongemm.cn>`_.

Write a "mdin.txt" file of SPONGE::

    basic test of Xponge
        mode = NVT
        thermostat = Andersen_thermostat
        dt = 2e-3
        constrain_mode = SHAKE
        cutoff = 8.0
        step_limit = 5000
        default_in_file_prefix = protein

Run SPONGE::

    SPONGE -mdin mdin.txt
    #Or you can use according to your environment setting
    #Xponge.mdrun SPONGE -mdin mdin.txt

Then we can get::

    ---------------------------------------------------------------------------------------
            step =         1000,         time =        2.000,  temperature =       305.35,
       potential =        47.21,           LJ =        -3.73,          PME =      -232.74,
         nb14_LJ =        26.62,      nb14_EE =       196.30,         bond =        12.74,
           angle =        16.63,     dihedral =        32.51,
    ---------------------------------------------------------------------------------------
            step =         2000,         time =        4.000,  temperature =       194.94,
       potential =        51.01,           LJ =        -2.02,          PME =      -229.18,
         nb14_LJ =        30.27,      nb14_EE =       188.48,         bond =        12.04,
           angle =        22.31,     dihedral =        30.75,
    ---------------------------------------------------------------------------------------
            step =         3000,         time =        6.000,  temperature =       246.57,
       potential =        49.60,           LJ =        -4.42,          PME =      -222.75,
         nb14_LJ =        31.28,      nb14_EE =       189.09,         bond =        11.10,
           angle =        16.47,     dihedral =        30.80,
    ---------------------------------------------------------------------------------------
            step =         4000,         time =        8.000,  temperature =       273.81,
       potential =        55.53,           LJ =        -4.54,          PME =      -231.54,
         nb14_LJ =        30.79,      nb14_EE =       196.70,         bond =         9.60,
           angle =        18.88,     dihedral =        29.40,
    ---------------------------------------------------------------------------------------
            step =         5000,         time =       10.000,  temperature =       344.12,
       potential =        48.44,           LJ =        -5.35,          PME =      -230.99,
         nb14_LJ =        29.09,      nb14_EE =       198.53,         bond =         6.32,
           angle =        19.48,     dihedral =        31.79,
    ---------------------------------------------------------------------------------------


.. TIP::

    Some contents are the same in the following examples, so some results will not be shown in the following examples.

The complete python script is::

    import Xponge
    import Xponge.forcefield.amber.ff14sb
    protein = ACE + ALA * 3 + NME
    protein.box_length = [50.0, 50.0, 50.0]
    Save_Mol2(protein, "protein.mol2")
    Save_PDB(protein, "protein.pdb")
    Save_SPONGE_Input(protein, "protein")

from a pdb file
################

First of all, you need to import the force field module::

    import Xponge
    import Xponge.forcefield.amber.ff14sb
    
Then, you can load your file. Here I will use `this file <https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/0.15_80_10_pH6.5_6lzg.result.pdb>`_ as an example::

    protein = load_pdb("0.15_80_10_pH6.5_6lzg.result.pdb")

The complete python script is::

    import Xponge
    import Xponge.forcefield.amber.ff14sb
    protein = loadpdb("0.15_80_10_pH6.5_6lzg.result.pdb")
    Save_PDB(protein, "protein.pdb")
    Save_SPONGE_Input(protein, "protein")

.. ATTENTION::

    1. Protein Data Bank (PDB) format is a standard for files containing atomic coordinates. If you want to modify a PDB file yourself, please make sure that you understand the format
    2. Every residue type in the pdb file should be in Xponge.ResidueType.get_all_types()

.. TIP::    
    
    The pdb files downloaded from the online PDB database, `RCSB <https://www.rcsb.org/>`_ for example, do not include hydrogen. The simplest way to do is to use ``protein.Add_Missing_Atoms()``, but some more complicated pre-processes by online services such as `H++ <http://newbiophysics.cs.vt.edu/H++>`_ may be more accurate.


building unknown residue types
==============================================

A new Xponge.ResidueType instance needs to be created, and then the following steps are the same as `building known residue types`_

There are two ways to create a new Xponge.ResidueType instance, and here I take a water residue under gaff as an example here.

from a mol2 file
#################

The first way to create a new Xponge.ResidueType instance is to load a mol2 file. Here we write a mol2 file named "WAT.mol2"::

    @<TRIPOS>MOLECULE
    TP3
        3     2     1     0     1 
    SMALL
    USER_CHARGES
    @<TRIPOS>ATOM
      1 O       0.000000    0.000000    0.000000 oh    1 WAT    -0.8340 ****
      2 H1      0.957200    0.000000    0.000000 ho    1 WAT     0.4170 ****
      3 H2     -0.239988    0.926627    0.000000 ho    1 WAT     0.4170 ****
    @<TRIPOS>BOND
        1     1     2 1
        2     1     3 1
    @<TRIPOS>SUBSTRUCTURE
          1  WAT              1 ****               0 ****  **** 

Then in python::

    import Xponge
    import Xponge.forcefield.amber.gaff
    TP3 = load_mol2("WAT.mol2")
    Save_SPONGE_Input(TP3, "tp3")


.. ATTENTION::

    1. Xponge will create a new Xponge.ResidueType according to the residue names in the mol2 file ("WAT" here), and the return value of `load_mol2` is a Xponge.Molecule ("TP3" here). Only the residue type name will be used when loading the PDB files.
    2. The atom type in the mol2 file should match the atom types in the force field, for the example here, the atom type for the oxygen atom and the hydrogen atoms should be the atom types in gaff "oh" and "ho", instead of the element name "O" and "H".

from a python script
#####################

It is also okay to direct use a python script::

    import Xponge
    import Xponge.forcefield.amber.gaff
    #Create a new residue type
    WAT = Xponge.ResidueType(name = "WAT")
    #Add atoms to the residue type
    WAT.Add_Atom(name = "O", atom_type = "oh", x = 0, y = 0, z = 0)
    WAT.O.Update(**{"charge[e]": -0.8340})
    WAT.Add_Atom(name = "H1", atom_type = "ho", x = 0.9572, y = 0, z = 0)
    WAT.H1.Update(**{"charge[e]": 0.4170})
    WAT.Add_Atom(name = "H2", atom_type = "ho", x = -0.239988, y = 0.926627, z = 0)
    WAT.H2.Update(**{"charge[e]": 0.4170})
    #Add connectivity to the residue type
    WAT.Add_Connectivity(WAT.O, WAT.H1)
    WAT.Add_Connectivity(WAT.O, WAT.H2)

    Save_SPONGE_Input(WAT, "TP3")


linkage information
#######################

If you do not want to use ``+`` or ``*`` to get a new molecule from the residue types, you do not need to do more now. Otherwise, you should give the linkage information. Here, we will make a "fake alanine" (FALA) residue type linkable as an example.

.. TIP::

    In fact, all the information of  an alanine dipeptide molecule is known in ff14SB. We just take an example here.
    

Write a mol2 file named dipeptide.mol2 first::

    @<TRIPOS>MOLECULE
    ACE
     22 21 3 0 1
    SMALL
    USER_CHARGES
    @<TRIPOS>ATOM
         1   H1    0.466   -8.051    1.242   HC     1      ACE   0.112300
         2  CH3    0.289   -7.339    0.436   CT     1      ACE  -0.366200
         3   H2   -0.703   -6.901    0.548   HC     1      ACE   0.112300
         4   H3    0.352   -7.853   -0.524   HC     1      ACE   0.112300
         5    C    1.325   -6.213    0.457    C     1      ACE   0.597200
         6    O    2.209   -6.195    1.311    O     1      ACE  -0.567900
         7    N    1.275   -5.267   -0.433    N     2      FALA  -0.415700
         8    H    0.563   -5.249   -1.150    H     2      FALA   0.271900
         9   CA    2.256   -4.201   -0.413   CX     2      FALA   0.033700
        10   HA    2.199   -3.671    0.538   H1     2      FALA   0.082300
        11   CB    3.668   -4.752   -0.580   CT     2      FALA  -0.182500
        12  HB1    3.744   -5.277   -1.532   HC     2      FALA   0.060300
        13  HB2    4.384   -3.930   -0.561   HC     2      FALA   0.060300
        14  HB3    3.888   -5.443    0.234   HC     2      FALA   0.060300
        15    C    2.009   -3.207   -1.539    C     2      FALA   0.597300
        16    O    1.072   -3.368   -2.318    O     2      FALA  -0.567900
        17    N    2.791   -2.179   -1.682    N     3      NME  -0.415700
        18    H    3.571   -2.011   -1.063    H     3      NME   0.271900
        19  CH3    2.555   -1.233   -2.754   CT     3      NME  -0.149000
        20 HH31    1.686   -1.546   -3.331   H1     3      NME   0.097600
        21 HH32    2.374   -0.244   -2.333   H1     3      NME   0.097600
        22 HH33    3.429   -1.195   -3.405   H1     3      NME   0.097600
    @<TRIPOS>BOND
         1      1      2 1
         2      2      3 1
         3      2      4 1
         4      2      5 1
         5      5      6 1
         6      5      7 1
         7      7      8 1
         8      7      9 1
         9      9     10 1
        10      9     11 1
        11      9     15 1
        12     11     12 1
        13     11     13 1
        14     11     14 1
        15     15     16 1
        16     15     17 1
        17     17     18 1
        18     17     19 1
        19     19     20 1
        20     19     21 1
        21     19     22 1
    @<TRIPOS>SUBSTRUCTURE
        1      ACE      1 ****               0 ****  **** 
        2      ALA      7 ****               0 ****  **** 
        3      NME     17 ****               0 ****  **** 

Then load the mol2 file and write the linkage information::

    import Xponge
    import Xponge.forcefield.amber.ff14sb

    protein = load_mol2("dipeptide.mol2")
    
    res = Xponge.ResidueType.get_type("FALA")
    #The head of the residue（the atom to link the previous residue）is the atom with a atom name "N"
    res.head = "N"
    #The bond length to link the previous residue is 1.3 angstrom
    res.head_length = 1.3
    #The head next atom of the residue（the second atom to link the previous residue）is the atom with a atom name "CA"
    res.head_next = "CA"
    #The tail of the residue（the atom to link the latter residue）is the atom with a atom name "C"
    res.tail = "C"
    #The bond length to link the latter residue is 1.3 angstrom
    res.tail_length = 1.3
    #The tail next atom of the residue（the second atom to link the latter residue）is the atom with a atom name "CA"
    res.tail_next = "CA"

    import numpy as np
    #When linking to the previous residue, the angle of the CA-N-tail of the preivous residue is 120/180 * np.pi
    res.head_link_conditions.append({"atoms":["CA", "N"], "parameter": 120/180 * np.pi})
    #When linking to the previous residue, the dihedral of the H-CA-N-tail of the preivous residue is -np.pi
    res.head_link_conditions.append({"atoms":["H", "CA", "N"], "parameter": -np.pi})
    #When linking to the latter residue, the angle of the CA-C-head of the latter residue is 120/180 * np.pi
    res.tail_link_conditions.append({"atoms":["CA", "C"], "parameter": 120/180 * np.pi})
    #When linking to the latter residue, the dihedral of the O-CA-C-head of the latter residue is -np.pi
    res.tail_link_conditions.append({"atoms":["O", "CA", "C"], "parameter": -np.pi})

.. TIP::

    There are 6 degrees of freedom when linking two residues, 1 bond, 2 angles and 3 dihedral, where the only bond is determined by the bond length, and one of the dihedrals is the main chain dihedral (tail_next - tail - head - head_next). The main chain dihedral is set to 180 degrees to prevent overlap. Define an angle and a dihedral condition for the previous residue and the next residue respectively to fix the remaining 4 degrees.

.. NOTE::

    In PDB file, the residue name should be not greater than 3 letters. It is all right to use a 4-letter residue name in python, but it will cause error if you want to load or write a pdb file.
    
If you want the terminal residue to use a different structure when loading the PDB file, then you need to do like this::

    Xponge.GlobalSetting.Add_PDB_Residue_Name_Mapping("head", "FALA", "NALA")
    Xponge.GlobalSetting.Add_PDB_Residue_Name_Mapping("tail", "FALA", "CALA")

histidine and cysteine 
#######################

In particular, additional settings are required to  discriminate different proton states for histidine and disulfide bonds for cysteine::

    # These codes have been written in Xponge. You do not need to do this for the existing histidine and cysteine.
    Xponge.GlobalSetting.HISMap["DeltaH"] = "HD1"
    Xponge.GlobalSetting.HISMap["EpsilonH"] = "HE2"
    Xponge.GlobalSetting.HISMap["HIS"].update({"HIS": {"HID":"HID", "HIE":"HIE", "HIP":"HIP"}, 
                                        "CHIS":{"HID":"CHID", "HIE":"CHIE", "HIP":"CHIP"},
                                        "NHIS":{"HID":"NHID", "HIE":"NHIE", "HIP":"NHIP"}})
    Xponge.ResidueType.types["CYX"].connect_atoms["ssbond"] = "SG"

