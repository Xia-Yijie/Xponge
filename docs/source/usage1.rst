building systems
--------------------
.. TIP::
    
    You may need to change the commands and codes according to your environmental variables.

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

building from existing residue type
==============================================

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

    print(Xponge.ResidueType.types)

Then you can get::

    {'ACE': Type of Residue: ACE, 'ASH': Type of Residue: ASH, 'CYM': Type of Residue: CYM, 'GLH': Type of Residue: GLH, 'LYN': Type of Residue: LYN, 'NME': Type of Residue: NME, 'NALA': Type of Residue: NALA, 'ALA': Type of Residue: ALA, 'CALA': Type of Residue: CALA, 'NARG': Type of Residue: NARG, 'ARG': Type of Residue: ARG, 'CARG': Type of Residue: CARG, 'NASN': Type of Residue: NASN, 'ASN': Type of Residue: ASN, 'CASN': Type of Residue: CASN, 'NASP': Type of Residue: NASP, 'ASP': Type of Residue: ASP, 'CASP': Type of Residue: CASP, 'NCYS': Type of Residue: NCYS, 'CYS': Type of Residue: CYS, 'CCYS': Type of Residue: CCYS, 'NCYX': Type of Residue: NCYX, 'CYX': Type of Residue: CYX, 'CCYX': Type of Residue: CCYX, 'NGLN': Type of Residue: NGLN, 'GLN': Type of Residue: GLN, 'CGLN': Type of Residue: CGLN, 'NGLU': Type of Residue: NGLU, 'GLU': Type of Residue: GLU, 'CGLU': Type of Residue: CGLU, 'NGLY': Type of Residue: NGLY, 'GLY': Type of Residue: GLY, 'CGLY': Type of Residue: CGLY, 'NHID': Type of Residue: NHID, 'HID': Type of Residue: HID, 'CHID': Type of Residue: CHID, 'NHIE': Type of Residue: NHIE, 'HIE': Type of Residue: HIE, 'CHIE': Type of Residue: CHIE, 'NHIP': Type of Residue: NHIP, 'HIP': Type of Residue: HIP, 'CHIP': Type of Residue: CHIP, 'NHE': Type of Residue: NHE, 'HYP': Type of Residue: HYP, 'CHYP': Type of Residue: CHYP, 'NILE': Type of Residue: NILE, 'ILE': Type of Residue: ILE, 'CILE': Type of Residue: CILE, 'NLEU': Type of Residue: NLEU, 'LEU': Type of Residue: LEU, 'CLEU': Type of Residue: CLEU, 'NLYS': Type of Residue: NLYS, 'LYS': Type of Residue: LYS, 'CLYS': Type of Residue: CLYS, 'NMET': Type of Residue: NMET, 'MET': Type of Residue: MET, 'CMET': Type of Residue: CMET, 'NPHE': Type of Residue: NPHE, 'PHE': Type of Residue: PHE, 'CPHE': Type of Residue: CPHE, 'NPRO': Type of Residue: NPRO, 'PRO': Type of Residue: PRO, 'CPRO': Type of Residue: CPRO, 'NSER': Type of Residue: NSER, 'SER': Type of Residue: SER, 'CSER': Type of Residue: CSER, 'NTHR': Type of Residue: NTHR, 'THR': Type of Residue: THR, 'CTHR': Type of Residue: CTHR, 'NTRP': Type of Residue: NTRP, 'TRP': Type of Residue: TRP, 'CTRP': Type of Residue: CTRP, 'NTYR': Type of Residue: NTYR, 'TYR': Type of Residue: TYR, 'CTYR': Type of Residue: CTYR, 'NVAL': Type of Residue: NVAL, 'VAL': Type of Residue: VAL, 'CVAL': Type of Residue: CVAL, 'HIS': Type of Residue: HIE, 'NHIS': Type of Residue: NHIE, 'CHIS': Type of Residue: CHIE}
