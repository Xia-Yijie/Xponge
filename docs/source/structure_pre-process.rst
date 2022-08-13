structure pre-process
---------------------
.. include:: namespace.rst

Inner Coordinate Imposing
==========================

You can use ``Impose_Bond``, ``Impose_Angle``, ``Impose_Dihedral`` to impose the inner coordinates.

Impose_Bond
############################

.. code-block:: python

    import Xponge
    import Xponge.forcefield.amber.ff14sb
    t = ACE + NME
    Save_Mol2(t, "imposing.mol2")
    Impose_Bond(t, t.residues[0].C, t.residues[1].N, 5)
    Save_Mol2(t, "imposed.mol2")
    Impose_Bond(t, t.residues[0].C, t.residues[1].CH3, 5)
    Save_Mol2(t, "imposed2.mol2")
    #Atoms that are not bonded in the same residue, and atoms that in the same ring cannot be "Impose_Bond", because it is not known to move which atoms to move
    #Impose_Bond(t, t.residues[0].C, t.residues[0].H2, 5)
    #The code above will raise AssertionError

Visualize the three mol2 files in VMD and we can get 

.. image:: https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/3.png

.. image:: https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/4.png

.. image:: https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/5.png

Impose_Angle
############################

.. code-block:: python

    import Xponge
    import Xponge.forcefield.amber.ff14sb
    t = ACE + NME
    #Impose_Angle rotate the atoms after the third atom and the third atom around the line from the first atom to the second atom 
    Impose_Angle(t, t.residues[0].C, t.residues[1].N, t.residues[1].CH3, 3.1415926 / 2)
    Save_Mol2(t, "imposed.mol2")

Visualize the file in VMD and we can get 

.. image:: https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/6.png

Impose_Dihedral
############################

.. code-block:: python

    import Xponge
    import Xponge.forcefield.AMBER.ff14SB
    t = ACE + ALA * 10 + NME
    Save_Mol2(t, "imposing.mol2")
    for i in range(1,len(t.residues)-1):
        head = t.residues[i-1]
        res = t.residues[i]
        tail = t.residues[i+1]
    #Impose_Dihedral rotate the atoms after the third atom and the third atom around the line from the third atom to the second atom 
        Impose_Dihedral(t, head.C, res.N, res.CA, res.C, -3.1415926/3)
        Impose_Dihedral(t, res.N, res.CA, res.C, tail.N, -3.1415926/3)
    Save_Mol2(t, "imposed.mol2")

Visualize the file "imposing.mol2" in VMD and we can get

.. image:: https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/7.png

.. image:: https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/8.png

Visualize the file "imposed.mol2" in VMD and we can get

.. image:: https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/9.png

.. image:: https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/10.png

Add Solvents and Ions
==========================

You can use ``add_solvent_box`` to add solvents, and ``solvent_replace`` to replace some solvents to ions or small molecules::

    import Xponge
    import Xponge.forcefield.amber.ff14sb
    import Xponge.forcefield.amber.tip3p
    t = NALA + ARG + NME
    c = int(round(t.charge))
    
    #Add solvent box, 5 Å for x<0, 10 Å for y<0, 15 Å for z< 0, 30 Å for x>0, 25 Å for y>0, 20 Å for z>0 from the boundary
    add_solvent_box(t, WAT, [5,10,15,30,25,20])
    
    #It is also possible to simply provide a number
    #add_solvent_box(t, WAT, 30)
    
    #Substitution, replacing the part of the residue named "WAT" with a potassium ion or a chloride atom
    solvent_replace(t, WAT, {CL:30 + c, K:30})
    
    #check the charge to make sure it's right
    assert int(round(t.charge)) == 0
    
    Save_PDB(t, "test.pdb")
    Save_SPONGE_Input(t, "test")

The Main Axis Rotation
=======================

You can use ``main_axis_rotate`` to rotate the main axis::

    import Xponge
    import Xponge.forcefield.amber.ff14sb

    t = ACE + ALA * 100 + NME

    Save_Mol2(t,"before.mol2")

    #This function rotate the main axis of the molecule to the disired direction. See the API document for detail 
    Molecule_Rotate(t)

    Save_Mol2(t,"after.mol2")

Visualize the file "before.mol2" in VMD and we can get

.. image:: https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/11.png

Visualize the file "after.mol2" in VMD (look along the z-axis) and we can get

.. image:: https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/12.png

Hydrogen Mass Repartition
===========================

You can use ``HMass_Repartition`` to repartition the mass of hydrogen::

    import Xponge
    import Xponge.forcefield.amber.ff14sb
    import Xponge.forcefield.amber.tip3p

    t = ACE + ALA + NME

    Process_Box(t, WAT, [5,10,15,30,25,20])

    HMass_Repartition(t)

    Save_SPONGE_Input(t, "test")

.. NOTE::

    Some functions guess the element from the mass, so please repartition only before your saving.
