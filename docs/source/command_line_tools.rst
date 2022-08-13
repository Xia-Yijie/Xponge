command line tools
--------------------
.. include:: namespace.rst

Xponge has two command line tools in general ``Xponge`` and ``Xponge.mdrun``.

The arguments which you can pass to the tools can be seen by ``Xponge --help`` and ``Xponge.mdrun --help`` 

Xponge
=========

The command line tool ``Xponge`` has multiple sub-tools to do some frequent used and complicated calculations

.. code-block::

    usage: Xponge [-h] {test,maskgen,exgen,name2name,mol2rfe,converter} ...

    optional arguments:
      -h, --help            show this help message and exit

    subcommands:
      Tools for SPONGE. Use Xponge XXX -h for the help of tool 'XXX'.

      {test,maskgen,exgen,name2name,mol2rfe,converter}
                            subcommands
        test                test the basic function of Xponge
        maskgen             use VMD to generate a file to record the atom indexes
                            of the corresponding mask
        exgen               process bond-like, angle-like, dihedral-like files to
                            get the atoms to exclude
        name2name           change the atom names of a residue from one file to
                            another file
        mol2rfe             calculate the relative free energy of a small molecule
                            using SPONGE
        converter           convert the format of coordinate file

test
##########

.. code-block::

    usage: Xponge test [-h] [-o test] [-v -1] [-d [todo [todo ...]]]

    optional arguments:
      -h, --help            show this help message and exit
      -o test               the prefix for the output files
      -v -1, --verbose -1   the verbose level for output, 1 or -1
      -d [todo [todo ...]], --do [todo [todo ...]]
                            the unit tests need to do, should be 'all', or one or
                            more of 'base', 'assign', 'charmm27'

Here are the list all the tests you can do:

.. list-table::
    :widths: 10 50
    :header-rows: 1
    :stub-columns: 1
    
    * - command
      - what to test
    * - base
      - the basic funtions of ``Xponge``
    * - assign
      - the assignment functions and dependency
    * - charmm27
      - the forcefield of charmm27
    * - all
      - all above

maskgen
###############

maskgen is used to generate the atom index list file according to VMD mask::

    usage: Xponge maskgen [-h] -p P [-c C] -o O [--vmd vmd]

    optional arguments:
      -h, --help  show this help message and exit
      -p P        the topology file
      -c C        the coordinate file
      -o O        the output file
      --vmd vmd   the command to start vmd

exgen
###############

exgen is used to generate the exclude file according to the bonded force files::

    usage: Xponge exgen [-h] -n N -o O [-b BOND [BOND ...]] [-a ANGLE [ANGLE ...]]
                        [-d DIHEDRAL [DIHEDRAL ...]] [-v VIRTUAL [VIRTUAL ...]]
                        [-e EXCLUDE [EXCLUDE ...]]

    optional arguments:
      -h, --help            show this help message and exit
      -n N                  the atom numbers
      -o O                  output exclude file name
      -b BOND [BOND ...], --bond BOND [BOND ...]
                            bond-like input files: skip the first line, and there
                            are 2 atoms in the head of following lines
      -a ANGLE [ANGLE ...], --angle ANGLE [ANGLE ...]
                            angle-like input files: skip the first line, and there
                            are 3 atoms in the head of following lines
      -d DIHEDRAL [DIHEDRAL ...], --dihedral DIHEDRAL [DIHEDRAL ...]
                            dihedral-like input files: skip the first line, and
                            there are 4 atoms in the head of following lines
      -v VIRTUAL [VIRTUAL ...], --virtual VIRTUAL [VIRTUAL ...]
                            virtual-atom-like input files: the first number
                            indicates the virtual type
      -e EXCLUDE [EXCLUDE ...], --exclude EXCLUDE [EXCLUDE ...]
                            exclude-like input files: add the information of
                            another exclude file

name2name
###############

name2name is used to match the atom names in ``ffile`` and ``tfile``::

    usage: Xponge name2name [-h] -fformat {mol2,pdb,gaff_mol2} -ffile FROM_FILE
                            [-fres FROM_RESIDUE] -tformat {mol2,pdb,gaff_mol2}
                            -tfile TO_FILE [-tres TO_RESIDUE] -oformat
                            {mol2,pdb,mcs_pdb} -ofile OUT_FILE [-ores OUT_RESIDUE]
                            [-tmcs TMCS]

    optional arguments:
      -h, --help            show this help message and exit
      -fformat {mol2,pdb,gaff_mol2}, -from_format {mol2,pdb,gaff_mol2}
                            the format of the file which is needed to change from
      -ffile FROM_FILE, -from_file FROM_FILE
                            the name of the file which is needed to change from
      -fres FROM_RESIDUE, -from_residue FROM_RESIDUE
                            the residue name in ffile if fformat == pdb
      -tformat {mol2,pdb,gaff_mol2}, -to_format {mol2,pdb,gaff_mol2}
                            the format of the file which is needed to change to
      -tfile TO_FILE, -to_file TO_FILE
                            the name of the file which is needed to change to
      -tres TO_RESIDUE, -to_residue TO_RESIDUE
                            the residue name in tfile if tformat == pdb
      -oformat {mol2,pdb,mcs_pdb}, -out_format {mol2,pdb,mcs_pdb}
                            the format of the output file
      -ofile OUT_FILE, -out_file OUT_FILE
                            the name of the output file
      -ores OUT_RESIDUE, -out_residue OUT_RESIDUE
                            the name of the output residue
      -tmcs TMCS            the time to find max common structure

The meaning of the choices of format:

.. list-table::
    :widths: 10 50
    :header-rows: 1
    :stub-columns: 1
    
    * - choice
      - meaning
    * - mol2
      - the original mol2 format
    * - pdb
      - the original pdb format
    * - gaff_mol2
      - the modified mol2 format: atom types should be the type in gaff.
    * - mcs_pdb
      - the original pdb format, but only the max common structure of ``ffile`` and ``tfile``

converter
###############

converter is used to convert the format of the coordinate file::

    usage: Xponge converter [-h] -p TOP [-c CRD] -o OUT [-cf GUESS] [-of GUESS]

    optional arguments:
      -h, --help  show this help message and exit
      -p TOP      the name of the topology file
      -c CRD      the name of the coordinate file
      -o OUT      the name of the output file
      -cf GUESS   the format of the topology file, can be "guess", "sponge_crd" or
                  "sponge_traj"
      -of GUESS   the format of the output file, can be "guess", "sponge_crd" or
                  "sponge_traj"

The meaning of the choices of format:

.. list-table::
    :widths: 10 50
    :header-rows: 1
    :stub-columns: 1
    
    * - choice
      - meaning
    * - guess
      - the file format except the SPONGE coordinate file and the SPONGE trajectory file. The file format will be guessed by its suffix.
    * - sponge_crd
      - the SPONGE coordinate file
    * - sponge_traj
      - the SPONGE trajectory file

mol2rfe
####################

mol2rfe is used to calculate the relative binding free energy of molecules::

    usage: Xponge mol2rfe [-h] [-do [todo [todo ...]]] -pdb PDB -r2 R2 -r1 R1
                          [-r0 [R0 [R0 ...]]] [-ri 0] [-nl 20] [-dohmr] [-ff FF]
                          [-mi [MI [MI ...]]] [-pi PI] [-ei EI] [-ai AI]
                          [-method {TI}] [-temp TMP] [-tmcs 10] [-dt dt]
                          [-msteps MSTEPS MSTEPS MSTEPS MSTEPS MSTEPS MSTEPS]
                          [-pstep pre_equilibrium_step] [-estep 500000]
                          [-thermostat middle_langevin]
                          [-barostat andersen_barostat]

    optional arguments:
      -h, --help            show this help message and exit
      -do [todo [todo ...]]
                            the things need to do, should be one or more of
                            'build', 'min', 'pre_equilibrium', 'equilibrium',
                            'analysis'
      -pdb PDB              the initial conformation given by the pdb file
      -r2 R2, -residuetype2 R2
                            molecule mutated to by an Xponge ResidueType mol2 file
      -r1 R1, -residuetype1 R1
                            molecule mutated from by an Xponge ResidueType mol2
                            file
      -r0 [R0 [R0 ...]], -residuetype0 [R0 [R0 ...]]
                            small molecules that do not mutate
      -ri 0, -residue_index 0
                            the residue index of the molecule to mutate
      -nl 20, -lambda_numbers 20
                            the number of lambda groups - 1, default 20 for 0,
                            0.05, 0.10, 0.15..., 1.0
      -dohmr, -do_hydrogen_mass_repartition
                            use the hydrogen mass repartition method
      -ff FF, -forcefield FF
                            Use this force field file instead of the default
                            ff14SB and gaff
      -mi [MI [MI ...]], -min_mdin [MI [MI ...]]
                            Use the minimization mdin file(s) here instead of the
                            default ones
      -pi PI, -pre_equilibrium_mdin PI
                            Use this pre-equilibrium mdin file instead of the
                            default one
      -ei EI, -equilibrium_mdin EI
                            Use this equilibrium mdin file instead of the default
                            one
      -ai AI, -analysis_mdin AI
                            Use this analysis mdin file instead of the default one
      -method {TI}          the method to calculate the free energy
      -temp TMP             the temporary file name prefix
      -tmcs 10              the timeout parameter for max common structure in unit
                            of second
      -dt dt                the dt used for simulation when mdin is not provided
      -msteps MSTEPS MSTEPS MSTEPS MSTEPS MSTEPS MSTEPS
                            the minimization steps for all the lambda. Default
                            5000 for each minimization simulation. There are 6
                            minimization simulations.
      -pstep pre_equilibrium_step, -pre_equilibrium_step pre_equilibrium_step
                            the pre-equilibrium step used for simulation when mdin
                            is not provided
      -estep 500000, -equilibrium_step 500000
                            the equilibrium step used for simulation when mdin is
                            not provided
      -thermostat middle_langevin
                            the thermostat used for simulation when mdin is not
                            provided
      -barostat andersen_barostat
                            the barostat used for simulation when mdin is not
                            provided

Xponge.mdrun
==============

The command line tool ``Xponge.mdrun`` is used to call the MD program SPONGE, which is convinient to set the environmental variables if you build SPONGE and Xponge together.

.. code-block::

     mdrun: run the SPONGE md simulation
        Usage:
            mdrun, mdrun -h, mdrun --help: see this help
            mdrun -set BIN_PATH: set the SPONGE bin direction path to BIN_PATH
                                 BIN_PATH can be an absolute path or a relative path to this module file
            mdrun SPONGE*:  run SPONGE

There are 4 subtools now, and here is a list of the subtools:

.. list-table::
    :widths: 10 50
    :header-rows: 1
    :stub-columns: 1
    
    * - tool
      - description
    * - SPONGE
      - the normal MD simulation with periodic box conditions
    * - SPONGE_NOPBC
      - the normal MD simulation without periodic box conditions
    * - SPONGE_TI
      - the MD calculation for thermodynamic integration
    * - SPONGE_FEP
      - the MD calculation for free energy perturbation

.. toctree::
    :maxdepth: 3
    :caption: input arguments

    sponge_api
    