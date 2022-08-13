commands for SPONGE
--------------------------

the way to pass the commands
=======================================

The input arguments can be passed to SPONGE in 2 ways. The first way is to directly pass the command in the terminal, for example::

    SPONGE -cutoff 8.0 -skin 2.0

The second way is to write an ``mdin`` file. SPONGE will search the "mdin.txt" file in the current working path and read the commands in it by default. The file name to search can be modified by::

    SPONGE -mdin your_mdin_file_name

Except the command "mdin", you can combine the two ways together, for example::

    SPONGE -mdin your_mdin_file_name -coordinate_in_file example_coordinate.txt

the format of the mdin file
===========================

The first line should always be the title of the md task.

The commands look like ``command = value``::

    mode = NPT
    dt = 1e-3

Some commands with the same prefix can be merged into one ``{}``::

    middle_langevin
    {
        seed = 1
        gamma = 1
        velocity_max = 20.0
    }

The example above is totally the same as the following one::

    middle_langevin_seed = 1
    middle_langevin_gamma = 1
    middle_langevin_velocity_max = 20.0

single line comments start with ``#``, and multi line comments are surrounded by ``## anything { }``::

    #This is a comment
    #这是一条注释
    
    ## anything 
    {
        this is a comment
        this is a comment
        this is a comment
        这是一条注释
        这是一条注释
        anything可以换为任何字符
    }

command documentation
========================

basic setting
##############

mode = NVE
    can be rerun, minimization, NVE, NVT or NPT

dt = 1e-3
    [ps] time step

step_limit = 1000
    the number of the simulation steps

write_information_interval = 1000 
    the number of interval steps between saving the information

write_restart_file_interval = 1000  
    the number of interval steps between saving the restart file

molecule_map_output = 0
    whether to move the atoms in one molecule to the same periodic box when saving

minimization
#############

minimization_max_move = 0.1
    [A] the largest move for one step

minimization_dynamic_dt = 0
    whether to use the dynamic dt stratege
    
minimization_dt_increasing_rate = 1.01
    the rate to increase if the energy of last step is bigger than this step

minimization_dt_decreasing_rate = 0.01
    the rate to decrease if the energy of last step is less than this step

minimization_momentum_keep = 0
    the rate to keep the momentum of last step

NVE
#####
nve_velocity_max
    [ A/(20.455 fs)] the largest velocity for one step for the NVE iteration

    If not set, the velocity will not be restricted.

neighbor list
##############
cutoff = 10.0
    [A] the cutoff radius for the non-bonded forces

skin = 2.0
    [A] the extra cutoff radius to build the neighbor list

    ``cutoff`` + ``skin`` is the cutoff radius to build the neighbor list.

neighbor_list_refresh_interval = 0
    the number of interval steps between the refreshing the neighbor list

    If a non-positive value is given, the neighbor list will be refreshed if there is an atom moved farther than 0.5 * ``neighbor_list_skin_permit`` * ``skin`` from its position on the last refreshed step.

neighbor_list_skin_permit = 0.5
    the parameter for the neighbor list automatical refreshing

neighbor_list_max_atom_in_grid_numbers = 64
    the number of max atoms in a neighbor list cell

neighbor_list_max_neighbor_numbers = 800
    the number of max neighbor atoms of one atom

thermostat
###########

target_temperature = 300.0
    [K] the target temperature for simulation

thermostat
    can be langevin(deprecated), middle_langevin, nose_hoover_chain, andersen_thermostat, berendsen_thermostat

middle_langevin
~~~~~~~~~~~~~~~~

middle_langevin_seed
    the random seed
    
    If not set, a random number will be used.
    
middle_langevin_gamma = 1
    [ps^-1] the collision frequency, also known as the friction factor

middle_langevin_velocity_max
    [ A/(20.455 fs)] the largest velocity for one step for the iteration

    If not set, the velocity will not be restricted.


nose_hoover_chain
~~~~~~~~~~~~~~~~~~~~~

nose_hoover_chain_length = 10
    the length of the nose hoover chain

nose_hoover_chain_tau = 1.0
    [ps] time constant for nose hoover chain

nose_hoover_chain_velocity_max
    [ A/(20.455 fs)] the largest velocity for one step for the iteration

    If not set, the velocity will not be restricted.

nose_hoover_chain_restart_input
    the input file of the restart file
    
    If not set, the coordinates and velocities of the nose hoover particles will all set to 0.

nose_hoover_chain_output
    the output file of the restart file
    
    If not set, the coordinates and velocities of the nose hoover particles will not be saved.

nose_hoover_chain_crd
    the coordinate trajectory file of the nose hoover particles
    
    If not set, no trajectory will be saved.

nose_hoover_chain_vel
    the velocity trajectory file of the nose hoover particles
    
    If not set, no trajectory will be saved.

andersen_thermostat
~~~~~~~~~~~~~~~~~~~~

andersen_thermostat_update_interval = 1000
    the number of interval steps

    ``update_interval`` * ``dt`` is like the time constant in other thermostats.

andersen_thermostat_seed
    the random seed

    If not set, a random number will be used.

andersen_thermostat_velocity_max
    [ A/(20.455 fs)] the largest velocity for one step for the iteration

    If not set, the velocity will not be restricted.
    
    .. ATTENTION::
        
        ``andersen_thermostat`` controls the temperature every ``update_interval`` step, and on other steps the NVE iterator are used, so if you want to control the max velocity, you need to set ``NVE_velocity_max``, too.

berendsen_thermostat
~~~~~~~~~~~~~~~~~~~~~

.. NOTE::

    1. Use ``NVE_velocity_max`` to restrict the velocity because the NVE iterator is used for berendsen thermostat.

    2. If the temperature is too far from ``target_temperature``, the system may crash, especially when the temperature = 0.

berendsen_thermostat_tau = 1.0 
    [ps] the time constant for the berendsen thermostat

berendsen_thermostat_stochastic_term = 0
    whether to add the stochastic term
    
    With the stochastic term, we can get a right canonical ensemble, also known as Bussi thermostat.

berendsen_thermostat_seed
    the random seed used for ``berendsen_thermostat_stochastic_term``

    If not set, a random number will be used.

barostat
##########

target_pressure = 1.0
    [bar] the target pressure for simulation

barostat
    can be andersen_barostat, monte_carlo_barostat, berendsen_barosta, andersen_barostat

berendsen_barostat
~~~~~~~~~~~~~~~~~~~~~

    berendsen_barostat_tau = 1.0
        [ps] the time constant for berendsen barostat

    berendsen_barostat_compressibility = 4.5e-5
        [bar^-1] the compressibility factor

    berendsen_barostat_update_interval = 10
        the numer of interval steps to update

    berendsen_barostat_stochastic_term = 0
        whether to add the stochastic term
        
        With the stochastic term, we can get a right NPT ensemble.

monte_carlo_barostat 
~~~~~~~~~~~~~~~~~~~~~~~

    monte_carlo_barostat_update_interval = 100
        the numer of interval steps to update

    monte_carlo_barostat_initial_ratio = 0.001 
        the factor of initial delta Vmax to the system volume

    monte_carlo_barostat_check_interval = 20
        the numer of interval steps to check the delta Vmax value

    monte_carlo_barostat_molecule_scale = 1
        whether to keep the rigidity of the molecule when controlling the pressure

    monte_carlo_barostat_accept_rate_low = 30
        the lowest acception percentage for monte carlo simulation

    monte_carlo_barostat_accept_rate_high = 40
        the highest acception percentage for monte carlo simulation

    monte_carlo_barostat_couple_dimension = xyz
        .. list-table::
            :widths: 10 50
            :header-rows: 1
            :stub-columns: 1
            
            * - choice
              - description
            * - XYZ
              - isotropic controlling, all directions will be scaled at the same time.
            * - NO
              - every direction will be scaled respectively.
            * - XY
              - the X and Y directions will be scaled together, and Z will be scaled itself
            * - XZ
              - the X and Z directions will be scaled together, and Y will be scaled itself
            * - YZ
              - the Y and Z directions will be scaled together, and X will be scaled itself 

andersen_barostat 
~~~~~~~~~~~~~~~~~~

andersen_barostat_tau = 1.0
    [ps] the time constant for berendsen barostat

andersen_barostat_compressibility = 4.5e-5
    [bar^-1] the compressibility factor

andersen_barostat_dV/dt = 0 
    [A^3/(20.455 fs)] the initial dV/dt

restraint
##############

positional restraint
~~~~~~~~~~~~~~~~~~~~

restrain_atom_id
    the atom index list file to restrain
    
    If not set, the ``restrain`` module will not be initialized.
    
restrain_weight = 20.0
    [kcal/mol/A^2] the force constant for the restraints

restrain_coordinate_in_file
    The reference coordinate file in SPONGE format

restrain_amber_rst7
    The reference coordinate file in AMBER format

    If neither of ``restrain_coordinate_in_file`` and ``restrain_amber_rst7`` is set, the reference coordinates will be the same as the initial coordinate of the MD task.

constraint
#############

constrain_mode
    can be SHAKE, simple_constrain
    
    If not set, no constraints will be used.

constrain_mass = 3.3
    [dalton] the mass parameter

    The bonds containing atoms whose mass is less than ``constrain_mass`` will be constrained.

constrain_angle = 0
    whether to constrain the angle

constrain_in_file
    the input file to give the parameters

settle_disable = 0
    whether to disable the SETTLE algorithm

simple_constrain
~~~~~~~~~~~~~~~~~~~~

simple_constrain_iteration_numbers = 25
    the number of iteration steps to constrain

simple_constrain_step_length = 1
    the step length to constrain

shake
~~~~~~~~~~~~~
shake_iteration_numbers = 25
    the number of iteration steps to constrain

shake_step_length = 1
    the step length to constrain

PME
~~~~~~~~~

PME_Direct_Tolerance = 1e-6
    the acceptable relative error for the direct energy

PME_fftx
    the number of the FFT grids in X direction

    If not set, a value will be guessed.

PME_ffty
    the number of the FFT grids in Y direction

    If not set, a value will be guessed.

PME_fftz
    the number of the FFT grids in Z direction

    If not set, a value will be guessed.

input files
############

default_in_file_prefix
    the default prefix for the input files
    
    If set, when a command named ``XXX_in_file`` is not set, the program will search the current working path to set ``XXX_in_file = {default_in_file_prefix}_XXX.txt``

coordinates and velocity
~~~~~~~~~~~~~~~~~~~~~~~~~~

coordinate_in_file
    the coordinate input file
    
    format:

    .. code-block::

        Line 1: Natom time
        From Line 2 to Line Natom+1: atom_x atom_y atom_z
        Line Natom+2: box_length_x box_length_y box_length_z box_angle_x box_angle_y box_angle_z

    .. code-block::

        Natom: the number of atoms
        time: [ps] optional, the time of the coordinate
        atom_*: [A] the coordinate of the atom
        box_length_*: [A] the box length
        box_angle_*: [degree] the box angle

velocity_in_file
    the velocity input file

    format:

    .. code-block::

        Line 1: Natom time
        From Line 2 to Line Natom+1: atom_velocity_x atom_velocity_y atom_velocity_z

    .. code-block::

        Natom: the number of atoms
        time: [ps] optional, the time of the coordinate
        atom_velocity_*: [A] the velocity of the atom

amber_rst7
    the coordinate and velocity input file in AMBER format

topology and force field parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mass_in_file
    the mass input file

    format:

    .. code-block::

        Line 1: Natom
        From Line 2 to Line Natom+1: atom_mass

    .. code-block::

        Natom: the number of atoms
        atom_mass: [dalton] the mass of the atom
    
charge_in_file
    the charge input file
    
    format:

    .. code-block::

        Line 1: Natom
        From Line 2 to Line Natom+1: atom_charge

    .. code-block::

        Natom: the number of atoms
        atom_mass: [1.0/18.2223 e] the charge of the atom

bond_in_file
    the bond input file

    format:

    .. math::

        E_{ij} = k_{ij} * (\Delta r_{ij} - b_{ij}) ^ 2

    .. code-block::

        Line 1: Nbond
        From Line 2 to Line Nbond+1: atom_i atom_j k_ij b_ij

    - Nbond: the number of bonds
    - atom_*: the atoms in the bond
    - :math:`k_{ij}`\: [kcal/mol/A^2] the force constant of the bond
    - :math:`b_{ij}`\: [A] the distance to keep the balance
    - :math:`\Delta r_{ij}`\: [A] the distance between the bonded atoms

angle_in_file
    the angle input file

    format:

    .. math::

        E_{ijk} = k_{ijk} * (\theta_{ijk} - b_{ijk}) ^ 2

    .. code-block::

        Line 1: Nangle
        From Line 2 to Line Nangle+1: atom_i atom_j atom_k k_ijk b_ijk

    - Nangle: the number of angles
    - atom_*: the atoms in the angle
    - :math:`k_{ijk}`\: [kcal/mol/rad^2] the force constant of the angle
    - :math:`b_{ijk}`\: [rad] the angle to keep the balance
    - :math:`\theta_{ijk}`\: [rad] the angle of the bonded atoms

dihedral_in_file
    the dihedral input file

    format:

    .. math::

        E_{ijkl} = V_{ijkl} * (1 + cos(n_{ijkl} * \phi_{ijkl} - b_{ijkl}))

    .. code-block::

        Line 1: Ndihedral
        From Line 2 to Line Ndihedral+1: atom_i atom_j atom_k atom_l n_ijkl V_ijkl b_ijkl

    - Ndihedral: the number of dihedrals
    - atom_*: the atoms in the dihedral
    - :math:`n_{ijkl}`\: the periodic number
    - :math:`V_{ijkl}`\: [kcal/mol] the minimum energy
    - :math:`b_{ijkl}`\: [rad] the dihedral angle to keep the balance
    - :math:`\phi_{ijkl}`\: [rad] the dihedral angle of the bonded atoms

LJ_in_file
    the LJ input file

    format:

    .. math::

        E_{ij} = A_{ij} * \Delta r_{ij}^{-12} - B_{ij} * \Delta r_{ij}^{-6}

    .. code-block::

        Line 1: Natom Ntype
        From Line 3 to Line Ntype+2: A_ij for j < i
        From Line Ntype+4 to Line 2Ntype+6: B_ij for j < i
        From Line 2Ntype+8 to Line 2Ntype+Natom+9: atom_type

    - Natom: the number of atoms
    - Ntype: the number of LJ types
    - :math:`A_{ij}`\: the A parameter of type i and j
    - :math:`B_{ij}`\: the B parameter of type i and j
    - :math:`\Delta r_{ij}`\: [A] the distance between the atoms

nb14_in_file
    the nb14 input file

    format:

    .. math::

        E_{ij} = k_{LJ} * (A_{ij} * \Delta r_{ij}^{-12} - B_{ij} * \Delta r_{ij}^{-6}) + k_{ee} * q_i * q_j * \Delta r_{ij} ^ {-1}

    .. code-block::

        Line 1: Nbond
        From Line 2 to Line Nbond+1: atom_i atom_j k_LJ  k_ee

    - Nbond: the number of nb14
    - atom_*: the atoms in the nb14
    - :math:`k_{LJ}`\: the scale factor for LJ
    - :math:`k_{ee}`\: the scale factor for electrostatic energy
    - :math:`\Delta r_{ij}`\: [A] the distance between the atoms

exclude_in_file
    the exclude input file

    format:

    .. code-block::

        Line 1: Natom Nexclude
        From Line 2 to Line Natom+1: Ni excluded_atoms_i

    - Natom: the number of atoms
    - Nexclude: the total number of excluded atoms
    - Ni: the number of excluded atoms for atom i
    - excluded_atoms_i: the excluded atoms for atom i

residue_in_file
    the residue input file

    format:

    .. code-block::

        Line 1: Nresidue
        From Line 2 to Line Ndihedral+1: Natom_i

    - Nresidue: the number of residues
    - Natom_i: the number of atoms for the i-th residue

amber_parm7
    the topology and force field parameters input file in AMBER format

standard output files
#########################

mdout = mdout.txt
    the output file for energies

mdinfo = mdinfo.txt
    the output file for parameters

rst = restart
    the output file for restart


box = mdbox.txt
    the output file for the box trajectory

crd = mdcrd.dat
    the output file for the coordinate trajectory

frc 
    the output file for the forces in every frame
    
    If not set, the forces will not be saved.

vel 
    the output file for the velocities in every frame
    
    If not set, the velocities will not be saved.

others
###########

end_pause = 0
    whether to ask press any key to end the program

device = 0
    the GPU device to use