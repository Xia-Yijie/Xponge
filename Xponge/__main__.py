def main():
    from . import tools
    import argparse

    parser = argparse.ArgumentParser(prog = "Xponge")

    subparsers = parser.add_subparsers(help = "subcommands", description = "Tools for SPONGE. Use Xponge XXX -h for the help of tool 'XXX'.")

    mytest = subparsers.add_parser("test", help = "test the basic function of Xponge")
    mytest.add_argument("-o", metavar = "test", default="test", help="the prefix for the output files")
    mytest.set_defaults(func = tools.test)
    
    dat2nc = subparsers.add_parser("dat2nc", help = "convert a traj file from .dat to .nc")

    dat2nc.add_argument("-n",required=True,type=int, help="the number of atoms in the traj file",metavar="atom_numbers")
    dat2nc.add_argument("-frame",type=int, help="the number of frames you want convert. If not set, a value will be guessed from the file size.",metavar="frame_numbers")
    dat2nc.add_argument("-x",required=True, help="the name of the .dat traj file.",metavar="dat_traj_file")
    group = dat2nc.add_mutually_exclusive_group(required=True)
    group.add_argument("-box", help="the name of the box traj file",metavar="box_traj_file")
    group.add_argument("-xyz", type=float, nargs=3, help="the three constant box length",metavar=("X","Y","Z"))
    dat2nc.add_argument("-nc", required=True, help="the name of the .nc traj file.",metavar="nc_traj_file")
    dat2nc.add_argument("-dt", type=float, default=1.0, help="the time difference between two frames in the unit of ps. 1.0 for default.")
    dat2nc.add_argument("--title",default="DEFAULT", help="the title of the traj.")
    dat2nc.set_defaults(func = tools.dat2nc)

    gro2crd = subparsers.add_parser("gro2crd", help = "convert a coordinate file from gromacs to SPONGE")
    gro2crd.add_argument("-i",required=True, help="the name of the .gro file.")
    gro2crd.add_argument("-o", required=True, help="the name of the SPONGE coordinate file.")
    gro2crd.set_defaults(func = tools.gro2crd)

    nc2rst7 = subparsers.add_parser("nc2rst7", help = "convert a rst7 file from .nc to .rst7")
    nc2rst7.add_argument("-rst7",required=True, help="the name of the .rst7 restart file.")
    nc2rst7.add_argument("-nc", required=True, help="the name of the .nc restart file.")
    nc2rst7.set_defaults(func = tools.nc2rst7)

    maskgen = subparsers.add_parser("maskgen", help = "use VMD to generate a file to record the atom indexes of the corresponding mask")
    maskgen.add_argument("-p",required=True, help="the topology file")
    maskgen.add_argument("-c", help="the coordinate file")
    maskgen.add_argument("-o",required=True, help="the output file")
    maskgen.add_argument("--vmd", metavar = "vmd", default="vmd", help="the command to start vmd")
    maskgen.set_defaults(func = tools.maskgen)

    exgen = subparsers.add_parser("exgen", help = 'process bond-like, angle-like, dihedral-like files to get the atoms to exclude')
    exgen.add_argument('-n', type = int, required = True, help='the atom numbers')
    exgen.add_argument('-o', required = True, help='output exclude file name')
    exgen.add_argument('-b', '--bond', nargs='+', help='bond-like input files: skip the first line , and there are 2 atoms in the head of following lines')
    exgen.add_argument('-a', '--angle', nargs='+', help='angle-like input files: skip the first line , and there are 3 atoms in the head of following lines')
    exgen.add_argument('-d', '--dihedral', nargs='+', help='dihedral-like input files: skip the first line , and there are 4 atoms in the head of following lines')
    exgen.add_argument('-v', '--virtual', nargs='+', help='virtual-atom-like input files: the first number indicates the virtual type')
    exgen.add_argument('-e', '--exclude', nargs='+', help='exclude-like input files: add the information of another exclude file')
    exgen.set_defaults(func = tools.exgen)

    dat1frame = subparsers.add_parser("dat1frame", help = 'extract 1 frame from .dat traj file')
    dat1frame.add_argument('-n', type = int, required = True, help='the atom numbers')
    dat1frame.add_argument('-frame', type = int, required = True, help='the frame to extract')
    dat1frame.add_argument('-o', required = True, help='output file name')
    dat1frame.add_argument('-box', required = True, help='box traj file name')
    dat1frame.add_argument('-dat', required = True, help='dat traj file name')
    dat1frame.set_defaults(func = tools.dat1frame)

    crd2rst7 = subparsers.add_parser("crd2rst7", help = 'convert a coordinate file from SPONGE .txt to .rst7')
    crd2rst7.add_argument("-title", default = "converted by Xponge", help='the title of the rst7 file')
    crd2rst7.add_argument("-crd", required = True, help = 'the SPONGE coordinate file')
    crd2rst7.add_argument("-vel", help = 'the SPONGE velocity file')
    crd2rst7.add_argument("-rst7", required = True, help = 'the output rst7 file')
    crd2rst7.set_defaults(func = tools.crd2rst7)

    trr2dat = subparsers.add_parser("trr2dat", help = 'convert a trajectory file from GROMACS .trr to SPONGE .dat and .box')
    trr2dat.add_argument("-i", required = True, help= 'the trr file name')
    trr2dat.add_argument("-o", required = True, help ='the output file name prefix')
    trr2dat.set_defaults(func = tools.trr2dat)

    name2name = subparsers.add_parser("name2name", help = "change the atom names of a residue from one file to another file")
    name2name.add_argument("-fformat", "-from_format", dest = "from_format", choices = ["mol2", "pdb", "gaff_mol2"], required = True, help = "the format of the file which is needed to change from")
    name2name.add_argument("-ffile", "-from_file", dest = "from_file", required = True, help = "the name of the file which is needed to change from")
    name2name.add_argument("-fres", "-from_residue", dest = "from_residue", default = "", help = "the residue name in ffile if fformat == pdb")

    name2name.add_argument("-tformat", "-to_format", dest = "to_format", choices = ["mol2", "pdb", "gaff_mol2"], required = True, help = "the format of the file which is needed to change to")
    name2name.add_argument("-tfile", "-to_file", dest = "to_file", required = True, help = "the name of the file which is needed to change to")
    name2name.add_argument("-tres", "-to_residue", dest = "to_residue", default = "", help = "the residue name in tfile if tformat == pdb")

    name2name.add_argument("-oformat", "-out_format", dest = "out_format", choices = ["mol2", "pdb", "mcs_pdb"], required = True, help = "the format of the output file")
    name2name.add_argument("-ofile", "-out_file", dest = "out_file", required = True, help = "the name of the output file")
    name2name.add_argument("-ores", "-out_residue", dest = "out_residue", default = "ASN", help = "the name of the output residue")
    name2name.set_defaults(func = tools.name2name)

    mol2opt = subparsers.add_parser("mol2opt", help = 'optimize a mol2 file by using minimization-mode SPONGE and GAFF')
    mol2opt.add_argument("-sponge", default = "SPONGE", help = "SPONGE program command")
    mol2opt.add_argument("-i", required = True, help = "the input mol2 file")
    mol2opt.add_argument("-temp", default = "TMP", help = "the temporary file name prefix")
    #mol2opt.add_argument("-deltemp", action = "store_true", help = "delete the temporary files" )
    mol2opt_ = mol2opt.add_mutually_exclusive_group(required = True)
    mol2opt_.add_argument("-o", help = "the output mol2 file")
    mol2opt_.add_argument("-nomin", action = "store_true", help = "only build the mol2 and do not minimization")
    mol2opt.add_argument("-step1", default = 2500, type = int, help = "the number of steps for first step minimization")
    mol2opt.add_argument("-step2", default = 2500, type = int, help = "the number of steps for first step minimization")
    mol2opt.set_defaults(func = tools.mol2opt)

    mol2hfe = subparsers.add_parser("mol2hfe", help = 'calculate the hydration free energy of a small molecule using SPONGE')
    mol2hfe.add_argument("-lambda_numbers", type = int, default = 20, help = "the number of lambda groups - 1, default 20 for 0, 0.05, 0.10, 0.15..., 1.0")
    mol2hfe.add_argument("-prebalance_step", type = int, default = 50000, help = "the numbers of step to do the prebalance simulation")
    mol2hfe.add_argument("-balance_step", type = int, default = 500000, help = "the numbers of step to do the balance simulation")
    mol2hfe.add_argument("-barostat", default = "andersen_barostat", help = "the barostat used for the simulation")
    mol2hfe.add_argument("-thermostat", default = "middle_langevin", help = "the thermostat used for the simulation")
    mol2hfe.add_argument("-method", default = "TI", choices = ["TI"], help = "the method to calculate the free energy")
    mol2hfe.add_argument("-sponge", default = "SPONGE", help = "SPONGE program command")
    mol2hfe.add_argument("-sponge_ti", default = "SPONGE_TI", help = "SPONGE_TI program command")
    mol2hfe.add_argument("-temp", default = "TMP", help = "the temporary file name prefix")
    mol2hfe_ = mol2hfe.add_mutually_exclusive_group(required = True)
    mol2hfe_.add_argument("-assignment", help = "input by an Xponge assignment mol2 file")
    mol2hfe_.add_argument("-name", help = "input by the name of the molecule")
    mol2hfe_.add_argument("-smiles", help = "input by the smiles of the molecule")
    mol2hfe_.add_argument("-residuetype", help = "input by an Xponge ResidueType mol2 file")
    mol2hfe.set_defaults(func = tools.mol2hfe)

    mol2rfe = subparsers.add_parser("mol2rfe", help = 'calculate the relative free energy of a small molecule using SPONGE')
    mol2rfe.add_argument("-do", metavar = "todo", nargs = "*", action = "append", help = "the things need to do, should be one or more of 'build', 'min', 'prebalance', 'balance', 'analysis'", choices = ["build", "min", "prebalance", "balance", "analysis"])


    mol2rfe.add_argument("-pdb", required = True, help = "the initial conformation given by the pdb file")
    mol2rfe.add_argument("-r2", "-residuetype2", required = True, help = "molecule mutated to by an Xponge ResidueType mol2 file")
    mol2rfe.add_argument("-r1", "-residuetype1", required = True, help = "molecule mutated from by an Xponge ResidueType mol2 file")
    mol2rfe.add_argument("-r0", "-residuetype0", nargs = "*", default = [], help = "small molecules that do not mutate")
    mol2rfe.add_argument("-ri", "-residue_index", type = int, metavar = 0, default = 0, help = "the residue index of the molecule to mutate")
    mol2rfe.add_argument("-nl", "-lambda_numbers", metavar = 20, type = int, default = 20, help = "the number of lambda groups - 1, default 20 for 0, 0.05, 0.10, 0.15..., 1.0")

    mol2rfe.add_argument("-dohmr", "-do_hydrogen_mass_repartition", action = "store_true", help = "use the hydrogen mass repartition method")
    mol2rfe.add_argument("-ff", "-forcefield", help = "Use this force field file instead of the default ff14SB and gaff")
    mol2rfe.add_argument("-pi", "-prebalance_mdin", help = "Use this prebalance mdin file instead of the default one")
    mol2rfe.add_argument("-bi", "-balance_mdin", help = "Use this balance mdin file instead of the default one")
    mol2rfe.add_argument("-ai", "-analysis_mdin", help = "Use this analysis mdin file instead of the default one")

    mol2rfe.add_argument("-method", default = "TI", choices = ["TI"], help = "the method to calculate the free energy")
    mol2rfe.add_argument("-sponge", default = "SPONGE", help = "SPONGE program command")
    mol2rfe.add_argument("-sponge_ti", default = "SPONGE_TI", help = "SPONGE_TI program command")
    mol2rfe.add_argument("-sponge_fep", default = "SPONGE_FEP", help = "SPONGE_FEP program command")
    mol2rfe.add_argument("-temp", default = "TMP", metavar = "TMP", help = "the temporary file name prefix")
    
    mol2rfe.add_argument("-tmcs", default = 10, type = int, metavar = "10", help = "the timeout parameter for max common structure in unit of second")
    mol2rfe.add_argument("-dt", default = 2e-3, type = float, metavar = "dt", help = "the dt used for simulation when mdin is not provided")
    mol2rfe.add_argument("-m1steps", type = int, nargs = 5, help = "the first-stage minimization steps for the 0th lambda. Default 5000 for each minimization simulation. There are 5 minimization simulations.", default = [5000, 5000, 5000, 5000, 5000])
    mol2rfe.add_argument("-m2steps", type = int, nargs = 5, help = "the second-stage minimization steps for the 0th lambda. Default 5000 for each minimization simulation. There are 5 minimization simulations.", default = [5000, 5000, 5000, 5000, 5000])
    mol2rfe.add_argument("-msteps", type = int, nargs = 2, help = "the minimization steps for all the lambda. Default 5000 for each minimization simulation. There are 2 minimization simulations.", default = [5000, 5000])
    mol2rfe.add_argument("-pstep", "-prebalance_step", dest = "prebalance_step", default = 50000, type = int, metavar = "prebalance_step", help = "the prebalance step used for simulation when mdin is not provided")
    mol2rfe.add_argument("-bstep", "-balance_step", dest = "balance_step", default = 500000, type = int, metavar = "balance_step", help = "the balance step used for simulation when mdin is not provided")
    mol2rfe.add_argument("-thermostat", default = "middle_langevin", metavar = "middle_langevin", help = "the thermostat used for simulation when mdin is not provided")
    mol2rfe.add_argument("-barostat", default = "andersen_barostat", metavar = "andersen_barostat", help = "the barostat used for simulation when mdin is not provided")

    mol2rfe.set_defaults(func = tools.mol2rfe)

    crd2pdb = subparsers.add_parser("crd2pdb", help = 'add the coordinary from SPONGE file to pdb')
    crd2pdb.add_argument("-pdb", required = True, help = "the input pdb file")
    crd2pdb.add_argument("-crd", required = True, help = "the SPONGE coordinary file")
    crd2pdb.add_argument("-o", required = True, help = "the output pdb file")
    crd2pdb.set_defaults(func = tools.crd2pdb)

    args = parser.parse_args()

    if hasattr(args, "func"):
      args.func(args)
    else:
      parser.print_help()
      
if __name__ == "__main__":
    main()
