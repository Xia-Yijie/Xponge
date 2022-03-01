def main():
    from . import tools
    import argparse

    parser = argparse.ArgumentParser(prog = "Xponge")

    subparsers = parser.add_subparsers(help = "subcommands", description = "Tools for SPONGE. Use Xponge XXX -h for the help of tool 'XXX'.")

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


    args = parser.parse_args()

    if hasattr(args, "func"):
      args.func(args)
    else:
      parser.print_help()
      
if __name__ == "__main__":
    main()
