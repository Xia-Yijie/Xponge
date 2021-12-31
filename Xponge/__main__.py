from . import tools
import argparse

parser = argparse.ArgumentParser(prog = "Xponge")

subparsers = parser.add_subparsers(help = "subcommands", description = "Tools for SPONGE")

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

nc2rst7 = subparsers.add_parser("nc2rst7", help = "convert a rst7 file from .nc to .rst7")
nc2rst7.add_argument("-rst7",required=True, help="the name of the .rst7 restart file.")
nc2rst7.add_argument("-nc", required=True, help="the name of the .nc restart file.")
nc2rst7.set_defaults(func = tools.nc2rst7)

maskgen = subparsers.add_parser("maskgen", help = "use VMD to generate a file to record the atom indexes of the corresponding mask")
maskgen.add_argument("-p",required=True, help="the parm file")
maskgen.add_argument("-c", help="the restart file")
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

args = parser.parse_args()
try:
  args.func(args)
except:
  parser.print_help()
