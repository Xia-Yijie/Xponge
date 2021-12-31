from . import tools
import argparse

parser = argparse.ArgumentParser(prog = "Xponge Tools")

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

args = parser.parse_args()
args.func(args)
