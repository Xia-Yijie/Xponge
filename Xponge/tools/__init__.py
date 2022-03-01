

def dat2nc(args):
    from struct import unpack
    from netCDF4 import Dataset
    import numpy as np
    import os
    from itertools import islice
    data=Dataset(args.nc,"w",format="NETCDF3_64BIT_OFFSET")
    data.createDimension("frame",0)
    data.createDimension("spatial",3)
    data.createDimension("atom",args.n)
    data.createDimension("cell_spatial",3)
    data.createDimension("label",5)
    data.createDimension("cell_angular",3)

    data.title=args.title
    data.application="SPONGE"
    data.program="SPONGE"
    data.programVersion="1.2"
    data.Conventions="AMBER"
    data.ConventionVersion="1.0"

    data.createVariable("time","f",("frame"))
    data.variables["time"].units="picosecond"
    data.createVariable("spatial","S1",("spatial"))
    data.createVariable("coordinates","f",("frame","atom","spatial"))
    data.variables["coordinates"].units="angstrom"
    data.createVariable("cell_spatial","S1",("cell_spatial"))
    data.createVariable("cell_angular","S1",("cell_angular","label"))
    data.createVariable("cell_lengths","f8",("frame","cell_spatial"))
    data.variables["cell_lengths"].units="angstrom"
    data.createVariable("cell_angles","f8",("frame","cell_angular"))
    data.variables["cell_angles"].units="degree"

    data.variables["cell_angular"][:]=[[b'a',b'l',b'p',b'h',b'a'],
                                   [b'b',b'e',b't',b'a',b' '],
                                   [b'g',b'a',b'm',b'm',b'a']]
    data.variables["cell_spatial"][:]=[b'a',b'b',b'c']
    data.variables["spatial"][:]=[b'x',b'y',b'z']


    if not args.frame:
        size = os.path.getsize(args.x)
        args.frame = size // 12 // args.n
        print("frame = %d"%args.frame)

    if args.box:
        with open(args.box) as fr:
            box = np.loadtxt(islice(fr,args.frame), dtype = np.float32)
        box = box.reshape((-1,6))
        data.variables["cell_lengths"][:] = box[:,:3]
        data.variables["cell_angles"][:] = box[:,3:]
        box = None
    else:
        data.variables["cell_lengths"][:] = np.array(args.xyz) 
        data.variables["cell_angles"][:]=[90,90,90]
    
    data.variables["time"][:]= np.arange(args.frame, dtype = np.float32) * args.dt

    data.variables["coordinates"][:] = np.fromfile(args.x, dtype=np.float32)[:args.frame*args.n*3].reshape((args.frame, args.n, 3))
    data.close()

def nc2rst7(args):
   from netCDF4 import Dataset

   data = Dataset(args.nc)

   crds = data.variables["coordinates"][:]
   crds = list(crds.reshape(-1,3))
   vels = data.variables["velocities"][:]
   vels = list(vels.reshape(-1,3))

   cell_lengths = data.variables["cell_lengths"][:]
   cell_angles = data.variables["cell_angles"][:]
   title = data.title
   time = data.variables["time"][0]
   atom_numbers = data.dimensions["atom"].size

   towrite="%s\n%d %f\n"%(title, atom_numbers, time)
   for i,ci in enumerate(crds):
     towrite += " ".join(map(lambda x: "%12.7f"%x,ci))
     if i%2 == 1:
        towrite += "\n"
   if i%2 == 0:
     towrite += "\n"
    
   for i,ci in enumerate(vels):
     towrite += " ".join(map(lambda x: "%12.7f"%x,ci))
     if i%2 == 1:
        towrite += "\n"
   if i%2 == 0:
     towrite += "\n"

   towrite+= " ".join(map(lambda x: "%12.7f"%x,cell_lengths))
   towrite+= " ".join(map(lambda x: "%12.7f"%x,cell_angles))

   fw=open(args.rst7,"w")
   fw.write(towrite)
   fw.close()

def maskgen(args):
  import os

  s = input("Please Enter Your Selection Mask:\n")

  p = args.p.split(os.path.sep)
  p = "/".join(p)

  c=""
  if args.c:
    c = args.c.split(os.path.sep)
    c = "/".join(c)
    c = "mol addfile " + c

  temp_write = """set f [open "{0}" "w"]
mol new {1}
{2}
atomselect top "{3}"
puts $f [atomselect0 list]
close $f
quit
""".format(args.o, p, c, s)

  temp = open("maskgen_temp_tcl_file","w")
  temp.write(temp_write)
  temp.close()

  os.system("{0} -dispdev none -e maskgen_temp_tcl_file".format(args.vmd))
  os.remove("maskgen_temp_tcl_file")

def exgen(args):

  partners = [set([]) for i in range(args.n)]

  def exclude_2_atoms(words):
    i, j = int(words[0]), int(words[1])
    partners[i].add(j)
    partners[j].add(i)

  def exclude_3_atoms(words):
    i, j, k = int(words[0]), int(words[1]), int(words[2])
    partners[i].add(k)
    partners[k].add(i)


  def exclude_4_atoms(words):
    i, j, k, l= int(words[0]), int(words[1]), int(words[2]), int(words[3])
    partners[i].add(l)
    partners[l].add(i)


  if args.bond:
    for bond in args.bond:
        with open(bond) as f:
            f.readline()
            for line in f:
                words = line.split()
                exclude_2_atoms(words)

  if args.angle:
    for angle in args.angle:
        with open(angle) as f:
            f.readline()
            for line in f:
                words = line.split()
                exclude_3_atoms(words)

  if args.dihedral:
    for dihedral in args.dihedral:
        with open(dihedral) as f:
            f.readline()
            for line in f:
                words = line.split()
                exclude_4_atoms(words)

  if args.virtual:
    for virtual in args.virtual:
        with open(virtual) as f:
            for line in f:
                words = line.split()
                t = int(words[0])
                if t == 0:
                    exclude_2_atoms(words[1:])
                elif t == 1:
                    exclude_3_atoms(words[1:])
                elif t in (2,3):
                    exclude_4_atoms(words[1:])
                else:
                    raise Exception("virtual atom type wrong: are you sure this is a SPONGE virtual atom file?")

  if args.exclude:
    for exclude in args.exclude:
        with open(exclude) as f:
            f.readline()
            count = 0
            for line in f:
                words = line.split()
                t = set(words[1:])
                partners[count] = partners[count].union(t)
                count += 1

  total = 0
  towrite = "{} {}\n"
  for i, p in enumerate(partners):
    p = list(filter(lambda x: x > i, p))
    towrite += "%d "%len(p)
    towrite += ("{} "*len(p)).format(*p) + "\n"
    total += len(p)
  towrite = towrite.format(args.n, total)

  f = open(args.o, "w")
  f.write(towrite)
  f.close()

def dat1frame(args):
    f = open(args.box)
    for i, BOX in enumerate(f):
        if i == args.frame:
            break
    f.close()
    import numpy as np
    crds = np.fromfile(args.dat, np.float32, count = args.n * 3, offset = args.n * 12 * args.frame).reshape(-1,3)
    with open(args.o, "w") as f:
        f.write("%d\n"%args.n + "".join(list(map(lambda crd: "%12.7f %12.7f %12.7f\n"%(crd[0], crd[1], crd[2]), crds))) + BOX)

def gro2crd(args):
    with open(args.i) as f:
        f.readline()
        towrite = f.readline().strip() + "\n"
        t = int(towrite.strip())
        for i in range(t):
            towrite += " ".join([ "%f"%(float(s) * 10) for s in f.readline().split()[3:]]) + "\n"
        t = f.readline()
        assert len(t.split()) == 3, "orthogonal box needed"
        towrite += " ".join([ "%f"%(float(s) * 10) for s in t.split()[:3]]) + " 90 90 90"  
    with open(args.o, "w") as fo:
        fo.write(towrite)
        #print(towrite)

def crd2rst7(args):
    towrite = args.title + "\n"
    atom_numbers = 0
    start_time = None
    crds = []
    vels = []
    box = [0, 0, 0, 90, 90, 90]
    with open(args.crd) as f:
        t = f.readline().split()
        atom_numbers = int(t[0])
        if len(t) > 1:
            start_time = float(t[1])
        for i in range(atom_numbers):
            t = list(map(float, f.readline().split()))
            crds.append(t)
        box = [float(i) for i in f.readline().split()]
        
    if args.vel:
        second_line = " %d %f\n"%(atom_numbers, start_time)
        with open(args.vel) as f:
            t = f.readline().split()
            for i in range(atom_numbers):
                t = list(map(float, f.readline().split()))
                vels.append(t)
    
    else:
        second_line = " %d\n"%atom_numbers
        

    towrite += second_line
    for i in range(0, atom_numbers):
        towrite += "%12.7f%12.7f%12.7f"%(crds[i][0], crds[i][1], crds[i][2])
        if i % 2 == 1:
            towrite += "\n"
    if atom_numbers % 2 == 1:
        towrite += "\n"

    if args.vel:
      for i in range(0, atom_numbers):
        towrite += "%12.7f%12.7f%12.7f"%(vels[i][0], vels[i][1], vels[i][2])
        if i % 2 == 1:
            towrite += "\n"
    if atom_numbers % 2 == 1:
        towrite += "\n"
    towrite += "%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f\n"%(box[0], box[1], box[2], box[3], box[4], box[5])

    with open(args.rst7, "w") as f:
        f.write(towrite)

def mol2opt(args):
    import Xponge
    import Xponge.forcefield.AMBER.gaff
    import os

    t = Xponge.assign.Get_Assignment_From_Mol2(args.i)
    t.Save_As_Mol2("%s.mol2"%args.temp)
    t.Determine_Atom_Type("GAFF")
    temp = t.To_ResidueType(args.temp)
    Xponge.BUILD.Save_Mol2(temp)
    Xponge.forcefield.AMBER.gaff.parmchk2_gaff("%s.mol2"%args.temp, "%s.frcmod"%args.temp)
    Xponge.GlobalSetting.boxspace = 10
    Xponge.BUILD.Save_SPONGE_Input(temp)
    if args.nomin:
        return
    os.system("%s -mode minimization -rst %s -mdinfo %s.mdinfo -mdout %s.mdout -minimization_dynamic_dt 1 -default_in_file_prefix %s -step_limit %d -write_mdout_interval %d -write_restart_file_interval %d -write_information_interval 0 "%(args.sponge, args.temp, args.temp, args.temp, args.temp, args.step1, args.step1, args.step1))
    os.system("%s -mode minimization -dt 1e-3 -rst %s -mdinfo %s.mdinfo -mdout %s.mdout -default_in_file_prefix %s -step_limit %d -write_mdout_interval %d -write_restart_file_interval %d -write_information_interval 0 -coordinate_in_file %s_coordinate.txt"%(args.sponge, args.temp, args.temp, args.temp, args.temp, args.step2, args.step2, args.step2, args.temp))
    with open("%s_coordinate.txt"%args.temp) as f:
        f.readline()
        tt = f.read()
        tt = map(float, tt.split()[:-6])
    for i, xyzi in enumerate(tt):
        t.coordinate[i // 3][i % 3] = xyzi
    t.Save_As_Mol2(args.o)

def trr2dat(args):
    import numpy as np
    from MDAnalysis.lib.formats.libmdaxdr import TRRFile
    box = []
    crd = []
    with TRRFile(args.i) as f:
        for frame in f:
            crd.append(frame.x)
            box.append([frame.box[0][0] * 10, frame.box[1][1] * 10, frame.box[2][2] * 10, 90, 90, 90])
    box = np.array(box, dtype=np.float32)
    crd = np.array(crd, dtype=np.float32) * 10
    np.savetxt("{}.box".format(args.o), box)
    crd.tofile("{}.dat".format(args.o))
