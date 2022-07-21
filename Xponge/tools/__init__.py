"""
This **module** implements the terminal commands
"""
from ..helper import source, GlobalSetting


def _basic_test(args):
    """

    :param args:
    :return:
    """
    source("..")
    source("..forcefield.amber.ff14sb")
    source("..forcefield.amber.tip3p")
    source("__main__")

    t = ACE + ALA * 10 + NME

    for i in range(1, len(t.residues) - 1):
        head = t.residues[i - 1]
        res = t.residues[i]
        tail = t.residues[i + 1]
        Impose_Dihedral(t, head.C, res.N, res.CA, res.C, -3.1415926 / 3)
        Impose_Dihedral(t, res.N, res.CA, res.C, tail.N, -3.1415926 / 3)

    Save_Mol2(t, f"{args.o}.mol2")
    c = int(round(t.charge))
    Add_Solvent_Box(t, WAT, 10)
    Solvent_Replace(t, lambda res: res.type.name == "WAT", {CL: 30 + c, K: 30})
    t.residues.sort(key=lambda residue: {"CL": 2, "K": 1, "WAT": 3}.get(residue.type.name, 0))
    Save_PDB(t, f"{args.o}.pdb")
    Save_SPONGE_Input(t, f"{args.o}")


def _assign_test(args):
    """

    :param args:
    :return:
    """
    source("..")
    source("..forcefield.amber.gaff")
    source("__main__")

    t = assign.Assign()
    t.add_atom("O", 0, 0, 0)
    t.add_atom("H", 1, 0, 0)
    t.add_atom("H", 0, 1, 0)
    t.add_bond(0, 1, 1)
    t.add_bond(0, 2, 1)
    t.determine_ring_and_bond_type()
    t.determine_atom_type("gaff")
    equal_atoms = t.Determine_Equal_Atoms()
    t.calculate_charge("resp", opt=True, extra_equivalence=equal_atoms)
    Save_PDB(t, f"{args.o}.pdb")
    Save_Mol2(t, f"{args.o}.mol2")
    wat = t.to_residuetype("WAT")
    Save_PDB(wat, f"{args.o}_Residue.pdb")
    Save_Mol2(wat, f"{args.o}_Residue.pdb")
    Save_SPONGE_Input(WAT, f"{args.o}")


def _charmm27_test(args):
    """

    :param args:
    :return:
    """
    source("..")
    source("..forcefield.charmm27.protein")
    source("__main__")

    t = ACE + ALA * 10 + NME
    Save_SPONGE_Input(t, f"{args.o}")


def test(args):
    """

    :param args:
    :return:
    """
    GlobalSetting.verbose = args.verbose
    if not args.do:
        args.do = [["base"]]
    args.do = args.do[0]
    if "base" in args.do:
        _basic_test(args)

    if "assign" in args.do:
        _assign_test(args)

    if "charmm27" in args.do:
        _charmm27_test(args)


def maskgen(args):
    """

    :param args:
    :return:
    """
    import os

    s = input("Please Enter Your Selection Mask:\n")

    p = args.p.split(os.path.sep)
    p = "/".join(p)

    c = ""
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

    with open("maskgen_temp_tcl_file", "w") as temp:
        temp.write(temp_write)

    os.system("{0} -dispdev none -e maskgen_temp_tcl_file".format(args.vmd))
    os.remove("maskgen_temp_tcl_file")


def exgen(args):
    """

    :param args:
    :return:
    """
    partners = [set([]) for i in range(args.n)]

    def exclude_2_atoms(words):
        i, j = int(words[0]), int(words[1])
        partners[i].add(j)
        partners[j].add(i)

    def exclude_3_atoms(words):
        i, k = int(words[0]), int(words[2])
        partners[i].add(k)
        partners[k].add(i)

    def exclude_4_atoms(words):
        i, l = int(words[0]), int(words[3])
        partners[i].add(l)
        partners[l].add(i)

    for bond in args.bond:
        with open(bond) as f:
            f.readline()
            for line in f:
                words = line.split()
                exclude_2_atoms(words)

    for angle in args.angle:
        with open(angle) as f:
            f.readline()
            for line in f:
                words = line.split()
                exclude_3_atoms(words)
    for dihedral in args.dihedral:
        with open(dihedral) as f:
            f.readline()
            for line in f:
                words = line.split()
                exclude_4_atoms(words)

    for virtual in args.virtual:
        with open(virtual) as f:
            for line in f:
                words = line.split()
                t = int(words[0])
                if t == 0:
                    exclude_2_atoms(words[1:])
                elif t == 1:
                    exclude_3_atoms(words[1:])
                elif t in (2, 3):
                    exclude_4_atoms(words[1:])
                else:
                    raise Exception("virtual atom type wrong: are you sure this is a SPONGE virtual atom file?")

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
        newp = []
        for pi in p:
            if pi > i:
                newp.append(pi)
        towrite += "%d " % len(newp)
        towrite += ("{} " * len(newp)).format(*newp) + "\n"
        total += len(newp)
        towrite = towrite.format(args.n, total)

    with open(args.o, "w") as f:
        f.write(towrite)


def name2name(args):
    """

    :param args:
    :return:
    """
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
    source("..")
    rdktool = source("..helper.rdkit")

    if args.to_format == "mol2":
        to_ = assign.Get_Assignment_From_Mol2(args.to_file)
    elif args.to_format == "gaff_mol2":
        source("..forcefield.amber.gaff")
        to_ = load_mol2(args.to_file).residues[0]
        to_ = assign.Get_Assignment_From_ResidueType(to_)
    elif args.to_format == "pdb":
        to_ = assign.Get_Assignment_From_PDB(args.to_file, determine_bond_order=False,
                                                    only_residue=args.to_residue)

    if args.from_format == "mol2":
        from_ = assign.Get_Assignment_From_Mol2(args.from_file)
    elif args.from_format == "gaff_mol2":
        source("..forcefield.amber.gaff")
        from_ = load_mol2(args.from_file).residues[0]
        from_ = assign.Get_Assignment_From_ResidueType(from_)
    elif args.from_format == "pdb":
        from_ = assign.Get_Assignment_From_PDB(args.from_file, determine_bond_order=False,
                                               only_residue=args.from_residue)

    rdmol_a = rdktool.Assign2RDKitMol(to_, True)
    rdmol_b = rdktool.Assign2RDKitMol(from_, True)

    result = rdFMCS.FindMCS([rdmol_a, rdmol_b], completeRingsOnly=True, timeout=args.tmcs)
    rdmol_mcs = Chem.MolFromSmarts(result.smartsString)

    match_a = rdmol_a.GetSubstructMatch(rdmol_mcs)
    match_b = rdmol_b.GetSubstructMatch(rdmol_mcs)
    matchmap = {from_.names[match_b[j]]: to_.names[match_a[j]] for j in range(len(match_a))}
    from_.names = [matchmap.get(name, name) for name in from_.names]
    from_.name = args.out_residue

    if args.out_format == "mol2":
        from_.Save_As_Mol2(args.out_file)
    elif args.out_format == "pdb":
        from_.Save_As_PDB(args.out_file)
    elif args.out_format == "mcs_pdb":
        towrite = towrite = "REMARK   Generated By Xponge (Max Common Structure)\n"
        for i, atom in enumerate(from_.atoms):
            if i in match_b:
                towrite += "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%17s%2s\n" % (i + 1, from_.names[i],
                                                                                     from_.name, " ", 1,
                                                                                     from_.coordinate[i][0],
                                                                                     from_.coordinate[i][1],
                                                                                     from_.coordinate[i][2], " ", atom)
        with open(args.out_file, "w") as f:
            f.write(towrite)


def _mol2rfe_build(args, merged_from, merged_to):
    """

    :param args:
    :param merged_from:
    :param merged_to:
    :return:
    """
    import os
    import Xponge
    import Xponge.forcefield.special.fep as FEP
    import Xponge.forcefield.special.min as MIN

    if "build" in args.do:
        print("\nBUILDING TOPOLOGY\n")
        FEP.Save_Soft_Core_LJ()

        for i in range(args.nl + 1):
            if os.path.exists("%d" % i):
                os.system("rm -rf %d" % i)
            os.mkdir("%d" % i)
            tt = FEP.Merge_Force_Field(merged_from, merged_to, i / args.nl)
            if i == 0:
                MIN.save_min_bonded_parameters()
            elif i == 1:
                MIN.do_not_save_min_bonded_parameters()
            Xponge.BUILD.Save_SPONGE_Input(tt, "%d/%s" % (i, args.temp))


def _mol2rfe_min(args):
    """

    :param args:
    :return:
    """
    source("..")

    if "min" in args.do:
        for i in range(args.nl + 1):
            if os.path.exists("%d/min" % i):
                os.system("rm -rf %d/min" % i)
            os.mkdir("%d/min" % i)
            if i != 0 and args.mlast:
                cif = "-coordinate_in_file {1}/min/{0}_coordinate.txt".format(args.temp, i - 1)
            else:
                cif = ""
                os.system(
                    "{0} -mode minimization -minimization_dynamic_dt 1 -default_in_file_prefix {2}/{1} -mdinfo {2}/min/{1}.mdinfo -mdout {2}/min/{1}.mdout -rst {2}/min/{1} -crd {2}/min/{1}.dat -box {2}/min/{1}.box -step_limit {5} {4} -mass_in_file 0/{1}_fake_mass.txt -lambda_lj {3}".format(
                        args.sponge, args.temp, i, i / args.nl, cif, args.m1steps[0]))
                cif = "-coordinate_in_file {1}/min/{0}_coordinate.txt".format(args.temp, i)
                os.system(
                    "{0} -mode minimization -dt 1e-7 -default_in_file_prefix {2}/{1} -mdinfo {2}/min/{1}.mdinfo -mdout {2}/min/{1}.mdout -rst {2}/min/{1} -crd {2}/min/{1}.dat -box {2}/min/{1}.box -step_limit {5}  {4} -mass_in_file 0/{1}_fake_mass.txt -lambda_lj {3}".format(
                        args.sponge, args.temp, i, i / args.nl, cif, args.m1steps[1]))
                os.system(
                    "{0} -mode minimization -dt 1e-6 -default_in_file_prefix {2}/{1} -mdinfo {2}/min/{1}.mdinfo -mdout {2}/min/{1}.mdout -rst {2}/min/{1} -crd {2}/min/{1}.dat -box {2}/min/{1}.box -step_limit {5}  {4} -mass_in_file 0/{1}_fake_mass.txt -lambda_lj {3}".format(
                        args.sponge, args.temp, i, i / args.nl, cif, args.m1steps[2]))
                os.system(
                    "{0} -mode minimization -dt 1e-5 -default_in_file_prefix {2}/{1} -mdinfo {2}/min/{1}.mdinfo -mdout {2}/min/{1}.mdout -rst {2}/min/{1} -crd {2}/min/{1}.dat -box {2}/min/{1}.box -step_limit {5}  {4} -mass_in_file 0/{1}_fake_mass.txt -lambda_lj {3}".format(
                        args.sponge, args.temp, i, i / args.nl, cif, args.m1steps[3]))
                os.system(
                    "{0} -mode minimization -dt 1e-4 -default_in_file_prefix {2}/{1} -mdinfo {2}/min/{1}.mdinfo -mdout {2}/min/{1}.mdout -rst {2}/min/{1} -crd {2}/min/{1}.dat -box {2}/min/{1}.box -step_limit {5}  {4} -mass_in_file 0/{1}_fake_mass.txt -lambda_lj {3}".format(
                        args.sponge, args.temp, i, i / args.nl, cif, args.m1steps[4]))
                os.system(
                    "{0} -mode minimization -minimization_dynamic_dt 1 -default_in_file_prefix {2}/{1} -mdinfo {2}/min/{1}.mdinfo -mdout {2}/min/{1}.mdout -rst {2}/min/{1} -crd {2}/min/{1}.dat -box {2}/min/{1}.box -lambda_lj {3} -step_limit {5}  {4}".format(
                        args.sponge, args.temp, i, i / args.nl, cif, args.m2steps[0]))
                os.system(
                    "{0} -mode minimization -dt 1e-7 -default_in_file_prefix {2}/{1} -mdinfo {2}/min/{1}.mdinfo -mdout {2}/min/{1}.mdout -rst {2}/min/{1} -crd {2}/min/{1}.dat -box {2}/min/{1}.box -lambda_lj {3} -step_limit {5} -constrain_mode SHAKE {4}".format(
                        args.sponge, args.temp, i, i / args.nl, cif, args.m2steps[1]))
                os.system(
                    "{0} -mode minimization -dt 1e-6 -default_in_file_prefix {2}/{1} -mdinfo {2}/min/{1}.mdinfo -mdout {2}/min/{1}.mdout -rst {2}/min/{1} -crd {2}/min/{1}.dat -box {2}/min/{1}.box -lambda_lj {3} -step_limit {5} -constrain_mode SHAKE  {4}".format(
                        args.sponge, args.temp, i, i / args.nl, cif, args.m2steps[2]))
                os.system(
                    "{0} -mode minimization -dt 1e-5 -default_in_file_prefix {2}/{1} -mdinfo {2}/min/{1}.mdinfo -mdout {2}/min/{1}.mdout -rst {2}/min/{1} -crd {2}/min/{1}.dat -box {2}/min/{1}.box -lambda_lj {3} -step_limit {5} -constrain_mode SHAKE  {4}".format(
                        args.sponge, args.temp, i, i / args.nl, cif, args.m2steps[3]))
                os.system(
                    "{0} -mode minimization -dt 1e-4 -default_in_file_prefix {2}/{1} -mdinfo {2}/min/{1}.mdinfo -mdout {2}/min/{1}.mdout -rst {2}/min/{1} -crd {2}/min/{1}.dat -box {2}/min/{1}.box -lambda_lj {3} -step_limit {5} -constrain_mode SHAKE  {4}".format(
                        args.sponge, args.temp, i, i / args.nl, cif, args.m2steps[4]))

            os.system(
                "{0} -mode minimization -minimization_dynamic_dt 1 -default_in_file_prefix {2}/{1} -mdinfo {2}/min/{1}.mdinfo -mdout {2}/min/{1}.mdout -rst {2}/min/{1} -crd {2}/min/{1}.dat -box {2}/min/{1}.box -lambda_lj {3} -step_limit {5} {4} constrain_mode SHAKE".format(
                    args.sponge, args.temp, i, i / args.nl, cif, args.msteps[0]))
            os.system(
                "{0} -mode minimization -dt 1e-3 -default_in_file_prefix {2}/{1} -mdinfo {2}/min/{1}.mdinfo -mdout {2}/min/{1}.mdout -rst {2}/min/{1} -crd {2}/min/{1}.dat -box {2}/min/{1}.box -lambda_lj {3} -step_limit {5} -constrain_mode SHAKE {4}".format(
                    args.sponge, args.temp, i, i / args.nl, cif, args.msteps[1]))


def _mol2rfe_prebalance(args):
    """

    :param args:
    :return:
    """
    source("..")

    if "prebalance" in args.do:
        for i in range(args.nl + 1):
            if os.path.exists("%d/prebalance" % i):
                os.system("rm -rf %d/prebalance" % i)
            os.mkdir("%d/prebalance" % i)
            if not args.pi:
                os.system(
                    "{0} -mode NPT -dt {7} -default_in_file_prefix {2}/{1} -mdinfo {2}/prebalance/{1}.mdinfo -mdout {2}/prebalance/{1}.mdout -rst {2}/prebalance/{1} -crd {2}/prebalance/{1}.dat -box {2}/prebalance/{1}.box -lambda_lj {3} -constrain_mode SHAKE -step_limit {4} -barostat {5} -thermostat {6} -coordinate_in_file {2}/min/{1}_coordinate.txt".format(
                        args.sponge, args.temp, i, i / args.nl, args.prebalance_step, args.barostat, args.thermostat,
                        args.dt))
            else:
                os.system(
                    "{0} -default_in_file_prefix {2}/{1} -mdinfo {2}/prebalance/{1}.mdinfo -mdout {2}/prebalance/{1}.mdout -rst {2}/prebalance/{1} -crd {2}/prebalance/{1}.dat -box {2}/prebalance/{1}.box -lambda_lj {3} -coordinate_in_file {2}/min/{1}_coordinate.txt -mdin {4}".format(
                        args.sponge, args.temp, i, i / args.nl, args.pi))


def _mol2rfe_balance(args):
    """

    :param args:
    :return:
    """
    source("..")

    if "balance" in args.do:
        for i in range(args.nl + 1):
            if os.path.exists("%d/balance" % i):
                os.system("rm -rf %d/balance" % i)
            os.mkdir("%d/balance" % i)
            if not args.bi:
                os.system(
                    "{0} -mode NPT -dt {7} -default_in_file_prefix {2}/{1} -mdinfo {2}/balance/{1}.mdinfo -mdout {2}/balance/{1}.mdout -rst {2}/balance/{1} -crd {2}/balance/{1}.dat -box {2}/balance/{1}.box -lambda_lj {3} -constrain_mode SHAKE -step_limit {4} -barostat {5} -thermostat {6} -coordinate_in_file {2}/prebalance/{1}_coordinate.txt -velocity_in_file {2}/prebalance/{1}_velocity.txt -write_information_interval 100 -write_mdout_interval 5000 -write_restart_file_interval {4}".format(
                        args.sponge, args.temp, i, i / args.nl, args.balance_step, args.barostat, args.thermostat,
                        args.dt))
            else:
                os.system(
                    "{0} -default_in_file_prefix {2}/{1} -mdinfo {2}/balance/{1}.mdinfo -mdout {2}/balance/{1}.mdout -rst {2}/balance/{1} -crd {2}/balance/{1}.dat -box {2}/balance/{1}.box -lambda_lj {3} -coordinate_in_file {2}/prebalance/{1}_coordinate.txt -velocity_in_file {2}/prebalance/{1}_velocity.txt -mdin {4}".format(
                        args.sponge, args.temp, i, i / args.nl, args.bi))


def _mol2rfe_analysis(args, merged_from):
    """

    :param args:
    :param merged_from:
    :return:
    """
    source("..")

    if "analysis" in args.do:
        with open("dh_dlambda.txt", "w") as f:
            f.write("")
        if args.method == "TI":
            for i in range(args.nl + 1):
                if os.path.exists("%d/ti" % i):
                    os.system("rm -rf %d/ti" % i)
                os.mkdir("%d/ti" % i)
                if not args.ai:
                    os.system(
                        "{0} -LJ_soft_core_in_file {2}/{1}_LJ_soft_core.txt -exclude_in_file {2}/{1}_exclude.txt -charge_in_file {2}/{1}_charge.txt -chargeA_in_file 0/{1}_charge.txt -chargeB_in_file {4}/{1}_charge.txt -mdinfo {2}/ti/{1}.mdinfo -mdout {2}/ti/{1}.mdout -crd {2}/balance/{1}.dat -box {2}/balance/{1}.box -lambda_lj {3} -subsys_division_in_file {2}/{1}_subsys_division.txt  -charge_pertubated 1 -atom_numbers {5} -frame_numbers {6} -TI dh_dlambda.txt".format(
                            args.sponge_ti, args.temp, i, i / args.nl, args.nl, len(merged_from.atoms),
                                                          args.balance_step // 100))
                else:
                    os.system(
                        "{0} -LJ_soft_core_in_file {2}/{1}_LJ_soft_core.txt -exclude_in_file {2}/{1}_exclude.txt -charge_in_file {2}/{1}_charge.txt -chargeA_in_file 0/{1}_charge.txt -chargeB_in_file {4}/{1}_charge.txt -mdinfo {2}/ti/{1}.mdinfo -mdout {2}/ti/{1}.mdout -crd {2}/balance/{1}.dat -box {2}/balance/{1}.box -lambda_lj {3} -subsys_division_in_file {2}/{1}_subsys_division.txt  -charge_pertubated 1 -TI dh_dlambda.txt -mdin {5}".format(
                            args.sponge_ti, args.temp, i, i / args.nl, args.nl, args.ai))
            dh_dlambda = np.loadtxt("dh_dlambda.txt")
            dh = []
            dh_int = []
            tempall = 0
            for i in range(args.nl):
                temp = dh_dlambda[i] * 0.5 / args.nl
                temp += dh_dlambda[i + 1] * 0.5 / args.nl
                dh.append(temp)
                tempall += temp
                dh_int.append(tempall)
            with open("free_energy.txt", "w") as f:
                f.write("lambda_state\tFE(i+1)-FE(i)[kcal/mol]\tFE(i+1)-FE(0)[kcal/mol]\n")
                f.write("\n".join(["%d\t\t%.2f\t\t\t%.2f" % (i, dh[i], dh_int[i]) for i in range(args.nl)]))
        elif args.method == "FEP_BAR":
            raise NotImplementedError


def mol2rfe(args):
    """

    :param args:
    :return:
    """
    source("..")
    source("..forcefield.special.fep")
    min =  source("..forcefield.special.min")

    if not args.ff:
        source("..forcefield.amber.gaff")
        source("..forcefield.amber.ff14sb")
        source("..forcefield.amber.tip3p")
    else:
        idic, ipy = os.path.split(args.ff)
        sys.path.append(idic)
        ipy, isuffix = os.path.splitext(ipy)
        assert isuffix == ".py", "the input force field file should be an Xponge file written by python"
        __import__(ipy)

    if not args.do:
        args.do = [["build", "min", "prebalance", "balance", "analysis"]]
    args.do = args.do[0]

    from_res_type_ = load_mol2(args.r1).residues[0]
    from_ = assign.Get_Assignment_From_ResidueType(from_res_type_)
    if not args.ff:
        parmchk2_gaff(args.r1, args.temp + "_TMP1.frcmod")

    to_res_type_ = Xponge.load_mol2(args.r2).residues[0]
    to_ = Xponge.assign.Get_Assignment_From_ResidueType(to_res_type_)
    if not args.ff:
        parmchk2_gaff(args.r2, args.temp + "_TMP2.frcmod")

    for mol2file in args.r0:
        load_mol2(mol2file)

    rmol = load_pdb(args.pdb)

    merged_from, merged_to = Merge_Dual_Topology(rmol, rmol.residues[args.ri], to_res_type_, from_, to_, args.tmcs)

    if args.dohmr:
        H_Mass_Repartition(merged_from)
        H_Mass_Repartition(merged_to)

    mol2rfe_build(args, merged_from, merged_to)

    mol2rfe_min(args)

    mol2rfe_prebalance(args)

    mol2rfe_balance(args)

    mol2rfe_analysis(args, merged_from)
