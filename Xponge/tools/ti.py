"""
This **module** contains the functions for ti analysis
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from ..mdrun import run
from ..helper import Xopen
from ..analysis import MdoutReader


def ti_analysis(args, merged_from):
    """
    This **function** is used to do the ti analysis

    :param args: the arguments from the command line
    :return: None
    """
    ses = []
    prefix_sum = []
    suffix_sum = []
    frame = args.equilibrium_step // 100
    for i in range(args.nl + 1):
        if os.path.exists("%d/ti" % i):
            os.system("rm -rf %d/ti" % i)
        os.mkdir("%d/ti" % i)
        inprefix = f"{i}/{args.temp}"
        command = f"SPONGE_TI -LJ_soft_core_in_file {inprefix}_LJ_soft_core.txt"
        command += " -exclude_in_file {0}_exclude.txt -charge_in_file {0}_charge.txt".format(inprefix)
        command += f" -chargeA_in_file 0/{args.temp}_charge.txt"
        command += f" -chargeB_in_file {args.nl}/{args.temp}_charge.txt"
        lambda_ = i / args.nl
        command += f" -lambda_lj {lambda_}"
        command += f" -subsys_division_in_file {inprefix}_subsys_division.txt  -charge_pertubated 1"
        inprefix = f"{i}/ti/{args.temp}"
        command += f" -mdinfo {inprefix}.mdinfo -mdout {inprefix}.mdout"
        inprefix = f"{i}/equilibrium/{args.temp}"
        command += f" -crd {inprefix}.dat -box {inprefix}.box -TI dh_dlambda.txt"
        command += f" -atom_numbers {len(merged_from.atoms)}"
        command += f" -frame_numbers {frame}"
        if not args.ai:
            run(command)
        else:
            command += f" -mdin {args.ai}"
            run(command)
        temp = MdoutReader(f"{i}/ti/{args.temp}.mdout").dH_dlambda
        prefix_sum.append(np.cumsum(temp[::]) / (np.arange(frame) + 1))
        suffix_sum.append(np.cumsum(temp[::-1]) / (np.arange(frame) + 1))
        ses.append(np.std(temp) / np.sqrt(frame))
    prefix_sum = np.array(prefix_sum)
    dh_dlambda = np.loadtxt("dh_dlambda.txt")
    ses = np.array(ses)
    ses *= ses
    dh = []
    dh_int = []
    dh_se = []
    dh_int_se = []
    tempall = 0
    temp_se_all = 0
    space = 0.5 / args.nl
    if os.path.exists("time_check"):
        os.system("rm -rf time_check")
    os.mkdir("time_check")
    for i in range(args.nl):
        temp = (prefix_sum[i] + prefix_sum[i + 1]) * space
        time = args.dt * 0.1 * (np.arange(frame) + 1)
        plt.plot(time, temp, label="forward")
        temp_ses = np.vstack((time, temp))
        temp = (suffix_sum[i] + suffix_sum[i + 1]) * space
        temp_ses = np.vstack((temp_ses, temp))
        plt.plot(time, temp, label="backward")
        plt.xlabel("time[ns]")
        plt.ylabel("free energy difference[kcal/mol]")
        plt.legend()
        plt.savefig(f"time_check/{i}-{i + 1}.png")
        np.savetxt(f"time_check/{i}-{i + 1}.csv", temp_ses.transpose(),
                   header="time[ps],forward DeltaG[kcal/mol],backward DeltaG[kcal/mol]", comments="", delimiter=",")
        plt.clf()
        temp = dh_dlambda[i] * space
        temp += dh_dlambda[i + 1] * space
        temp_ses = space * (ses[i] + ses[i + 1])
        dh.append(temp)
        dh_se.append(np.sqrt(temp_ses))
        tempall += temp
        temp_se_all += temp_ses
        dh_int.append(tempall)
        dh_int_se.append(np.sqrt(temp_se_all))

    temp_ses **= 0.5
    f = Xopen("free_energy.txt", "w")
    f.write("lambda_state\tFE(i+1)-FE(i)[kcal/mol]\tFE(i+1)-FE(0)[kcal/mol]\n")
    f.write("\n".join(
        [f"{i}\t\t{dh[i]: .2f} +- {dh_se[i]:.2f}\t\t{dh_int[i]: .2f} +- {dh_int_se[i]:.2f}" for i in range(args.nl)]))
    f.close()
