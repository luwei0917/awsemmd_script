#!/usr/bin/env python3
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
import imp
import glob
from time import sleep
import fileinput
import numpy as np
import pandas as pd
from small_script.variable_test import variable_test
from small_script.variable_test2 import variable_test2
import subprocess
from small_script.myFunctions import compute_theta_for_each_helix
from small_script.myFunctions import *

# Useful codes
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# tr " " "\n"
# sed 1d
# sort -u -k 3
# sed -e 's/+T//'
# import re
# numbers = re.compile(r'(\d+)')
# def numericalSort(value):
#     parts = numbers.split(value)
#     parts[1::2] = map(int, parts[1::2])
#     return parts
# mypath = os.environ["PATH"]
# os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
# my_env = os.environ.copy()

parser = argparse.ArgumentParser(description="This is my playground for current project")
parser.add_argument("-r", "--run", help="test mode",
                    action="store_true")
parser.add_argument("-s", "--see", help="test mode",
                    action="store_true")
# parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("-m", "--mode", type=int, default=0)
parser.add_argument("-d", "--day", type=str, default="someday")
parser.add_argument("-l", "--label", type=str, default="label")
parser.add_argument("-t", "--test", action="store_true", default=False)
args = parser.parse_args()

if args.test:
    do = print
else:
    do = os.system
cd = os.chdir

base_run_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST

module load GCC/4.9.3 OpenMPI/1.8.8
srun /home/wl45/build/awsem_lipid_fluctuations/src/lmp_mpi -p 12x1 -in 2xov_{}.in\n'''

base_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=01:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun {}\n'''


scavenge_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=commons
#SBATCH --partition=scavenge
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=04:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun {}\n'''


def replace(TARGET, FROM, TO):
    do("sed -i.bak 's/{}/{}/g' {}".format(FROM,TO,TARGET))



def continueRunConvertion(n=12, rerun=0, name="2xov", convert_read_data=False):
    rerun_plus_one = rerun + 1
    do(f"cp {name}_0.in {name}_{rerun_plus_one}.in")
    fileName = f"{name}_{rerun_plus_one}.in"
    replace(fileName, "variable r world "+ " ".join(str(i) for i in list(range(n))), "")
    replace(fileName, "# read_restart restart.25000000", "variable r world "+ " ".join(str(i) for i in list(range(n))))
    initial_steps = 20000000 * rerun_plus_one
    replace(fileName, "read_restart restart.extended", f"read_restart restart.$r.{initial_steps}")
    replace(fileName, "read_restart restart.native_topology", f"read_restart restart.$r.{initial_steps}")
    if convert_read_data:
        replace(fileName, "read_data data.2xov", f"read_restart restart.$r.{initial_steps}")
    replace(fileName, "0\/", f"{rerun_plus_one}\/")
    cmd = 'tail -n 1 log.lammps | cut -d" " -f2-'
    line = getFromTerminal(cmd).rstrip()
    replace(fileName, "reset_timestep	0", "variable w world " + line)
    replace(fileName, "fix     xbias all colvars colvars.x output x", "fix     xbias all colvars colvars.x output x.$r")
    cmd = f'grep "temper" {name}_{rerun_plus_one}.in'
    line = getFromTerminal(cmd).rstrip()
    replace(fileName, line, line + " $w")




def rerun(extra="", mode=0, useNewTable=False):
    print(mode)
    for i in range(12):
        name = f"{extra}_{i}"
        do(f"mkdir -p {name}")
        cd(str(name))
        do("cp ~/opt/2xov_eval/* .")
        if mode == 0:
            do("cp 2xov_eval_mar31.in 2xov_eval.in")
            do(f"cp 2xov_{extra}.seq 2xov.seq")
        elif mode == 1:
            do("cp 2xov_eval_apr13.in 2xov_eval.in")
        if mode == 2:
            do("cp 2xov_eval_n_side.in 2xov_eval.in")
        if mode == 3:
            do("cp 2xov_eval_n_side_3_helix.in 2xov_eval.in")
        if mode == 4:
            do("cp 2xov_eval_n_side_4_helix.in 2xov_eval.in")
        if useNewTable:
            create_zim("2xov.seq", isNew=False)
        else:
            create_zim("2xov.seq", isNew=True)
        fileName = "2xov_eval.in"
        if mode == 1:  # change lipid environment
            replace(fileName, "DMPC", extra)
        replace(fileName, "MY_RUN", str(i))
        # do("cp ../../zim .")
        replace("run.slurm", "ctbp-common", "commons")
        # replace("run.slurm", "ctbp-common", "interactive")
        do("sbatch run.slurm")
        cd("..")

# def extra(fileName, offset=0):
#     replace(fileName, "1 1 30 5", "1 1 30 {}".format(offset))
# def rerun(extra=extra, offset=0):
#     do("cp 2xov_0.in rerun_{}.in".format(offset))
#     fileName = "rerun_{}.in".format(offset)
#     replace(fileName, "fix               1 all nve", "")
#     replace(fileName, "fix               2 all langevin", "#")
#     replace(fileName, "0\/", "recompute_offset_{}\/".format(offset))
#     replace(fileName, "dump		1 all atom 4000", "#")
#     replace(fileName, "dump_modify	1 sort id", "")
#     replace(fileName, "run", "#")
#     replace(fileName, "restart         100000 restart", "rerun 0\/dump.lammpstrj dump x y z")

#     # replace(fileName, "1 1 30 5", "1 1 30 0")
#     # extra(fileName, offset)
#     slurm = "rerun_{}.slurm".format(offset)
#     do("cp ~/opt/2xov_eval/run.slurm " + slurm)
#     replace(slurm, "2xov_eval", "rerun_{}".format(offset))
#     # replace(slurm, "commons", "ctbp-common")
#     do("mkdir recompute_offset_{}".format(offset))
#     do("sbatch " + slurm)

def scancel_jobs_in_folder(folder):
    cd(bias)
    cmd = "find -name 'slurm-*' | rev | awk -F'[-.]' '{print $2}' | rev"
    lines = getFromTerminal(cmd).splitlines()
    for line in lines:
        # print(line)
        do("scancel " + line)
    cd("..")

localQ_slurm = '''#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:30:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun python3 ~/opt/gg_server.py -d feb28 -m 2
'''
quick_slurm = '''#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:30:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun python3 ~/opt/gg_server.py -d mar03 -m 3
'''

freeEnergy = """\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=23:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun python2 ~/opt/pulling_compute-pmf.py {}
"""
quick_template_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=01:30:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun {}
'''
quick_template_large_mem_slurm = '''#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=01:30:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun python3 ~/opt/gg_server.py {}
'''

def compute_quantity(cmd="", queue=1, sim_list=["0"], bias=""):
    """Compute some quantity.
    Queue 0 is ctbp-common, 1 is interactive, 2 is commons.
    bias is pre name of simulation folder
    """
    # print(cmd)
    simulation_list = glob.glob(f"{bias}*")
    print(simulation_list)
    dateAndTime = datetime.today().strftime('%d_%h_%H%M%S')
    # print(sim_list)
    for sim in sim_list:
        for folder in simulation_list:
            cd(folder)
            cd(sim)
            print(folder)
            quick = quick_template_slurm.format(cmd)
            with open(f"quick_{dateAndTime}.slurm", "w") as f:
                if queue == 1:
                    quick = quick.replace("--time=01:30:00", "--time=00:30:00")
                    quick = quick.replace("#SBATCH --account=ctbp-common", "")
                    quick = quick.replace("ctbp-common", "interactive")
                if queue == 2:
                    quick = quick.replace("ctbp-common", "commons")
                f.write(quick)
            do(f"sbatch quick_{dateAndTime}.slurm")
            cd("../..")

def compute_completeZ(temper=False, **kwargs):
    print("compute completeZ")
    if temper:
        cmd = "python3 ~/opt/gg_server.py -d mar10 -m 7"
    else:
        cmd = "python3 ~/opt/gg_server.py -d mar16 -m 3"
    compute_quantity(cmd=cmd, **kwargs)

def compute_disReal(temper=False, targetMode=0, name="2xov", **kwargs):
    print("compute DisReal")
    # if targetMode == 0:
    if temper:
        cmd = f"python3 ~/opt/small_script/temper_compute_HQ.py {name} -t {targetMode}"
    else:
        cmd = f"python3 ~/opt/small_script/find_distance.py {name} -t {targetMode}"
    # elif targetMode == 1:
    #     if temper:
    #         cmd = "python3 ~/opt/gg_server.py -d apr01 -m 3"
    #     else:
    #         cmd = f"python3 ~/opt/small_script/find_distance.py -t {targetMode}"
    compute_quantity(cmd=cmd, **kwargs)

# def compute_NsideEnergy(temper=False, targetMode=0, name="2xov", **kwargs):
#     print("compute DisReal")
#     # if targetMode == 0:
#     if temper:
#         cmd = f"python3 ~/opt/gg_server.py -d mar17 -m 1"
#     else:
#         cmd = f"python3 ~/opt/small_script/find_distance.py {name} -t {targetMode}"
#     # elif targetMode == 1:
#     #     if temper:
#     #         cmd = "python3 ~/opt/gg_server.py -d apr01 -m 3"
#     #     else:
#     #         cmd = f"python3 ~/opt/small_script/find_distance.py -t {targetMode}"
#     compute_quantity(cmd=cmd, **kwargs)

def let_compute_localQ(temper=False, **kwargs):
    print("compute LocalQ")
    if temper:
        cmd = "python3 ~/opt/gg_server.py -d mar28 -m 1"
    # if temper:
    #     cmd = "python3 ~/opt/gg_server.py -d feb28 -m 2"
    # else:
    #     pass
        # cmd = "python3 ~/opt/small_script/find_distance.py"
        #         native_contacts_table = compute_localQ_init()
        # for i in range(12):
        #     compute_localQ(native_contacts_table, pre=".", ii=i)
    compute_quantity(cmd=cmd, **kwargs)

if args.day == "common":
    if args.mode == 4:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_30.0']
        # simulation_list = ['dis_116.0', 'dis_332.0', 'dis_128.0', 'dis_266.0', 'dis_296.0', 'dis_290.0', 'dis_314.0', 'dis_176.0', 'dis_272.0', 'dis_284.0', 'dis_158.0', 'dis_338.0']
        # simulation_list = ['dis_326.0', 'dis_206.0', 'dis_254.0', 'dis_344.0', 'dis_308.0', 'dis_134.0', 'dis_152.0', 'dis_194.0', 'dis_320.0', 'dis_200.0', 'dis_212.0', 'dis_110.0', 'dis_248.0', 'dis_188.0', 'dis_242.0', 'dis_218.0', 'dis_350.0', 'dis_164.0', 'dis_236.0', 'dis_146.0', 'dis_182.0', 'dis_140.0', 'dis_122.0', 'dis_302.0', 'dis_224.0', 'dis_230.0', 'dis_278.0', 'dis_260.0', 'dis_170.0']
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            # cd("0")
            cd("1")
            # cd("2")
            # cd("3")
            # rerun(extra="Go", mode=2)
            # rerun(extra="Go_3helix", mode=3)
            rerun(extra="Go_4helix", mode=4)
            cd("../..")
    if args.mode == 3:
        goEnergy = False
        goEnergy3H = False
        goEnergy4H = True
        rerun = 3
        end = 2
        cwd = os.getcwd()
        print(cwd)
        pre = '/'.join(cwd.split("/")[:-2]) + "/"
        print(pre)
        # exit()
        # pre = "/scratch/wl45/apr_2018/sixth/"
        data_folder = "/scratch/wl45/aug_2018/02_week/freeEnergy/all_data_folder/"
        # folder_list = ["rg_0.15_lipid_1.0_mem_1_go_0.8_long"]
        folder_list = [cwd.split("/")[-2]]
        # label = "sixth_long"
        label = args.label
        # with open(label, "w") as f:
        #     f.write("test\n")
        # cd("simulation")
        # exit()
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]

        for attempt in range(2):
            try:
                process_complete_temper_data_3(pre, data_folder, folder_list, rerun=rerun, end=end, average_z=True, disReal=True, dis_h56=True, localQ=False, goEnergy=goEnergy, goEnergy3H=goEnergy3H, goEnergy4H=goEnergy4H, label=label)
                break
            except FileNotFoundError:
                bias = "dis"
                simulation_list = glob.glob(f"{bias}_*")
                # simulation_list = ['dis_30.0']
                # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
                print(simulation_list)
                for dis in simulation_list:
                    print(dis)
                    cd(dis)
                    i = rerun
                    i_plus_one = i +1
                    if os.path.exists(f"log{i}"):
                        do(f"mv log{i} back_log{i}")  # in case of override
                    do(f"mkdir -p log{i}")
                    do(f"cp log.lammps log{i}/")
                    cd("..")
    if args.mode == 2:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_30.0']
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 2
            i_plus_one = i +1
            # do(f"mv log{i} back_log{i}")  # in case of override
            # do(f"mkdir -p log{i}")
            # do(f"cp log.lammps log{i}/")

            # my_file = sys.Path(f"{i_plus_one}")
            # if my_file.is_dir():
            #     print("Attension")
            #     exit()
            continueRunConvertion(n=12, rerun=i)
            # continueRunConvertion(n=12, rerun=i, convert_read_data=True)
            do(f"mkdir {i_plus_one}")
            # do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            # replace(f"run_{i_plus_one}.slurm", "/home/ns24/lmp_mpi", "/home/wl45/build/awsem_lipid_fluctuations/src/lmp_mpi")
            # do(f"sbatch run_{i_plus_one}.slurm")
            run_slurm = base_run_slurm.format(i_plus_one)
            with open(f"run_{i_plus_one}.slurm", "w") as r:
                r.write(run_slurm)
            # replace(f"run_{i_plus_one}.slurm", "ctbp-common", "commons")
            do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
    if args.mode == 1:
        queue = 2
        # i = "0"
        # i = "1"
        i = "2"
        # i = "3"
        # let_compute_localQ(temper=True, bias="dis_", sim_list=[i], queue=1)

        # compute_disReal(temper=True, bias="dis_", sim_list=[i], queue=queue)
        # compute_disReal(temper=True, targetMode=1, bias="dis_", sim_list=[i], queue=queue)
        # compute_completeZ(temper=True, bias="dis_", sim_list=[i], queue=queue)

        compute_disReal(temper=True, targetMode=2, bias="dis_", sim_list=[i], queue=queue)
        compute_disReal(temper=True, targetMode=3, bias="dis_", sim_list=[i], queue=queue)

def removeResXfromlist():
    with open(f"database/cath-dataset-nonredundant-S20Clean.list", "w") as out2:
        with open(f"database/cath-dataset-nonredundant-S20Clean.atom.fa", "w") as out:
            with open("database/cath-dataset-nonredundant-S20.atom.fa", "r") as f:
                count = 0
                for l in f:
                    if count % 2 == 0:
                        #extract protein id
                        assert(l[0] == ">")
            #             print(l)
                        tmp = l
                        name = re.search('>cath\|(.*)\|(\w{7})\/(.*)', l).group(2)
            #             name = "test"
            #             print(name)
                    else:
                        assert(l[0] != ">")
            #             print(l)
                        if "X" in l:
                            pass
                        else:
                            out.write(tmp)
                            out.write(l)
                            out2.write(name+"\n")
                    count += 1

def removeExtraName():
    '''
    remove The 'B' or possible 'A' at position 16
    for example, chagne from
    ATOM    193  CB BMET A  30     -20.305 -21.245 -45.095  0.50 10.77
    to
    ATOM    193  CB  MET A  30     -20.305 -21.245 -45.095  0.50 10.77
    '''
    p_list = glob.glob("database/dompdb_origin/*")
    for name in p_list:
        toName = name.replace("dompdb_origin", "dompdb_cleaned")
        with open(toName, "w") as out:
            with open(name, "r") as f:
                for l in f:
                    tmp = list(l)
                    tmp[16] = " "
                    out.write("".join(tmp))


def isComplete(a):
    for model in a:
        for chain in model:
            for res in chain:
                try:
                    if res["CA"] is not None:
                        pass
                    if res.get_resname() == "GLY" or res["CB"] is not None:
                        pass
                    if res["N"] is not None:
                        pass
                    if res["C"] is not None:
                        pass
                except:
                    print(res)
                    return 0
    return 1


if args.day == "nov25":
    if args.mode == 6:
        for i in range(180, 256):
            with open(f"run_{i}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python2 ../top_fix_pick.py -i {i}"))
            do(f"sbatch run_{i}.slurm")
    if args.mode == 5:
        for i in range(256):
            with open(f"run_{i}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python2 ../top_fix_random.py -i {i}"))
            do(f"sbatch run_{i}.slurm")
    if args.mode == 1:
        for i in range(400):
            with open(f"run_{i}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python2 ../randomGen_nov25.py t_{i}.csv -n 1e7"))
            do(f"sbatch run_{i}.slurm")
    if args.mode == 2:
        for i in range(400):
            with open(f"run_{i}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python2 ../geneticAlgorithm_fast.py"))
            do(f"sbatch run_{i}.slurm")
    if args.mode == 3:
        a = pd.read_csv("../rr_nov25.csv", index_col=0)
        size = len(a)
        N = 200
        n = int(size/N)+1
        cc = 0
        for i in range(N):
            if cc+n >= size:
                end = size-1
            else:
                end = cc+n
            start = cc
            with open(f"run_{i}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python2 ../playAgainstEachOther.py -s {start} -e {end}"))
            do(f"sbatch run_{i}.slurm")
            cc += n
    if args.mode == 4:
        a = pd.read_csv("../rr_2_nov25.csv", index_col=0)
        size = len(a)
        N = 20
        n = int(size/N)+1
        cc = 0
        for i in range(N):
            if cc+n >= size:
                end = size-1
            else:
                end = cc+n
            start = cc
            with open(f"run_{i}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python2 ../playAgainstEachOther.py -s {start} -e {end}"))
            do(f"sbatch run_{i}.slurm")
            cc += n

        # for i in range(100):
        #     with open(f"run_{i}.slurm", "w") as out:
        #         out.write(scavenge_slurm.format(f"python2 ../geneticAlgorithm.py"))
        #     do(f"sbatch run_{i}.slurm")
if args.day == "nov20":
    if args.mode == 3:
        name_list = ["tr894", "tr882", "tr594", "tr869", "tr898", "tr862", "tr877", "tr872", "tr885", "tr866", "tr868", "tr884", "tr895", "tr896", "tr870", "tr921", "tr922", "tr891", "tr948", "tr947"]
        for name in name_list:
            do(f"mkdir {name}")
            # do(f"cp /scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/selection/{name}/PCselection/* {name}/")
            do(f"cp /scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/selection/PCselection/{name}/* {name}/")
    if args.mode == 2:
        name_list = ["tr894", "tr882", "tr594", "tr869", "tr898", "tr862", "tr877", "tr872", "tr885", "tr866", "tr868", "tr884", "tr895", "tr896", "tr870", "tr921", "tr922", "tr891", "tr948"]
        for name in name_list:
            do(f"mkdir {name}")
            cd(name)
            do(f"cp /scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/selection/{name}/towardsPC1/awsem_energy/bias.log .")
            do(f"cp /scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/selection/{name}/towardsPC1/awsem_energy/awsem.log .")
            do(f"cp /scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/selection/{name}/towardsPC1/awsem_energy/rwplusScore.txt .")
            do(f"cp /scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/selection/{name}/towardsPC1/awsem_energy/rmsd-angstrom.xvg .")
            cd("..")
    if args.mode == 1:
        name_list = ["tr894", "tr882", "tr594", "tr869", "tr898", "tr862", "tr877", "tr872", "tr885", "tr866", "tr868", "tr884", "tr895", "tr896", "tr870", "tr921", "tr922", "tr891", "tr948"]
        for name in name_list:
            do(f"mkdir {name}")
            # do(f"cp /scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/selection/{name}/PCselection/* {name}/")
            do(f"cp /scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/selection/PCselection/{name}/* {name}/")
if args.day == "nov09":
    pdb_list = "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", ")
    if args.mode == 1:
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            if name == "2fha" or name == "1mba":
                steps = 20
            else:
                steps = 40
            cd(name)
            do(f"run.py -n 20 {name} --commons 2 -s {steps} --runs 3")
            cd("..")

if args.day == "nov01":
    if args.mode == 1:
        with open("database/cath-dataset-nonredundant-S20Clean.list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list_{i}.txt", "w") as out:
                for ii in range(30):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i-1)
    if args.mode == 2:
        # clean up
        # do("rm database/S20_seq/*")
        # do("rm database/dompdb/*")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p phis")
    if args.mode == 3:
        for i in range(475):
            with open(f"run_{i}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py proteins_name_list_{i}.txt"))
            do(f"sbatch run_{i}.slurm")
    if args.mode == 4:
        # clean up
        # do("rm database/S20_seq/*")
        do("rm *.out")
        do("rm *.slurm")
        # do("rm database/dompdb/*")
        do("rm -r decoys/*")
        do("rm gammas/*")
        do("rm -r phis")
        do("mkdir -p phis")

if args.day == "oct26":
    if args.mode == 1:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        complete_proteins = "database/cath-dataset-nonredundant-S20Clean.list"
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_xl23(complete_proteins, "phi_list.txt", decoy_method='shuffle', num_decoys=1000, noise_filtering=True, jackhmmer=False)
        # individual_gammas_randomized_decoy=read_all_gammas("phi_list.txt", complete_proteins, training_decoy_method="shuffle", noise_filtering=True)
    if args.mode == 2:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        complete_proteins = "database/test.list"
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_xl23(complete_proteins, "phi_list.txt", decoy_method='shuffle', num_decoys=1000, noise_filtering=True, jackhmmer=False)
        # individual_gammas_randomized_decoy=read_all_gammas("phi_list.txt", complete_proteins, training_decoy_method="shuffle", noise_filtering=True)



if args.day == "oct25":
    if args.mode == 1:
        with open("database/cath-dataset-nonredundant-S20Clean.list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        while pos < n:
            with open(f"proteins_name_list_{i}.txt", "w") as out:
                for ii in range(50):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i-1)
    if args.mode == 2:
        removeExtraName()
    if args.mode == 3:
        p_list = glob.glob("database/dompdb/*")
        for p in p_list:
            from pyCodeLib import *
            import warnings
            warnings.filterwarnings('ignore')
            a = parse_pdb(p[:-4])
            if isComplete(a):
                pass
            else:
                print(p)
    if args.mode == 4:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # I may want to just extract good part, instead of remove completely.
        with open(f"database/cath-dataset-nonredundant-S20Clean.list", "w") as out2:
            with open(f"database/cath-dataset-nonredundant-S20Clean.atom.fa", "w") as out:
                with open("database/list_back/cath-dataset-nonredundant-S20Clean.atom.fa", "r") as f:
                    count = 0
                    for l in f:
                        if count % 2 == 0:
                            #extract protein id
                            assert(l[0] == ">")
                #             print(l)
                            tmp = l
                            name = re.search('>cath\|(.*)\|(\w{7})\/(.*)', l).group(2)
                #             name = "test"
                #             print(name)
                        else:
                            assert(l[0] != ">")
                #             print(l)
                            a = parse_pdb(f"database/dompdb/{name}")
                            if "X" in l:
                                pass
                            elif not isComplete(a):
                                pass
                            else:
                                out.write(tmp)
                                out.write(l)
                                out2.write(name+"\n")
                        count += 1

if args.day == "oct24":
    if args.mode == 1:
        # clean up
        # do("rm database/S20_seq/*")
        do("rm *.out")
        do("rm *.slurm")
        # do("rm database/dompdb/*")
        do("rm -r decoys/*")
        do("rm gammas/*")
        do("rm phis/*")
    if args.mode == 2:
        for i in range(285):
            with open(f"run_{i}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py proteins_name_list_{i}.txt"))
            do(f"sbatch run_{i}.slurm")
    if args.mode == 3:
        with open("proteins_name_list_56.txt") as f:
            content = f.readlines()
        # you may also want to remove whitespace characters like `\n` at the end of each line
        name_list = [x.strip() for x in content]
        for n in name_list:
            do(f"echo {n}")
            # do(f"grep 'CA' database/dompdb/{n}.pdb | wc")
            # do(f"grep 'CA' database/dompdb/{n}.pdb | wc")
            do(f"grep 'UNK' database/dompdb/{n}.pdb")
    if args.mode == 4:
        # clean up the database list
        removeResXfromlist()
    if args.mode == 5:
        # debug list 56
        with open("proteins_name_list_56.txt") as f:
            content = f.readlines()
        # you may also want to remove whitespace characters like `\n` at the end of each line
        pos = 0
        for i in range(40):
            with open(f"test_list_{i}.txt", "w") as out:
                for ii in range(5):
                    if pos < len(content):
                        out.write(content[pos])
                    pos += 1
    if args.mode == 6:
        for i in range(40):
            with open(f"run_{i}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python3 compute_phis.py test_list_{i}.txt"))
            do(f"sbatch run_{i}.slurm")


if args.day == "oct14":
    if args.mode == 1:
        for i in range(400):
            with open(f"run_{i}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python2 ../randomGen_nov19.py t_{i}.csv -n 2e6"))
            do(f"sbatch run_{i}.slurm")


if args.day == "sep21":
    if args.mode == 1:
        pdb_list = glob.glob(f"lowTstructure*.pdb")
        n = len(pdb_list)
        n_in_thousand = (n // 1000) + 1
        print(n)
        i = 0
        for i in range(10):
            for thousand in range(n_in_thousand):
                # computeMutualQ(i, j, thousand, n)
                do(f"cp run_template.slurm run_{i}_{thousand}.slurm")
                with fileinput.FileInput(f"run_{i}_{thousand}.slurm", inplace=True, backup='.bak') as file:
                    for line in file:
                        tmp = line.replace("FROM", str(i))
                        tmp = tmp.replace("NUMBER", str(n))
                        tmp = tmp.replace("THOUSAND", str(thousand))
                        # tmp = BETA_ATOMS
                        print(tmp, end='')
                do(f"sbatch run_{i}_{thousand}.slurm")
if args.day == "sep20":
    if args.mode == 1:
        name_list = ["tr894", "tr882", "tr594", "tr869", "tr898", "tr862", "tr877", "tr872", "tr885", "tr866", "tr868", "tr884", "tr895", "tr896", "tr870", "tr921", "tr922", "tr891", "tr948"]
        for name in name_list:
            do(f"mkdir {name}")
            # do(f"cp /scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/selection/{name}/PCselection/* {name}/")
            do(f"cp /scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/selection/PCselection/{name}/* {name}/")

if args.day == "sep16":
    if args.mode == 1:
        name_list = ["tr894", "tr882", "tr594", "tr869", "tr898", "tr862", "tr877", "tr872", "tr885", "tr866", "tr868", "tr884", "tr895", "tr896", "tr870", "tr921", "tr922", "tr891", "tr948"]
        for name in name_list:
            do(f"mkdir {name}")
            cd(name)
            do(f"cp /scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/selection/{name}/towardsPC1/awsem_energy/bias.log .")
            do(f"cp /scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/selection/{name}/towardsPC1/awsem_energy/awsem.log .")
            do(f"cp /scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/selection/{name}/towardsPC1/awsem_energy/rwplusScore.txt .")
            do(f"cp /scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/selection/{name}/towardsPC1/awsem_energy/rmsd-angstrom.xvg .")
            cd("..")

if args.day == "aug29":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_zAverage_dis":"17", "2d_z_qw":"13"}
        bias_list = {"2d_zAverage_dis":"17"}
        # bias_list = {"2d_dis_z56":"21", "2d_zAverage_dis":"22", "56_z_dis":"23", "2d_dis12_z14":"24", "2d_dis34_z14":"25", "2d_dis56_z14":"26"}
        # bias_list = {"2d_dis_z56":"21"}
        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"how_force_bias_added_5"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        # folder_list = [f"second_start_extended_combined_may19"]
        folder_list = [f"include_h14_aug18"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["300", "335", "373"]}
        temp_dic = {"_280-350":["335", "373", "417"]}
        temp_dic = {"_280-350":["335", "373"]}
        # temp_dic = {"_280-350":["373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=13)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    # nsample = len(folder_list)*2500
                    nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 32 --force 0")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "aug16":
    if args.mode == 1:
        for i in range(12):
            compute_theta_for_each_helix(output=f"angles_{i}.csv", dumpName=f"dump.lammpstrj.{i}")
    if args.mode == 2:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_zAverage_dis":"17", "2d_z_qw":"13"}
        bias_list = {"2d_zAverage_dis":"17"}
        bias_list = {"2d_dis_z56":"21", "2d_zAverage_dis":"22", "56_z_dis":"23", "2d_dis12_z14":"24", "2d_dis34_z14":"25", "2d_dis56_z14":"26"}
        # bias_list = {"2d_dis_z56":"21"}
        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"combined_more_rc"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"include_h14"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["300", "335", "373"]}
        temp_dic = {"_280-350":["335", "373", "417"]}
        temp_dic = {"_280-350":["335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=7, chosen_mode=9)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    # nsample = len(folder_list)*2500
                    nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 30 --force 10")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "aug15":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_zAverage_dis":"17", "2d_z_qw":"13"}
        bias_list = {"2d_zAverage_dis":"17"}
        bias_list = {"2d_zAverage_dis":"20", "56_z_dis":"18"}
        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"first_half"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"part1"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["300", "335", "373"]}
        temp_dic = {"_280-350":["335", "373", "417"]}
        temp_dic = {"_280-350":["335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=6, chosen_mode=9)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    # nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 29 --force 10")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 2:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_zAverage_dis":"17", "2d_z_qw":"13"}
        bias_list = {"2d_zAverage_dis":"17"}
        bias_list = {"2d_zAverage_dis":"20", "56_z_dis":"18"}
        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_half"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"part2"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["300", "335", "373"]}
        temp_dic = {"_280-350":["335", "373", "417"]}
        temp_dic = {"_280-350":["335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=6, chosen_mode=9)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    # nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 29 --force 10")
                    cd("..")
                cd("..")
        cd("..")



'''
if args.day == "aug11":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_zAverage_dis":"17", "2d_z_qw":"13"}
        bias_list = {"2d_zAverage_dis":"17"}
        bias_list = {"2d_zAverage_dis":"20", "56_z_dis":"18"}
        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"combined_more_force_fix_order"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_start_extended_combined_may19"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["300", "335", "373"]}
        temp_dic = {"_280-350":["335", "373", "417"]}
        temp_dic = {"_280-350":["335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=6, chosen_mode=9)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    # nsample = len(folder_list)*2500
                    nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 29 --force 10")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "aug09":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_zAverage_dis":"17", "2d_z_qw":"13"}
        bias_list = {"2d_zAverage_dis":"17"}
        bias_list = {"2d_zAverage_dis":"20", "56_z_dis":"18"}
        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"combined"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_start_extended_combined_may19"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["300", "335", "373"]}
        temp_dic = {"_280-350":["335", "373", "417"]}
        temp_dic = {"_280-350":["335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=6, chosen_mode=9)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    # nsample = len(folder_list)*2500
                    nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 29 --force 8")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 2:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_zAverage_dis":"17", "2d_z_qw":"13"}
        bias_list = {"2d_zAverage_dis":"17"}
        bias_list = {"2d_zAverage_dis":"20", "56_z_dis":"18"}
        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"combined_2"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_start_extended_combined_may19"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["300", "335", "373"]}
        temp_dic = {"_335-417":["335", "373", "417"]}
        # temp_dic = {"_280-350":["335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=6, chosen_mode=9)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    # nsample = len(folder_list)*2500
                    nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 29 --force 8")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "aug08":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_zAverage_dis":"17", "2d_z_qw":"13"}

        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_start_extended_combined"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_start_extended_combined_may19"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["300", "335", "373"]}
        temp_dic = {"_280-350":["335", "373", "417"]}
        temp_dic = {"_280-350":["335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=9)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    # nsample = len(folder_list)*2500
                    nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 21 --force 8")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 2:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_zAverage_dis":"17", "2d_z_qw":"13"}
        bias_list = {"2d_zAverage_dis":"17"}
        bias_list = {"2d_zAverage_dis":"20", "56_z_dis":"18"}
        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_start_extended_combined_2"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_start_extended_combined_may19"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["300", "335", "373"]}
        temp_dic = {"_280-350":["335", "373", "417"]}
        temp_dic = {"_280-350":["335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=6, chosen_mode=9)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    # nsample = len(folder_list)*2500
                    nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 28 --force 8")
                    cd("..")
                cd("..")
        cd("..")


if args.day == "aug04":
    if args.mode == 1:
        do("run.py -n 10 --restart 1 -s 9 -m 3 abeta42_12/")
if args.day == "jun21":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_zAverage_dis":"17", "2d_z_qw":"13"}
        bias_list = {"2d_zAverage_dis":"17"}
        # bias_list = {"1d_dis56":"9", "2d_dis56_z56":"14", "1d_z56":"11"}
        bias_list = {"1d_dis56":"91", "1d_h5":"92", "1d_h6":"93", "1d_h56":"94", "2d_dis56_h56":"18"}
        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_start_extended_dis56_z56"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_start_extended_combined_may19"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["300", "335", "373"]}
        temp_dic = {"_280-350":["335", "373", "417"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=6, chosen_mode=12)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    # nsample = len(folder_list)*2500
                    nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 27 --force 9")
                    cd("..")
                cd("..")
        cd("..")
def pick_structure_generate_show_script(n=2):
    with open("show.pml", "w") as f:
        for structure_index in range(0, n):
            f.write("load structure_%s.pdb\n" % structure_index)
            f.write("cealign structure_0, structure_%s\n" % structure_index)
            f.write("spectrum count, rainbow_rev, structure_%s, byres=1\n" % structure_index)
        # f.write("hide lines, all\n")
        # f.write("show cartoon, all\n")
        # f.write("hide nonbonded, all\n")

if args.day == "jun23":
    if args.mode == 2:
        cmd_pre = "python2 ~/opt/script/BuildAllAtomsFromLammps.py"
        # location_pre = "/Users/weilu/Research/server/apr_2018/sixth/rg_0.15_lipid_1.0_mem_1_go_0.8/simulation"

        # location_pre = "/Users/weilu/Research/server/may_2018/second_long/simulation"
        # cmd = cmd_pre + " " + location + " structure_2 4080 -seq ~/opt/pulling/2xov.seq"
        pick_list = ["low_e_jun23"]
        for picked in pick_list:
            do(f"mkdir {picked}")
            cd(picked)
            tt = pd.read_csv(f"/scratch/wl45/jun_2018/{picked}.csv", index_col=0)
            sample = tt.reset_index(drop=True)
            # sample["Frame"] = ((sample["Step"] - 2e7*rerun)/4000).astype("int")
            sample["rerun"] = (sample["Step"] // 2e7).astype(int)
            sample["Frame"] = ((sample["Step"] % 2e7)/4000).astype("int")
            for index, row in sample.iterrows():
                BiasTo = row["BiasTo"]
                Run = row["Run"]
                Frame = row["Frame"]
                rerun = row["rerun"]
                print(BiasTo, Run, Frame)

                # try:
                location_pre = "/scratch/wl45/may_2018/second/simulation"
                location = location_pre + f"/dis_{BiasTo}/{rerun}/dump.lammpstrj.{int(Run)}"
                cmd = cmd_pre + " " + location + f" structure_{index} {int(Frame)} -seq ~/opt/pulling/2xov.seq"
                hasProblem = do(cmd)
                if hasProblem == 0:
                    continue
                print(do(cmd))
                    # print(cmd)
                # except IOError:
                # print("-----------hi------------")
                location_pre = "/scratch/wl45/may_2018/second_long/simulation"
                location = location_pre + f"/dis_{BiasTo}/{rerun}/dump.lammpstrj.{int(Run)}"
                cmd2 = cmd_pre + " " + location + f" structure_{index} {int(Frame)} -seq ~/opt/pulling/2xov.seq"
                do(cmd2)
                    # print(cmd2)
            pick_structure_generate_show_script(n=len(sample))
            cd("..")

if args.day == "jun17":
    if args.mode == 2:
        cmd_pre = "python2 ~/opt/script/BuildAllAtomsFromLammps.py"
        # location_pre = "/Users/weilu/Research/server/apr_2018/sixth/rg_0.15_lipid_1.0_mem_1_go_0.8/simulation"

        # location_pre = "/Users/weilu/Research/server/may_2018/second_long/simulation"
        # cmd = cmd_pre + " " + location + " structure_2 4080 -seq ~/opt/pulling/2xov.seq"
        pick_list = ["low_e_jun01_h56", "low_e_jun01_h34", "low_e_jun01_h12", "low_e_jun01_out", "low_e_jun01_pre",
                        "low_e_jun01_transition", "low_e_jun01_post_transition"]
        pick_list = ["low_e_path1", "low_e_path2"]
        for picked in pick_list:
            do(f"mkdir {picked}")
            cd(picked)
            tt = pd.read_csv(f"/scratch/wl45/jun_2018/{picked}.csv", index_col=0)
            sample = tt.reset_index(drop=True)
            # sample["Frame"] = ((sample["Step"] - 2e7*rerun)/4000).astype("int")
            sample["rerun"] = (sample["Step"] // 2e7).astype(int)
            sample["Frame"] = ((sample["Step"] % 2e7)/4000).astype("int")
            for index, row in sample.iterrows():
                BiasTo = row["BiasTo"]
                Run = row["Run"]
                Frame = row["Frame"]
                rerun = row["rerun"]
                print(BiasTo, Run, Frame)

                # try:
                location_pre = "/scratch/wl45/may_2018/second/simulation"
                location = location_pre + f"/dis_{BiasTo}/{rerun}/dump.lammpstrj.{int(Run)}"
                cmd = cmd_pre + " " + location + f" structure_{index} {int(Frame)} -seq ~/opt/pulling/2xov.seq"
                hasProblem = do(cmd)
                if hasProblem == 0:
                    continue
                print(do(cmd))
                    # print(cmd)
                # except IOError:
                # print("-----------hi------------")
                location_pre = "/scratch/wl45/may_2018/second_long/simulation"
                location = location_pre + f"/dis_{BiasTo}/{rerun}/dump.lammpstrj.{int(Run)}"
                cmd2 = cmd_pre + " " + location + f" structure_{index} {int(Frame)} -seq ~/opt/pulling/2xov.seq"
                do(cmd2)
                    # print(cmd2)
            pick_structure_generate_show_script(n=len(sample))
            cd("..")

        # tt = pd.read_csv("/scratch/wl45/jun_2018/low_e_jun01_h56.csv", index_col=0)
        # tt = pd.read_csv("/scratch/wl45/jun_2018/low_e_jun01_h34.csv", index_col=0)
        # tt = pd.read_csv("/scratch/wl45/jun_2018/low_e_jun01_h12.csv", index_col=0)
        # tt = pd.read_csv("/scratch/wl45/jun_2018/low_e_jun01_out.csv", index_col=0)
        # tt = pd.read_csv("/scratch/wl45/jun_2018/low_e_jun01_pre.csv", index_col=0)
        # tt = pd.read_csv("/scratch/wl45/jun_2018/low_e_jun01_transition.csv", index_col=0)
        # tt = pd.read_csv("/scratch/wl45/jun_2018/low_e_jun01_post_transition.csv", index_col=0)
        # # rerun = 1
        # # sample = tt.sample(5).reset_index(drop=True)
        # ample = tt.reset_index(drop=True)
        # # sample["Frame"] = ((sample["Step"] - 2e7*rerun)/4000).astype("int")
        # sample["rerun"] = (sample["Step"] // 2e7).astype(int)
        # sample["Frame"] = ((sample["Step"] % 2e7)/4000).astype("int")
        # for index, row in sample.iterrows():
        #     BiasTo = row["BiasTo"]
        #     Run = row["Run"]
        #     Frame = row["Frame"]
        #     rerun = row["rerun"]
        #     print(BiasTo, Run, Frame)

        #     location = location_pre + f"/dis_{BiasTo}/{rerun}/dump.lammpstrj.{int(Run)}"
        #     cmd = cmd_pre + " " + location + f" structure_{index} {int(Frame)} -seq ~/opt/pulling/2xov.seq"
        #     print(cmd)
        #     do(cmd)
        # pick_structure_generate_show_script(n=len(sample))
    if args.mode == 1:
        data = pd.read_feather("/scratch/wl45/may_2018/03_week/all_data_folder/second_start_extended_combined_may19.feather")
        data = data.reset_index(drop=True)
        # data["BiasedEnergy"] = data["TotalE"] + 0.2*data["AMH_4H"]
        data["BiasedEnergy"] = data["Lipid"] + data["Rg"] + data["Membrane"] + data["AMH-Go"] + 0.2*data["AMH_4H"]
        data["BiasEnergy"] = 0.02 * (data["BiasTo"] - data["DisReal"])**2
        data["Energy_with_all_bias"] = data["BiasEnergy"] + data["BiasedEnergy"]

        t_pos = data.query("TempT == 373 and DisReal > 52 and DisReal < 57 and z_average > -4 and z_average < 0").reset_index(drop=True)
        chosen = t_pos.query("Lipid1 < -0.5").sort_values("Energy_with_all_bias").head(n=20)
        # chosen = t_pos.sort_values("Energy_with_all_bias").head(n=20)
        chosen.to_csv("low_e_jun01_pre.csv")

        t_pos = data.query("TempT == 373 and DisReal > 57 and DisReal <63 and z_average > -5 and z_average < -2").reset_index(drop=True)
        chosen = t_pos.query("Lipid1 < -0.5 and Lipid10 < -0.5").sort_values("Energy_with_all_bias").head(n=20)
        # chosen = t_pos.sort_values("Energy_with_all_bias").head(n=20)
        chosen.to_csv("low_e_jun01_transition.csv")

        t_pos = data.query("TempT == 373 and DisReal > 63 and DisReal <72 and z_average > -6 and z_average < -3").reset_index(drop=True)
        chosen = t_pos.query("Lipid1 < -0.5").sort_values("Energy_with_all_bias").head(n=20)
        # chosen = t_pos.sort_values("Energy_with_all_bias").head(n=20)
        chosen.to_csv("low_e_jun01_post_transition.csv")

        t_pos = data.query("TempT == 373 and DisReal > 80 and DisReal < 100 and z_average > -8 and z_average < -4").reset_index(drop=True)
        chosen = t_pos.query("Lipid1 < -0.5").sort_values("Energy_with_all_bias").head(n=20)
        chosen.to_csv("low_e_jun01_h56.csv")

        t_pos = data.query("TempT == 373 and DisReal > 140 and DisReal < 180 and z_average > -14 and z_average < -8").reset_index(drop=True)
        chosen = t_pos.sort_values("Energy_with_all_bias").head(n=20)
        chosen.to_csv("low_e_jun01_h34.csv")

        t_pos = data.query("TempT == 373 and DisReal > 220 and DisReal < 250 and z_average > -14 and z_average < -10").reset_index(drop=True)
        chosen = t_pos.sort_values("Energy_with_all_bias").head(n=20)
        chosen.to_csv("low_e_jun01_h12.csv")

        t_pos = data.query("TempT == 373 and DisReal > 260 and z_average < -16").reset_index(drop=True)
        chosen = t_pos.sort_values("Energy_with_all_bias").head(n=20)
        chosen.to_csv("low_e_jun01_out.csv")
if args.day == "jun12":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_zAverage_dis":"17", "2d_z_qw":"13"}
        bias_list = {"2d_zAverage_dis":"17"}
        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_start_extended"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_start_extended_combined_may19"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["300", "335", "373"]}
        temp_dic = {"_280-350":["335", "373", "417"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=9)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    # nsample = len(folder_list)*2500
                    nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 21 --force 8")
                    cd("..")
                cd("..")
        cd("..")


if args.day == "jun10":
    protein_list = ["T0833", "T0815", "T0803", "T0766"]
    for protein in protein_list:
        # do(f"mkdir {protein}")
        # cd(protein)
        do(f"scp wl45@davinci.rice.edu:/work/cms16/xl23/shared/IAAWSEM/AWSEM_HO_Results/protein_pool/06042018/{protein.lower()}/iter/post-processing/Qw.out {protein}/")
        do(f"scp wl45@davinci.rice.edu:/work/cms16/xl23/shared/IAAWSEM/AWSEM_HO_Results/protein_pool/06042018/{protein.lower()}/iter/rmsd/fromScwrl/totrmsd-angstrom.xvg {protein}/")
        # do(f"scp -r wl45@davinci.rice.edu:/work/cms16/xl23/shared/IAAWSEM/MC_DATA_28Feb2018/{protein}/AWSEM_energy_update AWSEM_energy")
        # do(f"scp -r wl45@davinci.rice.edu:/work/cms16/xl23/shared/IAAWSEM/MC_DATA_28Feb2018/{protein}/lowTstructure_update lowTstructure")
        # cd("..")

if args.day == "jun09":
    a = pd.read_csv("../selected.csv")
    for name, one in a.groupby("Name"):
        ii = 0
        for i, tmp in one.iterrows():
            print(i, name, tmp["Step"])
            do(f"cp ../{name}/lowTstructure/lowTstructure{tmp['Step']}.pdb {name}_{ii}.pdb")
            ii += 1

if args.day == "jun08":
    protein_list = ["T0833", "T0815", "T0803", "T0766"]
    for protein in protein_list:
        do(f"mkdir {protein}")
        cd(protein)
        do(f"scp -r wl45@davinci.rice.edu:/work/cms16/xl23/shared/IAAWSEM/MC_DATA_28Feb2018/{protein}/AWSEM_energy_update AWSEM_energy")
        do(f"scp -r wl45@davinci.rice.edu:/work/cms16/xl23/shared/IAAWSEM/MC_DATA_28Feb2018/{protein}/lowTstructure_update lowTstructure")
        cd("..")
if args.day == "may24":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_zAverage_dis":"17"}
        bias_list = {"2d_zAverage_dis":"17", "2d_z_qw":"13"}
        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        # freeEnergy_folder = f"second_start_topology"
        freeEnergy_folder = "enhance_go"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        # folder_list = [f"second_toplogy_may21"]
        folder_list = [f"second_start_extended_combined_may19"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["300", "335", "373"]}
        temp_dic = {"_280-350":["335", "373", "417"]}
        temp_dic = {"_280-350":["335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=11)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    # nsample = len(folder_list)*2500
                    nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 21 --force 7")
                    cd("..")
                cd("..")
        cd("..")

if args.day == "casp":
    # name = "T0953S2"
    # name = "T0954"
    # name = "T0955"
    # name = "T0956"
    # name = "T0957S1"
    # name = "T0957S2"
    name = "T0958"
    # name = "T0960"
    if args.mode == 1:
        # name = "T0953S1"
        do(f"scp -r wl45@davinci.rice.edu:/scratch/xl23/casp13/human/result/{name}/post-processing {name}")
        my_from = f"/scratch/mc70/CASP13/{name}/run1/"
        # my_from = f"wl45@davinci.rice.edu:/scratch/mc70/CASP13/{name}/run1/"
        my_to = f"{name}/awsem"
        cmd = "rsync -a --exclude='dump.lammpstrj' --exclude='slurm-*' --exclude='run.pdb' {} {}".format(my_from, my_to)
        print(cmd)
        os.system(cmd)
        do(f"cp casp.in {name}/awsem/")
        cd(f"{name}/awsem")
        cmd = "tail -n 3 ../model.1/lowTstructure/lowTstructure0.pdb | head -n 1"
        line = getFromTerminal(cmd)
        size = int(line.split()[4])
        print(line)
        print(size)
        alpha_carbons = " ".join([str(i) for i in list(range(1, size*3+1, 3))])
        beta_atoms = " ".join([str(i) for i in list(range(3, size*3+1, 3))])
        oxygens = " ".join([str(i) for i in list(range(2, size*3+1, 3))])
        with fileinput.FileInput("casp.in", inplace=True, backup='.bak') as file:
            for line in file:
                tmp = line.replace("ALPHA_CARBONS", alpha_carbons)
                tmp = tmp.replace("BETA_ATOMS", beta_atoms)
                tmp = tmp.replace("OXYGENS", oxygens)
                # tmp = BETA_ATOMS
                print(tmp, end='')
        cd("..")
    if args.mode == 2:
        cd(name)
        # i = 2
        queue = "interactive"
        # queue = "ctbp"
        # queue = "commons"
        for i in range(1,4):
            # time.sleep(100)
            cd(f"model.{i}")
            do("mkdir awsem_energy")
            cd("awsem_energy")
            simulation_list = glob.glob(f"../lowTstructure/lowTstructure*.pdb")
            print(len(simulation_list))
            for i in range(len(simulation_list)):
                do(f"cp -r ../../awsem awsem{i}")
                cd(f"awsem{i}")
                do(f"~/opt/script/PdbCoords2Lammps.sh ../../lowTstructure/lowTstructure{i} temp")
                quick = quick_template_slurm.format("/scratch/wl45/lmp_serial_nots < casp.in")
                if queue == "interactive":
                    quick = quick.replace("--time=01:30:00", "--time=00:30:00")
                    quick = quick.replace("#SBATCH --account=ctbp-common", "")
                    quick = quick.replace("ctbp-common", "interactive")
                if queue == "commons":
                    quick = quick.replace("ctbp-common", "commons")
                with open("run.slurm", "w") as f:
                    f.write(quick)
                do("sbatch run.slurm")
                cd("..")
                # do("/scratch/wl45/lmp_serial_nots < casp.in")
                # do(f"mv energy.log energy{i}.log")
                # do(f"tail -n 1 energy{i}.log >> awsem.log")
            cd("../..")
    if args.mode == 3:
        cd(name)
        # i = 1
        for i in range(1,4):
            cd(f"model.{i}/awsem_energy")
            do("rm awsem.log")
            simulation_list = glob.glob(f"../lowTstructure/lowTstructure*.pdb")
            n = len(simulation_list)
            print(len(simulation_list))
            for i in range(n):
                do(f"tail -n 1 awsem{i}/energy.log >> awsem.log")
            cd("../..")
if args.day == "may22":
    if args.mode == 5:
        cmd = "tail -n 3 lowTstructure/lowTstructure0.pdb | head -n 1"
        line = getFromTerminal(cmd)
        size = int(line.split()[4])
        print(line)
        print(size)
        alpha_carbons = " ".join([str(i) for i in list(range(1, size*3+1, 3))])
        beta_atoms = " ".join([str(i) for i in list(range(3, size*3+1, 3))])
        oxygens = " ".join([str(i) for i in list(range(2, size*3+1, 3))])
        # print("hi")
        # replace("casp.in", "data.casp", "data.temp")
        # replace("casp.in", "velocity", "# velocity")
    if args.mode == 4:
        do("scp -r wl45@davinci.rice.edu:/scratch/xl23/casp13/human/result/T0950 .")
    if args.mode == 3:
        cd("model.3/awsem_energy")
        do("rm awsem.log")
        simulation_list = glob.glob(f"../lowTstructure/lowTstructure*.pdb")
        n = len(simulation_list)
        print(len(simulation_list))
        for i in range(n):
            do(f"tail -n 1 awsem{i}/energy.log >> awsem.log")
    if args.mode == 2:
        do("mkdir awsem_energy")
        cd("awsem_energy")
        simulation_list = glob.glob(f"../lowTstructure/lowTstructure*.pdb")
        print(len(simulation_list))
        for i in range(len(simulation_list)):
            do(f"cp -r ../../awsem awsem{i}")
            cd(f"awsem{i}")
            do(f"~/opt/script/PdbCoords2Lammps.sh ../../lowTstructure/lowTstructure{i} temp")
            with open("run.slurm", "w") as f:
                f.write(quick_template_slurm.format("/scratch/wl45/lmp_serial_nots < casp.in"))
            do("sbatch run.slurm")
            cd("..")
            # do("/scratch/wl45/lmp_serial_nots < casp.in")
            # do(f"mv energy.log energy{i}.log")
            # do(f"tail -n 1 energy{i}.log >> awsem.log")
    if args.mode == 1:
        for i in range(153):
            do(f"python ~/opt/small_script/CalcQValueFromTwoPdb_2.py lowTstructure10.pdb lowTstructure{i}.pdb | tail -n 1 >> Q10_all")
if args.day == "may21":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_zAverage_dis":"17"}
        bias_list = {"2d_zAverage_dis":"17", "2d_z_qw":"13"}
        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        # freeEnergy_folder = f"second_start_topology"
        freeEnergy_folder = "quick"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        # folder_list = [f"second_toplogy_may21"]
        folder_list = [f"second_start_extended_combined_may19"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["300", "335", "373"]}
        temp_dic = {"_280-350":["335", "373", "417"]}
        temp_dic = {"_280-350":["335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=10)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    # nsample = len(folder_list)*2500
                    nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 21 --force 6")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "may20":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_zAverage_dis":"17", "2d_z_qw":"13"}
        bias_list = {"2d_zAverage_dis":"17"}
        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_start_extended_combined_2"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_start_extended_combined_may19"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["300", "335", "373"]}
        temp_dic = {"_280-350":["335", "373", "417"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=9)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    # nsample = len(folder_list)*2500
                    nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 21 --force 5")
                    cd("..")
                cd("..")
        cd("..")

if args.day == "may19":
    if args.mode == 2:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_zAverage_dis":"17", "2d_z_qw":"13"}
        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_start_extended_combined"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_start_extended_combined_may19"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["300", "335", "373"]}
        temp_dic = {"_280-350":["335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=9)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    # nsample = len(folder_list)*2500
                    nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 21 --force 4")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13", "2d_zAverage_dis":"17"}
        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_may19"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=8)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    # nsample = len(folder_list)*2500
                    nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 21 --force 4")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "may17":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = 7
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_quick"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_may17"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=7)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    # nsample = len(folder_list)*2500
                    nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 22 --force 0")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "may15":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = 7
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_rerun3_with_goEnergyrerun_7_14_May_124103"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=6)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 22 --force 3")
                    cd("..")
                cd("..")
        cd("..")

if args.day == "may14":
    if args.mode == 6:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_short_and_long"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_rerun3_with_goEnergyrerun_7_14_May_124103", "second_longrerun_5"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=6)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 5:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"16", "2d_z_h56":"15", "1d_qw":"10"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_combine"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_rerun3_with_goEnergyrerun_7_14_May_124103", "second_start_topology_rerun2"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=6)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 4:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"16", "2d_z_h56":"15", "1d_qw":"10"}  # z and Dis_h56
        i = 5
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_long"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_longrerun_5_14_May_155146"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=6)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 3:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"16"}  # z and Dis_h56
        i = 7
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_enhance_n2"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_rerun3_with_goEnergyrerun_7_14_May_124103"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=6)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 22 --force 0")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 2:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"16", "2d_z_h56":"15", "1d_qw":"10"}  # z and Dis_h56
        i = 7
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_enhance_n"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_rerun3_with_goEnergyrerun_7_14_May_124103"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=6)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 1:
        goEnergy = False
        rerun = 2
        end = 1
        cwd = os.getcwd()
        print(cwd)
        pre = '/'.join(cwd.split("/")[:-2]) + "/"
        print(pre)
        # exit()
        # pre = "/scratch/wl45/apr_2018/sixth/"
        data_folder = "/scratch/wl45/may_2018/02_week/all_data_folder/"
        # folder_list = ["rg_0.15_lipid_1.0_mem_1_go_0.8_long"]
        folder_list = [cwd.split("/")[-2]]
        # label = "sixth_long"
        label = args.label
        with open(label, "w") as f:
            f.write("test\n")
        # cd("simulation")
        # exit()
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]

        for attempt in range(2):
            try:
                process_complete_temper_data_3(pre, data_folder, folder_list, rerun=rerun, end=end, average_z=True, disReal=True, dis_h56=True, localQ=False, goEnergy=goEnergy, label=label)
                break
            except FileNotFoundError:
                bias = "dis"
                simulation_list = glob.glob(f"{bias}_*")
                # simulation_list = ['dis_30.0']
                # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
                print(simulation_list)
                for dis in simulation_list:
                    print(dis)
                    cd(dis)
                    i = rerun
                    i_plus_one = i +1
                    if os.path.exists(f"log{i}"):
                        do(f"mv log{i} back_log{i}")  # in case of override
                    do(f"mkdir -p log{i}")
                    do(f"cp log.lammps log{i}/")
                    cd("..")

if args.day == "may13":
    if args.mode == 5:
        cmd = "tail -n 3 lowTstructure/lowTstructure0.pdb | head -n 1"
        line = getFromTerminal(cmd)
        size = int(line.split()[4])
        print(line)
        print(size)
        alpha_carbons = " ".join([str(i) for i in list(range(1, size*3+1, 3))])
        beta_atoms = " ".join([str(i) for i in list(range(3, size*3+1, 3))])
        oxygens = " ".join([str(i) for i in list(range(2, size*3+1, 3))])
        # print("hi")
        # replace("casp.in", "data.casp", "data.temp")
        # replace("casp.in", "velocity", "# velocity")
    if args.mode == 4:
        do("scp -r wl45@davinci.rice.edu:/scratch/xl23/casp13/human/result/T0950 .")
    if args.mode == 3:
        cd("model.3/awsem_energy")
        do("rm awsem.log")
        simulation_list = glob.glob(f"../lowTstructure/lowTstructure*.pdb")
        n = len(simulation_list)
        print(len(simulation_list))
        for i in range(n):
            do(f"tail -n 1 awsem{i}/energy.log >> awsem.log")
    if args.mode == 2:
        do("mkdir awsem_energy")
        cd("awsem_energy")
        simulation_list = glob.glob(f"../lowTstructure/lowTstructure*.pdb")
        print(len(simulation_list))
        for i in range(len(simulation_list)):
            do(f"cp -r ../../awsem awsem{i}")
            cd(f"awsem{i}")
            do(f"~/opt/script/PdbCoords2Lammps.sh ../../lowTstructure/lowTstructure{i} temp")
            with open("run.slurm", "w") as f:
                f.write(quick_template_slurm.format("/scratch/wl45/lmp_serial_nots < casp.in"))
            do("sbatch run.slurm")
            cd("..")
            # do("/scratch/wl45/lmp_serial_nots < casp.in")
            # do(f"mv energy.log energy{i}.log")
            # do(f"tail -n 1 energy{i}.log >> awsem.log")
    if args.mode == 1:
        for i in range(153):
            do(f"python ~/opt/small_script/CalcQValueFromTwoPdb_2.py lowTstructure10.pdb lowTstructure{i}.pdb | tail -n 1 >> Q10_all")
if args.day == "may12":
    if args.mode == 1:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_30.0']
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            cd(dis)
            cd("2")
            # cd("3")
            rerun(extra="Go", mode=2)
            cd("../..")
if args.day == "may11":
    if args.mode == 2:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"13", "2d_z_h56":"15", "1d_qw":"10"}  # z and Dis_h56
        i = 3
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_long"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_longrerun_3_11_May_131422"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=4)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"13", "2d_z_h56":"15", "1d_qw":"10", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_combine"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_rerun3rerun_7_10_May_222655", "second_start_topologyrerun_1_11_May_133146"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=4)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")

if args.day == "may10":
    if args.mode == 4:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"13", "2d_z_h56":"15", "1d_qw":"10"}  # z and Dis_h56
        i = 7
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_rerun3"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_rerun3rerun_7_10_May_222655"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=4)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 3:
        do("rm awsem.log")
        simulation_list = glob.glob(f"../lowTstructure/lowTstructure*.pdb")
        n = len(simulation_list)
        print(len(simulation_list))
        for i in range(n):
            do(f"tail -n 1 awsem{i}/energy.log >> awsem.log")
    if args.mode == 2:
        simulation_list = glob.glob(f"../lowTstructure/lowTstructure*.pdb")
        print(len(simulation_list))
        for i in range(len(simulation_list)):
            do(f"cp -r ../../awsem awsem{i}")
            cd(f"awsem{i}")
            do(f"~/opt/script/PdbCoords2Lammps.sh ../../lowTstructure/lowTstructure{i} temp")
            with open("run.slurm", "w") as f:
                f.write(quick_template_slurm.format("/scratch/wl45/lmp_serial_nots < casp.in"))
            do("sbatch run.slurm")
            cd("..")
            # do("/scratch/wl45/lmp_serial_nots < casp.in")
            # do(f"mv energy.log energy{i}.log")
            # do(f"tail -n 1 energy{i}.log >> awsem.log")
    if args.mode == 1:
        for i in range(153):
            do(f"python ~/opt/small_script/CalcQValueFromTwoPdb_2.py lowTstructure10.pdb lowTstructure{i}.pdb | tail -n 1 >> Q10_all")
if args.day == "may09":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        # bias_list = {"2d_z_qw":"13", "2d_z_h56":"15", "1d_qw":"10"}  # z and Dis_h56
        bias_list = {"2d_z_qw":"16", "2d_z_h56":"15", "1d_qw":"10"}  # z and Dis_h56
        i = 5
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_rerun2"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_rerun2rerun_5_09_May_225324"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=4)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "may08":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"13", "2d_z_h56":"15", "1d_qw":"10"}  # z and Dis_h56
        i = 1
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"thrid_start_extended"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"third_start_extendedrerun_1_08_May_135309"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=4)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")

if args.day == "may07":
    if args.mode == 1:
        rerun(extra="DMPC", mode=1)
    if args.mode == 2:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"13", "2d_z_h56":"15", "1d_qw":"10"}  # z and Dis_h56
        i = 3
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_rerun1"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_rerun1rerun_3_07_May_154957"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=4)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 3:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"13", "2d_z_h56":"15", "1d_qw":"10"}  # z and Dis_h56
        i = 1
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"first_long"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"first_longrerun_1_07_May_170020"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=4)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.5, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "may06":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"13", "2d_z_h56":"15", "1d_qw":"10"}  # z and Dis_h56
        i = 1
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_start_native"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_start_nativererun_1_06_May_223049"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=4)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")

if args.day == "may05":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"13", "2d_z_h56":"15", "1d_qw":"10"}  # z and Dis_h56
        i = 1
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"secondrerun_1_05_May_155022"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["300", "335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=4)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "may03":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"13", "2d_z_h56":"15", "1d_qw":"10"} # z and Dis_h56
        i = 3
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        freeEnergy_folder = f"sixth_long/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"sixth_longrerun_3_18_Apr_205152"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "300", "325", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=1, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=3)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "may01":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"13", "2d_z_h56":"15", "1d_qw":"10"} # z and Dis_h56
        i = 1
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        freeEnergy_folder = f"seventh/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"seventhrerun_1_01_May_234419"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["280", "300", "325", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=5) #chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.5, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "apr30":
    if args.mode == 3:
        print("how Constant force refolding")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0.8, 1, 0.6, 0.4, 0.2, 0]
        pressure_list = [1, 0.5, 0]
        # force_ramp_rate_list=[0.25]
        force_ramp_rate_list=[0.5]
        temperature_list=[300]
        # memb_k_list = [0.8, 1]
        memb_k_list = [0.8]
        rg_list = [0.1]
        change_list = [4]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        # force_list = [0.22, 0.25, 0.28]
        # force_list = [0.0, 0.01]
        # force_list = [0.0143472]
        # force_list = ["ramp"]
        # force_list = [1, 2]
        # force_list = [3, 4]
        # force_list = [5, 6]
        # force_list = [0, 0.02, 0.05]
        force_list = [0.02]
        repeat = 20
        # change_list = [15]
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        change_list=change_list,
                        commons=0,simulation_base_steps=4e7)
    if args.mode == 2:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"13", "2d_z_h56":"15", "1d_qw":"10"} # z and Dis_h56
        i = 3
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        freeEnergy_folder = f"twelve/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"twelvererun_3_01_May_025941"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["280", "300", "325", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=4) #chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"13", "2d_z_h56":"15"} # z and Dis_h56
        i = 1
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        freeEnergy_folder = f"third/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"fifth_3rerun_1_30_Apr_150140"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["280", "300", "325", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=4) #chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.5, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")





if args.day == "apr27":
    if args.mode == 1:
        rerun(extra="DMPC", mode=1)
if args.day == "apr24":
    if args.mode == 5:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        bias_list = {"2d_z_qw":"13", "2d_z_h56":"15"} # z and Dis_h56
        i = 1
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        freeEnergy_folder = f"third/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"third_with_h56rerun_1_24_Apr_230220"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["280", "300", "325", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=4) #chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 4:
        print("how Constant force unfolding")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0.8, 1, 0.6, 0.4, 0.2, 0]
        # force_ramp_rate_list=[0.25]
        force_ramp_rate_list=[0.5]
        temperature_list=[300]
        # memb_k_list = [0.8, 1]
        memb_k_list = [0.8]
        rg_list = [0.1]
        change_list = [4]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        # force_list = [0.22, 0.25, 0.28]
        # force_list = [0.0, 0.01]
        # force_list = [0.0143472]
        # force_list = ["ramp"]
        # force_list = [1, 2]
        # force_list = [3, 4]
        # force_list = [5, 6]
        # force_list = [0, 0.02, 0.05]
        force_list = [0]
        repeat = 10
        # change_list = [15]
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        change_list=change_list,
                        commons=0,simulation_base_steps=4e7)
    if args.mode == 3:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        i = 1
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        freeEnergy_folder = f"eleventh/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"tenth_k0.2rerun_1_22_Apr_164742"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["280", "300", "325", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=3)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 2:
        read_variable_folder(".")
    if args.mode == 1:
        folder_list = glob.glob("*_")
        queue = 1
        i = "0"
        for folder in folder_list:
            cd(folder)
            cd("simulation")
            compute_disReal(temper=False, bias="", sim_list=[i], queue=queue)
            compute_completeZ(temper=False, bias="", sim_list=[i], queue=queue)
            cd("../../")
if args.day == "apr22":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        i = 1
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        freeEnergy_folder = f"third_expectedEnergy/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"third_force_0.02rerun_1_16_Apr_143519"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["280", "300", "325", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=3)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 26 --force 1")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "apr21":
    if args.mode == 1:
        i = 1
        temp_list  = ["280", "300", "325", "350"]
        temp_list  = ["300"]
        temp_list  = ["350"]
        make_metadata_3(temps_list=temp_list,k=0.2, i=i)

if args.day == "apr19":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        i = 3
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        freeEnergy_folder = f"sixth_expectedEnergy/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        # folder_list = [f"sixth_longrerun_3_18_Apr_205152"]
        folder_list = [f"sixth_DMPCrerun_3_14_Apr_223431"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "300", "325", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=1, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=3)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 2")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 2:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        i = 3
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        freeEnergy_folder = f"sixth_long/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"sixth_longrerun_3_18_Apr_205152"]
        # folder_list = [f"sixth_DMPCrerun_3_14_Apr_223431"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "300", "325", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=1, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 3:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        i = 1
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        freeEnergy_folder = f"fourth/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"fourthrerun_1_20_Apr_001224"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["280", "300", "325", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")




if args.day == "apr18":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        i = 3
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        freeEnergy_folder = f"sixth_long/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"sixth_longrerun_3_18_Apr_205152"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "300", "325", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=1, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 2")
                    cd("..")
                cd("..")
        cd("..")

if args.day == "apr16":
    if args.mode == 2:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        i = 1
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        freeEnergy_folder = f"third/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"third_force_0.02rerun_1_16_Apr_143519"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["280", "300", "325", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 1:
        rerun = 0
        pre = "/scratch/wl45/apr_2018/third/"
        data_folder = "/scratch/wl45/apr_2018/02_week/all_data_folder/"
        folder_list = ["force_0.02_rg_0.15_lipid_1.0_mem_1_go_0.8"]
        label = "third_force_0.02"
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        for attempt in range(2):
            try:
                process_complete_temper_data_3(pre, data_folder, folder_list, rerun=rerun, average_z=True, disReal=True, localQ=False, label=label)
                break
            except FileNotFoundError:
                bias = "dis"
                simulation_list = glob.glob(f"{bias}_*")
                # simulation_list = ['dis_30.0']
                # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
                print(simulation_list)
                for dis in simulation_list:
                    print(dis)
                    cd(dis)
                    i = rerun
                    i_plus_one = i +1
                    if os.path.exists(f"log{i}"):
                        do(f"mv log{i} back_log{i}")  # in case of override
                    do(f"mkdir -p log{i}")
                    do(f"cp log.lammps log{i}/")
                    cd("..")
if args.day == "apr15":
    if args.mode == 1:
        print("how Constant force unfolding")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        # force_ramp_rate_list=[0.25]
        force_ramp_rate_list=[0.5]
        temperature_list=[300]
        # memb_k_list = [0.8, 1]
        memb_k_list = [0.8]
        rg_list = [0.1, 0.15]
        change_list = [1, 2, 4, 8]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        # force_list = [0.22, 0.25, 0.28]
        # force_list = [0.0, 0.01]
        # force_list = [0.0143472]
        # force_list = ["ramp"]
        # force_list = [1, 2]
        # force_list = [3, 4]
        # force_list = [5, 6]
        # force_list = [0, 0.02, 0.05]
        force_list = [0]
        repeat = 10
        # change_list = [15]
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        change_list=change_list,
                        commons=0,simulation_base_steps=4e7)

if args.day == "apr14":
    if args.mode == 2:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        i = 3
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_DMPC/"
        # freeEnergy_folder = "sixth_DOPC"
        # freeEnergy_folder = "sixth_POPC"
        freeEnergy_folder = "sixth_DPPC"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        # folder_list = [f"sixth_DOPCrerun_3_18_Apr_212102"]
        # folder_list = [f"sixth_POPCrerun_3_18_Apr_213447"]
        folder_list = [f"sixth_DPPCrerun_3_18_Apr_220421"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "300", "325", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=1, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 2")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 1:
        rerun = 1
        pre = "/scratch/wl45/apr_2018/first/"
        data_folder = "/scratch/wl45/apr_2018/02_week/all_data_folder/"
        folder_list = ["force_0.04_rg_0.15_lipid_1.0_mem_1_go_0.8"]
        label = "first_force_0.04"
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        for attempt in range(2):
            try:
                process_complete_temper_data_3(pre, data_folder, folder_list, rerun=rerun, average_z=True, disReal=True, localQ=False, label=label)
                break
            except FileNotFoundError:
                bias = "dis"
                simulation_list = glob.glob(f"{bias}_*")
                # simulation_list = ['dis_30.0']
                # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
                print(simulation_list)
                for dis in simulation_list:
                    print(dis)
                    cd(dis)
                    # i = rerun
                    # i_plus_one = i +1
                    # if os.path.exists(f"log{i}"):
                    #     do(f"mv log{i} back_log{i}")  # in case of override
                    # do(f"mkdir -p log{i}")
                    # do(f"cp log.lammps log{i}/")
                    # cd("..")
if args.day == "apr13":
    if args.mode == 4:
        pre = "/scratch/wl45/apr_2018/sixth/"
        data_folder = "/scratch/wl45/apr_2018/02_week/all_data_folder/"
        folder_list = ["rg_0.15_lipid_1.0_mem_1_go_0.8"]
        # label = "sixth_i235d"
        # label = "sixth_i255d"
        # label = "sixth_original"
        # label = "sixth_DMPC"
        # label = "sixth_DOPC"
        # label = "sixth_POPC"
        label = "sixth_DPPC"
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=1, average_z=True, disReal=True, dis_h56=True, localQ=False, label=label)
    if args.mode == 3:
        # extra="i235d"
        # extra="i255d"
        # extra = "original"
        # extra = "DMPC"
        # extra = "DOPC"
        # extra = "POPC"
        extra = "DPPC"
        bias = "dis_"
        simulation_list = glob.glob(f"{bias}*")
        print(simulation_list)
        # print(sim_list)
        for folder in simulation_list:
            print(folder)
            cd(folder)
            cd("1")
            # do("mkdir backup")
            # do("mv wham.* backup/")
            # do("mv energy.* backup/")
            # do("mv lipid.* backup/")
            # do("mv rgs.* backup/")
            for i in range(12):
                do(f"cp {extra}_{i}/wham.dat wham.{i}.dat")
                do(f"cp {extra}_{i}/energy.dat energy.{i}.dat")
                do(f"cp {extra}_{i}/lipid.dat lipid.{i}.dat")
                do(f"cp {extra}_{i}/rgs.dat rgs.{i}.dat")
            # rerun(extra="i255d")
            cd("../..")

    if args.mode == 2:
        # extra = "DMPC"
        # extra = "DOPC"
        # extra = "POPC"
        extra = "DPPC"
        bias = "dis_"
        simulation_list = glob.glob(f"{bias}*")
        print(simulation_list)
        # print(sim_list)
        for folder in simulation_list:
            print(folder)
            cd(folder)
            cd("1")
            rerun(extra=extra, mode=1)
            # rerun(extra="original")
            # rerun(extra="i255d")
            # rerun(extra="i235d")
            cd("../..")

if args.day == "apr12":
    if args.mode == 1:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_30.0']
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 0
            i_plus_one = i +1
            # do(f"mv log{i} back_log{i}")  # in case of override
            # do(f"mkdir -p log{i}")
            # do(f"cp log.lammps log{i}/")

            continueRunConvertion(n=12, rerun=i)
            do(f"mkdir {i_plus_one}")
            do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
if args.day == "apr10":
    if args.mode == 1:
        rerun = 0
        pre = "/scratch/wl45/apr_2018/second/"
        data_folder = "/scratch/wl45/apr_2018/01_week/all_data_folder/"
        folder_list = ["force_0.02_rg_0.1_lipid_1.0_mem_1_go_0.8"]
        label = "second_force_0.02"
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        for attempt in range(2):
            try:
                process_complete_temper_data_3(pre, data_folder, folder_list, rerun=rerun, average_z=True, disReal=True, localQ=False, label=label)
            except FileNotFoundError:
                bias = "dis"
                simulation_list = glob.glob(f"{bias}_*")
                # simulation_list = ['dis_30.0']
                # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
                print(simulation_list)
                for dis in simulation_list:
                    print(dis)
                    cd(dis)
                    i = rerun
                    i_plus_one = i +1
                    if os.path.exists(f"log{i}"):
                        do(f"mv log{i} back_log{i}")  # in case of override
                    do(f"mkdir -p log{i}")
                    do(f"cp log.lammps log{i}/")
                    cd("..")
if args.day == "apr09":
    if args.mode == 3:
        print("how Constant force refolding")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        force_ramp_rate_list=[0.25]
        temperature_list=[300]
        # memb_k_list = [0.8, 1]
        memb_k_list = [0.8]
        rg_list = [0.1]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        # force_list = [0.22, 0.25, 0.28]
        # force_list = [0.0, 0.01]
        # force_list = [0.0143472]
        # force_list = ["ramp"]
        # force_list = [1, 2]
        # force_list = [3, 4]
        # force_list = [5, 6]
        force_list = [0, 0.02, 0.05]
        repeat = 30
        change_list = [15]
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        change_list=change_list,
                        commons=0,simulation_base_steps=4e7)
    if args.mode == 2:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        i = 1
        freeEnergy_folder = f"second_freeEnergy_{i}/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"tenth_force_0.04rerun_3_06_Apr_143741"]
        folder_list = [f"second_force_0.04rerun_1_09_Apr_154731"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "310"]}

        # temp_dic = {"_280-350":["290", "300", "310"]}

        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 1:
        rerun = 0
        pre = "/scratch/wl45/apr_2018/second/"
        data_folder = "/scratch/wl45/apr_2018/01_week/all_data_folder/"
        folder_list = ["force_0.04_rg_0.15_lipid_1.0_mem_1_go_0.8"]
        label = "second_force_0.04"
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        for attempt in range(2):
            try:
                process_complete_temper_data_3(pre, data_folder, folder_list, rerun=rerun, average_z=True, disReal=True, localQ=False, label=label)
            except FileNotFoundError:
                bias = "dis"
                simulation_list = glob.glob(f"{bias}_*")
                # simulation_list = ['dis_30.0']
                # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
                print(simulation_list)
                for dis in simulation_list:
                    print(dis)
                    cd(dis)
                    i = rerun
                    i_plus_one = i +1
                    if os.path.exists(f"log{i}"):
                        do(f"mv log{i} back_log{i}")  # in case of override
                    do(f"mkdir -p log{i}")
                    do(f"cp log.lammps log{i}/")
                    cd("..")
if args.day == "apr08":
    if args.mode == 2:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        i = 1
        freeEnergy_folder = f"first_freeEnergy_{i}/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"tenth_force_0.04rerun_3_06_Apr_143741"]
        folder_list = [f"first_force_0.04rerun_1_08_Apr_145204"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "310"]}

        # temp_dic = {"_280-350":["290", "300", "310"]}

        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 1:
        rerun = 0
        pre = "/scratch/wl45/apr_2018/first/"
        data_folder = "/scratch/wl45/apr_2018/01_week/all_data_folder/"
        folder_list = ["force_0.04_rg_0.15_lipid_1.0_mem_1_go_0.8"]
        label = "first_force_0.04"
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        for attempt in range(2):
            try:
                process_complete_temper_data_3(pre, data_folder, folder_list, rerun=rerun, average_z=True, disReal=True, localQ=False, label=label)
            except FileNotFoundError:
                bias = "dis"
                simulation_list = glob.glob(f"{bias}_*")
                # simulation_list = ['dis_30.0']
                # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
                print(simulation_list)
                for dis in simulation_list:
                    print(dis)
                    cd(dis)
                    i = rerun
                    i_plus_one = i +1
                    if os.path.exists(f"log{i}"):
                        do(f"mv log{i} back_log{i}")  # in case of override
                    do(f"mkdir -p log{i}")
                    do(f"cp log.lammps log{i}/")
                    cd("..")

if args.day == "apr06":
    if args.mode == 2:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_qw":"13"}
        i = 5
        freeEnergy_folder = f"tenth_freeEnergy_{i}/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"tenth_force_0.04rerun_3_06_Apr_143741"]
        folder_list = [f"tenth_force_0.04rerun_5_07_Apr_160809"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "310"]}

        # temp_dic = {"_280-350":["290", "300", "310"]}

        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 1:
        name = "tmhc2"
        # i = "0"
        i = "1"
        # i = "2"
        # i = "3"
        # let_compute_localQ(temper=True, bias="dis_", sim_list=[i], queue=1)
        compute_disReal(name="tmhc2", temper=True, bias="dis_", sim_list=[i], queue=0)
        # compute_disReal(temper=True, targetMode=1, bias="dis_", sim_list=[i], queue=1)
        compute_completeZ(name="tmhc2", temper=True, bias="dis_", sim_list=[i], queue=0)
if args.day == "apr05":
    if args.mode == 3:
        rerun = 0
        pre = "/scratch/wl45/apr_2018/TMHC2/"
        data_folder = "/scratch/wl45/apr_2018/01_week/all_data_folder/"
        folder_list = ["tmhc2_bias"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        for attempt in range(2):
            try:
                process_complete_temper_data_3(pre, data_folder, folder_list, rerun=rerun, average_z=True, disReal=True, localQ=False, label="tmhc2")
            except FileNotFoundError:
                bias = "dis"
                simulation_list = glob.glob(f"{bias}_*")
                # simulation_list = ['dis_30.0']
                # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
                print(simulation_list)
                for dis in simulation_list:
                    print(dis)
                    cd(dis)
                    i = rerun
                    i_plus_one = i +1
                    if os.path.exists(f"log{i}"):
                        do(f"mv log{i} back_log{i}")  # in case of override
                    do(f"mkdir -p log{i}")
                    do(f"cp log.lammps log{i}/")
                    cd("..")
    if args.mode == 2:
        print("how Constant force refolding")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        force_ramp_rate_list=[0.25]
        temperature_list=[300]
        # memb_k_list = [0.8, 1]
        memb_k_list = [0.8]
        rg_list = [0.15]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        # force_list = [0.22, 0.25, 0.28]
        # force_list = [0.0, 0.01]
        # force_list = [0.0143472]
        # force_list = ["ramp"]
        # force_list = [1, 2]
        # force_list = [3, 4]
        # force_list = [5, 6]
        force_list = [1, 2, 3, 4, 5, 8]
        repeat = 10
        change_list = [15]
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        change_list=change_list,
                        commons=0,simulation_base_steps=4e7)
    if args.mode == 1:
        rerun = 2
        pre = "/scratch/wl45/apr_2018/tenth/"
        data_folder = "/scratch/wl45/apr_2018/01_week/all_data_folder/"
        folder_list = ["force_0.04_rg_0.15_lipid_1.0_mem_1_go_0.8"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        for attempt in range(2):
            try:
                process_complete_temper_data_3(pre, data_folder, folder_list, rerun=rerun, average_z=True, disReal=True, localQ=False, label="tenth_force_0.04")
            except FileNotFoundError:
                bias = "dis"
                simulation_list = glob.glob(f"{bias}_*")
                # simulation_list = ['dis_30.0']
                # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
                print(simulation_list)
                for dis in simulation_list:
                    print(dis)
                    cd(dis)
                    i = rerun
                    i_plus_one = i +1
                    if os.path.exists(f"log{i}"):
                        do(f"mv log{i} back_log{i}")  # in case of override
                    do(f"mkdir -p log{i}")
                    do(f"cp log.lammps log{i}/")
                    cd("..")
if args.day == "apr04":
    if args.mode == 2:
        create_zim("tmhc2.seq", False)
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        i = 7
        freeEnergy_folder = f"ninth_freeEnergy_{i}_less_temp/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"ninth_force_0.06rerun_7_04_Apr_231330"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335"]}
        # temp_dic = {"_280-350":["280", "290", "300", "310"]}
        temp_dic = {"_280-350":["280", "290", "300"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "apr03":
    if args.mode == 6:
        temp_list = ["300"]
        i = 5
        make_metadata_3(temps_list=temp_list,k=0.02, i=i, biasHigh=85)
    if args.mode == 5:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        i = 3
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        freeEnergy_folder = f"sixth_orignal/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"sixth_originalrerun_3_04_Apr_153620"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "300", "325", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=1, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 2")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 1:
        create_zim("2xov.seq", False)
    if args.mode == 4:
        pre = "/scratch/wl45/apr_2018/"
        data_folder = "/scratch/wl45/apr_2018/01_week/all_data_folder/"
        folder_list = ["rg_0.15_lipid_1.0_mem_1_go_0.8"]
        # label = "sixth_i235d"
        # label = "sixth_i255d"
        label = "sixth_original"
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=1, average_z=True, disReal=True, dis_h56=True, localQ=False, label=label)
    if args.mode == 3:
        # extra="i235d"
        # extra="i255d"
        extra = "original"
        bias = "dis_"
        simulation_list = glob.glob(f"{bias}*")
        print(simulation_list)
        # print(sim_list)
        for folder in simulation_list:
            print(folder)
            cd(folder)
            cd("1")
            # do("mkdir backup")
            # do("mv wham.* backup/")
            # do("mv energy.* backup/")
            # do("mv lipid.* backup/")
            # do("mv rgs.* backup/")
            for i in range(12):
                do(f"cp {extra}_{i}/wham.dat wham.{i}.dat")
                do(f"cp {extra}_{i}/energy.dat energy.{i}.dat")
                do(f"cp {extra}_{i}/lipid.dat lipid.{i}.dat")
                do(f"cp {extra}_{i}/rgs.dat rgs.{i}.dat")
            # rerun(extra="i255d")
            cd("../..")


    if args.mode == 2:
        bias = "dis_"
        simulation_list = glob.glob(f"{bias}*")
        print(simulation_list)
        # print(sim_list)
        for folder in simulation_list:
            print(folder)
            cd(folder)
            cd("1")
            # rerun(extra="original")
            rerun(extra="i255d")
            # rerun(extra="i235d")
            cd("../..")

if args.day == "apr02":
    if args.mode == 3:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        i = 5
        freeEnergy_folder = f"ninth_freeEnergy_{i}_less_temp/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"ninth_force_0.06rerun_5_03_Apr_011519"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335"]}
        # temp_dic = {"_280-350":["280", "290", "300", "310"]}
        temp_dic = {"_280-350":["290", "300"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 2:
        rerun = 3
        pre = "/scratch/wl45/mar_2018/ninth/"
        data_folder = "/scratch/wl45/apr_2018/01_week/all_data_folder/"
        folder_list = ["force_0.06_rg_0.15_lipid_1.0_mem_1_go_0.8"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=rerun, average_z=True, disReal=True, localQ=False, label="ninth_force_0.06")
        # bias = "dis"
        # simulation_list = glob.glob(f"{bias}_*")
        # # simulation_list = ['dis_30.0']
        # # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        # print(simulation_list)
        # for dis in simulation_list:
        #     print(dis)
        #     cd(dis)
        #     i = rerun
        #     i_plus_one = i +1
        #     if os.path.exists(f"log{i}"):
        #         do(f"mv log{i} back_log{i}")  # in case of override
        #     do(f"mkdir -p log{i}")
        #     do(f"cp log.lammps log{i}/")
        #     cd("..")
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        i = 5
        freeEnergy_folder = f"ninth_freeEnergy_{i}/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"ninth_force_0.06rerun_5_03_Apr_011519"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335"]}
        # temp_dic = {"_280-350":["280", "290", "300", "310"]}
        # temp_dic = {"_280-350":["290", "300"]}
        temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "apr01":
    if args.mode == 6:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_30.0']
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 0
            i_plus_one = i +1
            # do(f"mv log{i} back_log{i}")  # in case of override
            # do(f"mkdir -p log{i}")
            # do(f"cp log.lammps log{i}/")

            continueRunConvertion(n=12, rerun=i, name="tmhc2")
            do(f"mkdir {i_plus_one}")
            name = "tmhc2"
            do(f"sed 's/{name}_{i}/{name}_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
    if args.mode == 5:
        temp_list = ["290"]
        i = 3
        make_metadata_3(temps_list=temp_list,k=0.02, i=i)
    if args.mode == 4:
        i = "0"
        # i = "1"
        # i = "2"
        # i = "3"
        # let_compute_localQ(temper=True, bias="dis_", sim_list=[i], queue=1)
        compute_disReal(temper=True, bias="dis_", sim_list=[i], queue=1)
        # compute_disReal(temper=True, targetMode=1, bias="dis_", sim_list=[i], queue=1)
        compute_completeZ(temper=True, bias="dis_", sim_list=[i], queue=1)
    if args.mode == 3:
        for i in range(12):
            cmd = f"python3 ~/opt/small_script/find_distance.py -m {i} -t 1"
            do(cmd)
    if args.mode == 2:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        i = 3
        freeEnergy_folder = f"ninth_freeEnergy_3/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"ninth_force_0.06rerun_3_01_Apr_144311"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 1:
        rerun = 3
        pre = "/scratch/wl45/mar_2018/ninth/"
        data_folder = "/scratch/wl45/apr_2018/01_week/all_data_folder/"
        folder_list = ["force_0.06_rg_0.15_lipid_1.0_mem_1_go_0.8"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        for attempt in range(2):
            try:
                process_complete_temper_data_3(pre, data_folder, folder_list, rerun=rerun, average_z=True, disReal=True, localQ=False, label="ninth_force_0.06")
            except FileNotFoundError:
                bias = "dis"
                simulation_list = glob.glob(f"{bias}_*")
                # simulation_list = ['dis_30.0']
                # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
                print(simulation_list)
                for dis in simulation_list:
                    print(dis)
                    cd(dis)
                    i = rerun
                    i_plus_one = i +1
                    if os.path.exists(f"log{i}"):
                        do(f"mv log{i} back_log{i}")  # in case of override
                    do(f"mkdir -p log{i}")
                    do(f"cp log.lammps log{i}/")
                    cd("..")

if args.day == "mar31":
    if args.mode == 10:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        i = 1
        freeEnergy_folder = f"ninth_freeEnergy/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"ninth_force_0.06rerun_1_31_Mar_182712"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 9:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_30.0']
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 1
            i_plus_one = i +1
            # do(f"mv log{i} back_log{i}")  # in case of override
            # do(f"mkdir -p log{i}")
            # do(f"cp log.lammps log{i}/")

            continueRunConvertion(n=12, rerun=i)
            do(f"mkdir {i_plus_one}")
            do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
    if args.mode == 8:
        rerun = 1
        pre = "/scratch/wl45/mar_2018/ninth/"
        data_folder = "/scratch/wl45/mar_2018/05_week/all_data_folder/"
        folder_list = ["force_0.06_rg_0.15_lipid_1.0_mem_1_go_0.8"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        for attempt in range(2):
            try:
                process_complete_temper_data_3(pre, data_folder, folder_list, rerun=rerun, average_z=True, disReal=True, localQ=False, label="ninth_force_0.06")
            except FileNotFoundError:
                bias = "dis"
                simulation_list = glob.glob(f"{bias}_*")
                # simulation_list = ['dis_30.0']
                # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
                print(simulation_list)
                for dis in simulation_list:
                    print(dis)
                    cd(dis)
                    i = rerun
                    i_plus_one = i +1
                    my_file = Path(f"log{i}")
                    if my_file.is_dir():
                        do(f"mv log{i} back_log{i}")  # in case of override
                    do(f"mkdir -p log{i}")
                    do(f"cp log.lammps log{i}/")
                    cd("..")

    if args.mode == 7:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        i = 3
        freeEnergy_folder = f"sixth_i235d/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"sixth_i235drerun_3_31_Mar_175942"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "300", "320", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=1, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 6 --force 2")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 4:
        pre = "/scratch/wl45/mar_2018/sixth/"
        data_folder = "/scratch/wl45/mar_2018/05_week/all_data_folder/"
        folder_list = ["rg_0.15_lipid_1.0_mem_1_go_0.8"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=1, average_z=True, disReal=True, localQ=False, label="sixth_i235d")
    if args.mode == 3:
        extra="i235d"
        bias = "dis_"
        simulation_list = glob.glob(f"{bias}*")
        print(simulation_list)
        # print(sim_list)
        for folder in simulation_list:
            print(folder)
            cd(folder)
            cd("1")
            do("mkdir backup")
            do("mv wham.* backup/")
            do("mv energy.* backup/")
            do("mv lipid.* backup/")
            do("mv rgs.* backup/")
            for i in range(12):
                do(f"cp {extra}_{i}/wham.dat wham.{i}.dat")
                do(f"cp {extra}_{i}/energy.dat energy.{i}.dat")
                do(f"cp {extra}_{i}/lipid.dat lipid.{i}.dat")
                do(f"cp {extra}_{i}/rgs.dat rgs.{i}.dat")
            # rerun(extra="i255d")
            cd("../..")


    if args.mode == 2:
        bias = "dis_"
        simulation_list = glob.glob(f"{bias}*")
        print(simulation_list)
        # print(sim_list)
        for folder in simulation_list:
            print(folder)
            cd(folder)
            cd("1")
            rerun(extra="i235d")
            # rerun(extra="i255d")
            cd("../..")
    if args.mode == 1:
        rerun(offset=0)
if args.day == "mar28":
    if args.mode == 8:
        pre = "/scratch/wl45/mar_2018/eighth/"
        data_folder = "/scratch/wl45/mar_2018/05_week/all_data_folder/"
        folder_list = ["force_0.03_rg_0.15_lipid_1.0_mem_1_go_0.8"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=1, average_z=True, disReal=True, localQ=False, label="eighth_force_0.03")
    if args.mode == 7:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        i = 3
        freeEnergy_folder = f"sixth_freeEnergy/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"sixth_new_localQrerun_3_28_Mar_162233"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "300", "320", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=1, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                # for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 6 --force 2")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 6:
        read_variable_folder(".", match="*_", average_z=True, disReal=True)
    if args.mode == 5:
        i = "0"
        # i = "1"
        # let_compute_localQ(temper=True, bias="dis_", sim_list=[i], queue=1)
        compute_disReal(temper=False, bias="", sim_list=[i], queue=0)
        compute_completeZ(temper=False, bias="", sim_list=[i], queue=0)
    if args.mode == 4:
        pre = "/scratch/wl45/mar_2018/sixth/"
        data_folder = "/scratch/wl45/mar_2018/05_week/all_data_folder/"
        folder_list = ["rg_0.15_lipid_1.0_mem_1_go_0.8"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=1, average_z=True, disReal=True, localQ=True, label="sixth_new_localQ")
    if args.mode == 3:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        i = 3
        freeEnergy_folder = f"sixth_new_localQ/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"sixth_new_localQrerun_3_28_Mar_162233"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "300", "320", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=1, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=1)

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                # for bias, mode in bias_list.items():
                for bias in range(36):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m 13 --commons 0 --nsample {nsample} --submode 24 --subsubmode {bias} --force 2")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 2:
        # i = "0"
        i = "1"
        # let_compute_localQ(temper=True, bias="dis_", sim_list=[i], queue=1)
        compute_disReal(temper=True, bias="dis_", sim_list=[i], queue=0)
        compute_completeZ(temper=True, bias="dis_", sim_list=[i], queue=0)
    if args.mode == 1:
        native_contacts_table = compute_localQ_init(MAX_OFFSET=6)
        for i in range(12):
            compute_localQ(native_contacts_table, MAX_OFFSET=6, pre=".", ii=i)
if args.day == "mar27":
    if args.mode == 7:
        print("how Constant force refolding")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        force_ramp_rate_list=[0.25]
        temperature_list=[300]
        # memb_k_list = [0.8, 1]
        memb_k_list = [1]
        rg_list = [0.15]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        # force_list = [0.22, 0.25, 0.28]
        # force_list = [0.0, 0.01]
        # force_list = [0.0143472]
        # force_list = ["ramp"]
        # force_list = [1, 2]
        # force_list = [3, 4]
        # force_list = [5, 6]
        force_list = [7, 8]
        repeat = 20
        change_list = [15]
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        change_list=change_list,
                        commons=0,simulation_base_steps=4e7)
    if args.mode == 6:
        temp_list = ["all"]
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        data_folder = "all_data_folder/"
        i = 3
        freeEnergy_folder = f"eighth_with_real_distance/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"eighthrerun_{i}_27_Mar_231139" for ii in range(3,4)]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["280", "300", "320"]}
        temp_dic = {"_280-350":["280", "290", "300", "310", "320"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=3, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)


        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 5:
        temp_list = ["all"]
        bias_list = {"first":"23", "second":"24", "thrid":"25"}
        data_folder = "all_data_folder/"
        i = 3
        freeEnergy_folder = f"sixth_localQ/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"sixthrerun_3_25_Mar_001526"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "300", "320", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=1, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=1)


        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                # for bias, mode in bias_list.items():
                for bias in range(18):
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m 13 --commons 0 --nsample {nsample} --submode 23 --subsubmode {bias} --force 2")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 4:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_30.0']
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 1
            i_plus_one = i +1
            do(f"mv log{i} back_log{i}")  # in case of override
            do(f"mkdir -p log{i}")
            do(f"cp log.lammps log{i}/")

            # continueRunConvertion(n=12, rerun=i)
            # do(f"mkdir {i_plus_one}")
            # do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            # do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
    if args.mode == 3:
        pre = "/scratch/wl45/mar_2018/eighth/"
        data_folder = "/scratch/wl45/mar_2018/05_week/all_data_folder/"
        folder_list = ["rg_0.15_lipid_1.0_mem_1_go_0.8"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=1, average_z=True, disReal=True, localQ=False, label="eighth")
    if args.mode == 2:
        i = "1"
        # i = "1"
        # let_compute_localQ(temper=True, bias="dis_", sim_list=[i], queue=1)
        compute_disReal(temper=True, bias="dis_", sim_list=[i], queue=0)
        compute_completeZ(temper=True, bias="dis_", sim_list=[i], queue=1)
    if args.mode == 1:
        cmd = f"python3 ~/opt/small_script/find_distance.py"
        do(cmd)
if args.day == "mar26":
    if args.mode == 1:
        print("how Constant force refolding")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        force_ramp_rate_list=[0.5]
        temperature_list=[260, 280]
        memb_k_list = [0.8]
        rg_list = [0.15]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        # force_list = [0.22, 0.25, 0.28]
        # force_list = [0.0, 0.01]
        force_list = [0.0143472]
        repeat = 30
        change_list = [15]
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        change_list=change_list,
                        commons=0,simulation_base_steps=4e7)
if args.day == "mar24":
    if args.mode == 8:
        temp_list = ["all"]
        bias_list = {"first":"23", "second":"24", "thrid":"25"}
        data_folder = "all_data_folder/"
        i = 3
        freeEnergy_folder = f"sixth_localQ/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"sixthrerun_3_25_Mar_001526"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "300", "320", "350"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=1, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=1)


        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m 13 --commons 0 --nsample {nsample} --submode {mode} --force 2")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 7:
        pre = "/scratch/wl45/mar_2018/sixth/"
        data_folder = "/scratch/wl45/mar_2018/05_week/all_data_folder/"
        folder_list = ["rg_0.15_lipid_1.0_mem_1_go_0.8"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=1, average_z=True, disReal=True, localQ=True, label="sixth")
    if args.mode == 6:
        # i = "0"
        i = "1"
        let_compute_localQ(temper=True, bias="dis_", sim_list=[i], queue=1)
        # compute_disReal(temper=True, bias="dis_", sim_list=[i], queue=0)
        # compute_completeZ(temper=True, bias="dis_", sim_list=[i], queue=1)
    if args.mode == 5:
        temp_list = ["all"]
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        data_folder = "all_data_folder/"
        i = 1
        freeEnergy_folder = f"seventh_with_real_distance/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"seventhrerun_{ii}_24_Mar_172933" for ii in range(1,2)]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["280", "300", "320"]}
        temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=2, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)


        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 4:
        temp_list = ["all"]
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        data_folder = "all_data_folder/"
        i = -2
        freeEnergy_folder = f"sixth_with_real_distance/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"sixth_without_directionrerun_{ii}_24_Mar_173616" for ii in range(6,8)]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "300", "320"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=1, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)


        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 3:
        pre = "/scratch/wl45/mar_2018/seventh/"
        data_folder = "/scratch/wl45/mar_2018/05_week/all_data_folder/"
        folder_list = ["rg_0.12_lipid_1.0_mem_1_go_0.8"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=0, average_z=True, disReal=True, localQ=False, label="seventh")
    if args.mode == 2:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_30.0']
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 0
            i_plus_one = i +1
            # do(f"mv log{i} back_log{i}")  # in case of override
            # do(f"mkdir -p log{i}")
            # do(f"cp log.lammps log{i}/")

            continueRunConvertion(n=12, rerun=i)
            do(f"mkdir {i_plus_one}")
            do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
    if args.mode == 1:
        pre = "/scratch/wl45/mar_2018/sixth/"
        data_folder = "/scratch/wl45/mar_2018/05_week/all_data_folder/"
        folder_list = ["rg_0.15_lipid_1.0_mem_1_go_0.8"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=3, average_z=True, disReal=True, localQ=False, label="sixth_without_direction")
if args.day == "mar22":
    if args.mode == 5:
        i = -2
        temp_list = ["280", "300"]
        make_metadata_3(temps_list=temp_list,k=0.02, i=i)
    if args.mode == 4:
        temp_list = ["all"]
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        data_folder = "all_data_folder/"
        i = -2
        freeEnergy_folder = f"second_sixth_with_real_distance/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"sixth_without_directionrerun_{ii}_23_Mar_134412" for ii in range(2,4)]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "300", "320", "350"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)


        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 3:
        pre = "/scratch/wl45/mar_2018/sixth/"
        data_folder = "/scratch/wl45/mar_2018/04_week/all_data_folder/"
        folder_list = ["rg_0.15_lipid_1.0_mem_1_go_0.8_without_direction"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=1, average_z=True, disReal=True, localQ=False, label="sixth_without_direction")
    if args.mode == 2:
        # i = "0"
        i = "3"
        compute_disReal(temper=True, bias="dis_", sim_list=[i], queue=0)
        compute_completeZ(temper=True, bias="dis_", sim_list=[i], queue=1)
    if args.mode == 1:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_30.0']
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 0
            i_plus_one = i +1
            do(f"mv log{i} back_log{i}")  # in case of override
            do(f"mkdir -p log{i}")
            do(f"cp log.lammps log{i}/")

            # continueRunConvertion(n=12, rerun=i)
            # do(f"mkdir {i_plus_one}")
            # do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            # do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
if args.day == "mar21":
    if args.mode == 2:
        i = 3
        temp_list = ["280"]
        make_metadata_3(temps_list=temp_list,k=0.00, i=i)
    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        i = -2
        freeEnergy_folder = f"sixth_with_real_distance_3/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"sixth_disRealrerun_{i}_22_Mar_011035" for i in range(4,6)]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "300", "320", "350"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)


        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "mar20":
    if args.mode == 4:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        i = 1
        freeEnergy_folder = f"sixth_with_real_distance/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"sixth_disRealrerun_{i}_20_Mar_200631" for i in range(1,2)]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_280-350":["280", "300", "320", "350"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)


        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 22 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 3:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_30.0']
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            if dis == "dis_30.0":
                continue
            cd(dis)
            for a in glob.glob("*.out"):
                if pd.read_table(a).iloc[-1].values[0].split(" ")[0] == '40000000':
                    # print(a)
                    os.system(f"tail -n 5004 {a} > log1/log.lammps")
                    os.system("cp log1/log.lammps log.lammps")
            cd("..")
    if args.mode == 2:
        pre = "/scratch/wl45/mar_2018/sixth/"
        data_folder = "/scratch/wl45/mar_2018/04_week/all_data_folder/"
        folder_list = ["rg_0.15_lipid_1.0_mem_1_go_0.8_without_direction"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=0, average_z=True, disReal=True, localQ=False, label="sixth_without_direction")
    if args.mode == 1:
        compute_disReal(temper=True, bias="dis_", sim_list=["2"], queue=1)
        compute_completeZ(temper=True, bias="dis_", sim_list=["2"], queue=1)
if args.day == "mar19":
    if args.mode == 5:
        print("start")
        time.sleep(30)
        print("hi")
    if args.mode == 4:
        time.sleep(4*60*60)
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_30.0']
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            if dis == "dis_30.0":
                continue
            print(dis)
            cd(dis)
            i = 0
            i_plus_one = i +1
            do(f"mv log{i} back_log{i}")  # in case of override
            do(f"mkdir -p log{i}")
            do(f"mv log.lammps.* log{i}/")
            do(f"cp log.lammps log{i}/")
            do(f"cp x.* log{i}/")

            continueRunConvertion(n=12, rerun=i)
            do(f"mkdir {i_plus_one}")
            do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
    if args.mode == 3:
        pre = "/scratch/wl45/mar_2018/fifth/"
        data_folder = "/scratch/wl45/mar_2018/04_week/all_data_folder/"
        folder_list = ["rg_0.2_lipid_1.0_mem_1_go_0.7"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=0, average_z=True, disReal=True, localQ=False, label="fifth_disReal")
    if args.mode == 2:
        # compute_disReal(temper=True, bias="dis_", sim_list=["0"], queue=0)
        compute_completeZ(temper=True, bias="dis_", sim_list=["0"], queue=0)
    if args.mode == 1:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 2
            i_plus_one = i +1
            do(f"mkdir -p log{i}")
            do(f"mv log.* log{i}/")
            do(f"cp log{i}/log.lammps .")
            do(f"cp x.* log{i}/")

            continueRunConvertion(n=12, rerun=i)
            do(f"mkdir {i_plus_one}")
            do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
if args.day == "mar18":
    if args.mode == 7:
        print("how Constant force refolding")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        force_ramp_rate_list=[1]
        temperature_list=[300]
        memb_k_list = [0.6, 0.7, 0.8]
        rg_list = [0.15]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        # force_list = [0.22, 0.25, 0.28]
        # force_list = [0.0, 0.01]
        force_list = [0.0143472]
        repeat = 30
        change_list = [15]
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        change_list=change_list,
                        commons=0,simulation_base_steps=4e7)
    if args.mode == 6:
        temp_list = ["all"]
        bias_list = {"1d_qw":"1"}
        data_folder = "all_data_folder/"
        i = 1
        freeEnergy_folder = f"freeEnergy_third_q_bias/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = ["third_q_bias"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"300":["300"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, sample_range_mode=i, sub_mode_name=temp_mode, average_z=0, chosen_mode=0)


        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=100, i=i)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 4 --force 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 5:
        data = read_folder(".", match="qbias")
        data["BiasTo"] = data["Run"].apply(lambda x: float(x.split("_")[1])*0.02)
        data["Temp"] = "T_defined"
        data["Step"] = data["Steps"]
        data.to_feather("../all_data_folder/third_q_bias.feather")
    if args.mode == 4:
        print("how Constant force refolding")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        force_ramp_rate_list=[1]
        temperature_list=[300]
        memb_k_list = [0.1, 0.6]
        rg_list = [0.1, 0.15, 0.2]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        # force_list = [0.22, 0.25, 0.28]
        # force_list = [0.0, 0.01]
        force_list = [0, 0.01]
        repeat = 10
        change_list = [15]
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        change_list=change_list,
                        commons=0,simulation_base_steps=2e7)
    if args.mode == 3:
        compute_disReal(bias="dis_", sim_list=["0"], queue=0)
        compute_completeZ(bias="dis_", sim_list=["0"], queue=1)
    if args.mode == 1:
        location = "."
        variables = glob.glob(os.path.join(location, "rg_*"))
        print(variables)
        data_list = []
        for variableFolder in variables:
            cd(variableFolder)
            cd("simulation")
            compute_disReal(bias="", sim_list=["0"], queue=0)
            compute_completeZ(bias="", sim_list=["0"], queue=1)
            cd("../..")
    if args.mode == 2:
        read_variable_folder(".", match="rg_*", average_z=True, disReal=True)
if args.day == "mar17":
    if args.mode == 5:
        data = read_simulation_2("/scratch/wl45/mar_2018/extra_03_week/refold/foldingTemp/0/")
        data.to_feather("test1.feather")
    if args.mode == 4:
        print("how Constant force refolding")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        force_ramp_rate_list=[1]
        temperature_list=[300, 320]
        memb_k_list = [0.1, 0.6]
        rg_list = [0.1]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        # force_list = [0.22, 0.25, 0.28]
        # force_list = [0.0, 0.01]
        force_list = [0.01, 0.02]
        repeat = 10
        change_list = [12, 15]
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        change_list=change_list,
                        commons=0,simulation_base_steps=2e7)
    if args.mode == 3:
        location = "."
        variables = glob.glob(os.path.join(location, "*_"))
        print(variables)
        data_list = []
        for variableFolder in variables:
            cd(variableFolder)
            cd("simulation")
            compute_disReal(bias="", sim_list=["0"], queue=0)
            cd("../..")
    if args.mode == 2:
        location = "."
        variables = glob.glob(os.path.join(location, "*_"))
        print(variables)
        data_list = []
        for variableFolder in variables:
            cd(variableFolder)
            cd("simulation")
            compute_completeZ()
            cd("../..")
    if args.mode == 1:
        read_variable_folder(".", average_z=True, disReal=True)
if args.day == "mar16":
    if args.mode == 8:
        print("how Constant force refolding")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        force_ramp_rate_list=[0.5]
        temperature_list=[300]
        memb_k_list = [0, 0.1, 0.5, 0.6, 1]
        rg_list = [0.1]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        # force_list = [0.22, 0.25, 0.28]
        # force_list = [0.0, 0.01]
        force_list = [0.0001, 0.01, 0.05]
        repeat = 10
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=0,simulation_base_steps=2e7)
    if args.mode == 7:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 1
            i_plus_one = i +1
            do(f"mkdir -p log{i}")
            do(f"mv log.* log{i}/")
            do(f"cp log{i}/log.lammps .")
            do(f"cp x.* log{i}/")

            # continueRunConvertion(n=12, rerun=i)
            # do(f"mkdir {i_plus_one}")
            # do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            # do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
    if args.mode == 6:
        do("freeEnergy_run.py -m 4 2xov/")
    if args.mode == 5:
        location = "."
        data = read_folder(location, average_z=True)
        data.to_feather("test.feather")
    if args.mode == 4:
        print("compute completeZ")
        # runFolders = os.listdir(location+"/simulation")
        # runFolders = [f for f in runFolders if re.match(r'[0-9]+', f)]
        # cd("simulation")
        bias = "dis"
        simulation_list = glob.glob("*")
        print(simulation_list)
        interactive = 1
        # sim_list = ["0", "1"]
        # sim_list = ["2"]
        sim_list = ["0"]
        for sim in sim_list:
            for folder in simulation_list:
                cd(folder)
                cd(sim)
                print(folder)
                quick = quick_template_slurm.format("-d mar16 -m 3")
                with open("computeZ.slurm", "w") as f:

                    if interactive:
                        quick = quick.replace("--time=01:30:00", "--time=00:30:00")
                        quick = quick.replace("#SBATCH --account=ctbp-common", "")
                        quick = quick.replace("ctbp-common", "interactive")
                    f.write(quick)
                    # f.write(quick.replace("ctbp-common", "commons"))
                do("sbatch computeZ.slurm")
                cd("../..")
    if args.mode == 3:
        compute_average_z_2("dump.lammpstrj", "z_complete.dat")
    if args.mode == 2:
        data = read_simulation_2("/scratch/wl45/mar_2018/extra_03_week/refold/foldingTemp_2/0/")
        data.to_feather("test.feather")
    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        i = 1
        freeEnergy_folder = f"with_force_fourth_with_real_distance/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"fourth_disRealrerun_{i}_16_Mar_220433" for i in range(2,4)]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, sample_range_mode=-2, sub_mode_name=temp_mode, average_z=5, chosen_mode=0)


        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=-2)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 2 --force 1")
                    cd("..")
                cd("..")
        cd("..")
if args.day == "mar15":
    if args.mode == 5:
        pre = "/scratch/wl45/mar_2018/fourth/"
        data_folder = "/scratch/wl45/mar_2018/extra_03_week/all_data_folder/"
        folder_list = ["rg_0.1_lipid_1.0_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=2, average_z=True, disReal=True, localQ=False, label="fourth_disReal")
    if args.mode == 4:
        print("compute DisReal")
        # print(native_contacts_table)
        # cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        sim_list = ["0"]
        sim_list = ["1"]
        # sim_list = ["2", "0", "1"]
        for sim in sim_list:
            for folder in simulation_list:
                cd(folder)
                cd(sim)
                print(folder)
                with open("quick.slurm", "w") as f:
                    quick_template_slurm = quick_template_slurm.format("-d mar14 -m 3")
                    f.write(quick_template_slurm)
                    # f.write(quick_template_slurm.replace("ctbp-common", "commons"))
                do("sbatch quick.slurm")
                cd("../..")
    if args.mode == 3:
        print("compute completeZ")
        # cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        sim_list = ["0", "1"]
        sim_list = ["2"]
        sim_list = ["0"]
        sim_list = ["1"]
        for sim in sim_list:
            for folder in simulation_list:
                cd(folder)
                cd(sim)
                print(folder)
                quick = quick_template_slurm.format("-d mar10 -m 7")
                with open("computeZ.slurm", "w") as f:
                    f.write(quick)
                    # f.write(quick.replace("ctbp-common", "commons"))
                do("sbatch computeZ.slurm")
                cd("../..")
    if args.mode == 2:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 0
            i_plus_one = i +1
            # do(f"mkdir -p log{i}")
            # do(f"mv log.* log{i}/")
            # do(f"cp log{i}/log.lammps .")
            # do(f"cp x.* log{i}/")

            continueRunConvertion(n=12, rerun=i)
            do(f"mkdir {i_plus_one}")
            do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
    if args.mode == 1:
        print("how Constant force refolding")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        force_ramp_rate_list=[0.5]
        temperature_list=[300]
        memb_k_list = [1]
        rg_list = [0.15, 0.1]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        # force_list = [0.22, 0.25, 0.28]
        # force_list = [0.0, 0.01]
        force_list = [0.0]
        repeat = 20
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=0,simulation_base_steps=2e7)

if args.day == "mar14":
    if args.mode == 7:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"

        freeEnergy_folder = f"second_combined_expectedDistanceReal_h5_h6/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"second_disRealrerun_{i}_14_Mar_202230" for i in range(4,6)]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, sample_range_mode=-1, sub_mode_name=temp_mode, average_z=4, chosen_mode=0)


        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=-1)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 3")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 6:
        temp_list = ["350", "400", "450", "500", "550"]
        make_metadata_3(temps_list=temp_list,k=0.01, i=-1)
    if args.mode == 5:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"

        freeEnergy_folder = f"second_combined_expectedDistanceReal/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"second_disRealrerun_{i}_14_Mar_202230" for i in range(4,6)]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, sample_range_mode=-1, sub_mode_name=temp_mode, average_z=3, chosen_mode=0)


        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=-1)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 3")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 4:
        pre = "/scratch/wl45/mar_2018/second/"
        data_folder = "/scratch/wl45/mar_2018/03_week/all_data_folder/"
        folder_list = ["rg_0.1_lipid_1.0_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=2, average_z=True, disReal=True, localQ=False, label="second_disReal")
    if args.mode == 3:
        for i in range(12):
            cmd = f"python3 ~/opt/small_script/find_distance.py -m {i}"
            do(cmd)
    if args.mode == 2:
        print("compute DisReal")
        # print(native_contacts_table)
        # cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # sim_list = ["0"]
        sim_list = ["2", "0", "1"]
        for sim in sim_list:
            for folder in simulation_list:
                cd(folder)
                cd(sim)
                print(folder)
                with open("quick.slurm", "w") as f:
                    quick_template_slurm = quick_template_slurm.format("-d mar14 -m 3")
                    f.write(quick_template_slurm)
                    # f.write(quick_template_slurm.replace("ctbp-common", "commons"))
                do("sbatch quick.slurm")
                cd("../..")
    if args.mode == 1:
        temp_list = ["350", "400", "450", "500", "550"]
        make_metadata_3(temps_list=temp_list,k=0.02, i=-1, biasLow=50)
if args.day == "mar12":
    if args.mode == 14:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"

        freeEnergy_folder = f"second_combined_expectedLocalQ_freeEnergy/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"second_rerun_{i}_13_Mar_041016" for i in range(4,6)]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=-1)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 5")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 13:
        print("Set up simulation and go")
        with open("quick.slurm", "w") as f:
            # f.write(localQ_slurm)
            quick_template_slurm = quick_template_slurm.format("-d mar12 -m 12")
            f.write(quick_template_slurm.replace("ctbp-common", "commons"))
        do("sbatch quick.slurm")
        cd("../..")
    if args.mode == 12:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"

        freeEnergy_folder = f"second_combined_expectedLocalQ_freeEnergy/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"second_rerun_{i}_13_Mar_041016" for i in range(4,6)]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, sample_range_mode=-1, sub_mode_name=temp_mode, average_z=2, chosen_mode=1)


        # cd(freeEnergy_folder)
        # for temp_mode, temp_list in temp_dic.items():
        #         cd(temp_mode)
        #         for bias, mode in bias_list.items():
        #             # name = "low_t_" + bias
        #             name = bias
        #             print(name)
        #             do("rm -r "+name)
        #             do("mkdir -p " + name)
        #             cd(name)
        #             make_metadata_3(temps_list=temp_list,k=0.02, i=-1)
        #             nsample = len(folder_list)*2500
        #             do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 5")
        #             cd("..")
        #         cd("..")
        # cd("..")
    if args.mode == 11:
        print("Process temper")
        with open("quick.slurm", "w") as f:
            # f.write(localQ_slurm)
            quick_template_slurm = quick_template_slurm.format("-d mar12 -m 10")
            f.write(quick_template_slurm.replace("ctbp-common", "commons"))
        do("sbatch quick.slurm")
        cd("../..")
    if args.mode == 10:
        pre = "/scratch/wl45/feb_2018/week_of_feb19/"
        data_folder = "/scratch/wl45/mar_2018/03_week/all_data_folder/"
        folder_list = ["rg_0.1_lipid_1.0_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=2, average_z=True, localQ=True, label="second_")
    if args.mode == 9:
        print("how Constant force refolding")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        force_ramp_rate_list=[0.5]
        temperature_list=[500]
        memb_k_list = [1]
        rg_list = [0.1]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        # force_list = [0.22, 0.25, 0.28]
        force_list = [0.0, 0.01, 0.02, 0.03, 0.05]
        repeat = 20
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=0,simulation_base_steps=2e7)
    if args.mode == 8:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"

        freeEnergy_folder = f"second_combined_expectedDistance_freeEnergy/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"second_rerun_{i}_12_Mar_211030" for i in range(4,6)]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, sample_range_mode=-1, sub_mode_name=temp_mode, average_z=2, chosen_mode=0)


        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=-1)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 2")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 7:
        print("compute localQ")
        # print(native_contacts_table)
        # cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # sim_list = ["0"]
        sim_list = ["2", "0", "1"]
        for sim in sim_list:
            for folder in simulation_list:
                cd(folder)
                cd(sim)
                print(folder)
                with open("localQ.slurm", "w") as f:
                    # f.write(localQ_slurm)
                    f.write(localQ_slurm.replace("ctbp-common", "commons"))
                do("sbatch localQ.slurm")
                cd("../..")
    if args.mode == 6:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"

        freeEnergy_folder = f"second_combined_z_6_freeEnergy/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        folder_list = [f"second_rerun_{i}_12_Mar_211030" for i in range(4,6)]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, sample_range_mode=-1, sub_mode_name=temp_mode, average_z=2, chosen_mode=0)


        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=-1)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 1")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 5:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(1, 6):
            freeEnergy_folder = f"second_z_6_freeEnergy_{sample_range_mode}/"
            print(freeEnergy_folder)
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            # folder_list = ["rerun_1_08_Mar_154259"]
            # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630"]
            folder_list = [f"second_rerun_{sample_range_mode}_12_Mar_211030"]
            # submode_list = ["_no_energy"]
            # submode_list = ["", "only_500"]
            # submode_list = ["350", "400", "450", "500", "550"]

            temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    move_data3(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=temp_mode, average_z=2, chosen_mode=0)

            cd(freeEnergy_folder)
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    cd(temp_mode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        make_metadata_3(temps_list=temp_list,k=0.02, i=sample_range_mode)
                        do("python3 ~/opt/pulling_analysis_2.py -m {} --commons 1 --nsample 2500 --submode 1".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
    if args.mode == 4:
        pre = "/scratch/wl45/feb_2018/week_of_feb19/"
        data_folder = "/scratch/wl45/mar_2018/03_week/all_data_folder/"
        folder_list = ["rg_0.1_lipid_1.0_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=2, average_z=True, label="second_")
    if args.mode == 3:
        print("compute localQ")
        # print(native_contacts_table)
        # cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # sim_list = ["0"]
        sim_list = ["0", "1", "2"]
        for sim in sim_list:
            for folder in simulation_list:
                cd(folder)
                cd(sim)
                print(folder)
                with open("localQ.slurm", "w") as f:
                    # f.write(localQ_slurm)
                    f.write(localQ_slurm.replace("ctbp-common", "commons"))
                do("sbatch localQ.slurm")
                cd("../..")
    if args.mode == 2:
        pre = "/scratch/wl45/feb_2018/week_of_feb19/"
        data_folder = "/scratch/wl45/mar_2018/03_week/all_data_folder/"
        folder_list = ["rg_0.1_lipid_1.0_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=2, average_z=True, label="first_")
    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(1, 6):
            freeEnergy_folder = f"first_z_6_freeEnergy_{sample_range_mode}/"
            print(freeEnergy_folder)
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            # folder_list = ["rerun_1_08_Mar_154259"]
            folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630"]
            # submode_list = ["_no_energy"]
            # submode_list = ["", "only_500"]
            # submode_list = ["350", "400", "450", "500", "550"]

            temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    move_data3(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=temp_mode, average_z=2, chosen_mode=0)

            cd(freeEnergy_folder)
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    cd(temp_mode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        make_metadata_3(temps_list=temp_list,k=0.02, i=sample_range_mode)
                        do("python3 ~/opt/pulling_analysis_2.py -m {} --commons 1 --nsample 2500 --submode 1".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
if args.day == "mar11":
    if args.mode == 1:
        pre = "/scratch/wl45/feb_2018/week_of_feb19/"
        data_folder = "/scratch/wl45/feb_2018/week_of_feb19/all_data_folder/"
        folder_list = ["rg_0.1_lipid_1.0_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_2(pre, data_folder, folder_list, rerun=2, average_z=True)
    if args.mode == 2:
        temp_list = ["all"]
        bias_list = {"2d_z_qw":"13", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12"}
        # bias_list = {"2d_z_qw":"13", "1d_dis":"9", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_z_dis":"14", "2d_qw_dis":"11"}
        data_folder = "all_data_folder/"
        sample_range_mode = -1
        freeEnergy_folder = f"z_6_freeEnergy_{sample_range_mode}/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        folder_list = ["rerun_2_10_Mar_183455"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_400-500":["400", "450", "500"]}
        for temp_mode, temp_list in temp_dic.items():
            for folder in folder_list:
                move_data2(data_folder, freeEnergy_folder, folder, sample_range_mode=-1, sub_mode_name=temp_mode, average_z=2, chosen_mode=0)
        cd(freeEnergy_folder)

        for temp_mode, temp_list in temp_dic.items():
            for folder in folder_list:
                cd(folder+temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata(temps_list=temp_list,k=0.02)
                    do("python3 ~/opt/pulling_analysis_2.py -m {} --commons 1 --nsample 5000 --submode 1".format(mode))
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 3:
        temp_list = ["all"]
        bias_list = {"2d_z_qw":"13", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12"}
        # bias_list = {"2d_z_qw":"13", "1d_dis":"9", "1d_qw":"10", "1d_z":"12"}
        # bias_list = {"2d_z_dis":"14", "2d_qw_dis":"11"}
        data_folder = "all_data_folder/"
        sample_range_mode = -1
        freeEnergy_folder = f"second_z_6_freeEnergy_{sample_range_mode}/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        folder_list = ["rerun_2_11_Mar_160738"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_400-500":["400", "450", "500"]}
        for temp_mode, temp_list in temp_dic.items():
            for folder in folder_list:
                move_data2(data_folder, freeEnergy_folder, folder, sample_range_mode=-1, sub_mode_name=temp_mode, average_z=2, chosen_mode=0)
        cd(freeEnergy_folder)

        for temp_mode, temp_list in temp_dic.items():
            for folder in folder_list:
                cd(folder+temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata(temps_list=temp_list,k=0.02)
                    do("python3 ~/opt/pulling_analysis_2.py -m {} --commons 1 --nsample 5000 --submode 1".format(mode))
                    cd("..")
                cd("..")
        cd("..")
if args.day == "mar10":
    if args.mode == 1:
        commons = 1
        # temp_list = temp_list = ["350", "400", "450"]
        temp_list = temp_list = ["400"]
        folder = "temp_default"
        do("mkdir "+folder)
        cd(folder)
        make_metadata_2(temps_list=temp_list,k=0.02)
        cmd = "-b 3 -e 1 -d 2 -f 0.05 -nsamples 2500 -v1 4 -v1n 30 -v2 2 -v2n 30 -ti 10 -st 380 -et 450 -ss y"
        freeEnergy = freeEnergy.format(cmd)
        if commons:
            freeEnergy = freeEnergy.replace("ctbp-common", "commons")
            freeEnergy = freeEnergy.replace("--time=23:00:00", "--time=08:00:00")
        # create freeEnergy.slurm
        with open("freeEnergy.slurm", "w") as f:
            f.write(freeEnergy)
        do("sbatch freeEnergy.slurm")
        cd("..")
    if args.mode == 2:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(3, 5):
            freeEnergy_folder = f"first_freeEnergy_{sample_range_mode}/"
            print(freeEnergy_folder)
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            # folder_list = ["rerun_1_08_Mar_154259"]
            folder_list = ["rerun_2_09_Mar_154823"]
            # submode_list = ["_no_energy"]
            # submode_list = ["", "only_500"]
            # submode_list = ["350", "400", "450", "500", "550"]

            temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    move_data2(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=temp_mode, average_z=True, chosen_mode=0)
            cd(freeEnergy_folder)
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    cd(folder+temp_mode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        make_metadata(temps_list=temp_list,k=0.02)
                        do("python3 ~/opt/pulling_analysis_2.py -m {} --commons 0 --nsample 2500 --submode 1".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
    if args.mode == 3:
        compute_average_z("dump.lammpstrj.1", "z_test.dat")
    if args.mode == 4:
        # protocol_list = ["er", "awsemer", "frag", "raptor"]
        protocol_list = ["raptor"]
        protein_list = ["2xov"]
        for protein in protein_list:
            for protocol in protocol_list:
                print("Work on protein: {}, protocol: {}".format(protein, protocol))
                do("mkdir -p {}".format(protocol))
                do("cp -r {} {}/".format(protein, protocol))
                cd(protocol)
                cd(protein)
                fileName = "{}_multi.in".format(protein)
                if protocol == "raptor":
                    backbone_file = "fix_backbone_coeff_er.data"
                    do("cp ~/opt/gremlin/protein/{}/raptor/go_rnativeC* .".format(protein))
                else:
                    backbone_file = "fix_backbone_coeff_{}.data".format(protocol)
                    do("cp ~/opt/gremlin/protein/{}/gremlin/go_rnativeC* .".format(protein))
                with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
                    for line in file:
                        tmp = line
                        tmp = tmp.replace("fix_backbone_coeff_er.data", backbone_file)
                        print(tmp, end='')
                cd("..")
                do("run.py -m 0 -n 20 {}".format(protein))
                cd("..")
    if args.mode == 5:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        sample_range_mode = -1
        freeEnergy_folder = f"abs_first_freeEnergy_{sample_range_mode}/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        folder_list = ["rerun_2_10_Mar_183455"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        temp_dic = {"_400-500":["400", "450", "500"]}
        for temp_mode, temp_list in temp_dic.items():
            for folder in folder_list:
                move_data2(data_folder, freeEnergy_folder, folder, sample_range_mode=-1, sub_mode_name=temp_mode, average_z=True, chosen_mode=0)
        cd(freeEnergy_folder)

        for temp_mode, temp_list in temp_dic.items():
            for folder in folder_list:
                cd(folder+temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata(temps_list=temp_list,k=0.02)
                    do("python3 ~/opt/pulling_analysis_2.py -m {} --commons 1 --nsample 5000 --submode 1".format(mode))
                    cd("..")
                cd("..")
        cd("..")
    if args.mode == 6:
        compute_theta_for_each_helix("dump.lammpstrj.0")
    if args.mode == 7:
        for i in range(12):
            compute_average_z_2(f"dump.lammpstrj.{i}", f"z_complete_{i}.dat")
    if args.mode == 8:
        # cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        sim_list = ["0", "1"]
        sim_list = ["2"]
        for sim in sim_list:
            for folder in simulation_list:
                cd(folder)
                cd(sim)
                print(folder)
                quick = quick_template_slurm.format("-d mar10 -m 7")
                with open("computeZ.slurm", "w") as f:
                    # f.write(quick_slurm)
                    f.write(quick.replace("ctbp-common", "commons"))
                do("sbatch computeZ.slurm")
                cd("../..")
    if args.mode == 9:
        pre = "/scratch/wl45/mar_2018/02_week/"
        data_folder = "/scratch/wl45/mar_2018/02_week/all_data_folder/"
        folder_list = ["rg_0.1_lipid_1.0_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_2(pre, data_folder, folder_list, rerun=2, average_z=True)
if args.day == "mar08":
    if args.mode == 3:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(1, 3):
            freeEnergy_folder = f"first_freeEnergy_{sample_range_mode}/"
            print(freeEnergy_folder)
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["rerun_1_08_Mar_154259"]
            # submode_list = ["_no_energy"]
            # submode_list = ["", "only_500"]
            # submode_list = ["350", "400", "450", "500", "550"]

            temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    move_data2(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=temp_mode, average_z=True, chosen_mode=0)
            cd(freeEnergy_folder)
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    cd(folder+temp_mode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        make_metadata(temps_list=temp_list,k=0.02)
                        do("python3 ~/opt/pulling_analysis_2.py -m {} --commons 0 --nsample 2500 --submode 1".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
    if args.mode == 2:
        pre = "/scratch/wl45/mar_2018/02_week/"
        data_folder = "/scratch/wl45/mar_2018/02_week/all_data_folder/"
        folder_list = ["rg_0.1_lipid_1.0_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_2(pre, data_folder, folder_list, rerun=2, average_z=True)
    if args.mode == 1:
        # cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        sim_list = ["2"]
        # sim_list = ["2"]
        for sim in sim_list:
            for folder in simulation_list:
                cd(folder)
                cd(sim)
                print(folder)
                with open("computeZ.slurm", "w") as f:
                    # f.write(quick_slurm)
                    f.write(quick_slurm.replace("ctbp-common", "commons"))
                do("sbatch computeZ.slurm")
                cd("../..")
if args.day == "mar06":
    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(3, 5):
            freeEnergy_folder = f"freeEnergy_{sample_range_mode}/"
            print(freeEnergy_folder)
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["rerun_2_06_Mar_225059"]
            # submode_list = ["_no_energy"]
            # submode_list = ["", "only_500"]
            # submode_list = ["350", "400", "450", "500", "550"]

            temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    move_data2(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=temp_mode, average_z=True, chosen_mode=0)
            cd(freeEnergy_folder)
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    cd(folder+temp_mode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        make_metadata(temps_list=temp_list,k=0.02)
                        do("python3 ~/opt/pulling_analysis_2.py -m {} --commons 0 --nsample 2500 --submode 1".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
if args.day == "mar05":
    if args.mode == 5:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 2
            i_plus_one = i +1
            do(f"mkdir -p log{i}")
            do(f"mv log.* log{i}/")
            do(f"cp log{i}/log.lammps .")
            do(f"cp x.* log{i}/")
            # continueRunConvertion(n=12, rerun=i)
            # do(f"mkdir {i_plus_one}")

            # run_slurm = base_run_slurm.format(i_plus_one)
            # with open(f"run_{i_plus_one}.slurm", "w") as r:
            #     r.write(run_slurm)

            # do(f"sbatch run_{i_plus_one}.slurm")

            # do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            # do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
    if args.mode == 4:
        pre = "/scratch/wl45/feb_2018/week_of_feb19/"
        data_folder = "/scratch/wl45/feb_2018/week_of_feb19/all_data_folder/"
        folder_list = ["rg_0.1_lipid_1.0_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_2(pre, data_folder, folder_list, rerun=2, average_z=True)
    if args.mode == 3:
        bias = "dis"
        # simulation_list = glob.glob(f"{bias}_*")
        simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 1
            i_plus_one = i +1
            # do(f"mkdir -p log{i}")
            # do(f"mv log.* log{i}/")
            # do(f"cp log{i}/log.lammps .")
            # do(f"cp x.* log{i}/")
            # continueRunConvertion(n=12, rerun=i)
            # do(f"mkdir {i_plus_one}")

            # run_slurm = base_run_slurm.format(i_plus_one)
            # with open(f"run_{i_plus_one}.slurm", "w") as r:
            #     r.write(run_slurm)

            # do(f"sbatch run_{i_plus_one}.slurm")

            # do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
    if args.mode == 1:
        print("compute localQ")
        # print(native_contacts_table)
        # cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # sim_list = ["0"]
        sim_list = ["0", "1"]
        for sim in sim_list:
            for folder in simulation_list:
                cd(folder)
                cd(sim)
                print(folder)
                with open("localQ.slurm", "w") as f:
                    # f.write(localQ_slurm)
                    f.write(localQ_slurm.replace("ctbp-common", "commons"))
                do("sbatch localQ.slurm")
                cd("../..")
    if args.mode == 2:
        pre = "/scratch/wl45/mar_2018/01_week/"
        data_folder = "/scratch/wl45/mar_2018/01_week/all_data_folder/"
        folder_list = ["rg_0.1_lipid_1.0_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_2(pre, data_folder, folder_list, rerun=1, average_z=True)
if args.day == "mar03":
    if args.mode == 1:
        simulation_list = glob.glob("dis_*")
        # simulation_list = ['dis_62.0', 'dis_48.0', 'dis_34.0', 'dis_36.0', 'dis_60.0', 'dis_40.0', 'dis_78.0', 'dis_74.0', 'dis_66.0', 'dis_52.0', 'dis_94.0', 'dis_82.0', 'dis_98.0', 'dis_68.0', 'dis_92.0', 'dis_64.0', 'dis_42.0', 'dis_58.0', 'dis_90.0', 'dis_100.0', 'dis_32.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 1
            i_plus_one = i +1
            # do(f"mkdir -p log{i}")
            # do(f"mv log.* log{i}/")
            # do(f"cp log{i}/log.lammps .")
            # do(f"cp x.* log{i}/")
            continueRunConvertion(n=12, rerun=i)
            do(f"mkdir {i_plus_one}")

            run_slurm = base_run_slurm.format(i_plus_one)
            with open(f"run_{i_plus_one}.slurm", "w") as r:
                r.write(run_slurm)

            # do(f"sbatch run_{i_plus_one}.slurm")

            # do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
    if args.mode == 2:
        # cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # sim_list = ["1"]
        sim_list = ["2"]
        for sim in sim_list:
            for folder in simulation_list:
                cd(folder)
                cd(sim)
                print(folder)
                with open("computeZ.slurm", "w") as f:
                    # f.write(quick_slurm)
                    f.write(quick_slurm.replace("ctbp-common", "commons"))
                do("sbatch computeZ.slurm")
                cd("../..")
    if args.mode ==3:
        for i in range(12):
            compute_average_z(f"dump.lammpstrj.{i}", f"z_{i}.dat")
if args.day == "feb28":
    if args.mode == 1:
        # print(native_contacts_table)
        # cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # sim_list = ["0"]
        sim_list = ["0", "1"]
        for sim in sim_list:
            for folder in simulation_list:
                cd(folder)
                cd(sim)
                print(folder)
                with open("localQ.slurm", "w") as f:
                    f.write(localQ_slurm.replace("ctbp-common", "commons"))
                do("sbatch localQ.slurm")
                cd("../..")
    if args.mode == 2:
        native_contacts_table = compute_localQ_init()
        for i in range(12):
            compute_localQ(native_contacts_table, pre=".", ii=i)
        # compute_localQ(native_contacts_table, pre=location, ii=i)
    if args.mode == 3:
        pre = "/scratch/wl45/feb_2018/week_of_feb19/"
        data_folder = "/scratch/wl45/feb_2018/week_of_feb19/all_data_folder/"
        folder_list = ["rg_0.1_lipid_1.0_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data(pre, data_folder, folder_list, rerun=1, average_z=True, localQ=True)
    if args.mode == 4:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(2, 3):
            freeEnergy_folder = f"2_compute_expected_localQ_{sample_range_mode}/"
            print(freeEnergy_folder)
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["expected_localQ"]
            # submode_list = ["_no_energy"]
            # submode_list = ["", "only_500"]
            # submode_list = ["350", "400", "450", "500", "550"]

            temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    move_data2(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=temp_mode, average_z=True, chosen_mode=1)
    if args.mode == 5:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(2, 3):
            freeEnergy_folder = f"2_compute_expected_localQ_{sample_range_mode}/"
            print(freeEnergy_folder)
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["expected_localQ"]
            # submode_list = ["_no_energy"]
            # submode_list = ["", "only_500"]
            # submode_list = ["350", "400", "450", "500", "550"]

            temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
            cd(freeEnergy_folder)
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    cd(folder+temp_mode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        make_metadata(temps_list=temp_list,k=0.02)
                        do("python3 ~/opt/pulling_analysis_2.py -m {} --commons 0 --nsample 2500 --submode 5".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
    if args.mode == 6:
        # simulation_list = glob.glob("dis_*")
        simulation_list = ['dis_62.0', 'dis_48.0', 'dis_34.0', 'dis_36.0', 'dis_60.0', 'dis_40.0', 'dis_78.0', 'dis_74.0', 'dis_66.0', 'dis_52.0', 'dis_94.0', 'dis_82.0', 'dis_98.0', 'dis_68.0', 'dis_92.0', 'dis_64.0', 'dis_42.0', 'dis_58.0', 'dis_90.0', 'dis_100.0', 'dis_32.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 0
            i_plus_one = i +1
            # do(f"mkdir -p log{i}")
            # do(f"mv log.* log{i}/")
            # do(f"cp log{i}/log.lammps .")
            # do(f"cp x.* log{i}/")
            # continueRunConvertion(n=12, rerun=i)
            # do(f"mkdir {i_plus_one}")

            # run_slurm = base_run_slurm.format(i_plus_one)
            # with open(f"run_{i_plus_one}.slurm", "w") as r:
            #     r.write(run_slurm)

            # do(f"sbatch run_{i_plus_one}.slurm")

            # do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
if args.day == "feb26":
    if args.mode == 1:
        pre = "/scratch/wl45/feb_2018/week_of_feb19/"
        data_folder = "/scratch/wl45/feb_2018/week_of_feb19/all_data_folder/"
        folder_list = ["rg_0.1_lipid_1.0_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data(pre, data_folder, folder_list, rerun=1, average_z=True)
    if args.mode == 2:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(3):
            freeEnergy_folder = f"freeEnergy_rg_0.1_lipid_1.0_mem_1_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["rg_0.1_lipid_1.0_mem_1"]
            # submode_list = ["_no_energy"]
            # submode_list = ["", "only_500"]
            # submode_list = ["350", "400", "450", "500", "550"]
            temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    move_data2(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=temp_mode, average_z=True)

            cd(freeEnergy_folder)
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    cd(folder+temp_mode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        make_metadata(temps_list=temp_list,k=0.02)
                        do("pulling_analysis.py -m {} --commons 0 --nsample 2500 --submode 5".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
if args.day == "feb25":
    if args.mode == 1:
        cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        sim_list = ["0"]
        for sim in sim_list:
            for folder in simulation_list:
                cd(folder)
                cd(sim)
                print(folder)
                for i in range(12):
                    compute_average_z(f"dump.lammpstrj.{i}", f"z_{i}.dat")
                cd("../..")
if args.day == "feb24":
    if args.mode == 1:
        simulation_list = glob.glob("dis_*")
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 1
            i_plus_one = i +1
            do(f"mkdir -p log{i}")
            do(f"mv log.* log{i}/")
            do(f"cp log{i}/log.lammps .")
            do(f"cp x.* log{i}/")
            continueRunConvertion(n=12, rerun=i)
            do(f"mkdir {i_plus_one}")

            run_slurm = base_run_slurm.format(i_plus_one)
            with open(f"run_{i_plus_one}.slurm", "w") as r:
                r.write(run_slurm)
            do(f"sbatch run_{i_plus_one}.slurm")

            # do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            # do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
if args.day == "feb22":
    if args.mode == 1:
        simulation_list = glob.glob("dis_*")
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 0
            i_plus_one = i +1
            do(f"mkdir -p log{i}")
            do(f"mv log.* log{i}/")
            do(f"cp log{i}/log.lammps .")
            do(f"cp x.* log{i}/")
            continueRunConvertion(n=12, rerun=i)
            do(f"mkdir {i_plus_one}")

            run_slurm = base_run_slurm.format(i_plus_one)
            with open(f"run_{i_plus_one}.slurm", "w") as r:
                r.write(run_slurm)
            do(f"sbatch run_{i_plus_one}.slurm")

            # do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            # do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")

if args.day == "feb19":
    if args.mode == 1:
        simulation_list = glob.glob("dis_*")
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 1
            i_plus_one = i +1
            # do(f"mkdir -p log{i}")
            # do(f"mv log.* log{i}/")
            # do(f"cp log{i}/log.lammps .")
            # do(f"cp x.* log{i}/")
            # continueRunConvertion(n=12, rerun=i)
            # do(f"mkdir {i_plus_one}")

            run_slurm = base_run_slurm.format(i_plus_one)
            with open(f"run_{i_plus_one}.slurm", "w") as r:
                r.write(run_slurm)
            do(f"sbatch run_{i_plus_one}.slurm")

            # do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            # do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
if args.day == "feb15":
    if args.mode == 1:
        cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        for folder in simulation_list:
            cd(folder)
            cd("1")
            print(folder)
            for i in range(12):
                compute_average_z(f"dump.lammpstrj.{i}", f"z_{i}.dat")
            cd("../..")
    if args.mode == 2:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(3):
            freeEnergy_folder = f"freeEnergy_rg_0.1_lipid_1.0_mem_1_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["rg_0.1_lipid_1.0_mem_1"]
            # submode_list = ["_no_energy"]
            # submode_list = ["", "only_500"]
            # submode_list = ["350", "400", "450", "500", "550"]
            temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    move_data2(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=temp_mode, average_z=True)

            cd(freeEnergy_folder)
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    cd(folder+temp_mode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        make_metadata(temps_list=temp_list,k=0.02)
                        do("pulling_analysis.py -m {} --commons 0 --nsample 2500 --submode 5".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
    if args.mode == 3:
        pre = "/scratch/wl45/feb_2018/week_of_feb05/"
        data_folder = "/scratch/wl45/feb_2018/week_of_feb12/all_data_folder/"
        folder_list = ["rg_0.1_lipid_1.0_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data(pre, data_folder, folder_list, rerun=1, average_z=True)
    if args.mode == 4:
        simulation_list = glob.glob("dis_*")
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 1
            i_plus_one = i +1
            do(f"mkdir -p log{i}")
            do(f"mv log.* log{i}/")
            do(f"cp log{i}/log.lammps .")
            do(f"cp x.* log{i}/")
            continueRunConvertion(n=12, rerun=i)
            do(f"mkdir {i_plus_one}")

            # run_slurm = base_run_slurm.format(i_plus_one)
            # with open(f"run_{i_plus_one}.slurm", "w") as r:
            #     r.write(run_slurm)
            # do(f"sbatch run_{i_plus_one}.slurm")

            # do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            # do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
    if args.mode == 5:
        print("how Constant force refolding")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        force_ramp_rate_list=[0.5]
        temperature_list=[500]
        memb_k_list = [1]
        rg_list = [0.1]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        # force_list = [0.22, 0.25, 0.28]
        force_list = [0.0, 0.01, 0.05, 0.08, 0.1, 0.15]
        repeat = 20
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=0,simulation_base_steps=2e7)
if args.day == "feb14":
    if args.mode == 1:
        # simulation_list = glob.glob("rg_0.1_")
        simulation_list = ["force_0.2_"]
        for folder in simulation_list:
            cd(folder)
            cd("simulation")
            for i in range(20):
                cd(str(i))
                cd("0")
                print(folder)
                compute_average_z(f"dump.lammpstrj", "z.dat")
                cd("../..")
            cd("../..")
    if args.mode == 2:
        simulation_list = glob.glob("dis_*")
        print(simulation_list)
        for dis in simulation_list:
            cd(dis)
            i = 0
            i_plus_one = i +1
            # do(f"mkdir -p log{i}")
            # do(f"mv log.* log{i}/")
            # do(f"cp log{i}/log.lammps .")
            # do(f"cp x.* log{i}/")
            # continueRunConvertion(n=12)
            # do(f"mkdir {i_plus_one}")
            run_slurm = base_run_slurm.format(i_plus_one)
            with open(f"run_{i_plus_one}.slurm", "w") as r:
                r.write(run_slurm)
            do(f"sbatch run_{i_plus_one}.slurm")
            # do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            # do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
if args.day == "feb13":
    if args.mode == 1:
        cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        for folder in simulation_list:
            cd(folder)
            cd("0")
            print(folder)
            for i in range(12):
                compute_average_z(f"dump.lammpstrj.{i}", f"z_{i}.dat")
            cd("../..")
    if args.mode == 2:
        print("Read variables data")
        read_variable_folder("/scratch/wl45/feb_2018/week_of_feb12/constant_force")
    if args.mode == 3:
        pre = "/scratch/wl45/feb_2018/week_of_feb05/"
        data_folder = "/scratch/wl45/feb_2018/week_of_feb12/all_data_folder/"
        folder_list = ["rg_0.1_lipid_1.0_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data(pre, data_folder, folder_list, rerun=-1, average_z=True)
    if args.mode == 4:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(1):
            freeEnergy_folder = f"freeEnergy_rg_0.1_lipid_1.0_mem_1_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["rg_0.1_lipid_1.0_mem_1"]
            # submode_list = ["_no_energy"]
            # submode_list = ["", "only_500"]
            # submode_list = ["350", "400", "450", "500", "550"]
            temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    move_data2(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=temp_mode, average_z=True)

            cd(freeEnergy_folder)
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    cd(folder+temp_mode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        make_metadata(temps_list=temp_list,k=0.02)
                        do("pulling_analysis.py -m {} --commons 0 --nsample 2500 --submode 5".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
if args.day == "feb12":
    if args.mode == 1:
        print("how Constant force unfolding")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        force_ramp_rate_list=[0.5]
        temperature_list=[500]
        memb_k_list = [1]
        rg_list = [0.1]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        # force_list = [0.22, 0.25, 0.28]
        force_list = [0.18, 0.20, 0.21]
        repeat = 20
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=0,simulation_base_steps=2e7)

if args.day == "feb07":
    if args.mode == 1:
        print("Read variables data")
        read_variable_folder("/scratch/wl45/feb_2018/week_of_feb05/refolding_pressure_1.0")
    if args.mode == 2:
        simulation_list = glob.glob("dis_*")
        print(simulation_list)
        for dis in simulation_list:
            cd(dis)
            i = 0
            i_plus_one = i +1
            do(f"mkdir -p log{i}")
            do(f"mv log.* log{i}/")
            do(f"cp log{i}/log.lammps .")
            do(f"cp x.* log{i}/")
            continueRunConvertion(n=12)
            do(f"mkdir {i_plus_one}")
            do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
if args.day == "feb06":
    if args.mode == 1:
        print("Read variables data")
        read_variable_folder("/scratch/wl45/feb_2018/week_of_feb05/pulling_pressure_1.0")
    if args.mode == 2:
        print("Read variables data")
        read_variable_folder("/scratch/wl45/feb_2018/week_of_feb05/different_rg_refold")
    if args.mode == 3:
        print("how Refolding change")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        force_ramp_rate_list=[0.5]
        temperature_list=[500]
        memb_k_list = [1]
        rg_list = [0.1, 0.15, 0.2, 0.3]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        force_list = [0]
        repeat = 20
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=0,simulation_base_steps=2e7)
    if args.mode == 4:
        simulation_list = glob.glob("rg_0.1_")
        for folder in simulation_list:
            cd(folder)
            cd("simulation")
            for i in range(20):
                cd(str(i))
                cd("0")
                print(folder)
                compute_average_z(f"dump.lammpstrj", "z.dat")
                cd("../..")
            cd("../..")
if args.day == "feb05":
    if args.mode == 1:
        print("how Unfolding change")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1.0]
        # pressure_list = [0.6, 1]
        force_ramp_rate_list=[0.5]
        temperature_list=[500]
        memb_k_list = [1]
        # rg_list = [0.3]
        rg_list = [0.1, 0.15, 0.2, 0.3]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        force_list = ["ramp"]
        repeat = 20
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=0,simulation_base_steps=1e7)
if args.day == "feb04":
    if args.mode == 1:
        print("Read variables data")
        read_variable_folder("/scratch/wl45/jan_2018/week_of_jan29/pulling")
    if args.mode == 2:
        print("how Refolding change")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0.6]
        force_ramp_rate_list=[0.5]
        temperature_list=[500]
        memb_k_list = [1]
        rg_list = [0.1, 0.2, 0.3]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        force_list = [0]
        repeat = 20
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=0,simulation_base_steps=2e7)
if args.day == "feb02":
    if args.mode == 1:
        print("how Unfolding change")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0.6]
        # pressure_list = [0.6, 1]
        force_ramp_rate_list=[0.5]
        temperature_list=[500]
        memb_k_list = [1]
        # rg_list = [0.3]
        rg_list = [0.15, 0.2]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        force_list = ["ramp"]
        repeat = 20
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=0,simulation_base_steps=2e7)
if args.day == "feb01":
    if args.mode == 1:
        print("how Refolding change")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0.6]
        force_ramp_rate_list=[0.5]
        temperature_list=[500]
        memb_k_list = [1]
        rg_list = [0.3]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        force_list = [0, 0.05, 0.1]
        repeat = 20
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=0,simulation_base_steps=2e7)
    if args.mode == 2:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(1):
            freeEnergy_folder = f"freeEnergy_rg_0.3_lipid_0.6_mem_1_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["rg_0.3_lipid_0.6_mem_1"]
            # submode_list = ["_no_energy"]
            # submode_list = ["", "only_500"]
            # submode_list = ["350", "400", "450", "500", "550"]
            temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    move_data2(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=temp_mode, average_z=True)

            cd(freeEnergy_folder)
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    cd(folder+temp_mode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        make_metadata(temps_list=temp_list,k=0.02)
                        do("pulling_analysis.py -m {} --commons 0 --nsample 2500 --submode 5".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
if args.day == "jan31":
    if args.mode == 1:
        pre = "/scratch/wl45/jan_2018/week_of_jan29/"
        data_folder = "/scratch/wl45/jan_2018/week_of_jan29/all_data_folder/"
        folder_list = ["rg_0.3_lipid_0.6_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data(pre, data_folder, folder_list, rerun=-1, average_z=True)
    if args.mode == 2:
        cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        for folder in simulation_list:
            cd(folder)
            cd("0")
            print(folder)
            for i in range(12):
                compute_average_z(f"dump.lammpstrj.{i}", f"z_{i}.dat")
            cd("../..")
    if args.mode == 3:
        print("how Unfolding change")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0.6]
        force_ramp_rate_list=[0.5]
        temperature_list=[500]
        memb_k_list = [1]
        rg_list = [0.3]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        force_list = ["ramp"]
        repeat = 20
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=0,simulation_base_steps=2e7)
if args.day == "jan25":
    if args.mode == 1:
        simulation_list = glob.glob("dis_*")
        print(simulation_list)
        for dis in simulation_list:
            cd(dis)
            do("mkdir -p log1")
            do("mv log.* log1/")
            do("cp log1/log.lammps .")
            do(f"cp x.* log1/")
            cd("..")
    if args.mode == 2:
        cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        for folder in simulation_list:
            cd(folder)
            cd("1")
            print(folder)
            for i in range(12):
                compute_average_z(f"dump.lammpstrj.{i}", f"z_{i}.dat")
            cd("../..")
    if args.mode == 3:
        pre = "/scratch/wl45/jan_2018/week_of_jan22/"
        data_folder = "/scratch/wl45/jan_2018/week_of_jan22/all_data_folder/"
        folder_list = ["rg_0.3_lipid_0.6_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data(pre, data_folder, folder_list, rerun=2, average_z=True)
    if args.mode == 4:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(3):
            freeEnergy_folder = f"freeEnergy_rg_0.3_lipid_0.6_mem_1_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["rg_0.3_lipid_0.6_mem_1"]
            # submode_list = ["_no_energy"]
            # submode_list = ["", "only_500"]
            # submode_list = ["350", "400", "450", "500", "550"]
            temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    move_data2(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=temp_mode, average_z=True)

            cd(freeEnergy_folder)
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    cd(folder+temp_mode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        make_metadata(temps_list=temp_list,k=0.02)
                        do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 5".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
    if args.mode == 5:
        print("how Refolding change")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0.6]
        force_ramp_rate_list=[0.5]
        temperature_list=[500]
        memb_k_list = [1]
        rg_list = [0.3]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        force_list = ["ramp"]
        repeat = 20
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=0,simulation_base_steps=2e7)
    if args.mode == 6:
        print("how Unfolding change")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0.6]
        force_ramp_rate_list=[0.5]
        temperature_list=[500]
        memb_k_list = [1]
        rg_list = [0.3]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        force_list = ["ramp"]
        repeat = 20
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        pressure_list=pressure_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=0)
if args.day == "jan24":
    if args.mode == 1:
        print("how unfolding change")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        # pressure_list = [0.6]
        force_ramp_rate_list=[1]
        temperature_list=[500]
        memb_k_list = [1, 2]
        rg_list = [0.3]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        force_list = ["ramp"]
        repeat = 3
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=1)
    if args.mode == 2:
        print("how unfolding change")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        # pressure_list = [0.6]
        force_ramp_rate_list=[1]
        temperature_list=[400, 500]
        memb_k_list = [1]
        rg_list = [0.3]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        # qbias_list = list(np.linspace(0.2,0.9,36))
        force_list = ["ramp"]
        repeat = 20
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=0)
if args.day == "jan22":
    if args.mode == 1:
        simulation_list = glob.glob("dis_*")
        print(simulation_list)
        for dis in simulation_list:
            cd(dis)
            do("mkdir -p log0")
            do("mv log.* log0/")
            do("cp log0/log.lammps .")
            do(f"cp x.* log0/")
            continueRunConvertion(n=12)
            do("mkdir 1")
            do("sed 's/2xov_0/2xov_1/g' run_0.slurm > run_1.slurm")
            do("sbatch run_1.slurm")
            cd("..")
if args.day == "jan20":
    if args.mode == 1:
        pre = "/scratch/wl45/jan_2018/"
        data_folder = "/scratch/wl45/jan_2018/all_data_folder/"
        folder_list = ["rg_0.3_lipid_0.6_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list)
    if args.mode == 2:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(1):
            freeEnergy_folder = f"freeEnergy_rg_0.3_lipid_0.6_mem_1_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["rg_0.3_lipid_0.6_mem_1"]
            # submode_list = ["_no_energy"]
            # submode_list = ["", "only_500"]
            # submode_list = ["350", "400", "450", "500", "550"]
            temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    move_data(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=temp_mode)

            cd(freeEnergy_folder)
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    cd(folder+temp_mode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        make_metadata(temps_list=temp_list,k=0.02)
                        do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
    if args.mode == 3:
        cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        for folder in simulation_list:
            cd(folder)
            cd("0")
            print(folder)
            for i in range(12):
                compute_average_z(f"dump.lammpstrj.{i}", f"z_{i}.dat")
            cd("../..")
if args.day == "dec04":
    i = 12
    if args.mode == 2:
        os.system("cp ../2xov.pdb .")
        with open(f"qnqc.slurm", "w") as f:
            f.write(qnqc_slurm.format(i))
        os.system(f"sbatch qnqc.slurm")
    if args.mode == 1:
        simulation_list = glob.glob("dis_*")
        print(simulation_list)
        for dis in simulation_list:
            cd(dis)
            for ii in range(2):
                cd(str(ii))
                os.system("cp ../2xov.pdb .")
                with open(f"qnqc.slurm", "w") as f:
                    f.write(qnqc_slurm.format(i))
                os.system(f"sbatch qnqc.slurm")
                cd("..")
            cd("..")

if args.day == "nov29":
    if args.mode == 2:
        simulation_list = glob.glob("dis_*")
        print(simulation_list)
        for dis in simulation_list:
            cd(dis)
            do("mkdir -p log1")
            do("mv log.* log1/")
            cd("..")
        print("hello!")
        pre = "/scratch/wl45/nov_2017/13nov/"
        data_folder = "/scratch/wl45/nov_2017/13nov/all_data_folder/"
        folder_list = ["no_side_contraint_memb_3_rg_0.4_lipid_0.6_extended"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list, rerun=2)
    if args.mode == 1:
        simulation_list = glob.glob("qbias_*")
        pre = "/scratch/wl45/nov_2017/27nov/"
        data_folder = "/scratch/wl45/nov_2017/27nov/all_data_folder/"
        folder_list = ["q_bias_temper_new"]
        print(simulation_list)
        for dis in simulation_list:
            cd(dis)
            do("mkdir -p log1")
            do("mv log.* log1/")
            cd("..")
        process_temper_data(pre, data_folder, folder_list, n=10, bias="qbias", rerun=2)
if args.day == "nov28":
    if args.mode == 1:
        print("hello!")
        pre = "/scratch/wl45/nov_2017/13nov/"
        data_folder = "/scratch/wl45/nov_2017/13nov/all_data_folder/"
        folder_list = ["bias_0.05_memb_3_rg_0.4_lipid_0.6_extended"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list)
if args.day == "nov27":
    if args.mode == 1:
        simulation_list = glob.glob("qbias_*")
        print(simulation_list)

        for dis in simulation_list:
            cd(dis)
            do("mkdir -p log0")
            do("mv log.* log0/")
            do("cp log0/log.lammps .")
            do(f"cp x.* log0/")
            continueRunConvertion(n=10)
            do("mkdir 1")

            do("sed 's/2xov_0/2xov_1/g' run_0.slurm > run_1.slurm")
            do("sbatch run_1.slurm")
            cd("..")
if args.day == "nov25":
    if args.mode == 1:
        print("hello!")
        pre = "/scratch/wl45/nov_2017/27nov/"
        data_folder = "/scratch/wl45/nov_2017/27nov/all_data_folder/"
        folder_list = ["q_bias_temper_new"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list, n=10, bias="qbias")

if args.day == "nov23":
    if args.mode == 4:
        print("how unfolding change")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        force_ramp_rate_list=[2]
        temperature_list=[500, 600, 650]
        memb_k_list = [3]
        rg_list = [0.3]
        # qbias_list = [0.25, 0.45, 0.65, 0.85]
        qbias_list = list(np.linspace(0.2,0.9,36))
        force_list = ["ramp"]
        repeat = 3
        variable_test2(temperature_list=temperature_list,
                        qbias_list=qbias_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=1)
    if args.mode == 3:
        simulation_list = glob.glob("dis_*")
        # print(simulation_list)
        for dis in simulation_list:
            do(f"mkdir -p {dis}/log{1}")
            do(f"cp {dis}/log.* {dis}/log{1}/")
    if args.mode == 2:
        print("hello!")
        pre = "/scratch/wl45/nov_2017/13nov/"
        data_folder = "/scratch/wl45/nov_2017/13nov/all_data_folder/"
        folder_list = ["no_side_contraint_memb_3_rg_0.4_lipid_0.6_extended"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list, rerun=2)
    if args.mode == 1:
        print("how unfolding change")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        force_ramp_rate_list=[2]
        temperature_list=[500]
        memb_k_list = [0,1,2,4,8]
        rg_list = [0, 0.1, 0.2, 0.3, 0.4, 0.8]
        force_list = ["ramp"]
        repeat = 10
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        mem_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=1)

if args.day == "nov21":
    if args.mode == 3:
        print("hello!")
        pre = "/scratch/wl45/nov_2017/13nov/"
        data_folder = "/scratch/wl45/nov_2017/13nov/all_data_folder/"
        folder_list = ["no_side_contraint_memb_3_rg_0.4_lipid_0.6_extended"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list, rerun=1)
    if args.mode == 2:
        dis_list = glob.glob("dis_*")
        print(dis_list)
        for dis in dis_list:
            num = float(dis.split("_")[-1])
            if num <= 90:
                print(num, "large")
                continue
            cd(dis)
            # do("mkdir log0")
            # do("mv log.* log0/")
            do("cp log0/log.lammps .")
            continueRunConvertion()
            do("mkdir 1")

            do("sed 's/2xov_0/2xov_1/g' run_0.slurm > run_1.slurm")
            do("sbatch run_1.slurm")
            cd("..")
    if args.mode == 1:
        simulation_list = glob.glob("dis_*")
        print(simulation_list)

        for dis in simulation_list:
            num = float(dis.split("_")[-1])
            if num <= 90:
                print(num, "large")
                continue
            do(f"mkdir {dis}/log{0}")
            do(f"cp {dis}/x.* {dis}/log{0}")
            do(f"mv {dis}/log.* {dis}/log{0}/")

if args.day == "nov19":
    if args.mode ==1:
        rerun()
if args.day == "nov18":
    if args.mode == 1:
        print("hello!")
        pre = "/scratch/wl45/nov_2017/06nov/"
        data_folder = "/scratch/wl45/nov_2017/06nov/all_data_folder_nov15/"
        folder_list = ["next_gen_native_based_memb_3_rg_0.4_lipid_0.6_extended"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list, rerun=2)
if args.day == "nov15":
    if args.mode == 6:
        dis_list = glob.glob("dis_*")
        print(dis_list)
        for dis in dis_list:
            cd(dis)
            # do("mkdir log0")
            # do("mv log.* log0/")
            do("cp log0/log.lammps .")
            continueRunConvertion()
            do("mkdir 1")

            do("sed 's/2xov_0/2xov_1/g' run_0.slurm > run_1.slurm")
            do("sbatch run_1.slurm")
            cd("..")
    if args.mode == 1:
        simulation_list = glob.glob("dis_*")
        # print(simulation_list)
        for dis in simulation_list:
            do(f"mkdir {dis}/log{1}")
            do(f"mv {dis}/log.* {dis}/log{1}/")
    if args.mode == 2:
        print("hello!")
        pre = "/scratch/wl45/nov_2017/06nov/"
        data_folder = "/scratch/wl45/nov_2017/06nov/all_data_folder_nov15/"
        folder_list = ["23oct/memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list, rerun=2)
    if args.mode == 3:
        print("hello!")
        pre = "/scratch/wl45/nov_2017/06nov/"
        data_folder = "/scratch/wl45/nov_2017/06nov/all_data_folder_nov15/"
        folder_list = ["23oct/memb_4_rg_0.1_lipid_1_extended",
                        "23oct/memb_4_rg_0.1_lipid_1_topology"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list)
    if args.mode == 4:
        print("hello!")
        pre = "/scratch/wl45/nov_2017/06nov/"
        data_folder = "/scratch/wl45/nov_2017/06nov/all_data_folder_nov15/"
        folder_list = ["23oct/rgWidth_memb_3_rg_0.1_lipid_1_extended",
                        "23oct/rgWidth_memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list)
    if args.mode == 5:
        print("hello!")
        pre = "/scratch/wl45/nov_2017/06nov/"
        data_folder = "/scratch/wl45/nov_2017/06nov/all_data_folder_nov15/"
        folder_list = ["expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended",
                        "next_gen_native_based_memb_3_rg_0.2_lipid_0.6_extended",
                        "next_gen_native_based_memb_3_rg_0.4_lipid_0.6_extended",
                        "next_gen_native_based_memb_3_rg_0.4_lipid_0.6_topology",
                        "stronger_bias_for_expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list)


if args.day == "nov14":
    if args.mode == 1:
        print("hello!")
        pre = "/scratch/wl45/nov_2017/06nov/"
        data_folder = "/scratch/wl45/nov_2017/06nov/all_data_folder/"
        folder_list = ["next_gen_native_based_memb_3_rg_0.4_lipid_0.6_topology"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list)
    if args.mode == 2:
        print("hello!")
        pre = "/scratch/wl45/nov_2017/06nov/"
        data_folder = "/scratch/wl45/nov_2017/06nov/all_data_folder/"
        folder_list = ["stronger_bias_for_expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list)
if args.day == "nov12":
    if args.mode == 2:
        rerun(offset=1)
    if args.mode == 1:
        print("hello")
        pre = "/scratch/wl45/nov_2017/06nov/"
        data_folder = "/scratch/wl45/nov_2017/06nov/all_data_folder/"
        folder_list = ["next_gen_native_based_memb_3_rg_0.4_lipid_0.6_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list)
if args.day == "nov11":
    if args.mode == 2:
        rerun(offset=1)
    if args.mode == 1:
        print("hello")
        pre = "/scratch/wl45/nov_2017/06nov/"
        data_folder = "/scratch/wl45/nov_2017/06nov/all_data_folder/"
        folder_list = ["next_gen_native_based_memb_3_rg_0.2_lipid_0.6_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list)
if args.day == "nov10":
    if args.mode == 1:
        print("how unfolding change")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0.8, 2]
        force_ramp_rate_list=[1]
        temperature_list=[500]
        memb_k_list = [3]
        rg_list = [0.3, 0.4]
        force_list = [0.6]
        repeat = 5
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=1)

if args.day == "nov08":
    if args.mode == 5:
        rerun()
    if args.mode == 4:
        print("how unfolding change")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [1, 3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0]
        force_ramp_rate_list=[10]
        temperature_list=[500]
        memb_k_list = [3]
        rg_list = [0.1]
        force_list = [0.4, 0.3, 0.2]
        repeat = 10
        variable_test2(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=1)
    if args.mode == 3:
        print("how unfolding change")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0]
        force_ramp_rate_list=[20]
        temperature_list=[500]
        memb_k_list = [3]
        rg_list = [0.1]
        # force_list = [0.4, 0.5, 0.6, 0.7, 0.8]
        force_list = [0.3, 0.2, 0.1, 0.0]
        repeat = 2
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=1)
    if args.mode == 2:
        print("how unfolding change")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0]
        force_ramp_rate_list=[10]
        temperature_list=[500]
        memb_k_list = [2, 4]
        rg_list = [0.1, 0.4]
        force_list = [0.0]
        repeat = 2
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=1)
    if args.mode == 1:
        print("hello")
        pre = "/scratch/wl45/nov_2017/06nov/"
        data_folder = "/scratch/wl45/nov_2017/06nov/all_data_folder/"
        folder_list = ["23oct/rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list)

if args.day == "nov05":
    if args.mode == 1:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended", "memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended", "rgWidth_memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        temp_list = ["all"]
        # bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        bias_list = {"2d_qw_dis":"11"}
        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    do("make_metadata.py -m 12")
                    do("mkdir t_" + temp)
                    do("mv metadatafile t_" + temp)
                    cd("t_" + temp)
                    do("pulling_analysis.py -m {} --commons 0 --nsample 2500 --submode 1".format(mode))
                    cd("..")
                cd("..")
            cd("..")

if args.day == "nov02":
    if args.mode == 1:
        for i in range(-10, 15, 1):
            rerun(offset=i)

if args.day == "nov01":
    if args.mode == 3:
        # location_list = ["strengthen_helix_1", "strengthen_helix_1_baseline_without_strengthen"]
        location_list = ["strengthen_helix_1_baseline_without_strengthen"]
        pre = os.getcwd() + "/"
        for location in location_list:
            folder_list = glob.glob(pathname=pre + location + "/*_")
            for folder in folder_list:
                print(folder)
                for i in range(10):
                    cd(folder + "/simulation/{}".format(i))
                    rerun()
    if args.mode == 2:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        folder_list = ["memb_3_rg_0.1_lipid_1_extended", "memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended", "rgWidth_memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "500"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                cd(bias)
                cmd = "find -name 'slurm-*' | rev | awk -F'[-.]' '{print $2}' | rev"
                lines = getFromTerminal(cmd).splitlines()
                for line in lines:
                    # print(line)
                    do("scancel " + line)
                cd("..")
                do("rm -r " + bias)
            cd("..")
        # cmd = "find -name 'slurm-*' | rev | awk -F'[-.]' '{print $2}' | rev"
        # lines = getFromTerminal(cmd).splitlines()
        # for line in lines:
        #     print(line)
        #     do("scancel " + line)
    if args.mode == 1:
        print("how unfolding change")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1]
        force_ramp_rate_list=[1]
        temperature_list=[350, 500]
        memb_k_list = [2, 4]
        rg_list = [0.1, 0.4]
        force_list = [0.0]
        repeat = 10
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=1)
if args.day == "oct31":
    if args.mode == 4:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        folder_list = ["memb_3_rg_0.1_lipid_1_extended", "memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended", "rgWidth_memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "500"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    if temp != "all":
                        do("make_metadata.py -m 10 --submode 3 -t " + temp)
                    else:
                        do("make_metadata.py -m 9 --submode 3")
                    do("mkdir t_" + temp)
                    do("mv metadatafile t_" + temp)
                    cd("t_" + temp)
                    do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 1".format(mode))
                    cd("..")
                cd("..")
            cd("..")
    if args.mode == 3:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        # folder_list = ["memb_2_rg_0.1_lipid_1_extended", "memb_2_rg_0.1_lipid_1_topology"]
        folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended", "rgWidth_memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "500"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    if temp != "all":
                        do("make_metadata.py -m 10 --submode 3 -t " + temp)
                    else:
                        do("make_metadata.py -m 9 --submode 3")
                    do("mkdir t_" + temp)
                    do("mv metadatafile t_" + temp)
                    cd("t_" + temp)
                    do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 1".format(mode))
                    cd("..")
                cd("..")
            cd("..")
    if args.mode == 2:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        # folder_list = ["memb_2_rg_0.1_lipid_1_extended", "memb_2_rg_0.1_lipid_1_topology"]
        folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended", "rgWidth_memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "500"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("rm -r " + bias)
            cd("..")
    if args.mode == 1:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        # folder_list = ["memb_2_rg_0.1_lipid_1_extended", "memb_2_rg_0.1_lipid_1_topology"]
        folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended", "rgWidth_memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "500"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    if temp != "all":
                        do("make_metadata.py -m 10 --submode 1 -t " + temp)
                    else:
                        do("make_metadata.py -m 9 --submode 1")
                    do("mkdir t_" + temp)
                    do("mv metadatafile t_" + temp)
                    cd("t_" + temp)
                    do("pulling_analysis.py -m {} --commons 0 --nsample 2500 --submode 1".format(mode))
                    cd("..")
                cd("..")
            cd("..")

if args.day == "oct28":
    if args.mode == 1:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        # folder_list = ["memb_2_rg_0.1_lipid_1_extended", "memb_2_rg_0.1_lipid_1_topology"]
        folder_list = ["memb_3_rg_0.1_lipid_1_extended", "memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "500"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    if temp != "all":
                        do("make_metadata.py -m 10 --submode 1 -t " + temp)
                    else:
                        do("make_metadata.py -m 9 --submode 1")
                    do("mkdir t_" + temp)
                    do("mv metadatafile t_" + temp)
                    cd("t_" + temp)
                    do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 1".format(mode))
                    cd("..")
                cd("..")
            cd("..")

if args.day == "oct25":
    if args.mode == 1:
        compute_theta_for_each_helix()
        print("Done")
    if args.mode == 2:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        # folder_list = ["memb_2_rg_0.1_lipid_1_extended", "memb_2_rg_0.1_lipid_1_topology"]
        folder_list = ["memb_3_rg_0.1_lipid_1_extended", "memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "500"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    if temp != "all":
                        do("make_metadata.py -m 10 -t " + temp)
                    else:
                        do("make_metadata.py -m 9")
                    do("mkdir t_" + temp)
                    do("mv metadatafile t_" + temp)
                    cd("t_" + temp)
                    do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 1".format(mode))
                    cd("..")
                cd("..")
            cd("..")
if args.day == "oct24":
    if args.mode == 2:
        cmd = 'tail -n 1 log.lammps | cut -d" " -f2-'
        line = getFromTerminal(cmd).rstrip()
        print(line)
    if args.mode == 1:
        dis_list = glob.glob("dis_*")
        print(dis_list)
        for dis in dis_list:
            cd(dis)
            do("cp log0/log.lammps .")
            continueRunConvertion()
            # do("mkdir log0")
            # do("mv log.* log0/")
            # do("mkdir 1")

            do("sed 's/2xov_0/2xov_1/g' run_0.slurm > run_1.slurm")
            do("sbatch run_1.slurm")
            cd("..")
        # continueRunConvertion()

if args.day == "oct21":
    if args.mode == 1:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        # folder_list = ["memb_2_rg_0.1_lipid_1_extended", "memb_2_rg_0.1_lipid_1_topology"]
        folder_list = ["memb_3_rg_0.1_lipid_1_extended", "memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "450"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    if temp != "all":
                        do("make_metadata.py -m 10 -t " + temp)
                    else:
                        do("make_metadata.py -m 9")
                    do("mkdir t_" + temp)
                    do("mv metadatafile t_" + temp)
                    cd("t_" + temp)
                    do("pulling_analysis.py -m {} --commons 0 --nsample 2500".format(mode))
                    cd("..")
                cd("..")
            cd("..")
    if args.mode == 2:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        folder_list = ["memb_2_rg_0.1_lipid_1_extended", "memb_2_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["important_all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    do("make_metadata.py -m 9")
                    do("rm t_t_" +temp)
                    do("mkdir -p t_t_" + temp)
                    do("mv metadatafile t_t_" + temp)
                    cd("t_t_" + temp)
                    do("pulling_analysis.py -m {} --commons 0 --nsample 2500".format(mode))
                    cd("..")
                cd("..")
            cd("..")
if args.day == "oct17":
    if args.mode == 1:
        print("can it fold")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 1, 2]
        force_ramp_rate_list=[1]
        temperature_list=[350, 500]
        memb_k_list = [1, 2, 4]
        rg_list = [0, 0.1, 0.4]
        force_list = [0.0]
        repeat = 2
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=1)
if args.day == "oct16":
    if args.mode == 1:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        folder_list = ["more_higher_temp"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "350", "400", "450", "500"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    if temp != "all":
                        do("make_metadata.py -m 10 -t " + temp)
                    else:
                        do("make_metadata.py -m 9")
                    do("mkdir t" + temp)
                    do("mv metadatafile t" + temp)
                    cd("t" + temp)
                    do("pulling_analysis.py -m {} --commons 1 --nsample 4000".format(mode))
                    cd("..")
                cd("..")
            cd("..")
if args.day == "oct13":
    if args.mode == 2:
        print("strong membrane fore ramp")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1, 2]
        force_ramp_rate_list=[1]
        temperature_list=[350]
        memb_k_list = [1, 2, 4, 8, 16]
        rg_list = [0, 0.4]
        force_list = [0.0]
        repeat = 10
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=0)
    if args.mode == 1:
        folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        temp_list = ["all", "350", "400", "450", "500"]
        for folder in folder_list:
            cd(folder)

            do("mkdir -p 2d_qw_dis")
            cd("2d_qw_dis")
            for temp in temp_list:
                if temp != "all":
                    do("make_metadata.py -m 10 -t " + temp)
                else:
                    do("make_metadata.py -m 9")
                do("mkdir t" + temp)
                do("mv metadatafile t" + temp)
                cd("t" + temp)
                do("pulling_analysis.py -m 11 --commons 1")
                cd("..")
            cd("..")
            cd("..")
if args.day == "oct12":
    if args.mode == 1:
        temp_list = ["400", "450", "500", "all"]
        for temp in temp_list:
            if temp != "all":
                do("make_metadata.py -m 10 -t " + temp)
            else:
                do("make_metadata.py -m 9")
            do("mkdir t" + temp)
            do("mv metadatafile t" + temp)
            cd("t" + temp)
            do("pulling_analysis.py -m 9")
            cd("..")
    if args.mode == 2:
        folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        temp_list = ["350", "400", "450", "500", "all"]
        for folder in folder_list:
            cd(folder)
            # do("mkdir -p 1d_dis")
            # cd("1d_dis")
            # for temp in temp_list:
            #     if temp != "all":
            #         do("make_metadata.py -m 10 -t " + temp)
            #     else:
            #         do("make_metadata.py -m 9")
            #     do("mkdir t" + temp)
            #     do("mv metadatafile t" + temp)
            #     cd("t" + temp)
            #     do("pulling_analysis.py -m 9")
            #     cd("..")
            # cd("..")
            do("mkdir -p 1d_qw")
            cd("1d_qw")
            for temp in temp_list:
                if temp != "all":
                    do("make_metadata.py -m 10 -t " + temp)
                else:
                    do("make_metadata.py -m 9")
                do("mkdir t" + temp)
                do("mv metadatafile t" + temp)
                cd("t" + temp)
                do("pulling_analysis.py -m 10")
                cd("..")
            cd("..")
            cd("..")
    if args.mode == 3:
        print("high temp refold")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [2]
        force_ramp_rate_list=[1000]
        temperature_list=[100]
        memb_k_list = [1, 2, 4, 8, 16, 32]
        rg_list = [0.4]
        force_list = [0.0]
        repeat = 1
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=1)

if args.day == "oct10":
    if args.mode == 5:
        do("mkdir combined")
        cd("combined")
        do("make_metadata.py -m 9")
        do("pulling_analysis.py -m 9 --commons 1")
    if args.mode == 4:
        temp_list = [350, 400, 450, 500]
        for temp in temp_list:
            do("mkdir t" + str(temp))
            cd("t" + str(temp))
            do("make_metadata.py -m 10 -t {}".format(temp))
            do("pulling_analysis.py -m 9 --commons 1")
            cd("..")
    if args.mode == 3:
        do("make_metadata.py -m 10 -t 400")
        do("pulling_analysis.py -m 9")
    if args.mode == 2:
        print("high temp refold")
        start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [2]
        force_ramp_rate_list=[1]
        temperature_list=[500, 550, 600]
        memb_k_list = [1]
        rg_list = [0.4]
        force_list = [0.0]
        repeat = 50
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat)
    if args.mode == 1:
        fold_list = ["rg_0.4_lipid_2_temp_300_extended", "rg_0.4_lipid_2_temp_300_topology",
                        "rg_0.4_lipid_2_temp_400_extended", "rg_0.4_lipid_2_temp_400_topology"]
        for folder in fold_list:
            cd(folder)
            cd("simulation")
            do("pulling_prepare.py -m 8 --submode 1")
            cd("..")
            do("mkdir freeEnergy_2")
            cd("freeEnergy_2")
            do("make_metadata.py -m 3 -k 0.02 -t 300 --submode 1")
            do("pulling_analysis.py -m 6")
            cd("../..")

if args.day == "oct05":
    if args.mode == 5:
        print("no go, constant force, refold")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [2]
        force_ramp_rate_list=[0.5]
        temperature_list=[400]
        memb_k_list = [0, 1, 2, 4]
        rg_list = [0.4]
        force_list = [0.0]
        repeat = 10
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat)
        # # start_from_list=["native", "extended", "topology"]
        # # start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        # mode_list = [3]  # lipid mediated interaction
        # # pressure_list = [0, 0.1, 1.0]
        # pressure_list = [0, 2]
        # force_ramp_rate_list=[0.5]
        # temperature_list=[400]
        # memb_k_list = [0, 1, 2, 4]
        # rg_list = [0, 0.4]
        # force_list = [0.0]
        # repeat = 10
        # variable_test(temperature_list=temperature_list,
        #                 start_from_list=start_from_list,
        #                 rg_list=rg_list,
        #                 memb_k_list=memb_k_list,
        #                 mode_list=mode_list,
        #                 pressure_list=pressure_list,
        #                 force_ramp_rate_list=force_ramp_rate_list,
        #                 force_list=force_list,
        #                 repeat=repeat)

    if args.mode == 4:
        print("membrane effect")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 0.1, 1, 2]
        force_ramp_rate_list=[100]
        temperature_list=[200, 300, 400, 500]
        memb_k_list = [1, 2, 4, 8, 16]
        rg_list = [0, 0.08, 0.4]
        force_list =["ramp"]
        repeat = 2
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat)
    if args.mode == 1:
        fold_list = ["rg_0.4_lipid_2_temp_300_extended", "rg_0.4_lipid_2_temp_300_topology",
                        "rg_0.4_lipid_2_temp_400_extended", "rg_0.4_lipid_2_temp_400_topology"]
        for folder in fold_list:
            cd(folder)
            cd("simulation")
            do("pulling_prepare.py -m 8 --submode 1")
            cd("..")
            do("mkdir freeEnergy_2")
            cd("freeEnergy_2")
            do("make_metadata.py -m 3 -k 0.02 -t 300 --submode 1")
            do("pulling_analysis.py -m 6")
            cd("../..")
    if args.mode == 2:
        print("2d wham")
        # cd("simulation")
        # do("pulling_prepare.py -m 8")
        # cd("..")
        fold_list = ["rg_0.4_lipid_2_temp_300_extended", "rg_0.4_lipid_2_temp_300_topology",
                        "rg_0.4_lipid_2_temp_400_extended", "rg_0.4_lipid_2_temp_400_topology"]
        for folder in fold_list:
            cd(folder)
            name = "wham2d"
            do("mkdir {}".format(name))
            cd(name)
            do("make_metadata.py -m 3 -k 0.02 -t 300 --submode 1")
            do("pulling_analysis.py -m 5")
            cd("../..")
    if args.mode == 3:
        print("qw wham")
        # cd("simulation")
        # do("pulling_prepare.py -m 8")
        # cd("..")
        fold_list = ["rg_0.4_lipid_2_temp_300_extended", "rg_0.4_lipid_2_temp_300_topology",
                        "rg_0.4_lipid_2_temp_400_extended", "rg_0.4_lipid_2_temp_400_topology"]
        for folder in fold_list:
            cd(folder)
            name = "qw_2"
            do("mkdir {}".format(name))
            cd(name)
            do("make_metadata.py -m 3 -k 0.02 -t 300 --submode 1")
            do("pulling_analysis.py -m 8")
            cd("../..")
if args.day == "oct04":
    if args.mode == 1:
        # cd("simulation")
        # do("pulling_prepare.py -m 8")
        # cd("..")
        do("mv freeEnergy old_freeEnergy")
        do("mkdir freeEnergy")
        cd("freeEnergy")
        do("make_metadata.py -m 3 -k 0.02 -t 300")
        do("pulling_analysis.py -m 6")
    if args.mode == 2:
        # cd("simulation")
        # do("pulling_prepare.py -m 8")
        # cd("..")
        fold_list = ["rg_0.4_lipid_2_temp_300_extended", "rg_0.4_lipid_2_temp_300_topology",
                        "rg_0.4_lipid_2_temp_400_extended", "rg_0.4_lipid_2_temp_400_topology"]
        for folder in fold_list:
            cd(folder)
            name = "qw"
            do("mkdir {}".format(name))
            cd(name)
            do("make_metadata.py -m 3 -k 0.02 -t 300")
            do("pulling_analysis.py -m 8")
            cd("../..")
if args.day == "oct03":
    if args.mode == 1:
        print("p3 force ramp")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0]
        force_ramp_rate_list=[0.5]
        temperature_list=[200, 300, 400]
        memb_k_list = [1]
        rg_list = [0, 0.08, 0.4]
        force_list =["ramp"]
        repeat = 5
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat)
    if args.mode == 2:
        print("p3 folding temperature")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 2]
        force_ramp_rate_list=[0.5]
        temperature_list=["ramp"]
        memb_k_list = [1, 2, 4]
        rg_list = [0, 0.08, 0.4]
        force_list =["ramp"]
        repeat = 2
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat)
if args.day == "sep26":
    if args.mode == 1:
        print("constant force")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0]
        force_ramp_rate_list=[0.5]
        temperature_list=[200, 300, 400]
        memb_k_list = [1]
        rg_list = [0, 0.08, 0.4]
        force_list =[0.1, 0.2, 0.25, 0.3, 0.4]
        repeat = 5
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat)
if args.day == "sep25":
    if args.mode == 5:
        print("constant force")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [2]
        force_ramp_rate_list=[0.5]
        temperature_list=[300]
        memb_k_list = [1]
        rg_list = [0.4]
        force_list =[0.6, 0.7, 0.75, 0.8, 0.85]
        repeat = 20
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat)
    if args.mode == 4:
        print("start with basic")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [2]
        force_ramp_rate_list=[0.5]
        temperature_list=[300]
        memb_k_list = [1]
        rg_list = [0.4]
        force_list =[0.0, 0.1, 0.2, 0.3]
        repeat = 50
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat)

    if args.mode == 3:
        print("start with basic")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0]
        force_ramp_rate_list=[1]
        temperature_list=[200, 300]
        memb_k_list = [1]
        rg_list = [0.08]
        repeat = 5
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)
    if args.mode == 2:
        print("Force ramp with high rg")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 1, 2, 3]
        force_ramp_rate_list=[1]
        temperature_list=[300]
        memb_k_list = [1, 2]
        rg_list = [0, 0.2, 0.4, 0.8, 1.6]
        repeat = 3
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)
    if args.mode == 1:
        print("Force ramp without go")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 1, 2, 3]
        force_ramp_rate_list=[1]
        temperature_list=[300]
        memb_k_list = [0, 1, 2, 4]
        rg_list = [0, 0.08, 0.2, 0.4, 1, 2]
        repeat = 3
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)
if args.day == "sep24":
    if args.mode == 1:
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 1, 2, 3]
        force_ramp_rate_list=[1]
        temperature_list=[200, 300]
        memb_k_list = [1, 2, 4]
        rg_list = [0, 0.02, 0.08, 0.2]
        repeat = 3
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)

if args.day == "sep23":
    if args.mode == 1:
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 2, 4, 6, 8]
        force_ramp_rate_list=[1]
        temperature_list=[200, 300, 400]
        memb_k_list = [1, 2, 4, 8, 16]
        repeat = 5
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)

if args.day == "sep20":
    if args.mode == 3:
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 0.1, 1, 2, 4]
        force_ramp_rate_list=[1]
        temperature_list=[300]
        memb_k_list = [1, 2, 4, 8, 16]
        repeat = 5
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)
    if args.mode == 1:
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0]
        force_ramp_rate_list=[10]
        temperature_list=["ramp"]
        memb_k_list = [1, 2, 4, 8, 16]
        repeat = 2
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)
    if args.mode == 2:
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0]
        force_ramp_rate_list=[10]
        temperature_list=["ramp"]
        repeat = 5
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)
if args.day == "sep17":
    if args.mode == 1:
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 0.1, 1, 2, 4]
        force_ramp_rate_list=[1, 10]
        temperature_list=[230, 300]
        repeat = 5
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)
if args.day == "sep11":
    if args.mode == 1:
        # folding tempearture
        pressure_list = [0.1]
        rg_list = [0.08]
        force_list =["ramp"]
        temperature_list=[230, 300]
        memb_k_list=[1]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go"]
        force_ramp_rate_list=[1]
        variable_test(temperature_list=temperature_list,
                        memb_k_list=memb_k_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=100,force_ramp_rate_list=force_ramp_rate_list)

if args.day == "sep10":
    if args.mode == 1:
        # folding tempearture
        pressure_list = [0.4, 0.8]
        rg_list = [0.4]
        force_list =[0.3, 0.35, 0.4]
        temperature_list=[230, 300]
        memb_k_list=[1, 4]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go"]
        force_ramp_rate_list=[2]
        variable_test(temperature_list=temperature_list,
                        memb_k_list=memb_k_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=10,force_ramp_rate_list=force_ramp_rate_list)
    if args.mode == 2:
        # folding tempearture
        pressure_list = [0.1]
        rg_list = [0.08, 0.2]
        memb_k_list=[1, 2, 4]
        force_list =[0.01, 0.02, 0.04, 0.08]
        temperature_list=[230, 300]
        start_from_list=["extended", "topology"]
        # start_from_list=["native"]
        simulation_model_list=["go"]
        force_ramp_rate_list=[2]
        variable_test(temperature_list=temperature_list,
                        memb_k_list=memb_k_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=4, force_ramp_rate_list=force_ramp_rate_list)
    if args.mode == 3:
        # folding tempearture
        pressure_list = [0.1]
        rg_list = [0.08]
        memb_k_list=[1]
        force_list =[0.2, 0.3, 0.35, 0.4]
        temperature_list=[230]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go"]
        force_ramp_rate_list=[2]
        variable_test(temperature_list=temperature_list,
                        memb_k_list=memb_k_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=20, force_ramp_rate_list=force_ramp_rate_list)
    if args.mode == 4:              # pending
        # folding tempearture
        pressure_list = [0.1]
        rg_list = [0.08]
        memb_k_list=[1, 2, 4]
        force_list =[0.2, 0.3, 0.35, 0.4]
        temperature_list=[230, 300]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go"]
        force_ramp_rate_list=[2]
        variable_test(temperature_list=temperature_list,
                        memb_k_list=memb_k_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=10, force_ramp_rate_list=force_ramp_rate_list)
if args.day == "sep09":
    if args.mode == 1:
        cd("simulation")
        do("pulling_prepare.py -m 7")
        cd("..")
        do("mkdir freeEnergy")
        cd("freeEnergy")
        do("make_metadata.py -m 3 -k 0.02 -t 230")
        do("pulling_analysis.py -m 6")
    if args.mode == 2:
        cmd = "gg.py -m 8"
        run_slurm = base_slurm.format(cmd)
        # folder_list = glob.glob("force_*")
        # folder_list = ['force_0.08', 'force_0.03', 'force_0.0']
        folder_list = ['force_0.055']
        # folder_list = ['force_0.07', 'force_0.02', 'force_0.045']
        # folder_list = ['force_0.06', 'force_0.04']
        print(folder_list)
        for folder in folder_list:
            cd(folder)
            cd("simulation")
            run_list = glob.glob("*")
            for run in run_list:
                cd(run)
                cd("0")
                with open("compute_angle.slurm", "w") as r:
                    r.write(run_slurm)
                do("sbatch compute_angle.slurm")
                cd("../..")
            cd("../..")
    if args.mode == 3:
        # folding tempearture
        pressure_list = [0.0, 0.1, 0.2]
        force_list =[0.0]
        rg_list = [0, 0.08, 0.1, 0.2]
        temperature_list=["ramp"]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go", "single"]
        variable_test(temperature_list=temperature_list,
                        simulation_model_list=simulation_model_list,
                        rg_list=rg_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=3,force_ramp_rate_list=[10, 1])
    if args.mode == 4:
        # folding tempearture
        pressure_list = [0, 0.1, 0.2]
        rg_list = [0, 0.08, 0.1, 0.2]
        force_list =[0.0]
        temperature_list=[230]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go", "single"]
        force_ramp_rate_list=[1, 10]
        variable_test(temperature_list=temperature_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=5,force_ramp_rate_list=force_ramp_rate_list)
    if(args.mode == 5):
        print("Extract qw and distance info.")
        for i in range(100):
            cd(str(i))
            cd("0")
            do("awk '{print $2}' wham.dat |  sed 's/,$//' | sed 1d > qw.dat")
            do("awk '{print $2}' addforce.dat |  sed 's/,$//' | sed 1d > distance.dat")
            cd("../..")
        print("create directory_list")
        with open("directory_list", "w") as f:
            for i in range(40):
                # print(os.getcwd())
                location = os.getcwd() + "/../"
                f.write(location+str(i)+"/0\n")
        do("cp ../../2xov/2xov.pdb .")
        do("python2 ~/opt/small_script/CalcLocalDistanceStats.py 2xov directory_list out")
    if args.mode == 6:
        # folding tempearture
        pressure_list = [0.1]
        rg_list = [0.1]
        memb_k_list=[1, 4]
        force_list =[0.1, 0.15, 0.05]
        temperature_list=[230, 300]
        start_from_list=["extended", "topology"]
        # start_from_list=["native"]
        simulation_model_list=["go", "single"]
        force_ramp_rate_list=[10]
        variable_test(temperature_list=temperature_list,
                        memb_k_list=memb_k_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=7,force_ramp_rate_list=force_ramp_rate_list)
    if args.mode == 7:
        # folding tempearture
        pressure_list = [0.4]
        rg_list = [0.4]
        memb_k_list=[4]
        force_list =[0.0]
        temperature_list=[230]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go", "single"]
        force_ramp_rate_list=[1]
        variable_test(temperature_list=temperature_list,
                        memb_k_list=memb_k_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=20,force_ramp_rate_list=force_ramp_rate_list)
    if args.mode == 8:
        # folding tempearture
        pressure_list = [0.4]
        rg_list = [0.4]
        memb_k_list=[4]
        force_list =[0.1, 0.15, 0.05]
        temperature_list=[230, 300]
        start_from_list=["extended", "topology"]
        # start_from_list=["native"]
        simulation_model_list=["go", "single"]
        force_ramp_rate_list=[2]
        variable_test(temperature_list=temperature_list,
                        memb_k_list=memb_k_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=6,force_ramp_rate_list=force_ramp_rate_list)
if args.day == "sep08":
    if args.mode == 1:
        pressure_list = [0.1, 0.2, 0.4]
        force_list =[0.01, 0.02, 0.04, 0.08]
        temperature_list=[230]
        start_from_list=["extended", "topology"]
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=20,force_ramp_rate_list=[1])
    if args.mode == 2:
        pressure_list = [0.1]
        force_list =[0.01, 0.02]
        temperature_list=[230, 240]
        # start_from_list=["extended", "topology"]
        start_from_list=["extended"]
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=50,force_ramp_rate_list=[0.5])
    if args.mode == 3:
        # folding tempearture
        pressure_list = [0.1, 0.2]
        force_list =[0.0]
        temperature_list=["ramp"]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go", "single"]
        variable_test(temperature_list=temperature_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=5,force_ramp_rate_list=[10])
    if args.mode == 4:
        # folding tempearture
        pressure_list = [0.1, 0.2]
        force_list =[0.0]
        temperature_list=[200]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go", "single"]
        force_ramp_rate_list=[1, 10]
        variable_test(temperature_list=temperature_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=20,force_ramp_rate_list=force_ramp_rate_list)
if args.day == "sep07":
    if args.mode == 9:
        # simulation_model_list=["go", "single"]
        force_list =[0.55, 0.5, 0.56, 0.57, 0.58, 0.59]
        variable_test(force_list=force_list, repeat=20)
    if args.mode == 8:
        # simulation_model_list=["go", "single"]
        # force_list =[0.6, 0.5, 0.4]
        # memb_k_list = [1, 2, 4]
        pressure_list = [0.1, 0.2, 0.4]
        # pressure_list = [0.8, 1.6]
        rg_list = [0.08, 0.04, 0.02]
        k_list = [0.1, 0.2, 0.3, 0.4]
        variable_test(k_list=k_list, pressure_list=pressure_list,
                        force_ramp_rate_list=[10], rg_list=rg_list, repeat=3)
    if args.mode == 7:
        cd("simulation")
        do("pulling_prepare.py -m 7")
        cd("..")
        do("mkdir freeEnergy")
        cd("freeEnergy")
        do("make_metadata.py -m 3 -k 0.02 -t 230")
        do("pulling_analysis.py -m 6")
    if args.mode == 1:
        memb_k_list = [1, 2, 4, 8, 16]
        variable_test(memb_k_list=memb_k_list, repeat=5)
    if args.mode == 2:
        force_ramp_rate_list = [1, 2, 4, 8, 16, 32]
        variable_test(force_ramp_rate_list=force_ramp_rate_list, repeat=5)
    if args.mode == 3:
        force_ramp_rate_list = [1, 2, 4, 8, 16, 32]
        variable_test(temperature_list=[200], force_ramp_rate_list=force_ramp_rate_list, repeat=5)
    if args.mode == 4:
        # simulation_model_list=["go", "single"]
        # force_list =[0.6, 0.5, 0.4]
        force_list =[0.8, 0.6, 0.4]
        memb_k_list = [1, 2, 4]
        variable_test(memb_k_list=memb_k_list, force_list=force_list, repeat=5)
    if args.mode == 6:
        # simulation_model_list=["go", "single"]
        # force_list =[0.6, 0.5, 0.4]
        force_list =[0.3]
        memb_k_list = [1, 2, 4]
        # pressure_list = [0.1, 0.2, 0.4]
        pressure_list = [0.8, 1.6]
        k_list = [10, 11, 12, 13]
        variable_test(k_list=k_list, pressure_list=pressure_list,
                        force_ramp_rate_list=[10], memb_k_list=memb_k_list, force_list=force_list, repeat=3)
if args.day == "sep06":
    if args.mode == 1:
        start_from_list=["extended", "topology"]
        simulation_model_list=["go", "single"]
        temperature_list = [200, 250, 300]
        pressure_list = [0, 0.1, 0.5, 1]
        variable_test(pressure_list=pressure_list,
                        start_from_list=start_from_list,
                        simulation_model_list=simulation_model_list,
                        repeat=5,
                        temperature_list=temperature_list,
                        commons=0)
    if args.mode == 2:
        start_from_list=["extended", "topology"]
        simulation_model_list=["go", "single"]
        temperature_list = [250, 300]
        memb_k_list = [0, 1, 2, 4]
        rg_list = [0, 0.1]
        variable_test(rg_list=rg_list, memb_k_list=memb_k_list,
                        start_from_list=start_from_list,
                        simulation_model_list=simulation_model_list,
                        repeat=3,
                        temperature_list=temperature_list,
                        commons=0)
# if args.mode == 20:
#     rg_list = [0]
#     temperature_list = [200]
#     variable_test(rg_list=rg_list, repeat=40, temperature_list=temperature_list, commons=True)


# if args.mode == 19:
#     rg_list = [0]
#     temperature_list = [175, 200, 225, 250]
#     variable_test(rg_list=rg_list, repeat=20, temperature_list=temperature_list, commons=True)

# if args.mode == 18:
#     rg_list = [0, 0.1, 0.2, 1]
#     memb_k_list = [0, 1, 2, 4]
#     pressure_list = [0, 0.1, 0.2, 0.4, 0.8, 1, 2]
#     # rg_list = [0.1]
#     # memb_k_list = [1]
#     # pressure_list = [0.1, 1]
#     variable_test(rg_list=rg_list, memb_k_list=memb_k_list, pressure_list=pressure_list, repeat=2)

# if args.mode == 17:
#     # protocol_list = ["er", "awsemer", "frag", "raptor"]
#     protocol_list = ["awsemer", "frag"]
#     protein_list = ["1occ"]
#     for protein in protein_list:
#         for protocol in protocol_list:
#             print("Work on protein: {}, protocol: {}".format(protein, protocol))
#             if protocol == "raptor":
#                 do("cp ~/opt/gremlin/protein/1occ/raptor/go_rnativeC* {}/".format(protein))
#             else:
#                 do("cp ~/opt/gremlin/protein/1occ/gremlin/go_rnativeC* {}/".format(protein))
#             do("mkdir -p {}".format(protocol))
#             do("cp -r {} {}/".format(protein, protocol))
#             cd(protocol)
#             cd(protein)
#             fileName = "{}_multi.in".format(protein)
#             if protocol == "raptor":
#                 backbone_file = "fix_backbone_coeff_er.data"
#                 do("cp ~/opt/gremlin/protein/{}/raptor/go_rnativeC* .".format(protein))
#             else:
#                 backbone_file = "fix_backbone_coeff_{}.data".format(protocol)
#                 do("cp ~/opt/gremlin/protein/{}/gremlin/go_rnativeC* .".format(protein))
#             with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
#                 for line in file:
#                     tmp = line
#                     tmp = tmp.replace("fix_backbone_coeff_er.data", backbone_file)
#                     print(tmp, end='')
#             cd("..")
#             do("run.py -m 0 -n 20 {}".format(protein))
#             cd("..")

# if args.mode == 16:
#     rg_list = [0, 0.1, 0.2, 0.4, 0.5, 1, 2, 4]
#     variable_test(rg_list=rg_list, repeat=1, commons=True)

# if(args.mode == 15):
#     print("create directory_list")
#     with open("directory_list", "w") as f:
#         for i in range(40):
#             # print(os.getcwd())
#             location = os.getcwd() + "/../"
#             f.write(location+str(i)+"/0\n")
#     do("cp ../../2xov/2xov.pdb .")
#     do("python2 ~/opt/small_script/CalcLocalDistanceStats.py 2xov directory_list out")
# if(args.mode == 14):
#     print("Extract qw and distance info.")
#     for i in range(100):
#         cd(str(i))
#         cd("0")
#         do("awk '{print $2}' wham.dat |  sed 's/,$//' | sed 1d > qw.dat")
#         do("awk '{print $2}' addforce.dat |  sed 's/,$//' | sed 1d > distance.dat")
#         cd("../..")

# if args.mode == 13:
#     rg_list = [0, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2]
#     memb_k_list = [0, 1, 2, 4, 8]
#     variable_test(rg_list=rg_list, memb_k_list=memb_k_list)
# if args.mode == 12:
#     rg_list = [0.1, 0.2, 0.4, 0.8, 1.6, 3.2]
#     variable_test(rg_list=rg_list)
# if args.mode == 11:
#     zim_type_list = ["aug04", "aug26"]
#     membrane_width_list = [30, 28.8]
#     for zim in zim_type_list:
#         for width in membrane_width_list:
#             folder = "zim_{}_width_{}".format(zim, width)
#             do("mkdir -p {}".format(folder))
#             cd(folder)
#             do("cp -r ../2xov .")
#             cd("2xov")
#             fixFile = "fix_backbone_coeff_single.data"
#             with fileinput.FileInput(fixFile, inplace=True, backup='.bak') as file:
#                 for line in file:
#                     print(line.replace("WIDTH", str(width)), end='')
#             do("cp zim_{} zim".format(zim))
#             cd("..")
#             do("run.py -n 2 2xov")
#             cd("..")
# if args.mode == 10:
#     distance_list = np.linspace(166, 180, 15)
#     for distance in distance_list:
#         folder = "dis_{}".format(distance)
#         cd(folder)
#         do("sbatch run_0.slurm")
#         cd("..")
# if args.mode == 9:
#     cmd = "python3 ~/opt/small_script/find_distance.py"
#     run_slurm = base_slurm.format(cmd)
#     folder_list = ['force_0.045']
#     print(folder_list)
#     for folder in folder_list:
#         cd(folder)
#         cd("simulation")
#         run_list = glob.glob("*")
#         for run in run_list:
#             cd(run)
#             cd("0")
#             with open("find_distance.slurm", "w") as r:
#                 r.write(run_slurm)
#             do("sbatch find_distance.slurm")
#             cd("../..")
#         cd("../..")
#
# if args.mode == 8:
#     cmd = "gg.py -m 8"
#     run_slurm = base_slurm.format(cmd)
#
#     # folder_list = glob.glob("force_*")
#     # folder_list = ['force_0.08', 'force_0.03', 'force_0.0']
#     folder_list = ['force_0.055']
#     # folder_list = ['force_0.07', 'force_0.02', 'force_0.045']
#     # folder_list = ['force_0.06', 'force_0.04']
#     print(folder_list)
#     for folder in folder_list:
#         cd(folder)
#         cd("simulation")
#         run_list = glob.glob("*")
#         for run in run_list:
#             cd(run)
#             cd("0")
#             with open("compute_angle.slurm", "w") as r:
#                 r.write(run_slurm)
#             do("sbatch compute_angle.slurm")
#             cd("../..")
#         cd("../..")

# if args.mode == 7:
#     for i in range(80):
#         do("mv {0} ../../../new_force_ramp/memb_0_force_ramp_rg_0_new/simulation/{1}".format(i,i+90))
# if args.mode == 6:
#     force_list = [0.55, 0.6, 0.65]
#     # force_list = [0.25, 0.35, 0.4, 0.45]
#     # force_list = [0.15, 0.2]
#     for force in force_list:
#         do("mkdir force_{}".format(force))
#         do("cp -r 2xov force_{}/".format(force))
#         cd("force_{}".format(force))
#         with fileinput.FileInput("2xov/2xov_multi.in", inplace=True, backup='.bak') as file:
#             for line in file:
#                 print(line.replace("MY_FORCE", str(force)), end='')
#         do("run.py -n 10 2xov/")
#         cd("..")
#
#
# if args.mode == 5:
#     # cd("start_misfolded")
#     distance_list = np.linspace(0, 30, 16)
#     for dis in distance_list:
#         do("mkdir -p dis_{}".format(dis))
#         do("cp -r ../2xov/ dis_{}".format(dis))
#         do("cp ../../freeEnergy/go_model_start_unfolded/simulation/dis_{0}/restart.25000000 dis_{0}/2xov/".format(dis))
#         cd("dis_{}".format(dis))
#         do("run.py -n 10 2xov/")
#         cd("..")

# if args.mode == 4:
#     do("rm data")
#     for i in range(100):
#         do("cat dis_{}/0/data >> data.dat".format(i))
#     do("awk '{print $1}' data.dat  > e.dat")
#     do("awk '{print $2}' data.dat  > p.dat")
#     do("awk '{print $3}' data.dat  > qw.dat")
#
# if args.mode == 1:
#     cd("simulation")
#     do("pulling_prepare.py")
#     cd("..")
#     do("mkdir freeEnergy")
#     cd("freeEnergy")
#     do("make_metadata.py -k 0.05 -t 300")
#     do("pulling_analysis.py -m 3 -p 2")
#
# if args.mode == 2:
#     print("80 bins.")
#     # cd("simulation")
#     # do("pulling_prepare.py")
#     # cd("..")
#     do("mkdir more_bin")
#     cd("more_bin")
#     do("make_metadata.py -k 0.05 -t 600")
#     do("pulling_analysis.py -m 3 -p 1")
#
# if args.mode == 3:
#     # cd("simulation")
#     # do("pulling_prepare.py")
#     # cd("..")
#     do("mkdir -p only_less_than_100")
#     cd("only_less_than_100")
#     do("make_metadata.py -k 0.05 -t 300 -m 2")
#     do("pulling_analysis.py -m 3 -p 1")
#     # for i in range(52, 70):
#     #     do("mv {}/{} .".format(i, i-40))
#     #     do("mv {} {}".format(i, i-20))
#     # for i in range(50):
#     #     do("mv force_0.8_2/{} force_0.8/{}".format(i, i+50))
#         # do("mv half_contact_force_0.8_memb1_rg1_2/{} half_contact_force_0.8_memb1_rg1/{}".format(i, i+20))
#

if(args.run):
    print("Hello World")

    name = "T0833"
    n = 21
    do("mkdir "+name)
    cd(name)
    for i in range(1, n):
        do("mkdir -p job.{}".format(i))
        do("cp ../preparation_files/myjob_nots.slurm job.{}".format(i))
        do("cp ../preparation_files/loopsubmit.bash job.{}".format(i))
        do("cp -r ../preparation_files/{0}_runpackage job.{1}/runpackage".format(name, i))
        do("cp ../preparation_files/{1}_tpr/run.{0}.tpr job.{0}/runpackage/run.tpr".format(i, name))
    for i in range(1, n):
        cd("job.{}".format(i))
        fileName = "myjob_nots.slurm"
        with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace("T0766", name), end='')

        do("bash loopsubmit.bash")
        cd("..")

    # for i in range(1, 6):
    #     do("mkdir job.{}".format(i))
    #     do("cp myjob_nots.slurm job.{}".format(i))
    #     do("cp loopsubmit.bash job.{}".format(i))
    #     do("cp -r runpackage job.{}".format(i))
    #     do("cp run.{0}.tpr job.{0}/runpackage/run.tpr".format(i))
    # for i in range(1, 6):
    #     cd("job.{}".format(i))
    #     fileName = "myjob_nots.slurm"
    #     name = "T0833"
    #     with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
    #         for line in file:
    #             print(line.replace("T0766", name), end='')
    #     do("bash loopsubmit.bash")
    #     cd("..")

    # force_list = [0.3, 0.5, 0.7, 0.9, 1.1, 1.3]
    # for force in force_list:
    #     folder = "1d_force_" + str(force)
    #     do("mkdir -p " + folder)
    #     cd(folder)
    #     do("cp ../metadatafile .")
    #     do("~/bin/python3/bin/python3 ~/opt/pulling_analysis.py -f -m 5 --force {}".format(force))
    #     do("sbatch freeEnergy.slurm")
    #     cd("..")
    # do("rm -r test")
    # do("cp -r 2xov test")
    # cd("test")
    # do("test_run.py test.in")

if(args.see):
    do("head test/0/addforce.dat")
def fix_error_run():
    n = args.number
    for i in range(n):
        os.chdir(str(i))
        os.system("cp ~/opt/pulling/qnqc.slurm .")
        os.system("sbatch qnqc.slurm")
        os.chdir("..")
    # os.system("grep 'srun: error: Application launch failed: Socket timed out on send/recv operation' . -r | cut -d':' -f1 | rev | cut -d"/" -f2- | rev > list")
    # array = []
    # cwd = os.getcwd()
    # print(cwd)
    # with open('list', 'r') as ins:
    #     for line in ins:
    #         target = line.strip('\n')
    #         array.append(target)
    # for i in array:
    #     os.chdir(i)
    #     os.system("rm slurm-*")
    #     os.system("sbatch rerun.slurm")
    #     sleep(0.5)  # Time in seconds.
    #     os.system("sbatch qnqc.slurm")
    #     sleep(0.5)  # Time in seconds.
    #     # os.system("pwd")
    #     os.chdir(cwd)

    # os.system("cut -d'/' -f2 list >")
# if(args.fix):
#     fix_error_run()


# def rerun():
#     n = args.number
#     for i in range(n):
#         os.system("cp -r {0} rerun_{0}".format(str(i)))
#         # source = "~/opt/gagb/gagb_constant200_rerun.in"
#         # target = " 2lhc.in"
#         source = "~/opt/pulling/2xov_force_load_dis.in"
#         target = " 2xov.in"
#         os.system("cp "+source + target)
#         os.system("cp "+target+" rerun_{}/".format(str(i)))
#         os.chdir("rerun_"+str(i))
#         os.system("rm slurm*")
#         os.system("sbatch run.slurm")
#         os.chdir("..")
# if(args.rerun):
#     rerun()


def continue_run():
    n = 10
    os.chdir("simulation")
    for i in range(n):
        os.system("cp -r {0} continue_{0}".format(str(i)))
        # source = "~/opt/gagb/gagb_constant200_rerun.in"
        # target = " 2lhc.in"
        source = "~/opt/pulling/2xov_continue_run.in"
        target = " 2xov.in"
        os.system("cp "+source + target)
        os.system("cp "+target+" continue_{}/".format(str(i)))
        os.chdir("continue_"+str(i))
        os.system("rm slurm*")
        os.system(  # replace RANDOM with a radnom number
            "sed -i.bak 's/RANDOM/'" +
            str(randint(1, 10**6)) +
            "'/g' 2xov.in")
        os.system("sbatch run.slurm")
        os.chdir("..")
# if(args.go):
#     continue_run()
# parser = argparse.ArgumentParser(
#         description="This is my playground for current project")
# parser.add_argument("protein", help="the name of protein")
# # parser.add_argument("template", help="the name of template file")
# args = parser.parse_args()
# # protein_name = args.template.split('_', 1)[-1].strip('/')
# protein_name = args.protein.strip('/')
    # name = "ga_2m"
## -------------Pulling--------
# os.system("cp ~/opt/small_script/springForce.plt .")
# os.system("cp ~/opt/small_script/springForce_smooth.plt .")
# os.system("gnuplot springForce.plt")
# os.system("gnuplot springForce_smooth.plt")
# os.system("open springForce.pdf")
# os.system("open springForce_smooth.pdf")
# SpringConstant_list = [3e-05, 5e-05, 1e-06, 3e-06, 5e-06, 1e-05, 1e-07]
# # SpringConstant_list = [3e-06, 5e-06]
# for SpringConstant in SpringConstant_list:
#     name = "spring"+str(SpringConstant)
#     os.system("mkdir "+name)
#     os.chdir(name)
#     os.system("cp -r ../2xov/ .")
#     os.system("cp ../variables.dat .")
#     os.chdir("2xov")
#     os.system(  # replace SIMULATION_STEPS with specific steps
#         "sed -i.bak 's/SpringForce/'" +
#         str(SpringConstant) +
#         "'/g' "+protein_name+".in")
#     os.chdir("..")
#     os.system("run.py 2xov/ -s 8 -n 2")
#     os.chdir("..")

# number_of_run_list = [2, 4, 8, 16]
# for n in number_of_run_list:
#     name = "ga_"+str(n)+"m"
#     # os.system("mkdir "+name)
#     os.system("cp -r 2lhd.pdb "+name)
#
#     # os.system("cp -r 2lhc variables.dat "+name)
#     os.chdir(name)
#     for i in range(20):
#         os.chdir("analysis/"+str(i))
#         os.system("cp ../../2lhd.pdb .")
#         os.system("python2 ~/opt/script/CalcQValue.py 2lhd.pdb dump.lammpstrj q_gb.dat")
#         os.system("python2 ~/opt/script/CalcQValue.py 2lhc.pdb dump.lammpstrj q_ga.dat")
#         os.system("cp ~/opt/small_script/qw_gagb.plt .")
#         os.system("gnuplot qw_gagb.plt")
#         os.system("mv qw_gagb.pdf ../../results/qw_gagb_{0}.pdf".format(str(i)))
#         os.chdir("../..")
#     os.chdir("..")
#
# for n in number_of_run_list:
#     name = "ga_"+str(n)+"m"
#     # os.system("mkdir "+name)
#     os.system("cp -r 2lhd.pdb "+name)
#
#     # os.system("cp -r 2lhc variables.dat "+name)
#     os.chdir(name)
#     for i in range(20):
#         os.chdir("analysis/"+str(i))
#         os.system("paste q_ga.dat q_gb.dat > q_gagb.dat")
#         os.system("cp ~/opt/small_script/qw_ga-gb.plt .")
#         os.system("gnuplot qw_ga-gb.plt")
#         os.system("mv qw_ga-gb.pdf ../../results/qw_ga-gb_{0}.pdf".format(str(i)))
#         os.chdir("../..")
#     os.system("cp ~/opt/small_script/qw_ga_all.plt .")
#     os.system("gnuplot qw_ga_all.plt")
#     os.system("cp ~/opt/small_script/qw_gb_all.plt .")
#     os.system("gnuplot qw_gb_all.plt")
#     os.system("cp ~/opt/small_script/qw_diff_all.plt .")
#     os.system("gnuplot qw_diff_all.plt")
#     os.chdir("..")

# simulation_steps = 4 * 10**6
# warm_up_steps = 10 * 10**5
#
# seed(datetime.now())
# n= 20
# vmd = "/Applications/VMD\ 1.9.2.app/Contents/MacOS/startup.command"
#
# os.system("BuildAllAtomsFromLammps.py dump.lammpstrj movie")
# os.system("cp ~/opt/plot_scripts/2xov_movie.tcl .")
# os.system(vmd+" -e 2xov_movie.tcl ")
# os.system("mkdir -p MyResults")
# for i in range(n):
#     print(i)
#     os.chdir("analysis/"+str(i))
#     os.system("cp ~/opt/plot_scripts/2xov_movie_screenshot.tcl .")
#     os.system(vmd+" -e 2xov_movie_screenshot.tcl")
#     os.system("cp frame1000.tga ../../MyResults/frame"+str(i)+"_1000.tga")
#     #os.system("cp frame450.tga ../Results/frame"+folder_name+"_450.tga")
#     # os.system("movie.py "+protein_name)
#     os.chdir("../..")
#     # analysis

# folder_name = ""
# result_folder = "WeiLu_Aug_07"

# protein_list = ['T089', 'T120', 'T251', 'TOP7', '1UBQ']
# sublist = ['']
# # sublist = ['_ha', '_he']
# # sublist = ['_lp', '_he_lp']
# # folder_list = []
# for protein in protein_list:
#     for sub in sublist:
#         folder_name = protein+sub
#         os.chdir(folder_name)
#         os.chdir("best_2nd")
#         os.system("pymol ~/opt/plot_scripts/align.pml > matrix.dat")
#         os.system("head -n 70 matrix.dat | tail -n 20 > cealign_matrix.dat")
#         # for i in range(19, -1, -1):
#         #     os.system("mv {}.pdb {}.pdb".format(i, i+1))
#         os.chdir("../..")
    # os.chdir(protein)
    # os.chdir("best_1st")
    # os.system("python3 ~/opt/small_script/cross_q.py")
    # os.chdir("..")
    # os.chdir("best_2nd")
    # os.system("python3 ~/opt/small_script/cross_q.py")
    # os.chdir("..")
    # os.chdir("..")
# n = 3
# for i in range(n):
#     # simulation set up
#     folder_name = str(i)
#     os.system("mkdir -p "+folder_name)
#     os.system("cp -r "+args.protein+"* "+folder_name)
#     os.chdir(folder_name)
#     os.system("cp ../../helix_less/simulation/"+str(i)+"/restart.4000000 .")
#     os.system(  # replace SIMULATION_STEPS with specific steps
#         "sed -i.bak 's/WARM_UP_STEPS/'" +
#         str(warm_up_steps) +
#         "'/g' "+protein_name+".in")
#     os.system(  # replace RANDOM with a radnom number
#             "sed -i.bak 's/RANDOM/'" +
#             str(randint(1, 10**6)) +
#             "'/g' "+protein_name+".in")
#     os.system(  # replace SIMULATION_STEPS with specific steps
#             "sed -i.bak 's/SIMULATION_STEPS/'" +
#             str(simulation_steps) +
#             "'/g' "+protein_name+".in")
# # if(platform.system() == 'Darwin'):
# #     os.system("/Users/weilu/Documents/lammps-9Oct12_modified/src/lmp_serial \
# #     < "+protein_name+".in")
#     if(platform.system() == 'Darwin'):
#         os.system("/Users/weilu/Documents/lammps-9Oct12_modified/src/lmp_serial \
#         < "+protein_name+".in")
#     elif(platform.system() == 'Linux'):
#         os.system("cp ~/opt/run.slurm .")
#         os.system(  # replace PROTEIN with pdb name
#                 "sed -i.bak 's/PROTEIN/'" +
#                 protein_name +
#                 "'/g' run.slurm")
#         os.system("sbatch run.slurm")
#     else:
#         print("system unkown")
#     os.chdir("..")
# exit(1)

# w_helix_list = [0.1, 0.5, 1, 1.5]
# m_helix_list = [0.1, 0.5, 1, 1.5]
#
# for i in range(len(w_helix_list)):
#     w = w_helix_list[i]
#     for j in range(len(m_helix_list)):
#
#         # m = m_helix_list[j]
#         folder_name = str(i)+"_"+str(j)
#         # os.system("cd "folder_name)
#         os.chdir(folder_name)
#         # os.system("analysis.py 2xov/")
#         # os.system("echo "+folder_name+" >> ../all")
#         os.system("sort -k 3 analysis/list_of_max_q > ../data/"+folder_name)
#         os.chdir("..")
#         # os.system("mkdir "+folder_name)
#         # os.chdir(folder_name)
#         # os.system("cp -r ../2xov .")
#         # os.chdir("2xov")
#         # os.system(
#         #         "sed -i.bak 's/W_HELIX/'" +
#         #         str(w) +
#         #         "'/g' fix_backbone_coeff.data")
#         # os.system(
#         #         "sed -i.bak 's/M_HELIX/'" +
#         #         str(m) +
#         #         "'/g' fix_backbone_coeff.data")
#         # os.chdir("..")
#         # os.system("run.py 2xov/ -n 5")

# os.system("cp ~/opt/gg.py this_gg.py")
# for i in range(5):
#     os.system("mkdir "+str(i))
#     os.chdir(str(i))
#     os.system("cp -r ../2xov/ .")
#     os.system("cp ../../2xov_strong_single_memory_600to500/simulation/"+str(i)+"/restart.2000000 2xov/")
#     os.system("run.py -s 4 -n 2 2xov/")
#     os.chdir("..")

# # rama_list = [6, 8, 16]
# # rama_list = [4]
# melt_t_list = [400, 500, 600]
# for variable in melt_t_list:
#     folder_name = str(variable)
#     os.system("mkdir "+folder_name)
#     os.chdir(folder_name)
#     os.system("cp -r ../1qjp .")
#     os.chdir("1qjp")
#     os.system(
#             "sed -i.bak 's/MELTT/'" +
#             str(variable) +
#             "'/g' 1qjp.in")
#     os.chdir("..")
#     # os.system("pwd")
#     os.system("run.py 1qjp/ -n 5 -s 5")
#     os.chdir("..")
# os.system("cp ~/opt/gg.py this_gg.py")
#
# exec (open("config.py").read())
# n = number_of_run
# steps = simulation_steps
#
# protein_name = args.protein.strip('/')
#
# temp = 400
# folder_name = "{}_t{}_q100_test11".format(protein_name, str(temp))
# print("all going to "+folder_name)
# os.system("mkdir -p "+folder_name)
# os.system("rm -f "+folder_name + "/*")
# command = 'cat simulation/{}/%d/wham11 \
# >> {}/all_wham.dat'.format(temp, folder_name)
# # cal rmsd
# os.chdir("simulation/"+str(temp))
# for i in range(n):
#     os.chdir(str(i))
#     os.system("awk '{print>\"file1\"(NR>(n/2)?2:1)}' n=\"$(wc -l <file1)\" file1")
#     os.system("cat file11 >> ../../../"+folder_name+"/rmsd_total")
#     # os.system("sed 1d wham.dat > wham1d.dat")
#     os.system("awk '{print>\"wham1\"(NR>(n/2)?2:1)}' n=\"$(wc -l <wham1)\" wham1")
#     os.chdir("..")
# os.chdir("../..")
# for i in range(n):
#     cmd = command % i
#     os.system(cmd)
# os.chdir(folder_name)
# os.system("awk '{print $2}' all_wham.dat > Qw_total")
# os.system("awk '{print $3}' all_wham.dat > rg_total")
# os.system("awk '{print $4}' all_wham.dat > p_total")
# os.system("awk '{print $5}' all_wham.dat > tc_total")
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# os.system("cp ~/opt/wham_analysis/*.m .")
# os.chdir("..")
# os.system("~/opt/script/wham/fused_calc_cv.sc {} top7 50 400 350 450 5 50 100 0 0.98".format(folder_name))
#
#
# folder_name = "{}_t{}_q100_test12".format(protein_name, str(temp))
# print("all going to "+folder_name)
# os.system("mkdir -p "+folder_name)
# os.system("rm -f "+folder_name + "/*")
# command = 'cat simulation/{}/%d/wham12 \
# >> {}/all_wham.dat'.format(temp, folder_name)
# # cal rmsd
# os.chdir("simulation/"+str(temp))
# for i in range(n):
#     os.chdir(str(i))
#     os.system("cat file12 >> ../../../"+folder_name+"/rmsd_total")
#     os.chdir("..")
# os.chdir("../..")
# for i in range(n):
#     cmd = command % i
#     os.system(cmd)
# os.chdir(folder_name)
# os.system("awk '{print $2}' all_wham.dat > Qw_total")
# os.system("awk '{print $3}' all_wham.dat > rg_total")
# os.system("awk '{print $4}' all_wham.dat > p_total")
# os.system("awk '{print $5}' all_wham.dat > tc_total")
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# os.system("cp ~/opt/wham_analysis/*.m .")
# os.chdir("..")
#
#
#
# os.system("~/opt/script/wham/fused_calc_cv.sc {} top7 50 400 350 450 5 50 100 0 0.98".format(folder_name))

#
# result_folder = "WeiLu_Aug_07"
# os.system("mkdir -p "+result_folder)
# protein_list = ['T089', 'T120', 'T251', 'top7', '1UBQ']
# # sublist = ['_ha', '_he']
# sublist = ['_lp', '_he_lp']
# folder_list = []
# for protein in protein_list:
#     for sub in sublist:
#         folder_list += [protein+sub]
# print(folder_list)
# # exit(1)
# # awk '{print>'file'(NR>(n/2)?2:1)}' n='$(wc -l <test)' test
# for folder in folder_list:
#     print(folder)
#     os.chdir(folder)
#     exec (open("config.py").read())
#     n = number_of_run
#     steps = simulation_steps
#     os.system("mkdir -p ../{}/".format(result_folder)+folder+"/best_q")
#     os.system("sort analysis/list_of_max_q > ../{}/q_".format(result_folder)+folder+".dat")
#     for i in range(n):
#         # move
#         os.chdir("analysis/"+str(i))
#         os.system("cp chosen.pdb ../../../{}/".format(result_folder) + folder+"/best_q/"+str(i)+".pdb")
#         os.chdir("../..")
#     os.chdir("..")


# result_folder = "WeiLu_Aug_07"
# os.system("mkdir -p "+result_folder)
# protein_list = ['T089', 'T120', 'T251', 'top7', '1UBQ']
# # sublist = ['_ha', '_he']
# sublist = ['_lp', '_he_lp']
# folder_list = []
# for protein in protein_list:
#     for sub in sublist:
#         folder_list += [protein+sub]
# print(folder_list)
# # exit(1)
#
# for folder in folder_list:
#     print(folder)
#     os.chdir(folder)
#     exec (open("config.py").read())
#     n = number_of_run
#     steps = simulation_steps
#     os.system("mkdir -p ../{}/".format(result_folder)+folder+"/best_q")
#     os.system("sort analysis/list_of_max_q > ../{}/q_".format(result_folder)+folder+".dat")
#     for i in range(n):
#         # move
#         os.chdir("analysis/"+str(i))
#         os.system("cp chosen.pdb ../../../{}/".format(result_folder) + folder+"/best_q/"+str(i)+".pdb")
#         os.chdir("../..")
#     os.chdir("..")



'''