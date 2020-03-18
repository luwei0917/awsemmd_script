#!/usr/bin/env python3
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
# import imp
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
from collections import defaultdict

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

# if args.test:
#     do = print
# else:
#     do = os.system
with open('cmd_gg_server.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')


def do(cmd, get=False, show=True):
    if get:
        out = subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True).communicate()[0].decode()
        if show:
            print(out, end="")
        return out
    else:
        return subprocess.Popen(cmd, shell=True).wait()
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
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH -o outs/slurm-%j.out
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
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=04:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH -o outs/slurm-%j.out
echo "My job ran on:"
echo $SLURM_NODELIST
srun {}\n'''

gpu_base_slurm = '''\
#!/bin/bash

#SBATCH --account=commons
#SBATCH --partition=scavenge
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:volta:1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --export=ALL
#SBATCH -o outs/slurm-%j.out
#SBATCH --mail-user=luwei0917@gmail.com

module load GCC/8.3.0 CUDA/10.1.168

srun {}
'''

gpu_commons_slurm = '''\
#!/bin/bash

#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:volta:1
#SBATCH --time=23:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --export=ALL
#SBATCH -o outs/slurm-%j.out
#SBATCH --mail-user=luwei0917@gmail.com

module load GCC/8.3.0 CUDA/10.1.168

srun {}
'''

# def replace(TARGET, FROM, TO):
#     do("sed -i.bak 's/{}/{}/g' {}".format(FROM,TO,TARGET))



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
quick_template_large_mem_slurm = '''\
#!/bin/bash
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
                        # extract protein id
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

def generate_multiShuffle(fullName, alignmentLocation="./", location="./", num_decoys=1000, nameMode=0):
    if nameMode == 0:
        with open(alignmentLocation+f"alignments/{fullName}_filtered_0.05.seqs") as f:
            a = f.readlines()
    elif nameMode == 1:
        with open(alignmentLocation+f"alignments/{fullName}.seqs") as f:
            a = f.readlines()

    # with open(location+f"../database/S20_seq/{fullName}.seq") as f:
    #     b = f.readlines()

    size = len(a)
    if size < num_decoys:
        a = a * (num_decoys//size+1)
        print(fullName, len(a), size)
    with open(location+f"decoys/multiShuffle/{fullName}.decoys", "w") as out:
        for seq in random.sample(a, num_decoys):
            s = seq.strip()
            shuffled_seq = ''.join(random.sample(s,len(s)))
            out.write(shuffled_seq+"\n")
        # print(shuffled_seq)

def waitForJobs(jobIdList, sleepInterval=30):
    from datetime import datetime as dt
    if len(jobIdList) == 0:
        return
    previousJobNotFinished = True
    while previousJobNotFinished:
        print(f"Waiting for previous jobs {jobIdList}", dt.now())
        time.sleep(sleepInterval)
        previousJobNotFinished = False
        a = getFromTerminal("squeue -u wl45")
        for jobId in jobIdList:
            if jobId in a:
                previousJobNotFinished = True
    print("Continue Next Script")

# dataset = {"old":("1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "), 40),
#             "new":("1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "), 80),
#             "test":(["t089", "t120", "t251", "top7", "1ubq", "t0766", "t0778", "t0782", "t0792", "t0803", "t0815", "t0833", "t0842", "t0844"], 40)}

# dataset = {"old":("1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR".split(", "), 40),
#             "new":("1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "), 80),
#             "test":(["t089", "t120", "t251", "top7", "1ubq", "t0766", "t0778", "t0782", "t0792", "t0803", "t0815", "t0833", "t0842", "t0844"], 40),
#             }

dataset = {"old":"1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "),
            "new":"1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "),
            "test":["t089", "t120", "t251", "top7", "1ubq", "t0766", "t0778", "t0782", "t0792", "t0803", "t0815", "t0833", "t0842", "t0844"]}
dataset["combined"] = dataset["old"] + dataset["new"]

dataset["may13"] = ['1r69', '3icb', '256b', '4cpv', '2mhr', '1mba', '2fha', '1fc2', '1enh', '2gb1', '2cro', '1ctf', '4icb']
dataset["membrane"] = ["2bg9", "1j4n", "1py6_SD", "2bl2", "1rhz", "1iwg", "2ic8", "1pv6", "1occ", "1kpl", "2bs2", "1py6", "1u19"]
dataset["hybrid"] = ["2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91", "4nv6", "4p79", "5dsg", "6g7o", "6a93", "2jo1", "1py6", "1pv6", "1u19"]
dataset["optimization"] = ['1e0m', '1w4e', '1e0g', '2wqg', '1jo8', '1fex', '2l6r', '1c8c', '1g6p', '1mjc', '2jmc', '1hdn', '1st7', '1n88', '1d6o', '1hcd', '2ga5', '1j5u', '3o4d', '1k0s']
dataset["optimization_cath"] = ['1a75A00', '1bekA01', '1bqbA02', '1cpcB00', '1cscA02', '1cy5A00', '1dv5A00', '1e8yA05', '1evyA02', '1in4A03', '1l1fA03', '1vq8P01', '1xmkA00', '1zcaA02', '2grhA00', '2ii2A04', '2q6fB03', '2wh6A00', '3g0vA00', '3geuA00', '3h99A03', '3hrdD02', '3ju5A01', '3p1wA03', '4cxfA01', '4i2aA01', '4i4tB03', '4i6uB00', '5kn9A02']
# '1hcd', '1k0s' have problem in beta sheets. (the crystal structure is more round while prediction is more like a straight sheet)
dataset["optimization_v2"] = ['1e0m', '1w4e', '1e0g', '2wqg', '1jo8', '1fex', '2l6r', '1c8c', '1g6p', '1mjc', '2jmc', '1hdn', '1st7', '1n88', '1d6o', '2ga5', '1j5u', '3o4d']
dataset["optimization_cbd"] = ['1hoe', '1hyp', '1tif', '1vcc', '1by9', '1bdo', '451c', '1cc5', '1bb9', '1pht', '1opd', '1a32', '1ptf', '1cyo', '1tig', '1ctj', '1fna', '1rzl', '1who', '2cbp', '2acy', '1plc', '1bm8', '1opc', '3vub', '1tul', '1kte', '1erv', '1btn', '1a1x', '3cyr', '1bkf', '1ycc', '1sfp', '1kpf', '2mcm', '2pii', '1a6f', '1by2', '1bea', '1rmd', '1poa', '1tmy', '2a0b', '1mai', '1neu', '1dun', '1b6e', '2sak', '1dhn', '1cxc', '1bgf', '7rsa', '1bqk', '3pyp', '1bfg', '1opy', '1rlw', '1rie', '3chy', '1rcb', '1cpq', '1pdo', '3lzt', '1hmt', '1htp', '1c52', '1kuh', '1crb', '1poc', '1aqt', '2end', '5nul', '1pne', '1lcl', '2sns', '1flp', '1tfe', '1ax8', '1pkp', '1rss', '1jon', '1vls', '1lba', '1aly', '1mba', '2hbg', '1akr', '1osa', '1div']
# def returnSteps(p):
#     if p in "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "):
#         steps = 80
#     elif p in "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR".split(", "):
#         steps = 40
#     elif p in ["1MBA", "2FHA"]:
#         steps = 30
#     return steps


def get_aligned_info(p1, p2):
    cmd = f"/home/wl45/build/TMalign/TMalign {p1} {p2} | grep -A 1 'Aligned'"
    output = getFromTerminal(cmd)
    # print(output)
    line = output.split("\n")[0]
    line2 = output.split("\n")[1]
    aligned_length,rmsd,seqid = line.split(",")
    aligned_length = int(aligned_length.split("=")[1])
    rmsd = float(rmsd.split("=")[1])
    tmscore = float(line2.split(",")[0].split("=")[-1].strip().split(" ")[0])
    seqid = float(seqid.split("=")[-1])
    # print("aligned_length, rmsd, tmscore, seqid")
    # print(aligned_length, rmsd, tmscore, seqid)
    return aligned_length, rmsd, tmscore, seqid

def readList(fileName):
    with open(fileName) as f:
        a = f.readlines()
    theList = [b.strip() for b in a]
    return theList

def slurmRun(slurmFileName, cmd, template=scavenge_slurm, memory=1, thread=1):
    os.system("mkdir -p outs")
    with open(slurmFileName, "w") as out:
        out.write(template.format(cmd))
        # out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 1 proteins_name_list/proteins_name_list_{name}.txt"))
    replace(slurmFileName, "#SBATCH --mem-per-cpu=1G", f"#SBATCH --mem-per-cpu={memory}G")
    replace(slurmFileName, "#SBATCH --cpus-per-task=1", f"#SBATCH --cpus-per-task={thread}")

    a = getFromTerminal(f"sbatch {slurmFileName}")
    jobId = a.split(" ")[-1].strip()
    return jobId

if args.day == "tmpred":
    pdb_list = ["4nv6", "4p79", "5dsg", "6g7o", "6a93", "5xpd", "3kp9", "4a2n", "5d91", "2jo1"]
    pdb_list += ["1py6", "1pv6", "1u19"]
    pdb_list += ["2xov_complete", "6e67A"]
    if args.mode == 1:
        do("mkdir -p TM_pred")
        cd("TM_pred")
        # pdb = "4rws"
        # pdb = args.label
        # do(f"/projects/pw8/wl45/topology_prediction/PureseqTM_Package/PureseqTM.sh -i ../setup/{pdb}/{pdb}.fasta")
        for pdb in pdb_list:
            do(f"/projects/pw8/wl45/topology_prediction/PureseqTM_Package/PureseqTM_proteome.sh -i ../setup/{pdb}/{pdb}.fasta")
        cd("..")
    if args.mode == 2:
        for pdb in pdb_list:
            PredictedZimAndForceSetupFile(pdb, forceLocation="forces")
    if args.mode == 3:
        for pdb in pdb_list:
            print(pdb)
            do(f"cat setup/{pdb}/{pdb}.fasta")

if args.day == "dmp":
    dmp_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=commons
#SBATCH --partition=scavenge
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --time=04:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH -o outs/slurm-%j.out
echo "My job ran on:"
echo $SLURM_NODELIST
srun {}\n'''

    if args.label[-6:] == ".fasta":
        pdb = args.label[:-6]
    else:
        pdb = args.label
    # folder = pdb
    # folder = "DMP"
    # do(f"mkdir -p {folder}")
    # os.chdir(folder)
    # os.system(f"cp ../{pdb}.fasta .")
    cmd = f"bash /projects/pw8/wl45/DeepMetaPSICOV/run_DMP.sh -i {pdb}.fasta"
    do("mkdir -p outs")
    out = slurmRun(f"{pdb}.slurm", cmd, template=dmp_slurm)
    print(out)
    # os.chdir("..")

if args.day == "plot":
    import matplotlib.pyplot as plt
    cutoff_list = [100, 200, 300, 400, 500, 600]
    # cutoff_list += [10, 20, 30, 40, 50, 80]
    save_gamma_pre = "saved_gammas"
    # trial_name = "iter1"
    trial_name = args.label
    os.system(f"mkdir -p {save_gamma_pre}/figures")
    for cutoff_i in cutoff_list:
        # cutoff_i = 400
        name = f"{save_gamma_pre}/{trial_name}_cutoff{cutoff_i}_impose_Aprime_constraint"
        filtered_gamma = np.loadtxt(name)
        figureName = f"{save_gamma_pre}/figures/{trial_name}_cutoff{cutoff_i}_equal_legend"
        title = f"{trial_name}_cutoff{cutoff_i}"
        show_together_v2(filtered_gamma, figureName, title=title, inferBound=2)

def read_simulation_info(folder_list, pdb_list, simulationType, base_path, run_n):
    all_data = []
    for folder in folder_list:
        for pdb in pdb_list:
            for i in range(run_n):
                    pre = f"{base_path}/{simulationType}/{folder}/{pdb}/{i}"
                    info_file = "info.dat"
                    location = f"{pre}/{info_file}"
                    try:
                        tmp = pd.read_csv(location, sep="\s+")
                        tmp = tmp.assign(Run=i, Protein=pdb, Folder=folder)
                        all_data.append(tmp)
                    except:
                        print(pdb, i, folder)
                        pass
    data = pd.concat(all_data)
    today = datetime.today().strftime('%m-%d')
    outFile = f"{simulationType}.csv"
    data.reset_index(drop=True).to_csv(outFile)
    print(outFile)
    return outFile

if args.day == "mar16":
    # 3dRobot decoys.
    # mass specific decoys iteration.
    a = glob.glob("/scratch/wl45/mar_2020/3DROBOT_DATA/3DRobot_set/*")
    pdb_list = []
    for pdb in a:
        if pdb[-4:] == ".bz2":
            continue
        name = pdb.split("/")[-1]
        pdb_list.append(name)
    print(pdb_list)
    print(len(pdb_list))

    folder_list = ["3DRobot_set"]
    # folder_list = ["iteration_1"]
    # folder_list = ["iteration_2"]
    run_n = 2
    base_path = "/scratch/wl45/mar_2020"
    simulationType = "mass_iterative_run"
    lastN_frame = 50
    if args.mode == 1:
        do("mkdir -p database/dompdb")
        do("mkdir -p database/S20_seq")
        for pdb in pdb_list:
            do(f"cp /scratch/wl45/mar_2020/3DROBOT_DATA/3DRobot_set/{pdb}/native.pdb database/dompdb/{pdb}.pdb")
        def getSeq(fileLocation):
            p = PDBParser()
            s = p.get_structure("test", fileLocation)
            seq = ""
            residues = list(s.get_residues())
            for residue in residues:
                res_id = residue.get_id()[0]
                if res_id==' ':
                    residue_name = residue.get_resname()
                    seq += three_to_one(residue_name)
            return seq
        for pdb in pdb_list:
            fileLocation = f"database/dompdb/{pdb}.pdb"
            seq = getSeq(fileLocation)
            opt_pdbName = pdb
            fileLocation = f"database/S20_seq/{opt_pdbName}.seq"
            with open(fileLocation, "w") as out:
                out.write(seq+"\n")

    if args.mode == 2:
        # do("cp /home/wl45/opt/optimization/phi_list_contact.txt phi_list.txt")
        do("cp ~/opt/optimization/phi_com.txt phi_list.txt")
        do("mkdir -p proteins_name_list")
        do("mkdir -p slurms")
        do('mkdir -p outs')
        # do("mkdir -p decoys/multiShuffle")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p ../phis")
        # with open("protein_list", "w") as out:
        #     for pdb in pdb_list:
        #         out.write(f"{pdb}\n")
        with open("protein_list", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
    if args.mode == 321:
        pdb = "3ZOOA"
        pre = f"{base_path}/3DROBOT_DATA/3DRobot_set/{pdb}/"
        nativePdb = f"{pre}/native.pdb"
        decoyPdb = f"{pre}/decoy73_1.pdb"
        compute_Q_from_two_pdb(nativePdb, decoyPdb)
    if args.mode == 3:
        to_folder = "."
        for folder in folder_list:
            os.system(f"mkdir -p {to_folder}/decoys/3DRobot")
            for i, pdb in enumerate(pdb_list):
                print(pdb, i)
                pre = f"{base_path}/3DROBOT_DATA/3DRobot_set/{pdb}/"
                nativePdb = f"{pre}/native.pdb"
                decoy_list = glob.glob(f"{pre}/decoy*.pdb")
                info_ = []
                p = PDBParser(QUIET=True)
                for decoyPdb in decoy_list:
                    decoyName = decoyPdb.split("/")[-1].split(".")[0]
                    # decoyPdb = decoy_list[0]
                    s = p.get_structure("decoy", decoyPdb)
                    decoy_s = s[0]
                    # q = do(f"python ~/opt/small_script/CalcQValueFromTwoPdb_2.py {nativePdb} {decoyPdb}", get=True)
                    q = compute_Q_from_two_pdb(nativePdb, decoyPdb)
                    try:
                        q = float(q)
                    except:
                        print("??", q)
                        print(decoyPdb)
                        continue
                    info_.append([decoyName, q, decoy_s])
                data = pd.DataFrame(info_, columns=["Protein", "Qw", "structure"]).sort_values("Qw").reset_index(drop=True)
                data.to_pickle(f"{to_folder}/decoys/3DRobot/{pdb}_{folder}.pkl")
    if args.mode == 33:
        # time.sleep(3600)
        dataFile = f"{simulationType}.csv"
        simulation_folder = simulationType
        data = pd.read_csv(dataFile, index_col=0)
        parser = PDBParser()
        # folder = "first"
        # for folder in ["first", "first_cpu2"]:
        for folder in folder_list:
            pre = f"{base_path}/{simulation_folder}/{folder}"
            to_folder = "."
            os.system(f"mkdir -p {to_folder}/decoys/openMM")
            # pdb_list = ['1j5u']
            for pdb in pdb_list:
                complete_models = []
                print(pdb)
                for i in range(run_n):
                    movieFile = f"{pre}/{pdb}/{i}/movie.pdb"
                    allFrames, n, size = getAllFrames(movieFile)
                    num_of_frames = int(n/size)
                    first_chosen_frame = num_of_frames - lastN_frame
                    # last_chosen_frame = num_of_frames
                    oneFrame = allFrames[size*first_chosen_frame:size*(num_of_frames)]

                    p = PDBParser()
                    f = io.StringIO("".join(oneFrame))
                    s = p.get_structure("test", f)
                    complete_models += list(s.get_models())
                t = data.query(f"Protein == '{pdb}' and Folder == '{folder}' and Steps > 1").reset_index(drop=True)
                t = t.groupby("Run").tail(50).reset_index(drop=True)
                t["structure"] = complete_models
                t = t.rename(columns={"Q":"Qw"})
                # last50 = t.groupby("Run").tail(50).reset_index(drop=True)
                last50 = t
                # print(last50.head())
                # print(last50.tail())
                # print(last50.shape)
                # print(len(complete_models))
                to_folder = "."
                last50.to_pickle(f"{to_folder}/decoys/openMM/{pdb}_{folder}.pkl")

    if args.mode == 23:
        # split proteins_name_list
        n = split_proteins_name_list(pdb_per_txt=2)
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d mar05 -m 22 -l {proteins}", template=base_slurm, memory=10)
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
    if args.mode == 22:
        protein = args.label
        do(f"python3 ~/opt/compute_phis.py -m 7 {protein}")
        # do("python3 ~/opt/compute_phis.py -m 0 test_protein")
    if args.mode == 23:
        # split proteins_name_list
        n = split_proteins_name_list(pdb_per_txt=2)
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d mar05 -m 22 -l {proteins}", template=base_slurm, memory=10)
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
        waitForJobs(jobIdList, sleepInterval=60)
        # folder_list = ["iteration_0_cbd", "iteration_1_cbd", "iteration_2_cbd"]
        # folder_list = ["iteration_start_native"]
        # folder_list = ["iteration_start_native", "iteration_start_native_iter2"]
        folder_list = ["iteration_start_native", "iteration_start_native_iter2",        "iteration_start_native_iter3", "iteration_new_1", "iteration_start_native_iter4"]
        folder_list += ["iteration_start_native_iter5", "iteration_new_2", "iteration_new_3"]
        # folder_list += ["iteration_new_3"]
        with open("protein_list_complete", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d mar05 -m 33"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")

    if args.mode == 33:
        # folder_list = ["iteration_0", "iteration_1", "iteration_2"]
        # folder_list = ["iteration_0_cbd", "iteration_1_cbd", "iteration_2_cbd"]
        # folder_list = ["iteration_start_native"]
        # folder_list = ["iteration_start_native", "iteration_start_native_iter2"]
        # folder_list = ["iteration_start_native", "iteration_start_native_iter2",        "iteration_start_native_iter3", "iteration_new_1", "iteration_start_native_iter4"]
        # folder_list += ["iteration_start_native_iter5", "iteration_new_2", "iteration_new_3"]
        folder_list = ["iteration_0_cbd", "iteration_1_cbd", "iteration_2_cbd"]
        folder_list += ["iteration_new_1", "iteration_new_2", "iteration_new_3"]
        with open("protein_list_complete", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d mar05 -m 33"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 33:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list_complete"
        n_decoys = 50 * run_n
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='openMM',
                                        withBiased=True, oneMinus=True, decoyBiasName='decoysQ', num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, )
        # next: optimization_analysze.py iter{i}
    if args.mode == 4:
        # folder_list = ["iteration_0", "iteration_1", "iteration_2"]
        # folder_list = ["iteration_0_cbd", "iteration_1_cbd", "iteration_2_cbd"]
        # with open("protein_list_complete", "w") as out:
        #     for folder in folder_list:
        #         for pdb in pdb_list:
        #             out.write(f"{pdb}_{folder}\n")
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d mar05 -m 44"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 44:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list_complete"
        n_decoys = 50 * run_n
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='openMM',
                                        withBiased=True, oneMinus=False, decoyBiasName='deocyQZBias', num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, )
        # next: optimization_analysze.py iter{i}
    if args.mode == 5:
        # i = 2
        # i = 2
        # i = 3
        # i = 0
        # i = 1
        # i = 2
        i = "native_iter2"
        i = "native_iter3"
        i = "native_iter4"
        i = "new_3"
        i = "new_4"
        fromFile = "/scratch/wl45/mar_2020/cath_dataset_shuffle_optimization/optimization_iter0_com_density_rmin_2p5/saved_gammas/iter0_cutoff600_impose_Aprime_constraint"

        do(f"optimization_analyze.py {i} --proteinList protein_list_complete --gammaFile {fromFile}")

    if args.mode == 6:
        # i = 0
        # i = 1
        # alpha = 0.1
        # i = 2
        i = "native_iter2"
        i = "native_iter3"
        i = "native_iter4"
        alpha = 0.9
        cutoff = 400
        # i = 2
        fromFile = "/scratch/wl45/mar_2020/cath_dataset_shuffle_optimization/optimization_iter0_com_density_rmin_2p5/saved_gammas/iter0_cutoff600_impose_Aprime_constraint"
        # ../cath_dataset_shuffle_optimization/optimization_iter0_com_density_rmin_2p5/Oct26_saved_gammas/iter0_cutoff600_impose_Aprime_constraint/gamma.dat
        # do("mix_gamma.py ~/opt/parameters/original_gamma saved_gammas/iter0_cutoff100_impose_Aprime_constraint -i 0")
        # do(f"mix_gamma.py ~/opt/parameters/original_gamma saved_gammas/iter{i}_cutoff300_impose_Aprime_constraint -i {i}")
        # do(f"mix_gamma.py ~/opt/parameters/original_gamma saved_gammas/iter{i}_cutoff400_impose_Aprime_constraint -i {i}")
        do(f"mix_gamma.py {fromFile} saved_gammas/{i}_cutoff{cutoff}_impose_Aprime_constraint -i {i} -a {alpha}")
    if args.mode == 66:
        import matplotlib.pyplot as plt
        # i = 2
        save_gamma_pre = "saved_gammas"
        trial_name = f"iter{i}"
        # trial_name = args.label
        # os.system(f"mkdir -p {save_gamma_pre}/figures")

        cutoff_i = 400
        name = f"iter_{i}_30"
        filtered_gamma = np.loadtxt(name)
        figureName = f"{trial_name}_cutoff{cutoff_i}_equal_legend"
        title = f"{trial_name}_cutoff{cutoff_i}"
        show_together_v2(filtered_gamma, figureName, title=title, inferBound=2)


if args.day == "mar17":
    if args.mode == 1:
        pdb_list = ["1r69"]
        for pdb in pdb_list:
            pre = f"./"
            # pdb = "1baj"
            original_openAWSEM_input = f"{pre}/{pdb}-openmmawsem.pdb"
            new_openAWSEM_input = f"{pre}/cbd-openmmawsem.pdb"
            all_atom_pdb_file = f"{pre}/crystal_structure-cleaned.pdb"
            replace_CB_coord_with_CBD_for_openAWSEM_input(original_openAWSEM_input, new_openAWSEM_input, all_atom_pdb_file)
    if args.mode == 2:
        pdb = "1r69"
        pre = f"."
        fromFile = f"{pre}/crystal_structure.pdb"
        toFile = f"{pre}/cbd_{pdb}.pdb"
        convert_all_atom_pdb_to_cbd_representation(fromFile, toFile)
        cmd = f"python /projects/pw8/wl45/openawsem/helperFunctions/Pdb2Gro.py {pre}/cbd_{pdb}.pdb {pre}/cbd_{pdb}.gro"
        do(cmd)
        do(f"cp {pre}/single_frags.mem {pre}/cbd_single_frags.mem")
        replace(f"{pre}/cbd_single_frags.mem", pdb, f"cbd_{pdb}")
if args.day == "mar15":
    # generate MSA
    # mass specific decoys iteration.
    # should remove 1puc, 1skz
    # pdb_list = dataset["optimization_cbd"]
    print("mar15")
    data = pd.read_csv("training_set.csv")
    specific_decoys = data.query("Length < 150 and Length > 70").reset_index(drop=True)
    pdb_list = specific_decoys["Protein"].to_list()
    pdb_list = [a.lower() for a in pdb_list]
    skip_pdb_list = ["1puc", "1skz"]
    skip_pdb_list += ["1msc", "1fmb", "1gvp", "2tgi", "1whi", "1baj", "1rmd", "1div"]  #dimer.
    skip_pdb_list += ["1aqe"]  # lots of ligand
    filtered_pdb_list = [x for x in pdb_list if x not in skip_pdb_list]
    pdb_list = filtered_pdb_list
    # pdb_list = dataset["optimization_v2"]
    # pdb_location = "cleaned_pdbs/first_test_set"
    pdb_location = "database/dompdb"
    if args.mode == 1:
        def pdbToFasta(pdb, pdbLocation, fastaFile, chains="A"):
            import textwrap
            from Bio.PDB.PDBParser import PDBParser
            from Bio.PDB.Polypeptide import three_to_one
            # pdb = "1r69"
            # pdbLocation = "/Users/weilu/Research/server/may_2019/family_fold/1r69.pdb"
            # chains = "A"
            # fastaFile = "/Users/weilu/Research/server/may_2019/family_fold/1r69.fasta"
            s = PDBParser(QUIET=True).get_structure("X", pdbLocation)
            # m = s[0]  # model 0
            seq = ""
            with open(fastaFile, "w") as out:
                chain = "A"
                out.write(f">{pdb.upper()}:{chain.upper()}|PDBID|CHAIN|SEQUENCE\n")
                seq = ""
                for res in s.get_residues():
                    resName = three_to_one(res.resname)
                    seq += resName
                out.write("\n".join(textwrap.wrap(seq, width=80))+"\n")
            return seq

        # a = glob.glob("cleaned_pdbs/*.pdb")
        jobIdList = []
        # a = glob.glob("../membrane_only_contact_optimization/database/dompdb/*.pdb")
        do("mkdir -p fasta")
        for pdb in pdb_list:
            print(pdb)
            # pdb = line.split("/")[-1].split(".")[0]
            pdbToFasta(pdb, f"{pdb_location}/{pdb}.pdb", f"fasta/{pdb}.fasta")
            cmd = f"/projects/pw8/wl45/hh-suite/build/bin/hhblits -i fasta/{pdb}.fasta -d ../../uniclust30_2018_08/uniclust30_2018_08 -o {pdb} -oa3m {pdb}.a3m -n 2"
            # do(cmd)
            jobId = slurmRun(f"slurms/{pdb}.slurm", cmd, memory=10)
            jobIdList.append(jobId)
    if args.mode == 2:
        do("mkdir MSA")
        for pdb in pdb_list:
            do(f"mv {pdb} MSA/")
            do(f"mv {pdb}.a3m MSA/")
    if args.mode == 3:
        do("mkdir -p alignments")
        for pdb in pdb_list:
            a = pdb
            print(a)
            data = get_MSA_data(f"MSA/{a}.a3m")
            with open(f"alignments/{a}_filtered_0.05.seqs", "w") as out:
                for line in data:
                    out.write(line+"\n")




if args.day == "mar13":
    # mass specific decoys iteration.
    # should remove 1puc, 1skz
    # pdb_list = dataset["optimization_cbd"]
    data = pd.read_csv("training_set.csv")
    specific_decoys = data.query("Length < 150 and Length > 70").reset_index(drop=True)
    pdb_list = specific_decoys["Protein"].to_list()
    pdb_list = [a.lower() for a in pdb_list]
    skip_pdb_list = ["1puc", "1skz"]
    skip_pdb_list += ["1msc", "1fmb", "1gvp", "2tgi", "1whi", "1baj", "1rmd", "1div"]  #dimer.
    skip_pdb_list += ["1aqe"]  # lots of ligand
    filtered_pdb_list = [x for x in pdb_list if x not in skip_pdb_list]
    pdb_list = filtered_pdb_list
    # print(pdb_list)
    if args.mode == 1:
        # pre = "/Users/weilu/Research/server/feb_2020/casp13_targets/setups/T0953s2-D1"
        for pdb in pdb_list:
            pre = f"setups/{pdb}/"
            # pdb = "1baj"
            original_openAWSEM_input = f"{pre}/{pdb}-openmmawsem.pdb"
            new_openAWSEM_input = f"{pre}/cbd-openmmawsem.pdb"
            all_atom_pdb_file = f"{pre}/crystal_structure-cleaned.pdb"
            replace_CB_coord_with_CBD_for_openAWSEM_input(original_openAWSEM_input, new_openAWSEM_input, all_atom_pdb_file)

if args.day == "mar10":
    if args.mode == 1:
        folder = "optimization_new_3"
        do(f"mkdir -p {folder}")
        cd(folder)
        do("cp ../training_set.csv .")
        do("gg_server.py -d mar05 -m 1")
        do("gg_server.py -d mar05 -m 12")
        do("gg_server.py -d mar05 -m 13")
        do("gg_server.py -d mar05 -m 23")
        # do("gg_server.py -d mar05 -m 5")
if args.day == "mar09":
    # mass specific decoys iteration.
    # should remove 1puc, 1skz
    data = pd.read_csv("training_set.csv")
    specific_decoys = data.query("Length < 150 and Length > 70").reset_index(drop=True)
    pdb_list = specific_decoys["Protein"].to_list()
    pdb_list = [a.lower() for a in pdb_list]
    skip_pdb_list = ["1puc", "1skz"]
    skip_pdb_list += ["1msc", "1fmb", "1gvp", "2tgi", "1whi", "1baj", "1rmd", "1div"]  #dimer.
    skip_pdb_list += ["1aqe"]  # lots of ligand
    filtered_pdb_list = [x for x in pdb_list if x not in skip_pdb_list]
    pdb_list = filtered_pdb_list
    print(f"number of pdb: {len(pdb_list)}")
    if args.mode == 1:
        do("mkdir -p outs")
        do("mkdir -p slurms")
        # i = 1
        # i = 3
        # iteration = f"iteration_1_cbd"
        # subMode = 2
        # iteration = f"iteration_start_native"
        # subMode = 1
        # iteration = f"iteration_start_native_iter2"
        # subMode = 4

        iteration = f"iteration_start_native_iter3"
        subMode = 5
        force_setup = "forces_setup.py"

        for pdb in pdb_list:
            for i in range(2):
                cmd = f"python mm_run.py setups/{pdb}/{pdb} --to {iteration}/{pdb}/{i} -s 1e6 --tempStart 300 --tempEnd 300 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 2:
        do("mkdir -p outs")
        do("mkdir -p slurms")
        # i = 1
        # i = 3
        # iteration = f"iteration_1_cbd"
        # subMode = 2
        # iteration = f"iteration_start_native"
        # subMode = 1
        # iteration = f"iteration_start_native_iter2"
        # subMode = 4

        iteration = f"iteration_new_1"
        subMode = 5
        force_setup = "forces_setup.py"

        for pdb in pdb_list:
            for i in range(2):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 3:
        do("mkdir -p outs")
        do("mkdir -p slurms")

        iteration = f"iteration_start_native_iter4"
        subMode = 6
        force_setup = "forces_setup.py"

        for pdb in pdb_list:
            for i in range(2):
                cmd = f"python mm_run.py setups/{pdb}/{pdb} --to {iteration}/{pdb}/{i} -s 1e6 --tempStart 300 --tempEnd 300 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 4:
        do("mkdir -p outs")
        do("mkdir -p slurms")

        iteration = f"iteration_start_native_iter5"
        subMode = 7
        force_setup = "forces_setup.py"

        for pdb in pdb_list:
            for i in range(2):
                cmd = f"python mm_run.py setups/{pdb}/{pdb} --to {iteration}/{pdb}/{i} -s 1e6 --tempStart 300 --tempEnd 300 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 5:
        do("mkdir -p outs")
        do("mkdir -p slurms")
        # i = 1
        # i = 3
        # iteration = f"iteration_1_cbd"
        # subMode = 2
        # iteration = f"iteration_start_native"
        # subMode = 1
        # iteration = f"iteration_start_native_iter2"
        # subMode = 4

        iteration = f"iteration_new_2"
        subMode = 7
        force_setup = "forces_setup.py"

        for pdb in pdb_list:
            for i in range(2):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 6:
        do("mkdir -p outs")
        do("mkdir -p slurms")
        # i = 1
        # i = 3
        # iteration = f"iteration_1_cbd"
        # subMode = 2
        # iteration = f"iteration_start_native"
        # subMode = 1
        # iteration = f"iteration_start_native_iter2"
        # subMode = 4

        iteration = f"iteration_new_3"
        subMode = 8
        force_setup = "forces_setup.py"

        for pdb in pdb_list:
            for i in range(2):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 7:
        do("mkdir -p outs")
        do("mkdir -p slurms")
        # i = 1
        # i = 3
        # iteration = f"iteration_1_cbd"
        # subMode = 2
        # iteration = f"iteration_start_native"
        # subMode = 1
        # iteration = f"iteration_start_native_iter2"
        # subMode = 4

        iteration = f"iteration_new_4"
        subMode = 9
        force_setup = "forces_setup.py"

        for pdb in pdb_list:
            for i in range(2):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 8:
        do("mkdir -p outs")
        do("mkdir -p slurms")

        iteration = f"iteration_native_new_4"
        subMode = 9
        force_setup = "forces_setup.py"

        for pdb in pdb_list:
            for i in range(2):
                cmd = f"python mm_run.py setups/{pdb}/cbd --to {iteration}/{pdb}/{i} -s 1e6 --tempStart 300 --tempEnd 300 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 9:
        # mar09
        do("mkdir -p outs")
        do("mkdir -p slurms")

        iteration = f"iteration_new_4_without_burial"
        subMode = 10
        force_setup = "forces_setup.py"

        for pdb in pdb_list:
            for i in range(2):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 10:
        # mar09
        do("mkdir -p outs")
        do("mkdir -p slurms")

        iteration = f"iteration_new_4_without_burial_shift_well"
        subMode = 11
        force_setup = "forces_setup.py"

        for pdb in pdb_list:
            for i in range(2):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)

if args.day == "mar08":
    print("mar08")
    from pyCodeLib import *
    import warnings
    warnings.filterwarnings('ignore')
    if args.mode == 1:
        # cutoff = 400
        cutoff = 600
        pre = "."
        gamma_pre = f"{pre}/saved_gammas"
        # trial_name = "iter1"
        trial_name = "iter2"
        trial_name = "native_iter2"
        trial_name = "native_iter3"
        trial_name = "native_iter4"
        trial_name = "new_3"
        trial_name = "new_4"
        gamma_file_name = f"{gamma_pre}/{trial_name}_cutoff{cutoff}_impose_Aprime_constraint"
        # gamma_file_name = f"{pre}/iter_2_30"
        data = validate_hamiltonian_wei("phi_list.txt", "protein_list_complete", gamma_file_name, "openMM", 100, mode=0)
        dataFile = f"optimization_2020.csv"
        data.to_csv(dataFile)
    if args.mode == 2:
        phi_list_file_name = "phi_list.txt"
        training_set_file = "protein_list_complete"
        decoyBiasName = "deocyQZBias"
        decoy_method = "openMM"
        dataFile = f"optimization_2020.csv"
        data = pd.read_csv(dataFile, index_col=0)
        data["normalized_Z_weight"] = (10/(data["Z_scores"] - data["Z_scores"].min() + 10))
        data["normalized_Z_weight_sq"] = (10/(data["Z_scores"] - data["Z_scores"].min() + 10))**2
        data["normalized_Z_weight_half"] = (10/(data["Z_scores"] - data["Z_scores"].min() + 10))**0.5
        phi_list = read_phi_list(phi_list_file_name)
        training_set = read_column_from_file(training_set_file, 1)

        protein = training_set[0]
        for i, line in data.iterrows():
            protein = line["Protein"]
            normalized_Z_weight = line["normalized_Z_weight"]
            normalized_Z_weight_sq = line["normalized_Z_weight_sq"]
            normalized_Z_weight_half = line["normalized_Z_weight_half"]
            for i_phi_function, phi_and_parameters in enumerate(phi_list):
                # print(protein, i_phi_function, phi_and_parameters)
                phi = phi_and_parameters[0]
                parameters = phi_and_parameters[1]
                i_phi = phi_list.index(phi_and_parameters)
                parameters_string = get_parameters_string(parameters)
                # try:
                fromFile = f"../phis/{phi}_{protein}_decoysQ_{decoy_method}_{parameters_string}"
                a = np.loadtxt(fromFile)
                # print(a)
                # b = (1-a) * normalized_Z_weight
                # toFile = f"../phis/{phi}_{protein}_{decoyBiasName}_{decoy_method}_{parameters_string}"
                # np.savetxt(toFile, b, fmt='%.4f')
                # b = (1-a) * normalized_Z_weight
                # b = (1-a) * normalized_Z_weight_sq
                b = (1-a) * normalized_Z_weight_half
                toFile = f"../phis/{phi}_{protein}_{decoyBiasName}_{decoy_method}_{parameters_string}"
                np.savetxt(toFile, b, fmt='%.4f')

if args.day == "mar06":
    if args.mode == 1:
        import matplotlib.pyplot as plt
        i = 2
        # save_gamma_pre = "saved_gammas"
        # trial_name = f"iter{i}"
        # trial_name = args.label
        # os.system(f"mkdir -p {save_gamma_pre}/figures")

        cutoff_i = 400
        # name = f"iter_{i}_30"
        name = args.label
        trial_name = name
        filtered_gamma = np.loadtxt(name)
        figureName = f"{trial_name}_cutoff{cutoff_i}_equal_legend"
        title = f"{trial_name}_cutoff{cutoff_i}"
        show_together_v2(filtered_gamma, figureName, title=title, inferBound=2)

if args.day == "mar05":
    # mass specific decoys iteration.
    data = pd.read_csv("training_set.csv")
    specific_decoys = data.query("Length < 150 and Length > 70").reset_index(drop=True)
    pdb_list = specific_decoys["Protein"].to_list()
    pdb_list = [a.lower() for a in pdb_list]
    skip_pdb_list = ["1puc", "1skz"]
    skip_pdb_list += ["1msc", "1fmb", "1gvp", "2tgi", "1whi", "1baj", "1rmd", "1div"]  #dimer.
    skip_pdb_list += ["1aqe"]  # lots of ligand
    filtered_pdb_list = [x for x in pdb_list if x not in skip_pdb_list]
    pdb_list = filtered_pdb_list

    # folder_list = ["iteration_0"]
    # folder_list = ["iteration_0_cbd"]
    # folder_list = ["iteration_1_cbd"]
    # folder_list = ["iteration_2_cbd"]
    # folder_list = ["iteration_start_native"]
    # folder_list = ["iteration_start_native_iter2"]
    # folder_list = ["iteration_start_native_iter3", "iteration_new_1"]
    folder_list = ["iteration_start_native_iter4"]
    folder_list = ["iteration_start_native_iter5", "iteration_new_2"]
    folder_list = ["iteration_new_3"]
    folder_list = ["iteration_new_4_without_burial"]
    # folder_list = ["iteration_1"]
    # folder_list = ["iteration_2"]
    run_n = 2
    base_path = "/scratch/wl45/mar_2020"
    simulationType = "mass_iterative_run"
    lastN_frame = 50
    if args.mode == 111:
        for folder in folder_list:
            cd(folder)
            for pdb in pdb_list:
                for i in range(run_n):
                    cd(f"{pdb}/{i}")
                    # do()
                    cmd = "python /projects/pw8/wl45/openawsem/helperFunctions/convertOpenmmTrajectoryToStandardMovie.py movie.pdb"
                    out = slurmRun(f"convert.slurm", cmd, template=scavenge_slurm)
                    print(out)
                    cd("../..")
            cd("..")

    if args.mode == 1:
        # do("cp /home/wl45/opt/optimization/phi_list_contact.txt phi_list.txt")
        # do("cp ../phi_com.txt phi_list.txt")
        do("cp ../phi_com_withoutBurial.txt phi_list.txt")
        do("mkdir proteins_name_list")
        do("mkdir slurms")
        do('mkdir -p outs')
        # do("mkdir -p decoys/multiShuffle")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p ../phis")
        # with open("protein_list", "w") as out:
        #     for pdb in pdb_list:
        #         out.write(f"{pdb}\n")
        with open("protein_list", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
    if args.mode == 11:
        do("mkdir -p database/dompdb")
        do("mkdir -p database/S20_seq")
        for pdb in pdb_list:
            do(f"cp ../mass_iterative_run/setups/{pdb}/{pdb}.pdb database/dompdb/")
        def getSeq(fileLocation):
            p = PDBParser()
            s = p.get_structure("test", fileLocation)
            seq = ""
            residues = list(s.get_residues())
            for residue in residues:
                res_id = residue.get_id()[0]
                if res_id==' ':
                    residue_name = residue.get_resname()
                    seq += three_to_one(residue_name)
            return seq
        for pdb in pdb_list:
            fileLocation = f"database/dompdb/{pdb}.pdb"
            seq = getSeq(fileLocation)
            opt_pdbName = pdb
            fileLocation = f"database/S20_seq/{opt_pdbName}.seq"
            with open(fileLocation, "w") as out:
                out.write(seq+"\n")
    if args.mode == 12:
        outFile = read_simulation_info(folder_list, pdb_list, simulationType, base_path, run_n)
        # today = datetime.today().strftime('%m-%d')
        outFile = f"{simulationType}.csv"
    if args.mode == 13:
        # time.sleep(3600)
        dataFile = f"{simulationType}.csv"
        simulation_folder = simulationType
        data = pd.read_csv(dataFile, index_col=0)
        parser = PDBParser()
        # folder = "first"
        # for folder in ["first", "first_cpu2"]:
        for folder in folder_list:
            pre = f"{base_path}/{simulation_folder}/{folder}"
            to_folder = "."
            os.system(f"mkdir -p {to_folder}/decoys/openMM")
            # pdb_list = ['1j5u']
            for pdb in pdb_list:
                complete_models = []
                print(pdb)
                for i in range(run_n):
                    movieFile = f"{pre}/{pdb}/{i}/movie.pdb"
                    allFrames, n, size = getAllFrames(movieFile)
                    num_of_frames = int(n/size)
                    first_chosen_frame = num_of_frames - lastN_frame
                    # last_chosen_frame = num_of_frames
                    oneFrame = allFrames[size*first_chosen_frame:size*(num_of_frames)]

                    p = PDBParser()
                    f = io.StringIO("".join(oneFrame))
                    s = p.get_structure("test", f)
                    complete_models += list(s.get_models())
                t = data.query(f"Protein == '{pdb}' and Folder == '{folder}' and Steps > 1").reset_index(drop=True)
                t = t.groupby("Run").tail(50).reset_index(drop=True)
                t["structure"] = complete_models
                t = t.rename(columns={"Q":"Qw"})
                # last50 = t.groupby("Run").tail(50).reset_index(drop=True)
                last50 = t
                # print(last50.head())
                # print(last50.tail())
                # print(last50.shape)
                # print(len(complete_models))
                to_folder = "."
                last50.to_pickle(f"{to_folder}/decoys/openMM/{pdb}_{folder}.pkl")

    if args.mode == 2:
        # split proteins_name_list
        n = split_proteins_name_list(pdb_per_txt=2)
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d mar05 -m 22 -l {proteins}", template=base_slurm, memory=10)
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
    if args.mode == 22:
        protein = args.label
        do(f"python3 ~/opt/compute_phis.py -m 7 {protein}")
        # do("python3 ~/opt/compute_phis.py -m 0 test_protein")
    if args.mode == 23:
        # split proteins_name_list
        n = split_proteins_name_list(pdb_per_txt=2)
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d mar05 -m 22 -l {proteins}", template=base_slurm, memory=10)
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
        waitForJobs(jobIdList, sleepInterval=60)
        # folder_list = ["iteration_0_cbd", "iteration_1_cbd", "iteration_2_cbd"]
        # folder_list = ["iteration_start_native"]
        # folder_list = ["iteration_start_native", "iteration_start_native_iter2"]
        folder_list = ["iteration_start_native", "iteration_start_native_iter2",        "iteration_start_native_iter3", "iteration_new_1", "iteration_start_native_iter4"]
        folder_list += ["iteration_start_native_iter5", "iteration_new_2", "iteration_new_3"]
        # folder_list += ["iteration_new_3"]
        with open("protein_list_complete", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d mar05 -m 33"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")

    if args.mode == 3:
        # folder_list = ["iteration_0", "iteration_1", "iteration_2"]
        # folder_list = ["iteration_0_cbd", "iteration_1_cbd", "iteration_2_cbd"]
        # folder_list = ["iteration_start_native"]
        # folder_list = ["iteration_start_native", "iteration_start_native_iter2"]
        # folder_list = ["iteration_start_native", "iteration_start_native_iter2",        "iteration_start_native_iter3", "iteration_new_1", "iteration_start_native_iter4"]
        # folder_list += ["iteration_start_native_iter5", "iteration_new_2", "iteration_new_3"]
        folder_list = ["iteration_0_cbd", "iteration_1_cbd", "iteration_2_cbd"]
        folder_list += ["iteration_new_1", "iteration_new_2", "iteration_new_3"]
        with open("protein_list_complete", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d mar05 -m 33"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 33:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list_complete"
        n_decoys = 50 * run_n
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='openMM',
                                        withBiased=True, oneMinus=True, decoyBiasName='decoysQ', num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, )
        # next: optimization_analysze.py iter{i}
    if args.mode == 4:
        # folder_list = ["iteration_0", "iteration_1", "iteration_2"]
        # folder_list = ["iteration_0_cbd", "iteration_1_cbd", "iteration_2_cbd"]
        # with open("protein_list_complete", "w") as out:
        #     for folder in folder_list:
        #         for pdb in pdb_list:
        #             out.write(f"{pdb}_{folder}\n")
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d mar05 -m 44"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 44:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list_complete"
        n_decoys = 50 * run_n
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='openMM',
                                        withBiased=True, oneMinus=False, decoyBiasName='deocyQZBias', num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, )
        # next: optimization_analysze.py iter{i}
    if args.mode == 5:
        # i = 2
        # i = 2
        # i = 3
        # i = 0
        # i = 1
        # i = 2
        i = "native_iter2"
        i = "native_iter3"
        i = "native_iter4"
        i = "new_3"
        i = "new_4"
        do(f"optimization_analyze.py {i} --proteinList protein_list_complete -c -100")
        # do(f"optimization_analyze.py {i} --proteinList protein_list_complete")
        # fromFile = "/scratch/wl45/mar_2020/cath_dataset_shuffle_optimization/optimization_iter0_com_density_rmin_2p5/saved_gammas/iter0_cutoff600_impose_Aprime_constraint"

        # do(f"optimization_analyze.py {i} --proteinList protein_list_complete --gammaFile {fromFile}")

    if args.mode == 6:
        # i = 0
        # i = 1
        # alpha = 0.1
        # i = 2
        i = "native_iter2"
        i = "native_iter3"
        i = "native_iter4"
        alpha = 0.9
        cutoff = 400
        # i = 2
        fromFile = "/scratch/wl45/mar_2020/cath_dataset_shuffle_optimization/optimization_iter0_com_density_rmin_2p5/saved_gammas/iter0_cutoff600_impose_Aprime_constraint"
        # ../cath_dataset_shuffle_optimization/optimization_iter0_com_density_rmin_2p5/Oct26_saved_gammas/iter0_cutoff600_impose_Aprime_constraint/gamma.dat
        # do("mix_gamma.py ~/opt/parameters/original_gamma saved_gammas/iter0_cutoff100_impose_Aprime_constraint -i 0")
        # do(f"mix_gamma.py ~/opt/parameters/original_gamma saved_gammas/iter{i}_cutoff300_impose_Aprime_constraint -i {i}")
        # do(f"mix_gamma.py ~/opt/parameters/original_gamma saved_gammas/iter{i}_cutoff400_impose_Aprime_constraint -i {i}")
        do(f"mix_gamma.py {fromFile} saved_gammas/{i}_cutoff{cutoff}_impose_Aprime_constraint -i {i} -a {alpha}")
    if args.mode == 66:
        import matplotlib.pyplot as plt
        # i = 2
        save_gamma_pre = "saved_gammas"
        trial_name = f"iter{i}"
        # trial_name = args.label
        # os.system(f"mkdir -p {save_gamma_pre}/figures")

        cutoff_i = 400
        name = f"iter_{i}_30"
        filtered_gamma = np.loadtxt(name)
        figureName = f"{trial_name}_cutoff{cutoff_i}_equal_legend"
        title = f"{trial_name}_cutoff{cutoff_i}"
        show_together_v2(filtered_gamma, figureName, title=title, inferBound=2)


    if args.mode == 7:
        # i = 1
        i = "native_iter3"
        i = "native_iter4"
        i = "new_3"
        i = "new_4"
        cutoff = 600
        for pdb in pdb_list:
            do(f"cp Oct26_saved_gammas/{i}_cutoff{cutoff}_impose_Aprime_constraint/gamma.dat {base_path}/{simulationType}/setups/{pdb}/{i}_withoutBurial_gamma.dat")
            # do(f"cp Oct26_saved_gammas/{i}_cutoff{cutoff}_impose_Aprime_constraint/gamma.dat {base_path}/{simulationType}/setups/{pdb}/{i}_gamma.dat")
            # do(f"cp Oct26_saved_gammas/{i}_cutoff{cutoff}_impose_Aprime_constraint/burial_gamma.dat  {base_path}/{simulationType}/setups/{pdb}/{i}_burial_gamma.dat")
    if args.mode == 77:
        # i = 1
        i = "native"
        alpha = 10

        # i = 0

        for pdb in pdb_list:
            do(f"cp for_simulation/iteration_iter_{i}_gamma_{alpha}.dat {base_path}/{simulationType}/setups/{pdb}/")
            do(f"cp for_simulation/iteration_iter_{i}_burial_gamma_{alpha}.dat {base_path}/{simulationType}/setups/{pdb}/")

    if args.mode == 8:
        do("mkdir -p outs")
        do("mkdir -p slurms")
        # i = 1
        # i = 3
        # iteration = f"iteration_1_cbd"
        # subMode = 2
        iteration = f"iteration_2_cbd"
        subMode = 3
        force_setup = "forces_setup.py"

        for pdb in pdb_list:
            for i in range(2):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)

if args.day == "mar02":
    if args.mode == 1:
        # Cornichon
        # start from specific topology, 3TH
        # run = "simulation_chainE_positive_inside_submode_4"

        subModeList = [4]
        # force_setup = "forces_setup_chainE_3TH.py"
        # for i in range(5):
        #     for subMode in subModeList:
        #         cmd = f"python mm_run_with_pulling_start.py chain_E/extended -f {force_setup} --subMode {subMode} --to {run}/3TH_{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
        #         out = slurmRun(f"slurms/gpu_feb05_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
        #         print(out)

        # Cornichon
        # start from specific topology, 4TH
        run = "simulation_chainE_positive_inside_submode_4"
        force_setup = "forces_setup_chainE_4TH_hand_modify.py"
        for i in range(5):
            for subMode in subModeList:
                cmd = f"python mm_run_with_pulling_start.py chain_E/extended -f {force_setup} --subMode {subMode} --to {run}/4TH_{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                out = slurmRun(f"slurms/gpu_feb05_4th_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)


if args.day == "mar03":
    data = pd.read_csv("training_set.csv")
    specific_decoys = data.query("Length < 150 and Length > 70").reset_index(drop=True)
    pdb_list = specific_decoys["Protein"].to_list()
    pdb_list = [a.lower() for a in pdb_list]
    if args.mode == 1:
        iteration = 0
        for pdb in pdb_list:
            fromPath = "/scratch/wl45/mar_2020/cath_dataset_shuffle_optimization/optimization_iter0_com_density_rmin_2p5/Oct26_saved_gammas/iter0_cutoff600_impose_Aprime_constraint"
            toPath = f"/scratch/wl45/mar_2020/mass_iterative_run/setups/{pdb}"
            do(f"cp {fromPath}/gamma.dat {toPath}/iter{iteration}_gamma.dat")
            do(f"cp {fromPath}/burial_gamma.dat {toPath}/iter{iteration}_burial_gamma.dat")
    if args.mode == 2:
        do("mkdir -p outs")
        do("mkdir -p slurms")
        iteration = "iteration_0"
        force_setup = "forces_setup.py"
        subMode = 1
        for pdb in pdb_list:
            for i in range(2):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 3:
        do("mkdir -p outs")
        do("mkdir -p slurms")
        iteration = "iteration_0"
        force_setup = "forces_setup.py"
        subMode = 1
        for pdb in pdb_list:
            for i in range(1):
                cmd = f"python mm_run.py setups/{pdb}/{pdb} -f forces_setup.py --to native_4e4/{pdb} -s 4e4 --reportFrequency 1000 --tempStart 300 --tempEnd 300 --subMode 1 --platform CPU"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/quick_gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                print(out)

if args.day == "mar01":
    pdb_list = dataset["may13"]
    pdb_list += ["1uzc", "1ccr", "1jwe", "T0172_2"]
    # 	T0172 is 1n2x
    # T0129 is 1izm
    # T170 is 1H40 superseded by 1uzc
    if args.mode == 1:
        run = "run2"
        do("mkdir -p slurms")
        do("mkdir -p outs")
        subModeList = [11]
        force_setup = "forces_setup.py"
        for i in range(10):
            for pdb in pdb_list:
                for subMode in subModeList:
                    cmd = f"python mm_run.py setups/{pdb}/extended -f {force_setup} --subMode {subMode} --to {run}/{pdb}/{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                    out = slurmRun(f"slurms/gpu_feb16_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                    print(out)
    if args.mode == 2:
        run = "run2"
        do("mkdir -p slurms")
        do("mkdir -p outs")
        subModeList = [12]
        force_setup = "forces_setup.py"
        for i in range(10):
            for pdb in pdb_list:
                for subMode in subModeList:
                    cmd = f"python mm_run.py setups/{pdb}/extended -f {force_setup} --subMode {subMode} --to {run}/{pdb}/{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                    out = slurmRun(f"slurms/gpu_feb16_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                    print(out)
    if args.mode == 3:
        run = "server_native"
        do("mkdir -p slurms")
        do("mkdir -p outs")
        subModeList = [7, 11, 12]
        force_setup = "forces_setup.py"
        for i in range(3):
            for pdb in pdb_list:
                for subMode in subModeList:
                    cmd = f"python mm_run.py setups/{pdb}/{pdb} -f {force_setup} --subMode {subMode} --to {run}/{pdb}/{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 300 --tempEnd 300 --platform CUDA"
                    out = slurmRun(f"slurms/gpu_feb16_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                    print(out)

if args.day == "feb18":
    # casp13_targets
    pdb_list = ["T0951-D1", "T0953s2-D1", "T0955-D1", "T0957s1-D1", "T0957s1-D2", "T0958-D1", "T0960-D5", "T0963-D3", "T0968s1-D1", "T1008-D1"]

    if args.mode == 1:
        run = "side_chain_run1"
        do("mkdir -p slurms")
        do("mkdir -p outs")
        subModeList = [4, 5]
        force_setup = "forces_setup.py"
        for i in range(1, 10):
            for pdb in pdb_list:
                for subMode in subModeList:
                    cmd = f"python mm_run.py setups/{pdb}/extended -f {force_setup} --subMode {subMode} --to {run}/{pdb}/{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                    out = slurmRun(f"slurms/feb18_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                    print(out)
    if args.mode == 2:
        run = "native"
        do("mkdir -p slurms")
        do("mkdir -p outs")
        subModeList = [4, 5]
        # subModeList = [4]
        subModeList = [6]
        force_setup = "forces_setup.py"
        for i in range(1):
            for pdb in pdb_list:
                for subMode in subModeList:
                    cmd = f"python mm_run.py setups/{pdb}/{pdb} -f {force_setup} --subMode {subMode} --to {run}/{pdb}/{subMode}_{i} -s 2e4 --reportFrequency 1000 --tempStart 300 --tempEnd 300 --platform CPU"
                    out = slurmRun(f"slurms/native_feb18_{run}_{subMode}_{i}.slurm", cmd, template=base_slurm, memory=10)
                    print(out)
    if args.mode == 3:
        iteration = 0
        for pdb in pdb_list:
            do(f"cp iter{iteration}_cutoff500_impose_Aprime_constraint/gamma.dat /scratch/wl45/feb_2020/casp13_targets/setups/{pdb}/iter{iteration}_gamma.dat")
            do(f"cp iter{iteration}_cutoff500_impose_Aprime_constraint/burial_gamma.dat /scratch/wl45/feb_2020/casp13_targets/setups/{pdb}/iter{iteration}_burial_gamma.dat")
    if args.mode == 4:
        iteration = 0
        for pdb in pdb_list:
            do(f"cp iter{iteration}_cutoff500_impose_Aprime_constraint/gamma.dat /scratch/wl45/feb_2020/casp13_targets/setups/{pdb}/rmin_35_iter{iteration}_gamma.dat")
            do(f"cp iter{iteration}_cutoff500_impose_Aprime_constraint/burial_gamma.dat /scratch/wl45/feb_2020/casp13_targets/setups/{pdb}/rmin_35_iter{iteration}_burial_gamma.dat")
if args.day == "feb16":
    pdb_list = dataset["may13"]
    pdb_list += ["1uzc", "1ccr", "1jwe", "T0172_2"]
    # 	T0172 is 1n2x
    # T0129 is 1izm
    # T170 is 1H40 superseded by 1uzc

    if args.mode == 1:
        # center of mass. side chain
        run = "side_chain_run1"
        do("mkdir -p slurms")
        do("mkdir -p outs")
        subModeList = [0, 1]
        force_setup = "forces_setup.py"
        for i in range(5):
            for subMode in subModeList:
                cmd = f"python mm_run.py extended -f {force_setup} --subMode {subMode} --to {run}/{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                out = slurmRun(f"slurms/gpu_feb16_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 2:
        time.sleep(3600)
        run = "side_chain_run1"
        do("mkdir -p slurms")
        do("mkdir -p outs")
        subModeList = [0, 1]
        force_setup = "forces_setup.py"
        for i in range(1):
            for pdb in pdb_list:
                for subMode in subModeList:
                    cmd = f"python mm_run.py setups/{pdb}/extended -f {force_setup} --subMode {subMode} --to {run}/{pdb}/{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                    out = slurmRun(f"slurms/gpu_feb16_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                    print(out)

        time.sleep(3600)
        run = "side_chain_run1"
        do("mkdir -p slurms")
        do("mkdir -p outs")
        subModeList = [0, 1]
        force_setup = "forces_setup.py"
        for i in range(1, 10):
            for pdb in pdb_list:
                for subMode in subModeList:
                    cmd = f"python mm_run.py setups/{pdb}/extended -f {force_setup} --subMode {subMode} --to {run}/{pdb}/{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                    out = slurmRun(f"slurms/gpu_feb16_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                    print(out)
    if args.mode == 3:
        run = "side_chain_run1"
        do("mkdir -p slurms")
        do("mkdir -p outs")
        subModeList = [2, 3]
        force_setup = "forces_setup.py"
        for i in range(10):
            for pdb in pdb_list:
                for subMode in subModeList:
                    cmd = f"python mm_run.py setups/{pdb}/extended -f {force_setup} --subMode {subMode} --to {run}/{pdb}/{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                    out = slurmRun(f"slurms/gpu_feb16_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                    print(out)
    if args.mode == 4:
        run = "side_chain_run1"
        do("mkdir -p slurms")
        do("mkdir -p outs")
        subModeList = [4, 5]
        force_setup = "forces_setup.py"
        for i in range(10):
            for pdb in pdb_list:
                for subMode in subModeList:
                    cmd = f"python mm_run.py setups/{pdb}/extended -f {force_setup} --subMode {subMode} --to {run}/{pdb}/{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                    out = slurmRun(f"slurms/gpu_feb16_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                    print(out)
    if args.mode == 5:
        iteration = 0
        for pdb in pdb_list:
            do(f"cp iter{iteration}_cutoff500_impose_Aprime_constraint/gamma.dat /scratch/wl45/feb_2020/compare_side_chain_with_and_without/setups/{pdb}/iter{iteration}_gamma.dat")
            do(f"cp iter{iteration}_cutoff500_impose_Aprime_constraint/burial_gamma.dat /scratch/wl45/feb_2020/compare_side_chain_with_and_without/setups/{pdb}/iter{iteration}_burial_gamma.dat")
    if args.mode == 6:
        iteration = 0
        for pdb in pdb_list:
            do(f"cp iter{iteration}_cutoff500_impose_Aprime_constraint/gamma.dat /scratch/wl45/feb_2020/compare_side_chain_with_and_without/setups/{pdb}/rmin_35_iter{iteration}_gamma.dat")
            do(f"cp iter{iteration}_cutoff500_impose_Aprime_constraint/burial_gamma.dat /scratch/wl45/feb_2020/compare_side_chain_with_and_without/setups/{pdb}/rmin_35_iter{iteration}_burial_gamma.dat")
    if args.mode == 7:
        run = "side_chain_run1"
        do("mkdir -p slurms")
        do("mkdir -p outs")
        subModeList = [6]
        force_setup = "forces_setup.py"
        for i in range(10):
            for pdb in pdb_list:
                for subMode in subModeList:
                    cmd = f"python mm_run.py setups/{pdb}/extended -f {force_setup} --subMode {subMode} --to {run}/{pdb}/{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                    out = slurmRun(f"slurms/gpu_feb16_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                    print(out)
    if args.mode == 8:
        run = "side_chain_run1"
        do("mkdir -p slurms")
        do("mkdir -p outs")
        subModeList = [7]
        force_setup = "forces_setup.py"
        for i in range(10):
            for pdb in pdb_list:
                for subMode in subModeList:
                    cmd = f"python mm_run.py setups/{pdb}/extended -f {force_setup} --subMode {subMode} --to {run}/{pdb}/{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                    out = slurmRun(f"slurms/gpu_feb16_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                    print(out)

if args.day == "feb12":
    if args.mode == 1:
        # Cornichon
        # start from specific topology, 3TH
        run = "simulation_chainE_positive_inside_submode_4"

        subModeList = [4]
        force_setup = "forces_setup_chainE_3TH.py"
        for i in range(5):
            for subMode in subModeList:
                cmd = f"python mm_run_with_pulling_start.py chain_E/extended -f {force_setup} --subMode {subMode} --to {run}/3TH_{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                out = slurmRun(f"slurms/gpu_feb05_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)

        # Cornichon
        # start from specific topology, 4TH
        run = "simulation_chainE_positive_inside"
        force_setup = "forces_setup_chainE_4TH_hand_modify.py"
        for i in range(5):
            for subMode in subModeList:
                cmd = f"python mm_run_with_pulling_start.py chain_E/extended -f {force_setup} --subMode {subMode} --to {run}/4TH_{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                out = slurmRun(f"slurms/gpu_feb05_4th_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
if args.day == "feb11":
    if args.mode == 1:
        i = 3
        do("gg_server.py -d feb08 -m 111 -l ")
        cd("../mass_specific_decoys/")
        do(f"mkdir optimization_iter{i}")
        cd(f"optimization_iter{i}")
        do("cp ../training_set.csv .")
        do("gg_server.py -d feb08 -m 12 -l ")
        do("gg_server.py -d feb08 -m 13 -l ")
        do("gg_server.py -d feb08 -m 2 -l ")
if args.day == "feb10":
    if args.mode == 1:
        name = f"iter_0_30"
        filtered_gamma = np.loadtxt(name)
        figureName = f"iter_0_30"
        title = figureName
        # inferBound=2, equal vmin, vmax
        show_together_v2(filtered_gamma, figureName, title=title, inferBound=2, invert_sign=True)
    if args.mode == 2:
        name = f"/home/wl45/opt/parameters/original_gamma"
        filtered_gamma = np.loadtxt(name)
        figureName = f"original"
        title = figureName
        # inferBound=2, equal vmin, vmax
        show_together_v2(filtered_gamma, figureName, title=title, inferBound=2, invert_sign=True)
    if args.mode == 3:
        i = 1
        name = f"iter_{i}_30"
        filtered_gamma = np.loadtxt(name)
        figureName = f"iter_{i}_30"
        title = figureName
        # inferBound=2, equal vmin, vmax
        show_together_v2(filtered_gamma, figureName, title=title, inferBound=2, invert_sign=True)
    if args.mode == 4:
        i = 2
        name = f"iter_{i}_30"
        filtered_gamma = np.loadtxt(name)
        figureName = f"iter_{i}_30"
        title = figureName
        # inferBound=2, equal vmin, vmax
        show_together_v2(filtered_gamma, figureName, title=title, inferBound=2, invert_sign=True)
if args.day == "feb08":
    # mass specific decoys iteration.
    data = pd.read_csv("training_set.csv")
    specific_decoys = data.query("Length < 150 and Length > 70").reset_index(drop=True)
    pdb_list = specific_decoys["Protein"].to_list()

    # folder_list = ["iteration_0"]
    # folder_list = ["iteration_1"]
    folder_list = ["iteration_2"]
    run_n = 1
    base_path = "/scratch/wl45/feb_2020"
    simulationType = "mass_iterative_run"
    lastN_frame = 50
    if args.mode == 111:
        for folder in folder_list:
            cd(folder)
            for pdb in pdb_list:
                for i in range(run_n):
                    cd(f"{pdb}/{i}")
                    do("python /projects/pw8/wl45/openawsem/helperFunctions/convertOpenmmTrajectoryToStandardMovie.py movie.pdb")
                    cd("../..")
            cd("..")
    if args.mode == 1:
        do("cp /home/wl45/opt/optimization/phi_list_contact.txt phi_list.txt")
        do("mkdir proteins_name_list")
        do("mkdir slurms")
        do('mkdir -p outs')
        # do("mkdir -p decoys/multiShuffle")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p ../phis")
        # with open("protein_list", "w") as out:
        #     for pdb in pdb_list:
        #         out.write(f"{pdb}\n")
        with open("protein_list", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
    if args.mode == 11:
        do("mkdir -p database/dompdb")
        do("mkdir -p database/S20_seq")
        for pdb in pdb_list:
            do(f"cp ../mass_iterative_run/setups/{pdb}/{pdb}.pdb database/dompdb/")
        def getSeq(fileLocation):
            p = PDBParser()
            s = p.get_structure("test", fileLocation)
            seq = ""
            residues = list(s.get_residues())
            for residue in residues:
                res_id = residue.get_id()[0]
                if res_id==' ':
                    residue_name = residue.get_resname()
                    seq += three_to_one(residue_name)
            return seq
        for pdb in pdb_list:
            fileLocation = f"database/dompdb/{pdb}.pdb"
            seq = getSeq(fileLocation)
            opt_pdbName = pdb
            fileLocation = f"database/S20_seq/{opt_pdbName}.seq"
            with open(fileLocation, "w") as out:
                out.write(seq+"\n")
    if args.mode == 12:
        outFile = read_simulation_info(folder_list, pdb_list, simulationType, base_path, run_n)
        # today = datetime.today().strftime('%m-%d')
        outFile = f"{simulationType}.csv"
    if args.mode == 13:
        # time.sleep(3600)
        dataFile = f"{simulationType}.csv"
        simulation_folder = simulationType
        data = pd.read_csv(dataFile, index_col=0)
        parser = PDBParser()
        # folder = "first"
        # for folder in ["first", "first_cpu2"]:
        for folder in folder_list:
            pre = f"{base_path}/{simulation_folder}/{folder}"
            to_folder = "."
            os.system(f"mkdir -p {to_folder}/decoys/openMM")
            # pdb_list = ['1j5u']
            for pdb in pdb_list:
                complete_models = []
                print(pdb)
                for i in range(run_n):
                    movieFile = f"{pre}/{pdb}/{i}/movie.pdb"
                    allFrames, n, size = getAllFrames(movieFile)
                    num_of_frames = int(n/size)
                    first_chosen_frame = num_of_frames - lastN_frame
                    # last_chosen_frame = num_of_frames
                    oneFrame = allFrames[size*first_chosen_frame:size*(num_of_frames)]

                    p = PDBParser()
                    f = io.StringIO("".join(oneFrame))
                    s = p.get_structure("test", f)
                    complete_models += list(s.get_models())
                t = data.query(f"Protein == '{pdb}' and Folder == '{folder}' and Steps > 1").reset_index(drop=True)
                t = t.groupby("Run").tail(50).reset_index(drop=True)
                t["structure"] = complete_models
                t = t.rename(columns={"Q":"Qw"})
                # last50 = t.groupby("Run").tail(50).reset_index(drop=True)
                last50 = t
                # print(last50.head())
                # print(last50.tail())
                # print(last50.shape)
                # print(len(complete_models))
                to_folder = "."
                last50.to_pickle(f"{to_folder}/decoys/openMM/{pdb}_{folder}.pkl")

    if args.mode == 2:
        # split proteins_name_list
        n = split_proteins_name_list(pdb_per_txt=2)
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d feb08 -m 22 -l {proteins}", template=base_slurm, memory=10)
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
    if args.mode == 22:
        protein = args.label
        do(f"python3 ~/opt/compute_phis.py -m 7 {protein}")
        # do("python3 ~/opt/compute_phis.py -m 0 test_protein")
    if args.mode == 3:
        folder_list = ["iteration_0", "iteration_1", "iteration_2"]
        with open("protein_list_complete", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d feb08 -m 33"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 33:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list_complete"
        n_decoys = 50
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='openMM',
                                        withBiased=True, num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, )
        # next: optimization_analysze.py iter{i}
    if args.mode == 5:
        i = 2
        # i = 2
        # i = 3
        # i = 0
        do(f"optimization_analyze.py iter{i} --proteinList protein_list_complete")
    if args.mode == 6:
        # i = 1
        i = 2
        # do("mix_gamma.py ~/opt/parameters/original_gamma saved_gammas/iter0_cutoff100_impose_Aprime_constraint -i 0")
        # do(f"mix_gamma.py ~/opt/parameters/original_gamma saved_gammas/iter{i}_cutoff300_impose_Aprime_constraint -i {i}")
        do(f"mix_gamma.py ~/opt/parameters/original_gamma saved_gammas/iter{i}_cutoff400_impose_Aprime_constraint -i {i}")
    if args.mode == 7:
        # i = 1
        i = 2
        for pdb in pdb_list:
            do(f"cp ../mass_specific_decoys/optimization_iter{i}/for_simulation/iteration_iter_{i}_gamma_30.dat setups/{pdb}/")
            do(f"cp ../mass_specific_decoys/optimization_iter{i}/for_simulation/iteration_iter_{i}_burial_gamma_30.dat setups/{pdb}/")
    if args.mode == 8:
        do("mkdir -p outs")
        do("mkdir -p slurms")
        # i = 2
        i = 3
        iteration = f"iteration_{i}"
        force_setup = "forces_setup.py"
        subMode = i
        for pdb in pdb_list:
            for i in range(1):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
if args.day == "feb07":
    pdb_list = ["6f45"]
    if args.mode == 1:
        do("cp /home/wl45/opt/optimization/phi_list_contact.txt phi_list.txt")
        do("mkdir proteins_name_list")
        do("mkdir slurms")
        do('mkdir -p outs')
        # do("mkdir -p decoys/multiShuffle")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p ../phis")
        with open("protein_list", "w") as out:
            for pdb in pdb_list:
                out.write(f"{pdb}\n")
    if args.mode ==222:
        n = split_proteins_name_list(pdb_per_txt=5)
    if args.mode == 2:
        # split proteins_name_list
        n = split_proteins_name_list(pdb_per_txt=5)
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d feb07 -m 22 -l {proteins}", template=base_slurm, memory=10)
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
    if args.mode == 22:
        from pyCodeLib import *
        proteins = args.label
        n_decoys = 1000
        generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
        do(f"python3 ~/opt/compute_phis.py -m 0 {proteins}")
    if args.mode == 3:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d feb07 -m 33"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 33:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list"
        n_decoys = 1000
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='shuffle',
                                        withBiased=False, num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, )
        # next: optimization_analysze.py iter{i}
    if args.mode == 5:
        # i = 1
        # i = 2
        # i = 3
        i = 0
        do(f"optimization_analyze.py iter{i}")
    if args.mode == 6:
        # i = 1
        # i = 2
        # i = 3
        i = 0
        do(f"optimization_analyze.py iter{i} -c -87")

if args.day == "feb06":
    if args.mode == 1:
        import matplotlib.pyplot as plt
        cutoff_list = [100, 200, 300, 400, 500, 600, 700, 800]
        # cutoff_list += [10, 20, 30, 40, 50, 80]
        save_gamma_pre = "saved_gammas"
        # trial_name = "iter1"
        # trial_name = args.label
        trial_name = "iter0"
        os.system(f"mkdir -p {save_gamma_pre}/figures")
        for cutoff_i in cutoff_list:
            name = f"{save_gamma_pre}/{trial_name}_cutoff{cutoff_i}_impose_Aprime_constraint"
            filtered_gamma = np.loadtxt(name)
            figureName = f"{save_gamma_pre}/figures/{trial_name}_cutoff{cutoff_i}_equal_legend_4_well"
            title = f"{trial_name}_cutoff{cutoff_i}_4well"
            show_together_v2(filtered_gamma, figureName, title=title, inferBound=2, n=4,
            ax_title_list = ["Direct", "High density(protein)", "Partial", "Low density(water)"])
            # filtered_gamma = filtered_gamma[210:]
            # show_together_v2(filtered_gamma, figureName, title=title, inferBound=2)
            # figureName = f"{save_gamma_pre}/figures/{trial_name}_cutoff{cutoff_i}_equal_legend"
            # show_together(filtered_gamma, figureName, title=title, inferBound=2)
    if args.mode == 11:
        pdb_list = ["6f45"]
        do("mkdir -p database/dompdb")
        do("mkdir -p database/S20_seq")
        for pdb in pdb_list:
            do(f"cp ../multimer_folding/cleaned_pdbs/{pdb}.pdb database/dompdb/")
    if args.mode == 12:
        pdb_list = ["6f45"]
        def getSeq(fileLocation):
            p = PDBParser()
            s = p.get_structure("test", fileLocation)
            seq = ""
            residues = list(s.get_residues())
            for residue in residues:
                res_id = residue.get_id()[0]
                if res_id==' ':
                    residue_name = residue.get_resname()
                    seq += three_to_one(residue_name)
            return seq
        for pdb in pdb_list:
            fileLocation = f"database/dompdb/{pdb}.pdb"
            seq = getSeq(fileLocation)
            opt_pdbName = pdb
            fileLocation = f"database/S20_seq/{opt_pdbName}.seq"
            with open(fileLocation, "w") as out:
                out.write(seq+"\n")

    if args.mode == 2:
        name = f"protein_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_phi_native_summary.txt"
        filtered_gamma = np.loadtxt(name)
        figureName = f"native_phi_change_limit"
        title = figureName
        show_together(filtered_gamma, figureName, title=title, inferBound=1, invert_sign=False)

        name = f"protein_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_phi_decoy_summary.txt"
        filtered_gamma = np.loadtxt(name)
        figureName = f"decoy_phi_change_limit"
        title = figureName
        show_together(filtered_gamma, figureName, title=title, inferBound=1, invert_sign=False)
    if args.mode == 3:
        name = f"log_native.txt"
        filtered_gamma = np.loadtxt(name)
        figureName = f"native_phi_log"
        title = figureName
        show_together(filtered_gamma, figureName, title=title, inferBound=1, invert_sign=False)

        name = f"log_decoy.txt"
        filtered_gamma = np.loadtxt(name)
        figureName = f"decoy_phi_log"
        title = figureName
        show_together(filtered_gamma, figureName, title=title, inferBound=1, invert_sign=False)
    if args.mode == 4:
        name = f"protein_list_phi_pairwise_contact_well3.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_phi_native_summary.txt"
        filtered_gamma = np.loadtxt(name)
        figureName = f"native_phi_change_limit_r3_5"
        title = figureName
        show_together(filtered_gamma, figureName, title=title, inferBound=1, invert_sign=False)

        name = f"protein_list_phi_pairwise_contact_well3.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_phi_decoy_summary.txt"
        filtered_gamma = np.loadtxt(name)
        figureName = f"decoy_phi_change_limit_r3_5"
        title = figureName
        show_together(filtered_gamma, figureName, title=title, inferBound=1, invert_sign=False)
if args.day == "feb05":
    if args.mode == 1:
        # Cornichon
        # start from specific topology, 3TH
        run = "simulation_chainE_positive_inside"

        subModeList = [2]
        force_setup = "forces_setup_chainE_3TH.py"
        for i in range(5):
            for subMode in subModeList:
                cmd = f"python mm_run_with_pulling_start.py chain_E/extended -f {force_setup} --subMode {subMode} --to {run}/3TH_{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                out = slurmRun(f"slurms/gpu_feb05_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 2:
        # Cornichon
        # start from specific topology, 4TH
        run = "simulation_chainE_positive_inside"

        subModeList = [2]
        force_setup = "forces_setup_chainE_4TH_hand_modify.py"
        for i in range(5):
            for subMode in subModeList:
                cmd = f"python mm_run_with_pulling_start.py chain_E/extended -f {force_setup} --subMode {subMode} --to {run}/4TH_{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                out = slurmRun(f"slurms/gpu_feb05_4th_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 3:
        # Cornichon
        # start from specific topology, 3TH
        run = "simulation_chainE_positive_inside_submode_3"

        subModeList = [3]
        force_setup = "forces_setup_chainE_3TH.py"
        for i in range(5):
            for subMode in subModeList:
                cmd = f"python mm_run_with_pulling_start.py chain_E/extended -f {force_setup} --subMode {subMode} --to {run}/3TH_{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                out = slurmRun(f"slurms/gpu_feb05_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)

        # Cornichon
        # start from specific topology, 4TH
        run = "simulation_chainE_positive_inside"
        force_setup = "forces_setup_chainE_4TH_hand_modify.py"
        for i in range(5):
            for subMode in subModeList:
                cmd = f"python mm_run_with_pulling_start.py chain_E/extended -f {force_setup} --subMode {subMode} --to {run}/4TH_{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                out = slurmRun(f"slurms/gpu_feb05_4th_{run}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)

if args.day == "feb03":
    data = pd.read_csv("training_set.csv")
    specific_decoys = data.query("Length < 150 and Length > 70").reset_index(drop=True)
    pdb_list = specific_decoys["Protein"].to_list()
    if args.mode == 1:
        do("mkdir -p outs")
        do("mkdir -p slurms")
        iteration = "iteration_0"
        force_setup = "forces_setup.py"
        subMode = -1
        for pdb in pdb_list:
            for i in range(1):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)

if args.day == "jan31":
    if args.mode == 1:
        with open("../../cath_dataset_shuffle_optimization/optimization_iter0/protein_list") as f:
            a = f.readlines()
        b = [i.strip() for i in a]
        # print(b)
        for name in b:
            do(f"cp -r ../../cath_dataset/dompdb/{name}.pdb .")
    if args.mode == 2:
        name = f"protein_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0_phi_native_summary.txt"
        filtered_gamma = np.loadtxt(name)
        figureName = f"native_phi_change_limit"
        title = figureName
        show_together(filtered_gamma, figureName, title=title, inferBound=0, invert_sign=False)

        name = f"protein_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0_phi_decoy_summary.txt"
        filtered_gamma = np.loadtxt(name)
        figureName = f"decoy_phi_change_limit"
        title = figureName
        show_together(filtered_gamma, figureName, title=title, inferBound=0, invert_sign=False)
    if args.mode == 3:
        name = f"log_native.txt"
        filtered_gamma = np.loadtxt(name)
        figureName = f"native_phi_log"
        title = figureName
        show_together(filtered_gamma, figureName, title=title, inferBound=1, invert_sign=False)

        name = f"log_decoy.txt"
        filtered_gamma = np.loadtxt(name)
        figureName = f"decoy_phi_log"
        title = figureName
        show_together(filtered_gamma, figureName, title=title, inferBound=1, invert_sign=False)

if args.day == "jan30":
    if args.mode == 1:
        # Cornichon
        # start from native
        rotation_all = np.arange(0, 180, 30)
        run = "run5_constant_temp_300"
        x_shift_all = [50]
        subModeList = [4]
        for i in range(2):
            for subMode in subModeList:
                for x_shift in x_shift_all:
                    for degree in rotation_all:
                        cmd = f"python mm_run.py my_ABCDE_with_ligand/separate_4TH_centered_{x_shift}_{degree} -f forces_setup.py --subMode {subMode} --to {run}/separate_4TH_{x_shift}_{degree}_{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 300 --tempEnd 300 --platform CUDA --fasta my_ABCDE_with_ligand.fasta"
                        out = slurmRun(f"slurms/gpu_jan27_{run}_{x_shift}_{degree}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                        print(out)
    if args.mode == 2:
        # Cornichon
        # start from native
        rotation_all = np.arange(0, 180, 30)
        run = "run5_constant_temp_400"
        x_shift_all = [50]
        subModeList = [4]
        for i in range(2):
            for subMode in subModeList:
                for x_shift in x_shift_all:
                    for degree in rotation_all:
                        cmd = f"python mm_run.py my_ABCDE_with_ligand/separate_4TH_centered_{x_shift}_{degree} -f forces_setup.py --subMode {subMode} --to {run}/separate_4TH_{x_shift}_{degree}_{subMode}_{i} -s 2e7 --reportFrequency 4000 --tempStart 400 --tempEnd 400 --platform CUDA --fasta my_ABCDE_with_ligand.fasta"
                        out = slurmRun(f"slurms/gpu_jan27_{run}_{x_shift}_{degree}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                        print(out)
    if args.mode == 3:
        # Cornichon
        # start from native
        rotation_all = np.arange(0, 180, 30)
        run = "run5_constant_temp_350"
        x_shift_all = [50]
        subModeList = [4]
        for i in range(2):
            for subMode in subModeList:
                for x_shift in x_shift_all:
                    for degree in rotation_all:
                        cmd = f"python mm_run.py my_ABCDE_with_ligand/separate_4TH_centered_{x_shift}_{degree} -f forces_setup.py --subMode {subMode} --to {run}/separate_4TH_{x_shift}_{degree}_{subMode}_{i} -s 2e7 --reportFrequency 4000 --tempStart 350 --tempEnd 350 --platform CUDA --fasta my_ABCDE_with_ligand.fasta"
                        out = slurmRun(f"slurms/gpu_jan27_{run}_{x_shift}_{degree}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                        print(out)

if args.day == "jan27":
    if args.mode == 1:
        cmd = "python mm_run.py my_ABCDE_with_ligand/my_ABCDE_with_ligand --to run3/0 -s 8e6 --reportFrequency 4000 -f forces_setup.py --subMode 1 --fasta my_ABCDE_with_ligand.fasta --platform CUDA"
        out = slurmRun(f"slurms/gpu_jan27.slurm", cmd, template=gpu_commons_slurm)
    if args.mode == 2:
        # Cornichon
        # start from native
        rotation_all = np.arange(0, 180, 30)
        run = "run4"
        x_shift_all = [50]
        subModeList = [3, 4]
        subModeList = [4]
        for i in range(2):
            for subMode in subModeList:
                for x_shift in x_shift_all:
                    for degree in rotation_all:
                        cmd = f"python mm_run.py my_ABCDE_with_ligand/separate_4TH_centered_{x_shift}_{degree} -f forces_setup.py --subMode {subMode} --to {run}/separate_4TH_{x_shift}_{degree}_{subMode}_{i} -s 8e6 --reportFrequency 4000 --tempStart 600 --tempEnd 300 --platform CUDA --fasta my_ABCDE_with_ligand.fasta"
                        out = slurmRun(f"slurms/gpu_jan27_{run}_{x_shift}_{degree}_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                        print(out)
if args.day == "jan22":
    # Single Seq shuffle
    # data = pd.read_csv("data_info_3.csv", index_col=0).query("Problematic != 4")

    # data = pd.read_csv("pfam_selected.csv", index_col=0)
    # pdb_list = data["FullName"].to_list()

    data = pd.read_csv("has_structures_small_dataset_cleaned.csv")
    pdb_list = data["Protein"].to_list()

    if args.mode == 1:
        do("cp /home/wl45/opt/optimization/phi_list_contact.txt phi_list.txt")
        do("mkdir proteins_name_list")
        do("mkdir slurms")
        do('mkdir -p outs')
        # do("mkdir -p decoys/multiShuffle")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p ../phis")
        with open("protein_list", "w") as out:
            for pdb in pdb_list:
                out.write(f"{pdb}\n")
    if args.mode ==222:
        n = split_proteins_name_list(pdb_per_txt=5)
    if args.mode == 2:
        # split proteins_name_list
        n = split_proteins_name_list(pdb_per_txt=5)
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d jan22 -m 22 -l {proteins}", template=base_slurm, memory=10)
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
    if args.mode == 22:
        from pyCodeLib import *
        proteins = args.label
        n_decoys = 1000
        generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
        do(f"python3 ~/opt/compute_phis.py -m 0 {proteins}")
    if args.mode == 3:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d jan22 -m 33"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 33:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        # complete_proteins = "protein_list"
        complete_proteins = "protein_list"
        n_decoys = 1000
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='shuffle',
                                        withBiased=False, num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, )
        # next: optimization_analysze.py iter{i}
    if args.mode == 5:
        # i = 1
        # i = 2
        # i = 3
        i = 0
        do(f"optimization_analyze.py iter{i}")
    if args.mode == 6:
        # i = 1
        # i = 2
        # i = 3
        i = 0
        do(f"optimization_analyze.py iter{i} -c -87")

if args.day == "jan20":
    # generate MSA
    pdb_list = dataset["optimization_v2"]
    pdb_location = "cleaned_pdbs/first_test_set"
    if args.mode == 1:
        def pdbToFasta(pdb, pdbLocation, fastaFile, chains="A"):
            import textwrap
            from Bio.PDB.PDBParser import PDBParser
            from Bio.PDB.Polypeptide import three_to_one
            # pdb = "1r69"
            # pdbLocation = "/Users/weilu/Research/server/may_2019/family_fold/1r69.pdb"
            # chains = "A"
            # fastaFile = "/Users/weilu/Research/server/may_2019/family_fold/1r69.fasta"
            s = PDBParser(QUIET=True).get_structure("X", pdbLocation)
            # m = s[0]  # model 0
            seq = ""
            with open(fastaFile, "w") as out:
                chain = "A"
                out.write(f">{pdb.upper()}:{chain.upper()}|PDBID|CHAIN|SEQUENCE\n")
                seq = ""
                for res in s.get_residues():
                    resName = three_to_one(res.resname)
                    seq += resName
                out.write("\n".join(textwrap.wrap(seq, width=80))+"\n")
            return seq

        # a = glob.glob("cleaned_pdbs/*.pdb")
        jobIdList = []
        # a = glob.glob("../membrane_only_contact_optimization/database/dompdb/*.pdb")
        do("mkdir -p fasta")
        for pdb in pdb_list:
            print(pdb)
            # pdb = line.split("/")[-1].split(".")[0]
            pdbToFasta(pdb, f"{pdb_location}/{pdb}.pdb", f"fasta/{pdb}.fasta")
            cmd = f"/projects/pw8/wl45/hh-suite/build/bin/hhblits -i fasta/{pdb}.fasta -d ../../uniclust30_2018_08/uniclust30_2018_08 -o {pdb} -oa3m {pdb}.a3m -n 2"
            # do(cmd)
            jobId = slurmRun(f"slurms/{pdb}.slurm", cmd, memory=10)
            jobIdList.append(jobId)
    if args.mode == 2:
        do("mkdir MSA")
        for pdb in pdb_list:
            do(f"mv {pdb} MSA/")
            do(f"mv {pdb}.a3m MSA/")
    if args.mode == 3:
        do("mkdir -p alignments")
        for pdb in pdb_list:
            a = pdb
            print(a)
            data = get_MSA_data(f"MSA/{a}.a3m")
            with open(f"alignments/{a}_filtered_0.05.seqs", "w") as out:
                for line in data:
                    out.write(line+"\n")


if args.day == "jan19":
    pdb_list = dataset["optimization_v2"]
    # folder = "iter2_real_gpu"
    # folder = "iter3_gpu"
    # folder_list = [folder]
    folder_list = ["frag_mem_with_single_iter0_original", "frag_mem_with_single_multi_pfam_iter0_original"]
    run_n = 5
    base_path = "/scratch/wl45/jan_2020"
    simulationType = "iterative_optimization"
    if args.mode == 1:
        for folder in folder_list:
            cd(folder)
            for pdb in pdb_list:
                for i in range(run_n):
                    cd(f"{pdb}/{i}")
                    do("python /projects/pw8/wl45/openawsem/helperFunctions/convertOpenmmTrajectoryToStandardMovie.py movie.pdb")
                    cd("../..")
            cd("..")
    if args.mode == 11:
        do("mkdir -p database/dompdb")
        do("mkdir -p database/S20_seq")
        for pdb in pdb_list:
            do(f"cp ../../optimization_database/setups/{pdb}/cleaned_pdbs/{pdb}.pdb database/dompdb/")
    if args.mode == 12:
        def getSeq(fileLocation):
            p = PDBParser()
            s = p.get_structure("test", fileLocation)
            seq = ""
            residues = list(s.get_residues())
            for residue in residues:
                res_id = residue.get_id()[0]
                if res_id==' ':
                    residue_name = residue.get_resname()
                    seq += three_to_one(residue_name)
            return seq
        for pdb in pdb_list:
            fileLocation = f"database/dompdb/{pdb}.pdb"
            seq = getSeq(fileLocation)
            opt_pdbName = pdb
            fileLocation = f"database/S20_seq/{opt_pdbName}.seq"
            with open(fileLocation, "w") as out:
                out.write(seq+"\n")
    if args.mode == 13:

        outFile = read_simulation_info(folder_list, pdb_list, simulationType, base_path, run_n)
        today = datetime.today().strftime('%m-%d')
        outFile = f"{simulationType}.csv"
    if args.mode == 2:
        # change directory to optimization.
        # dataFile = "iterative_optimization_single_mem_iter7_01-13.csv"
        # dataFile = "iterative_optimization_frag_mem_iter1_01-14.csv"
        # dataFile = "iterative_optimization_frag_mem_iter2_01-15.csv"
        dataFile = f"{simulationType}.csv"
        simulation_folder = "iterative_optimization"
        data = pd.read_csv(dataFile, index_col=0)
        parser = PDBParser()
        # folder = "first"
        # for folder in ["first", "first_cpu2"]:
        for folder in folder_list:
            pre = f"{base_path}/{simulation_folder}/{folder}"
            to_folder = "."
            os.system(f"mkdir -p {to_folder}/decoys/openMM")
            # pdb_list = ['1j5u']
            for pdb in pdb_list:
                complete_models = []
                print(pdb)
                for i in range(run_n):
                    movieFile = f"{pre}/{pdb}/{i}/movie.pdb"
                    lastN_frame = 50
                    allFrames, n, size = getAllFrames(movieFile)
                    num_of_frames = int(n/size)
                    first_chosen_frame = num_of_frames - lastN_frame
                    # last_chosen_frame = num_of_frames
                    oneFrame = allFrames[size*first_chosen_frame:size*(num_of_frames)]

                    p = PDBParser()
                    f = io.StringIO("".join(oneFrame))
                    s = p.get_structure("test", f)
                    complete_models += list(s.get_models())
                t = data.query(f"Protein == '{pdb}' and Folder == '{folder}' and Steps > 1").reset_index(drop=True)
                t = t.groupby("Run").tail(50).reset_index(drop=True)
                t["structure"] = complete_models
                t = t.rename(columns={"Q":"Qw"})
                # last50 = t.groupby("Run").tail(50).reset_index(drop=True)
                last50 = t
                # print(last50.head())
                # print(last50.tail())
                # print(last50.shape)
                # print(len(complete_models))
                to_folder = "."
                last50.to_pickle(f"{to_folder}/decoys/openMM/{pdb}_{folder}.pkl")
    if args.mode == 3:
        do("mkdir -p proteins_name_list")
        do("mkdir -p slurms")
        do("mkdir -p gammas")
        do('mkdir -p outs')
        do("mkdir -p ../phis")
        # folder = "first"
        # folder_list = ["first", "first_cpu2"]
        # folder_list = [folder]
        with open("protein_list", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
        do("cp /home/wl45/opt/optimization/phi_list_contact.txt phi_list.txt")
    if args.mode == 4:
        do("gg_server.py -d jan12 -m 44")
    if args.mode == 44:
        # time.sleep(36000)
        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(1):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i)
        n = i
        i = 0
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d jan12 -m 444 -l {proteins}", template=base_slurm, memory=10)
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
    if args.mode == 444:
        protein = args.label
        do(f"python3 ~/opt/compute_phis.py -m 7 {protein}")
        # do("python3 ~/opt/compute_phis.py -m 0 test_protein")

    if args.mode == 55:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d jan12 -m 555"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 555:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list"
        n_decoys = 500
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        # simulation_location_list=["first", "first_cpu2"]
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='openMM',
                                        withBiased=True ,num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, )
        # next: optimization_analysze.py iter{i}
    if args.mode == 5:
        folder_list = ["iter0_gpu", "iter0_4e6_gpu", "iter0_8e6_gpu", "single_mem_iter7"]
        folder_list += ["frag_mem_iter1", "frag_mem_iter2"]
        # folder_list += ["iter3_gpu"]
        with open("protein_list", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
        do("gg_server.py -d dec14 -m 333")
    if args.mode == 6:
        # i = 1
        # i = 2
        i = 3
        do(f"optimization_analyze.py iter{i}")
    if args.mode == 7:
        # go to iteration folder
        # iteration = 3
        # iteration = 1
        # iteration = 2
        iteration = 3
        for pdb in pdb_list:
            do(f"cp ../optimization/optimization_iter{iteration}/Oct26_saved_gammas/iter{iteration}_cutoff500_impose_Aprime_constraint/gamma.dat setups/{pdb}/iter{iteration}_gamma.dat")
            do(f"cp ../optimization/optimization_iter{iteration}/Oct26_saved_gammas/iter{iteration}_cutoff500_impose_Aprime_constraint/burial_gamma.dat setups/{pdb}/iter{iteration}_burial_gamma.dat")
    if args.mode == 8:
        # modify force_steup.py
        force_setup = "forces_setup.py"
        # iteration = "frag_mem_iter1"
        # subMode = 3
        # iteration = "frag_mem_iter2"
        # subMode = 4
        iteration = "frag_mem_iter3"
        subMode = 5
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 4e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 9:
        # modify force_steup.py
        iteration = "native"
        force_setup = "forces_setup.py"
        subMode = 0
        for pdb in pdb_list:
            for i in range(1):
                cmd = f"python mm_run.py setups/{pdb}/{pdb} --to {iteration}/{pdb} -s 4e6 --tempStart 300 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)

if args.day == "jan18":
    pdb_list = dataset["optimization_v2"]
    if args.mode == 1:
        # modify force_steup.py
        force_setup = "forces_setup_single_and_frag.py"
        # iteration = "frag_mem_iter1"
        # subMode = 3
        # iteration = "frag_mem_iter2"
        # subMode = 4
        iteration = "frag_mem_with_single_iter0_original"
        subMode = 0
        for pdb in pdb_list:
            for i in range(5):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 2:
        # go to iteration folder

        iteration = 0
        cutoff = 600
        for pdb in pdb_list:
            do(f"cp ../multi_seq_Pfam/optimization_iter{iteration}/Oct26_saved_gammas/iter{iteration}_cutoff{cutoff}_impose_Aprime_constraint/gamma.dat setups/{pdb}/multi_pfam_iter{iteration}_gamma.dat")
            do(f"cp ../multi_seq_Pfam/optimization_iter{iteration}/Oct26_saved_gammas/iter{iteration}_cutoff{cutoff}_impose_Aprime_constraint/burial_gamma.dat setups/{pdb}/multi_pfam_iter{iteration}_burial_gamma.dat")
    if args.mode == 3:
        # modify force_steup.py
        force_setup = "forces_setup_single_and_frag.py"
        # iteration = "frag_mem_iter1"
        # subMode = 3
        # iteration = "frag_mem_iter2"
        # subMode = 4
        iteration = "frag_mem_with_single_multi_pfam_iter0_original"
        subMode = 1
        for pdb in pdb_list:
            for i in range(5):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)

if args.day == "jan16":
    # multi Seq shuffle
    # data = pd.read_csv("data_info_3.csv", index_col=0).query("Problematic != 4")
    data = pd.read_csv("pfam_selected.csv", index_col=0)
    pdb_list = data["FullName"].to_list()
    if args.mode == 1:
        do("cp /home/wl45/opt/optimization/phi_list_contact.txt phi_list.txt")
        do("mkdir proteins_name_list")
        do("mkdir slurms")
        do('mkdir -p outs')
        do("mkdir -p decoys/multiShuffle")
        do("mkdir -p gammas")
        do("mkdir -p ../phis")
        with open("protein_list", "w") as out:
            for pdb in pdb_list:
                out.write(f"{pdb}\n")
    if args.mode == 2:
        pdb = pdb_list[0]
        generate_multiShuffle(pdb, location="./", alignmentLocation="../", num_decoys=1000, nameMode=0)
    if args.mode == 3:
        # split proteins_name_list
        n = split_proteins_name_list(pdb_per_txt=5)
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d jan16 -m 33 -l {proteins}", template=base_slurm, memory=10)
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
    if args.mode == 33:
        from pyCodeLib import *
        protein = args.label
        pdb_list = read_column_from_file(protein, 1)
        for pdb in pdb_list:
            generate_multiShuffle(pdb, location="./", alignmentLocation="../", num_decoys=1000, nameMode=0)
        do(f"python3 ~/opt/compute_phis.py -m 4 {protein}")
    if args.mode == 4:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d jan16 -m 44"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 44:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list"
        n_decoys = 500
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='multiShuffle',
                                        withBiased=False, num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=True, )
        # next: optimization_analysze.py iter{i}
    if args.mode == 5:
        # i = 1
        # i = 2
        # i = 3
        i = 0
        do(f"optimization_analyze.py iter{i}")

if args.day == "jan14":
    if args.mode == 1:
        # Cornichon
        # start from native
        pdb = "my_ABCDF"
        rotation_all = np.arange(0, 180, 30)
        x_shift_all = [30]
        subModeList = [8]
        for i in range(2):
            for subMode in subModeList:
                for x_shift in x_shift_all:
                    for degree in rotation_all:
                        cmd = f"python mm_run.py setups/{pdb}/separate_4TH_{x_shift}_{degree} -f forces_setup_cut_out_LBD.py --subMode {subMode} --to run3/6ud8_my_ABCDF_separate_4TH_{x_shift}_{degree}_{subMode}_{i} -s 4e6 --reportFrequency 4000 --tempStart 600 --tempEnd 300 --platform CUDA --fasta crystal_structure.fasta"
                        out = slurmRun(f"slurms/gpu_separate_4TH_{x_shift}_{degree}_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)
                        print(out)
    if args.mode == 2:
        # Cornichon
        # 0.5 constraint, frag memeory, centered_
        pdb = "my_ABCDF"
        rotation_all = np.arange(0, 180, 30)
        x_shift_all = [50]
        subModeList = [8]
        for i in range(2):
            for subMode in subModeList:
                for x_shift in x_shift_all:
                    for degree in rotation_all:
                        cmd = f"python mm_run.py setups/{pdb}/separate_4TH_centered_{x_shift}_{degree} -f forces_setup_cut_out_LBD.py --subMode {subMode} --to run3/6ud8_my_ABCDF_separate_4TH_centered_{x_shift}_{degree}_{subMode}_{i} -s 4e6 --reportFrequency 4000 --tempStart 600 --tempEnd 300 --platform CUDA --fasta crystal_structure.fasta"
                        out = slurmRun(f"slurms/gpu_separate_4TH_centered_{x_shift}_{degree}_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)
                        print(out)

if args.day == "jan12":
    pdb_list = dataset["optimization_v2"]
    # folder = "iter2_real_gpu"
    # folder = "iter3_gpu"
    # folder_list = [folder]
    folder_list = ["iter0_gpu", "iter0_4e6_gpu", "iter0_8e6_gpu", "single_mem_iter7"]
    folder_list = ["frag_mem_iter1"]
    folder_list = ["frag_mem_iter2"]
    run_n = 10
    if args.mode == 1:
        for folder in folder_list:
            cd(folder)
            for pdb in pdb_list:
                for i in range(run_n):
                    cd(f"{pdb}/{i}")
                    do("python /projects/pw8/wl45/openawsem/helperFunctions/convertOpenmmTrajectoryToStandardMovie.py movie.pdb")
                    cd("../..")
            cd("..")
    if args.mode == 11:
        do("mkdir -p database/dompdb")
        do("mkdir -p database/S20_seq")
        for pdb in pdb_list:
            do(f"cp ../../optimization_database/setups/{pdb}/cleaned_pdbs/{pdb}.pdb database/dompdb/")
    if args.mode == 12:
        def getSeq(fileLocation):
            p = PDBParser()
            s = p.get_structure("test", fileLocation)
            seq = ""
            residues = list(s.get_residues())
            for residue in residues:
                res_id = residue.get_id()[0]
                if res_id==' ':
                    residue_name = residue.get_resname()
                    seq += three_to_one(residue_name)
            return seq
        for pdb in pdb_list:
            fileLocation = f"database/dompdb/{pdb}.pdb"
            seq = getSeq(fileLocation)
            opt_pdbName = pdb
            fileLocation = f"database/S20_seq/{opt_pdbName}.seq"
            with open(fileLocation, "w") as out:
                out.write(seq+"\n")
    if args.mode == 2:
        # change directory to optimization.
        # dataFile = "iterative_optimization_single_mem_iter7_01-13.csv"
        # dataFile = "iterative_optimization_frag_mem_iter1_01-14.csv"
        dataFile = "iterative_optimization_frag_mem_iter2_01-15.csv"
        simulation_folder = "iterative_optimization"
        data = pd.read_csv(dataFile, index_col=0)
        parser = PDBParser()
        # folder = "first"
        # for folder in ["first", "first_cpu2"]:
        for folder in folder_list:
            pre = f"/scratch/wl45/jan_2020/{simulation_folder}/{folder}"
            to_folder = "."
            os.system(f"mkdir -p {to_folder}/decoys/openMM")
            # pdb_list = ['1j5u']
            for pdb in pdb_list:
                complete_models = []
                print(pdb)
                for i in range(run_n):
                    movieFile = f"{pre}/{pdb}/{i}/movie.pdb"
                    lastN_frame = 50
                    allFrames, n, size = getAllFrames(movieFile)
                    num_of_frames = int(n/size)
                    first_chosen_frame = num_of_frames - lastN_frame
                    # last_chosen_frame = num_of_frames
                    oneFrame = allFrames[size*first_chosen_frame:size*(num_of_frames)]

                    p = PDBParser()
                    f = io.StringIO("".join(oneFrame))
                    s = p.get_structure("test", f)
                    complete_models += list(s.get_models())
                t = data.query(f"Protein == '{pdb}' and Folder == '{folder}' and Steps > 1").reset_index(drop=True)
                t = t.groupby("Run").tail(50).reset_index(drop=True)
                t["structure"] = complete_models
                t = t.rename(columns={"Q":"Qw"})
                # last50 = t.groupby("Run").tail(50).reset_index(drop=True)
                last50 = t
                # print(last50.head())
                # print(last50.tail())
                # print(last50.shape)
                # print(len(complete_models))
                to_folder = "."
                last50.to_pickle(f"{to_folder}/decoys/openMM/{pdb}_{folder}.pkl")
    if args.mode == 3:
        do("mkdir -p proteins_name_list")
        do("mkdir -p slurms")
        do("mkdir -p gammas")
        do('mkdir -p outs')
        do("mkdir -p ../phis")
        # folder = "first"
        # folder_list = ["first", "first_cpu2"]
        # folder_list = [folder]
        with open("protein_list", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
        do("cp /home/wl45/opt/optimization/phi_list_contact.txt phi_list.txt")
    if args.mode == 4:
        do("gg_server.py -d jan12 -m 44")
    if args.mode == 44:
        # time.sleep(36000)
        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(1):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i)
        n = i
        i = 0
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d jan12 -m 444 -l {proteins}", template=base_slurm, memory=10)
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
    if args.mode == 444:
        protein = args.label
        do(f"python3 ~/opt/compute_phis.py -m 7 {protein}")
        # do("python3 ~/opt/compute_phis.py -m 0 test_protein")

    if args.mode == 55:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d jan12 -m 555"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 555:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list"
        n_decoys = 500
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        # simulation_location_list=["first", "first_cpu2"]
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='openMM',
                                        withBiased=True ,num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, )
        # next: optimization_analysze.py iter{i}
    if args.mode == 5:
        folder_list = ["iter0_gpu", "iter0_4e6_gpu", "iter0_8e6_gpu", "single_mem_iter7"]
        folder_list += ["frag_mem_iter1", "frag_mem_iter2"]
        # folder_list += ["iter3_gpu"]
        with open("protein_list", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
        do("gg_server.py -d dec14 -m 333")
    if args.mode == 6:
        # i = 1
        # i = 2
        i = 3
        do(f"optimization_analyze.py iter{i}")
    if args.mode == 7:
        # go to iteration folder
        # iteration = 3
        # iteration = 1
        # iteration = 2
        iteration = 3
        for pdb in pdb_list:
            do(f"cp ../optimization/optimization_iter{iteration}/Oct26_saved_gammas/iter{iteration}_cutoff500_impose_Aprime_constraint/gamma.dat setups/{pdb}/iter{iteration}_gamma.dat")
            do(f"cp ../optimization/optimization_iter{iteration}/Oct26_saved_gammas/iter{iteration}_cutoff500_impose_Aprime_constraint/burial_gamma.dat setups/{pdb}/iter{iteration}_burial_gamma.dat")
    if args.mode == 8:
        # modify force_steup.py
        force_setup = "forces_setup.py"
        # iteration = "frag_mem_iter1"
        # subMode = 3
        # iteration = "frag_mem_iter2"
        # subMode = 4
        iteration = "frag_mem_iter3"
        subMode = 5
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 4e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 9:
        # modify force_steup.py
        iteration = "native"
        force_setup = "forces_setup.py"
        subMode = 0
        for pdb in pdb_list:
            for i in range(1):
                cmd = f"python mm_run.py setups/{pdb}/{pdb} --to {iteration}/{pdb} -s 4e6 --tempStart 300 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)

if args.day == "jan10":
    pdb_list = dataset["optimization_v2"]
    if args.mode == 1:
        # modify force_steup.py
        iteration = 0
        force_setup = "forces_setup.py"
        iteration = f"iter{iteration}_4e6_gpu"
        subMode = 0
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 4e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 2:
        # modify force_steup.py
        iteration = 0
        force_setup = "forces_setup.py"
        iteration = f"iter{iteration}_8e6_gpu"
        subMode = 0
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 3:
        # modify force_steup.py
        cmd = f"python mm_run.py 5a63 --to server_cpu1_step1e3 --platform CPU --thread 1 -s 1e3 --reportFrequency 10 --tempStart 300"
        # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
        out = slurmRun(f"slurms/cpu_1.slurm", cmd, template=scavenge_slurm, memory=15)
        print(out)
    if args.mode == 4:
        thread = 2
        subMode = 3
        i = 0
        cmd = f"python mm_run.py 5a63 --to server_cpu{thread}_{i}_step1e3 --subMode {subMode} --platform CPU --thread {thread} -s 1e3 --reportFrequency 10 --tempStart 300"
        # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
        out = slurmRun(f"slurms/cpu_{thread}_{subMode}_{i}.slurm", cmd, template=scavenge_slurm, memory=15, thread=thread)
        print(out)
    if args.mode == 5:
        subMode_list = [-1, 0, 1, 2, 3]
        thread_list = [1, 2, 4, 8]
        for i in range(3):
            for subMode in subMode_list:
                for thread in thread_list:
                    # modify force_steup.py
                    cmd = f"python mm_run.py 5a63 --to new_server_cpu{thread}_{subMode}_{i}_step1e3 --subMode {subMode} --platform CPU --thread {thread} -s 1e3 --reportFrequency 10 --tempStart 300"
                    # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                    out = slurmRun(f"slurms/cpu_{thread}_{subMode}_{i}.slurm", cmd, template=scavenge_slurm, memory=15, thread=thread)
                    print(out)
    if args.mode == 6:
        subMode_list = [-1, 0, 1, 2, 3]
        for i in range(3):
            for subMode in subMode_list:
                cmd = f"python mm_run.py 5a63 --to server_gpu_{subMode}_{i}_step1e4 --subMode {subMode} --platform CUDA -s 1e4 --reportFrequency 100 --tempStart 300"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/GPU_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 66:
        subMode_list = [4, 5, 6, 7, 8, 3]
        for i in range(2):
            for subMode in subMode_list:
                cmd = f"python mm_run.py 5a63 --to server_gpu_{subMode}_{i}_step1e5 --subMode {subMode} --platform CUDA -s 1e5 --reportFrequency 100 --tempStart 300"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/GPU_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 7:
        # go to iteration folder
        for pdb in pdb_list:
            do(f"cp saved_gammas/iter7_gamma.dat setups/{pdb}/single_mem_iter7_gamma.dat")
            do(f"cp saved_gammas/iter7_burial_gamma.dat setups/{pdb}/single_mem_iter7_burial_gamma.dat")
    if args.mode == 8:
        # modify force_steup.py
        iteration = "single_mem_iter7"
        force_setup = "forces_setup.py"
        subMode = 1
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 4e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)

if args.day == "jan05":
    pdb_list = dataset["optimization_v2"]
    if args.mode == 1:
        # modify force_steup.py
        iteration = 0
        force_setup = "forces_setup.py"
        iteration = f"iter{iteration}_gpu"
        subMode = 0
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 2e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
if args.day == "jan01":
    do("find . -maxdepth 1 -mindepth 1 -type d -exec tar cvf {}.tar {}  \;")
if args.day == "dec31":
    if args.mode == 1:
        # Cornichon
        # start from native
        pdb = "my_ABCDF"
        rotation_all = np.arange(0, 180, 30)
        x_shift_all = [30]
        subModeList = [7]
        for i in range(2):
            for subMode in subModeList:
                for x_shift in x_shift_all:
                    for degree in rotation_all:
                        cmd = f"python mm_run.py setups/{pdb}/separate_4TH_{x_shift}_{degree} -f forces_setup_cut_out_LBD.py --subMode {subMode} --to run3/6ud8_my_ABCDF_separate_4TH_{x_shift}_{degree}_{subMode}_{i} -s 4e6 --reportFrequency 4000 --tempStart 600 --tempEnd 300 --platform CUDA --fasta crystal_structure.fasta"
                        out = slurmRun(f"slurms/gpu_separate_4TH_{x_shift}_{degree}_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)
                        print(out)


if args.day == "dec29":
    pdb_list = dataset["optimization_cath"]
    # folder = "iter2_real_gpu"
    # folder = "iter3_gpu"
    # folder_list = [folder]
    folder_list = ["iter0_gpu"]
    run_n = 5
    if args.mode == 1:

        for folder in folder_list:
            cd(folder)
            for pdb in pdb_list:
                for i in range(run_n):
                    cd(f"{pdb}/{i}")
                    do("python /projects/pw8/wl45/openawsem/helperFunctions/convertOpenmmTrajectoryToStandardMovie.py movie.pdb")
                    cd("../..")
            cd("..")
    if args.mode == 11:
        do("mkdir -p database/dompdb")
        do("mkdir -p database/S20_seq")
        for pdb in pdb_list:
            do(f"cp ../../optimization_database/setups/{pdb}/cleaned_pdbs/{pdb}.pdb database/dompdb/")
    if args.mode == 12:
        def getSeq(fileLocation):
            p = PDBParser()
            s = p.get_structure("test", fileLocation)
            seq = ""
            residues = list(s.get_residues())
            for residue in residues:
                res_id = residue.get_id()[0]
                if res_id==' ':
                    residue_name = residue.get_resname()
                    seq += three_to_one(residue_name)
            return seq
        for pdb in pdb_list:
            fileLocation = f"database/dompdb/{pdb}.pdb"
            seq = getSeq(fileLocation)
            opt_pdbName = pdb
            fileLocation = f"database/S20_seq/{opt_pdbName}.seq"
            with open(fileLocation, "w") as out:
                out.write(seq+"\n")
    if args.mode == 2:
        # change directory to optimization.
        dataFile = "optimization_database_iter0_gpu_12-29.csv"
        simulation_folder = "optimization_database"
        data = pd.read_csv(dataFile, index_col=0)
        parser = PDBParser()
        # folder = "first"
        # for folder in ["first", "first_cpu2"]:
        for folder in folder_list:
            pre = f"/scratch/wl45/dec_2019/{simulation_folder}/{folder}"
            to_folder = "."
            os.system(f"mkdir -p {to_folder}/decoys/openMM")
            for pdb in pdb_list:
                complete_models = []
                print(pdb)
                for i in range(run_n):
                    movie = f"{pre}/{pdb}/{i}/movie.pdb"
                    s = parser.get_structure("X", movie)
                    complete_models += list(s.get_models())
                t = data.query(f"Protein == '{pdb}' and Folder == '{folder}' and Steps > 1").reset_index(drop=True)
                t["structure"] = complete_models
                t = t.rename(columns={"Q":"Qw"})
                last50 = t.groupby("Run").tail(50).reset_index(drop=True)
                to_folder = "."
                last50.to_pickle(f"{to_folder}/decoys/openMM/{pdb}_{folder}.pkl")
    if args.mode == 3:
        do("mkdir -p proteins_name_list")
        do("mkdir -p slurms")
        do("mkdir -p gammas")
        do('mkdir -p outs')
        do("mkdir -p ../phis")
        # folder = "first"
        # folder_list = ["first", "first_cpu2"]
        # folder_list = [folder]
        with open("protein_list", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
        do("cp /home/wl45/opt/optimization/phi_list_contact.txt phi_list.txt")
    if args.mode == 4:
        do("gg_server.py -d dec29 -m 44")
    if args.mode == 44:
        # time.sleep(36000)
        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(1):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i)
        n = i
        i = 0
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d dec29 -m 444 -l {proteins}", template=base_slurm, memory=10)
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
    if args.mode == 444:
        protein = args.label
        do(f"python3 ~/opt/compute_phis.py -m 7 {protein}")
        # do("python3 ~/opt/compute_phis.py -m 0 test_protein")

    if args.mode == 55:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d dec29 -m 555"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 555:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list"
        n_decoys = 500
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        # simulation_location_list=["first", "first_cpu2"]
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='openMM',
                                        withBiased=True ,num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, )
        # next: optimization_analysze.py iter{i}
    if args.mode == 5:
        folder_list = ["first_iter1_cpu4", "iter2_gpu", "first_gpu", "iter2_real_gpu"]
        folder_list += ["iter3_gpu"]
        with open("protein_list", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
        do("gg_server.py -d dec14 -m 333")
    if args.mode == 6:
        i = 4
        do(f"optimization_analyze.py iter{i}")
    if args.mode == 7:
        # go to iteration folder
        # iteration = 3
        iteration = 4
        for pdb in pdb_list:
            do(f"cp ../multiDensityOptimization/optimization_iteration1/optimization_iter{iteration}/Oct26_saved_gammas/iter{iteration}_cutoff500_impose_Aprime_constraint/gamma.dat setups/{pdb}/iter{iteration}_gamma.dat")
            do(f"cp ../multiDensityOptimization/optimization_iteration1/optimization_iter{iteration}/Oct26_saved_gammas/iter{iteration}_cutoff500_impose_Aprime_constraint/burial_gamma.dat setups/{pdb}/iter{iteration}_burial_gamma.dat")


if args.day == "dec28":
    pdb_list = dataset["optimization_cath"]
    if args.mode == 1:
        # modify force_steup.py
        iteration = 0
        force_setup = "forces_setup.py"
        iteration = f"iter{iteration}_gpu"
        subMode = 0
        for pdb in pdb_list:
            for i in range(5):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 2e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 2:
        # modify force_steup.py
        iteration = 0
        force_setup = "forces_setup.py"
        iteration = f"iter{iteration}_gpu_8millions_steps"
        subMode = 0
        for pdb in pdb_list:
            for i in range(5):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}_8millions.slurm", cmd, template=gpu_commons_slurm)
                print(out)
if args.day == "dec27":
    pdb_list = dataset["optimization_cath"]
    if args.mode == 1:
        def pdbToFasta(pdb, pdbLocation, fastaFile, chains="A"):
            import textwrap
            from Bio.PDB.PDBParser import PDBParser
            from Bio.PDB.Polypeptide import three_to_one
            # pdb = "1r69"
            # pdbLocation = "/Users/weilu/Research/server/may_2019/family_fold/1r69.pdb"
            # chains = "A"
            # fastaFile = "/Users/weilu/Research/server/may_2019/family_fold/1r69.fasta"
            s = PDBParser(QUIET=True).get_structure("X", pdbLocation)
            # m = s[0]  # model 0
            seq = ""
            with open(fastaFile, "w") as out:
                chain = "A"
                out.write(f">{pdb.upper()}:{chain.upper()}|PDBID|CHAIN|SEQUENCE\n")
                seq = ""
                for res in s.get_residues():
                    resName = three_to_one(res.resname)
                    seq += resName
                out.write("\n".join(textwrap.wrap(seq, width=80))+"\n")
            return seq

        # a = glob.glob("cleaned_pdbs/*.pdb")
        jobIdList = []
        # a = glob.glob("../membrane_only_contact_optimization/database/dompdb/*.pdb")
        do("mkdir -p fasta")
        for pdb in pdb_list:
            print(pdb)
            # pdb = line.split("/")[-1].split(".")[0]
            pdbToFasta(pdb, f"original_pdbs/{pdb}.pdb", f"fasta/{pdb}.fasta")
            cmd = f"/projects/pw8/wl45/hh-suite/build/bin/hhblits -i fasta/{pdb}.fasta -d ../../uniclust30_2018_08/uniclust30_2018_08 -o {pdb} -oa3m {pdb}.a3m -n 2"
            # do(cmd)
            jobId = slurmRun(f"slurms/{pdb}.slurm", cmd, memory=10)
            jobIdList.append(jobId)
    if args.mode == 2:
        for pdb in pdb_list:
            a = pdb
            print(a)
            data = get_MSA_data(f"MSA/{a}.a3m")
            with open(f"alignments/{a}_filtered_0.05.seqs", "w") as out:
                for line in data:
                    out.write(line+"\n")

    if args.mode == 3:
        # modify force_steup.py
        iteration = 0
        force_setup = "forces_setup.py"
        # iteration = f"iter{iteration}_gpu"
        iteration = "native"
        subMode = 0
        for pdb in pdb_list:
            for i in range(1):
                cmd = f"python mm_run.py setups/{pdb}/{pdb} --to {iteration}/{pdb}/{i} -s 2e6 --tempStart 300 --tempEnd 300 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)

    if args.mode == 4:
        # cornichon
        # start from native
        pdb = "my_ABCDF"
        i = 0
        # subModeList = [3, 1, 2]
        subModeList = [4, 5]
        for i in range(2):
            for subMode in subModeList:
                cmd = f"python mm_run.py setups/{pdb}/{pdb} -f forces_setup_cut_out_LBD.py --subMode {subMode} --to native/6ud8_my_ABCDF_{subMode}_{i}_no_beta_with_contact -s 2e6 --reportFrequency 4000 --tempStart 300 --tempEnd 300 --platform CUDA"
                out = slurmRun(f"slurms/gpu_native_my_ABCDF_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)
                print(out)


if args.day == "dec24":
    pdb_list = dataset["optimization"]
    # folder = "iter2_real_gpu"
    # folder = "iter3_gpu"
    # folder_list = [folder]
    # folder_list = ["iter0_gpu_less_beta", "iter4_gpu"]
    # folder_list = ["iter5_withBiased_gpu"]
    # folder_list = ["iter6_gpu"]
    folder_list = ["iter7_gpu"]
    if args.mode == 1:

        for folder in folder_list:
            cd(folder)
            for pdb in pdb_list:
                for i in range(10):
                    cd(f"{pdb}/{i}")
                    do("python /projects/pw8/wl45/openawsem/helperFunctions/convertOpenmmTrajectoryToStandardMovie.py movie.pdb")
                    cd("../..")
            cd("..")
    if args.mode == 2:
        # change directory to optimization.
        # dataFile = "iterative_optimization_first_cpu2_12-15.csv"
        # dataFile = "iterative_optimization_iter2_gpu_12-19.csv"
        # dataFile = "iterative_optimization_iter2_real_gpu_12-24.csv"
        # dataFile = "iterative_optimization_iter3_gpu_12-25.csv"
        # dataFile = "iterative_optimization_iter3_gpu_12-27.csv"
        # dataFile = "iterative_optimization_iter3_gpu_12-30.csv"
        # dataFile = "iterative_optimization_iter3_gpu_12-31.csv"
        dataFile = "iterative_optimization_iter3_gpu_01-02.csv"

        data = pd.read_csv(dataFile, index_col=0)
        parser = PDBParser()
        # folder = "first"
        # for folder in ["first", "first_cpu2"]:
        for folder in folder_list:
            pre = f"/scratch/wl45/dec_2019/iterative_optimization/{folder}"
            to_folder = "."
            os.system(f"mkdir -p {to_folder}/decoys/openMM")
            for pdb in pdb_list:
                complete_models = []
                print(pdb)
                for i in range(10):
                    movie = f"{pre}/{pdb}/{i}/movie.pdb"
                    s = parser.get_structure("X", movie)
                    complete_models += list(s.get_models())
                t = data.query(f"Protein == '{pdb}' and Folder == '{folder}' and Steps > 1").reset_index(drop=True)
                t["structure"] = complete_models
                t = t.rename(columns={"Q":"Qw"})
                last50 = t.groupby("Run").tail(50).reset_index(drop=True)
                to_folder = "."
                last50.to_pickle(f"{to_folder}/decoys/openMM/{pdb}_{folder}.pkl")
    if args.mode == 3:
        do("mkdir -p proteins_name_list")
        do("mkdir -p slurms")
        do("mkdir -p gammas")
        do('mkdir -p outs')
        do("mkdir -p ../phis")
        # folder = "first"
        # folder_list = ["first", "first_cpu2"]
        # folder_list = [folder]
        with open("protein_list", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
        do("cp /home/wl45/opt/optimization/phi_list_contact.txt phi_list.txt")
    if args.mode == 4:
        do("gg_server.py -d dec14 -m 222")
    if args.mode == 55:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d dec24 -m 555"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 555:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list"
        n_decoys = 500
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        # simulation_location_list=["first", "first_cpu2"]
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='openMM',
                                        withBiased=True ,num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, )
        # next: optimization_analysze.py iter{i}
    if args.mode == 5:
        folder_list = ["first_iter1_cpu4", "iter2_gpu", "first_gpu", "iter2_real_gpu"]
        folder_list += ["iter3_gpu"]
        folder_list += ["iter0_gpu_less_beta", "iter4_gpu"]
        folder_list += ["iter5_withBiased_gpu"]
        folder_list += ["iter6_gpu"]
        with open("protein_list", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
        do("gg_server.py -d dec24 -m 55")
    if args.mode == 6:
        i = 7
        do(f"optimization_analyze.py iter{i}")
    if args.mode == 7:
        # go to iteration folder
        # iteration = 3
        iteration = 7
        for pdb in pdb_list:
            do(f"cp ../multiDensityOptimization/optimization_iteration1/optimization_iter{iteration}/Oct26_saved_gammas/iter{iteration}_cutoff500_impose_Aprime_constraint/gamma.dat setups/{pdb}/iter{iteration}_gamma.dat")
            do(f"cp ../multiDensityOptimization/optimization_iteration1/optimization_iter{iteration}/Oct26_saved_gammas/iter{iteration}_cutoff500_impose_Aprime_constraint/burial_gamma.dat setups/{pdb}/iter{iteration}_burial_gamma.dat")

        # iteration = "5"
        # name = "5_withBiased"
        # for pdb in pdb_list:
        #     do(f"cp ../multiDensityOptimization/optimization_iteration1/optimization_iter{name}/Oct26_saved_gammas/iter{iteration}_cutoff500_impose_Aprime_constraint/gamma.dat setups/{pdb}/iter{name}_gamma.dat")
        #     do(f"cp ../multiDensityOptimization/optimization_iteration1/optimization_iter{name}/Oct26_saved_gammas/iter{iteration}_cutoff500_impose_Aprime_constraint/burial_gamma.dat setups/{pdb}/iter{name}_burial_gamma.dat")
    if args.mode == 8:
        # modify force_steup.py
        iteration = 7
        force_setup = "forces_setup.py"
        iteration = f"iter{iteration}_gpu"
        subMode = 9
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 2e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 888:
        # modify force_steup.py
        iteration = 7
        force_setup = "forces_setup.py"
        iteration = f"iter{iteration}_gpu_long"
        subMode = 9
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 88:
        # modify force_steup.py
        iteration = "5_withBiased"
        force_setup = "forces_setup.py"
        iteration = f"iter{iteration}_gpu"
        subMode = 7
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 2e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 9:
        # modify force_steup.py
        iteration = 0
        force_setup = "forces_setup.py"
        iteration = f"iter{iteration}_gpu_less_beta"
        subMode = 0
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 2e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)

if args.day == "dec23":
    pdb_list = dataset["optimization"]
    # folder_list = ["first_iter1_cpu4", "iter2_gpu", "first_gpu"]
    if args.mode == 1:
        for pdb in pdb_list:
            do(f"cp ../multiDensityOptimization/optimization_iteration1/optimization_iter2/Oct26_saved_gammas/iter2_cutoff500_impose_Aprime_constraint/gamma.dat setups/{pdb}/iter2_gamma.dat")
            do(f"cp ../multiDensityOptimization/optimization_iteration1/optimization_iter2/Oct26_saved_gammas/iter2_cutoff500_impose_Aprime_constraint/burial_gamma.dat setups/{pdb}/iter2_burial_gamma.dat")
    if args.mode == 2:
        pdb_list = dataset["optimization"]
        force_setup = "forces_setup.py"
        iteration = "iter2_real_gpu"
        subMode = 3
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 2e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 3:
        # start from native
        pdb = "cut_out_LBD"
        i = 0
        subModeList = [1, 2]
        for i in range(2):
            for subMode in subModeList:
                cmd = f"python mm_run.py setups/{pdb}/{pdb} -f forces_setup_cut_out_LBD.py --subMode {subMode} --to native/6ud8_cut_out_LBD_{subMode}_{i}_no_beta -s 2e6 --reportFrequency 4000 --tempStart 300 --tempEnd 300 --platform CUDA"
                out = slurmRun(f"slurms/gpu_native_LBD_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)
                print(out)

if args.day == "dec22":
    if args.mode == 1:
        # start from native
        pdb = "6ud8_ABCD_embeded"
        i = 0
        for subMode in [1, 2]:
            cmd = f"python mm_run.py setups/{pdb}/{pdb} -f forces_setup_ABCD.py --subMode {subMode} --to native/6ud8_ABCD_{subMode}_{i} -s 2e6 --reportFrequency 4000 --tempStart 300 --tempEnd 300 --platform CUDA"
            out = slurmRun(f"slurms/gpu_native_ABCD_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)
            print(out)
    if args.mode == 2:
        run = "run3"
        # subModeList = [5, 6]
        # subModeList = [7, 8, 9]
        subModeList = [5, 6, 10, 11]
        pdb = "6ud8_F_complete"
        for subMode in subModeList:
            for i in range(2):
                # i = 0
                # subMode = 2
                cmd = f"python mm_run_with_pulling_start.py setups/{pdb}/extended -f with_topology/forces_setup_6ud8_F_complete_3TH.py --subMode {subMode} --to {run}/3TH_{i}_subMode_{subMode} -s 5e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                out = slurmRun(f"slurms/gpu_3TH_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
                cmd = f"python mm_run_with_pulling_start.py setups/{pdb}/extended -f with_topology/forces_setup_6ud8_F_complete_4TH_hand_modify.py --subMode {subMode} --to {run}/4TH_{i}_subMode_{subMode} -s 5e6 --reportFrequency 4000  --tempStart 800 --tempEnd 300 --platform CUDA"
                out = slurmRun(f"slurms/gpu_4TH_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
            for i in range(2):
                # subMode = 2
                cmd = f"python mm_run.py setups/{pdb}/{pdb} -f with_topology/forces_setup_6ud8_F_complete_3TH.py --subMode {subMode} --to {run}/native_3TH_{i}_subMode_{subMode} -s 5e6 --reportFrequency 4000 --tempStart 300 --tempEnd 300 --platform CUDA"
                out = slurmRun(f"slurms/gpu_native_3TH_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
                cmd = f"python mm_run.py setups/{pdb}/{pdb} -f with_topology/forces_setup_6ud8_F_complete_4TH_hand_modify.py --subMode {subMode} --to {run}/native_4TH_{i}_subMode_{subMode} -s 5e6 --reportFrequency 4000 --tempStart 300 --tempEnd 300 --platform CUDA"
                out = slurmRun(f"slurms/gpu_native_4TH_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)
                print(out)

if args.day == "dec18":
    if args.mode == 1:
        # start from native
        do("mkdir -p outs")
        do("mkdir -p slurms")
        # i = 0
        # subMode = 3
        # subMode = 4
        subMode = 5
        pdb = "6ud8_F_complete"
        for i in range(3):
            # subMode = 2
            cmd = f"python mm_run.py setups/{pdb}/{pdb} -f with_topology/forces_setup_6ud8_F_complete_3TH.py --subMode {subMode} --to native/3TH_{i}_subMode_{subMode} -s 5e6 --reportFrequency 4000 --tempStart 300 --tempEnd 300 --platform CUDA"
            out = slurmRun(f"slurms/gpu_native_3TH_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)
            print(out)
            cmd = f"python mm_run.py setups/{pdb}/{pdb} -f with_topology/forces_setup_6ud8_F_complete_4TH.py --subMode {subMode} --to native/4TH_{i}_subMode_{subMode} -s 5e6 --reportFrequency 4000 --tempStart 300 --tempEnd 300 --platform CUDA"
            out = slurmRun(f"slurms/gpu_native_4TH_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)
            print(out)
    if args.mode == 2:
        # subMode = 4
        subMode = 6
        # subMode = 5
        run = "run1"
        subModeList = [5, 6]
        for subMode in subModeList:
            for i in range(2):
                # i = 0
                # subMode = 2
                pdb = "6ud8_F_complete"
                cmd = f"python mm_run_with_pulling_start.py setups/{pdb}/extended -f with_topology/forces_setup_6ud8_F_complete_3TH.py --subMode {subMode} --to {run}/3TH_{i}_subMode_{subMode} -s 5e6 --reportFrequency 4000 --tempStart 800 --tempEnd 300 --platform CUDA"
                out = slurmRun(f"slurms/gpu_3TH_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
                cmd = f"python mm_run_with_pulling_start.py setups/{pdb}/extended -f with_topology/forces_setup_6ud8_F_complete_4TH.py --subMode {subMode} --to {run}/4TH_{i}_subMode_{subMode} -s 5e6 --reportFrequency 4000  --tempStart 800 --tempEnd 300 --platform CUDA"
                out = slurmRun(f"slurms/gpu_4TH_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)
                print(out)

if args.day == "dec17":
    if args.mode == 1:
        do("mkdir -p outs")
        do("mkdir -p slurms")
        subMode = 4
        for i in range(3):
            cmd = f"python mm_run_with_pulling_start.py setups/6ud8_F/extended -f set_topology/forces_setup_6ud8_F_with_helical.py --subMode {subMode} --to run2/3TH_{subMode}_{i} -s 2e6 --platform CUDA"
            out = slurmRun(f"slurms/gpu_3TH_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
            print(out)
            cmd = f"python mm_run_with_pulling_start.py setups/6ud8_F/extended -f set_topology/forces_setup_6ud8_F_with_helical_4TH_topo.py --subMode {subMode} --to run2/4TH_{subMode}_{i} -s 2e6 --platform CUDA"
            out = slurmRun(f"slurms/gpu_4TH_{subMode}_{i}.slurm", cmd, template=gpu_commons_slurm)
            print(out)
    if args.mode == 2:
        subMode = 4
        for i in range(3):
            # i = 0
            cmd = f"python mm_run.py setups/6ud8_F/6ud8_F -f set_topology/forces_setup_6ud8_F_with_helical.py --subMode {subMode} --to native/3TH_{i}_subMode_{subMode} -s 5e6 --reportFrequency 4000 --tempStart 300 --tempEnd 300 --platform CUDA"
            out = slurmRun(f"slurms/gpu_native_3TH_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)
            print(out)
            cmd = f"python mm_run.py setups/6ud8_F/6ud8_F -f set_topology/forces_setup_6ud8_F_with_helical_4TH_topo.py --subMode {subMode} --to native/4TH_{i}_subMode_{subMode} -s 5e6 --reportFrequency 4000  --tempStart 300 --tempEnd 300 --platform CUDA"
            out = slurmRun(f"slurms/gpu_native_4TH_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)
            print(out)
if args.day == "dec16":
    if args.mode == 1:
        pdb = "6ud8_F"
        pdb = "6ud8_F_complete"
        do("mkdir -p TM_pred")
        cd("TM_pred")
        do(f"/projects/pw8/wl45/topology_prediction/PureseqTM_Package/PureseqTM_proteome.sh -i ../setups/{pdb}/{pdb}.fasta")
        cd("..")
    if args.mode == 2:
        pdb_list = dataset["optimization"]
        force_setup = "forces_setup.py"
        iteration = "first_gpu"
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 2e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 3:
        pdb_list = dataset["optimization"]
        force_setup = "forces_setup.py"
        iteration = "iter2_gpu"
        subMode = 2
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 2e6 --tempStart 800 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f {force_setup} --subMode {subMode}"
                # out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/gpu_{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)


if args.day == "dec15":
    # pdb_list = ['1e0m', '1w4e', '1e0g', '2wqg', '1jo8', '1fex', '2l6r', '1c8c', '1g6p', '1mjc', '2jmc', '1hdn', '1st7', '1n88', '1d6o', '1hcd', '2ga5', '1j5u', '3o4d', '1k0s']
    pdb_list = dataset["optimization"]
    folder_list = ["first_iter1_cpu4", "iter2_gpu", "first_gpu"]
    if args.mode == 1:
        for pdb in pdb_list:
            do(f"cp ../saved_gammas/Dec15_500/burial_gamma.dat setups/{pdb}/iter1_burial_gamma.dat")
            do(f"cp ../saved_gammas/Dec15_500/gamma.dat setups/{pdb}/iter1_gamma.dat")
    if args.mode == 2:
        force_setup = "forces_setup_iter1.py"
        iteration = "first_iter1_cpu4"
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 2e6 --tempStart 800 --tempEnd 200 --platform CPU --thread 4 --reportFrequency 4000 -f {force_setup}"
                out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                # out = slurmRun(f"slurms/{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 3:
        iteration = "first_iter1_cpu4"
        force_setup = "forces_setup_iter1.py"
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"mm_analysis.py setups/{pdb}/extended -t /scratch/wl45/dec_2019/iterative_optimization/{iteration}/{pdb}/{i}/movie.dcd --subMode -1 -f {force_setup}"
                out = slurmRun(f"slurms/analysis_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                # out = slurmRun(f"slurms/{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                print(out)
    if args.mode == 4:
        for folder in folder_list:
            cd(folder)
            for pdb in pdb_list:
                for i in range(10):
                    cd(f"{pdb}/{i}")
                    do("python /projects/pw8/wl45/openawsem/helperFunctions/convertOpenmmTrajectoryToStandardMovie.py movie.pdb")
                    cd("../..")
            cd("..")
if args.day == "dec14":
    # pdb_list = ['1e0m', '1w4e', '1e0g', '2wqg', '1jo8', '1fex', '2l6r', '1c8c', '1g6p', '1mjc', '2jmc', '1hdn', '1st7', '1n88', '1d6o', '1hcd', '2ga5', '1j5u', '3o4d', '1k0s']
    pdb_list = dataset["optimization"]
    folder_list = ["first_iter1_cpu4", "iter2_gpu", "first_gpu"]
    if args.mode == 1:
        # dataFile = "iterative_optimization_first_cpu2_12-15.csv"
        dataFile = "iterative_optimization_iter2_gpu_12-19.csv"
        data = pd.read_csv(dataFile, index_col=0)
        parser = PDBParser()
        # folder = "first"
        # for folder in ["first", "first_cpu2"]:
        for folder in folder_list:
            pre = f"/scratch/wl45/dec_2019/iterative_optimization/{folder}"
            to_folder = "."
            os.system(f"mkdir -p {to_folder}/decoys/openMM")
            for pdb in pdb_list:
                complete_models = []
                print(pdb)
                for i in range(10):
                    movie = f"{pre}/{pdb}/{i}/movie.pdb"
                    s = parser.get_structure("X", movie)
                    complete_models += list(s.get_models())
                t = data.query(f"Protein == '{pdb}' and Folder == '{folder}' and Steps > 1").reset_index(drop=True)
                t["structure"] = complete_models
                t = t.rename(columns={"Q":"Qw"})
                last50 = t.groupby("Run").tail(50).reset_index(drop=True)
                to_folder = "."
                last50.to_pickle(f"{to_folder}/decoys/openMM/{pdb}_{folder}.pkl")
    # if args.mode == 2:
    #     # folder = "first"
    #     for folder in folder_list:
    #         for pdb in pdb_list:
    #             # do(f"mv {folder}_{pdb} {folder}_{pdb}.pkl")
    #             do(f"mv {folder}_{pdb}.pkl {pdb}_{folder}.pkl")
    if args.mode == 3:
        # n = 10
        n = 500
        do("mkdir -p proteins_name_list")
        do("mkdir -p slurms")
        do("mkdir -p gammas")
        do('mkdir -p outs')
        do("mkdir -p ../phis")
        # folder = "first"
        # folder_list = ["first", "first_cpu2"]
        with open("protein_list", "w") as out:
            for folder in folder_list:
                for pdb in pdb_list:
                    out.write(f"{pdb}_{folder}\n")
    if args.mode == 4:
        do("mkdir -p database/dompdb")
        do("mkdir -p database/S20_seq")
        for pdb in pdb_list:
            do(f"cp ../../iterative_optimization/cleaned_pdbs/first_test_set/{pdb}.pdb database/dompdb/")
    if args.mode == 5:
        def getSeq(fileLocation):
            p = PDBParser()
            s = p.get_structure("test", fileLocation)
            seq = ""
            residues = list(s.get_residues())
            for residue in residues:
                res_id = residue.get_id()[0]
                if res_id==' ':
                    residue_name = residue.get_resname()
                    seq += three_to_one(residue_name)
            return seq
        for pdb in pdb_list:
            fileLocation = f"database/dompdb/{pdb}.pdb"
            seq = getSeq(fileLocation)
            opt_pdbName = pdb
            fileLocation = f"database/S20_seq/{opt_pdbName}.seq"
            with open(fileLocation, "w") as out:
                out.write(seq+"\n")
    if args.mode == 22:
        protein = args.label
        do(f"python3 ~/opt/compute_phis.py -m 7 {protein}")
        # do("python3 ~/opt/compute_phis.py -m 0 test_protein")
    if args.mode == 222:
        # time.sleep(36000)
        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(1):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i)
        n = i
        i = 0
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d dec14 -m 22 -l {proteins}", template=base_slurm, memory=10)
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
        exit()
        waitForJobs(jobIdList, sleepInterval=300)

        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d dec14 -m 33"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 333:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d dec14 -m 33"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 33:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list"
        n_decoys = 500
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        # simulation_location_list=["first", "first_cpu2"]
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='openMM',
                                        num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, )
        # next: optimization_analysze.py iter{i}
if args.day == "dec12":
    pdb_list = ["2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91"]
    pdb_list += ["4nv6", "4p79", "5dsg", "6g7o", "6a93", "2jo1", "1py6", "1pv6", "1u19"]
    if args.mode == 1:
        iteration = "third"
        # mode = 3  # hybrid
        # mode = 1  # water
        # mode =2   # membrane
        # for mode in [1, 2, 3]:
        # redo mode 1
        for mode in [1]:
            for i in range(0, 10):
                # batch prediction
                for pdb in pdb_list:
                    cmd = f"python mm_run_with_pulling_start.py setup/{pdb}/extended --to {iteration}_{mode}/{pdb}/{i} -s 8e6 --tempStart 1000 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f set_topology_force/forces_setup_{pdb}.py --subMode {mode}"
                    out = slurmRun(f"slurms/{iteration}_{mode}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                    print(out)
    if args.mode == 2:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        iteration = "first_gpu"
        d = pd.read_csv("original_pdbs/first_test_set.csv", index_col=0)
        pdb_list = d.PDB.str.lower().to_list()
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 2e6 --tempStart 800 --tempEnd 200 --reportFrequency 4000 -f forces_setup.py"
                # out = slurmRun(f"slurms/{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                out = slurmRun(f"slurms/{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 3:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        iteration = "first_gpu"
        d = pd.read_csv("original_pdbs/first_test_set.csv", index_col=0)
        pdb_list = d.PDB.str.lower().to_list()
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"mm_analysis.py setups/{pdb}/extended -t /scratch/wl45/dec_2019/iterative_optimization/first/{pdb}/{i}/movie.dcd --subMode -1 -f forces_setup.py"
                out = slurmRun(f"slurms/analysis_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                # out = slurmRun(f"slurms/{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                print(out)
    if args.mode == 4:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        iteration = "first_cpu2"
        d = pd.read_csv("original_pdbs/first_test_set.csv", index_col=0)
        pdb_list = d.PDB.str.lower().to_list()
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"python mm_run.py setups/{pdb}/extended --to {iteration}/{pdb}/{i} -s 2e6 --tempStart 800 --tempEnd 200 --platform CPU --thread 2 --reportFrequency 4000 -f forces_setup.py"
                out = slurmRun(f"slurms/cpu_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                # out = slurmRun(f"slurms/{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 5:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        iteration = "first_cpu2"
        d = pd.read_csv("original_pdbs/first_test_set.csv", index_col=0)
        pdb_list = d.PDB.str.lower().to_list()
        for pdb in pdb_list:
            for i in range(10):
                cmd = f"mm_analysis.py setups/{pdb}/extended -t /scratch/wl45/dec_2019/iterative_optimization/{iteration}/{pdb}/{i}/movie.dcd --subMode -1 -f forces_setup.py"
                out = slurmRun(f"slurms/analysis_{iteration}_{pdb}_{i}.slurm", cmd, template=scavenge_slurm)
                # out = slurmRun(f"slurms/{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                print(out)
if args.day == "dec10":
    if args.mode == 1:
        # get unfinished job.
        # pre = "/Users/weilu/Research/server/oct_2019/membrane_optimization"
        aa = glob.glob("outs/slurm-*.out")
        # print(aa)
        skip_pdbs = []
        for a in aa:
            cmd = f"grep 'CANCELLED' {a}"
            has_cancelled = getFromTerminal(cmd)
            if has_cancelled == "":
                pass
            else:
                print(has_cancelled)
                cmd = f"grep -A1 'phi_pairwise_contact_well' {a}"
                c = getFromTerminal(cmd)
                pdb = eval(c.split("\n")[1].strip())[0]
                skip_pdbs.append(pdb)
                print(pdb)
    if args.mode == 2:
        skip_pdbs = ['2FVMD_53-502', '5K1EA_21-567', '2W55H_135-705', '4EK7B_18-513', '6CTSA_44-423', '4S3RA_148-667', '5D51L_53-603', '1ULVA_281-956', '3VKHB_4007-4725', '2X51A_59-760', '5EEOA_157-822', '3RILD_78-427']
        with open("head500_protein_list") as f:
            a = f.readlines()
        a = [i.strip() for i in a]
        with open("filtered_protein_list", "w") as out:
            for pdb in a:
                if pdb in skip_pdbs:
                    pass
                else:
                    out.write(pdb+"\n")
        # print(aa)
        # skip_pdbs = []
        # for a in aa:
        #     cmd = f"grep -A1 'phi_pairwise_contact_well' {a}"
        #     c = getFromTerminal(cmd)
        #     if c == '':
        #         continue
        #     pdb = eval(c.split("\n")[1].strip())[0]
        #     skip_pdbs.append(pdb)
        #     # print(pdb)
        print(len(skip_pdbs))
        print(skip_pdbs)


if args.day == "dec09":
    # multi chain iter 0 globular.
    # based on jun10
    from pyCodeLib import *
    import warnings
    warnings.filterwarnings('ignore')
    n_decoys = 2000

    if args.mode == 44:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d dec09 -m 4"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 11:
        # n = 10
        n = 1000
        do("mkdir -p proteins_name_list")
        do("mkdir -p slurms")
        do("mkdir -p decoys/multiShuffle")
        do("mkdir -p gammas")
        do('mkdir -p outs')
        do("mkdir -p ../phis")
        with open("protein_list") as f:
            content = f.readlines()
        for line in content:
            name = line.strip()
            generate_multiShuffle(name, num_decoys=n, nameMode=0)
    if args.mode == 111:
        # filter out protein with small MSAs.
        with open("protein_list") as f:
            content = f.readlines()
        with open("protein_list_filtered", "w") as out:
            for line in content:
                name = line.strip()
                with open(f"alignments/{name}_filtered_0.05.seqs") as f:
                    a = f.readlines()
                if len(a) < 1000:
                    pass
                else:
                    out.write(line)
    if args.mode == 1111:
        # n = 10
        n = 1000
        do("mkdir -p proteins_name_list")
        do("mkdir -p slurms")
        do("mkdir -p decoys/multiShuffle")
        do("mkdir -p gammas")
        do('mkdir -p outs')
        do("mkdir -p ../phis")
        with open("protein_list") as f:
            content = f.readlines()
        for line in content:
            name = line.strip()
            generate_multiShuffle(name, num_decoys=n, nameMode=0)
    if args.mode == 1:
        # do("cp ../phi_list.txt .")
        # do("cp ~/opt/optimization/phi_list_contact.txt phi_list.txt")

        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(1):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i)
        n = i
        i = 0
        jobIdList = []
        # exit()
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d dec09 -m 2 -l {proteins}")
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
        waitForJobs(jobIdList, sleepInterval=300)
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d dec09 -m 4"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")

    if args.mode == 2:
        proteins = args.label
        # sampleK=2
        sampleK=1000
        evaluate_phis_over_training_set_for_native_structures_Wei(proteins, "phi_list.txt", decoy_method='multiShuffle', max_decoys=1e+10, tm_only=False, num_processors=1, multi_seq=True, sampleK=sampleK)

    if args.mode == 22:
        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        jobIdList = []
        for i in range(n):
            sub_list_file = f"proteins_name_list/proteins_name_list_{i}.txt"
            jobId = slurmRun(f"slurms/run_AB_{i}.slurm", f"python3 ~/opt/gg_server.py -d jun10 -m 3 -l {sub_list_file}", memory=5)
            jobIdList.append(jobId)
        # for testing toymodel
    if args.mode == 3:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # n = len(simulation_location_list)
        # n = 2
        # n = len(new_simulation_list)

        # import time
        # time.sleep(4000)
        # with open(f"proteins_name_list.txt", "w") as out:
        #     for data in ["new", "old"]:
        #         pdb_list, steps = dataset[data]
        #         for p in pdb_list:
        #             name = p.lower()[:4]
        #             out.write(f"{name}\n")
        # complete_proteins = "proteins_name_list/proteins_name_list_tiny.txt"
        # complete_proteins = "proteins_name_list.txt"
        # complete_proteins = "tiny_list.txt"
        complete_proteins = args.label
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='multiShuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=True)

        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        out = calculate_A_B_and_gamma_wl45_parallel(complete_proteins, "phi_list.txt", decoy_method='multiShuffle',
                                        num_decoys=1000, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=True)
        print("Success", out)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_xl23(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False, mode=1, simulation_location_list=simulation_location_list)
    # if args.mode == 4:
    #     with open(f"slurms/run_on_scavenge.slurm", "w") as out:
    #         out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr13 -m 3"))
    #     replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
    #     do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 4:
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list"
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='multiShuffle',
                                        num_decoys=1000, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=True)

if args.day == "dec07":
    pdb_list = ["2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91"]
    pdb_list += ["4nv6", "4p79", "5dsg", "6g7o", "6a93", "2jo1", "1py6", "1pv6", "1u19"]
    mode_dic = {"first_normalized":3, "second":1, "third":2}
    if args.mode == 1:
        iteration = "native_normalized"
        mode = 3  # hybrid
        # mode = 1  # water
        for mode in [1, 2, 3]:
            for i in range(1):
                # batch prediction
                for pdb in pdb_list:
                    cmd = f"python mm_run.py setup/{pdb}/{pdb} --to {iteration}_{mode}/{pdb}/{i} --platform CUDA -f set_topology_force/forces_setup_{pdb}.py --subMode {mode}  --reportFrequency 4000 -s 1e6 --tempStart 300 --tempEnd 300"
                    out = slurmRun(f"slurms/{iteration}_{mode}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                    print(out)
    if args.mode == 2:
        iteration = "second"
        # mode = 3  # hybrid
        # mode = 1  # water
        # mode =2   # membrane
        for mode in [1, 2, 3]:
            for i in range(0, 5):
                # batch prediction
                for pdb in pdb_list:
                    cmd = f"python mm_run_with_pulling_start.py setup/{pdb}/extended --to {iteration}_{mode}/{pdb}/{i} -s 8e6 --tempStart 1000 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f set_topology_force/forces_setup_{pdb}.py --subMode {mode}"
                    out = slurmRun(f"slurms/{iteration}_{mode}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                    print(out)
    if args.mode == 3:
        iteration = "third"
        # mode = 3  # hybrid
        # mode = 1  # water
        # mode =2   # membrane
        for mode in [1, 2, 3]:
            for i in range(0, 10):
                # batch prediction
                for pdb in pdb_list:
                    cmd = f"python mm_run_with_pulling_start.py setup/{pdb}/extended --to {iteration}_{mode}/{pdb}/{i} -s 8e6 --tempStart 1000 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f set_topology_force/forces_setup_{pdb}.py --subMode {mode}"
                    out = slurmRun(f"slurms/{iteration}_{mode}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                    print(out)

if args.day == "dec04":
    pdb_list = ["2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91"]
    pdb_list += ["4nv6", "4p79", "5dsg", "6g7o", "6a93", "2jo1", "1py6", "1pv6", "1u19"]

    if args.mode == 1:
        iteration = "first"
        mode = 3  # hybrid
        # mode = 1  # water
        # mode = 2  # membrane
        for i in range(1, 10):
            # batch prediction
            for pdb in pdb_list:
                cmd = f"python mm_run_with_pulling_start.py setup/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 1000 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f set_topology_force/forces_setup_{pdb}.py --subMode {mode}"
                out = slurmRun(f"slurms/{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 2:
        iteration = "first_beta_strength"
        mode = 3  # hybrid
        mode_dic = {"first_normalized":3, "second":1, "third":2}
        # mode = 1  # water
        for i in range(5, 20):
            # batch prediction
            for pdb in pdb_list:
                cmd = f"python mm_run_with_pulling_start.py setup/{pdb}/extended --to {iteration}/{pdb}/{i} -s 8e6 --tempStart 1000 --tempEnd 200 --platform CUDA --reportFrequency 4000 -f set_topology_force/forces_setup_{pdb}.py --subMode {mode}"
                out = slurmRun(f"slurms/{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
    if args.mode == 3:
        iteration = "native"
        mode = 3  # hybrid
        # mode = 1  # water
        for i in range(1):
            # batch prediction
            for pdb in pdb_list:
                cmd = f"python mm_run.py setup/{pdb}/{pdb} --to {iteration}/{pdb}/{i} --platform CUDA -f set_topology_force/forces_setup_{pdb}.py --subMode {mode}  --reportFrequency 4000 -s 1e6 --tempStart 300 --tempEnd 300"
                out = slurmRun(f"slurms/{iteration}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)
                print(out)
if args.day == "dec03":
    do("mkdir -p slurms")
    do("mkdir -p outs")
    jobIdList = []
    for i in range(10):
        for step in [1, 10, 100, 1000]:
            s = int(2e3*step)
            cmd = f"python mm_run.py setups/1r69/1r69 --to run1/step{step}/{i} --tempStart 300 --tempEnd 300 -s {s}  --reportFrequency {step} --platform CPU"
            # jobId = slurmRun(f"slurms/run_{step}_{i}.slurm", cmd, template=gpu_base_slurm)
            jobId = slurmRun(f"slurms/run_{step}_{i}.slurm", cmd, template=scavenge_slurm)

            # jobId = slurmRun(f"slurms/{simulation}_{pdb}_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)

            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
if args.day == "dec01":
    if args.mode == 1:
        # pdb_list = ["1fs3"]
        # pdb_list = ["1hn4"]
        pdb_list = ["1lmm", "1tcg", "1hn4", "1fs3", "1bpi"]
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        simulation = "run_without_er_gamma_yesCysCys"
        start_from ="extended"

        # force_file = "forces_setup_gamma_noCysCys.py"
        force_file = "forces_setup.py"
        subMode_list = []
        # subMode_list += [0, 1, 2, 3, 4, 5]
        # subMode_list += [10, 11, 12, 13, 14, 15]
        # subMode_list += [20, 21, 22, 23, 24, 25]
        # subMode_list += [30, 31, 32, 33, 34, 35]
        # subMode_list = [26]
        # subMode_list += [100, 101, 102, 103, 104, 105, 106]
        # subMode_list += [110, 111, 112, 113, 114, 115, 116]
        # subMode_list += [120, 121, 122, 123, 124, 125, 126]
        # subMode_list += [130, 131, 132, 133, 134, 135]
        # subMode_list += [0, 1, 2, 3]
        subMode_list += [30, 31, 32, 33]
        for i in range(20):
            for pdb in pdb_list:
                for subMode in subMode_list:
                    print(pdb)
                    cmd = f"python mm_run.py setups/{pdb}/{start_from} --to {simulation}/{pdb}/{subMode}_{i} -f {force_file} --subMode {subMode} --platform CUDA --reportFrequency 4000 -s 1e7 --tempStart 800 --tempEnd 200"
                    # jobId = slurmRun(f"slurms/run_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                    jobId = slurmRun(f"slurms/{simulation}_{pdb}_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)

                    # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
                    jobIdList.append(jobId)

if args.day == "nov26":
    if args.mode == 1:
        # pdb_list = ["1fs3"]
        # pdb_list = ["1hn4"]
        pdb_list = ["1lmm", "1tcg"]
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        simulation = "run_without_er_gamma_noCysCys_2"
        start_from ="extended"

        force_file = "forces_setup_gamma_noCysCys.py"
        subMode_list = []
        subMode_list += [0, 1, 2, 3, 4, 5]
        # subMode_list += [10, 11, 12, 13, 14, 15]
        # subMode_list += [20, 21, 22, 23, 24, 25]
        subMode_list += [30, 31, 32, 33, 34, 35]
        # subMode_list = [26]
        # subMode_list += [100, 101, 102, 103, 104, 105, 106]
        # subMode_list += [110, 111, 112, 113, 114, 115, 116]
        # subMode_list += [120, 121, 122, 123, 124, 125, 126]
        # subMode_list += [130, 131, 132, 133, 134, 135]
        for i in range(20):
            for pdb in pdb_list:
                for subMode in subMode_list:
                    print(pdb)
                    cmd = f"python mm_run.py setups/{pdb}/{start_from} --to {simulation}/{pdb}/{subMode}_{i} -f {force_file} --subMode {subMode} --platform CUDA --reportFrequency 4000 -s 1e7 --tempStart 800 --tempEnd 200"
                    # jobId = slurmRun(f"slurms/run_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                    jobId = slurmRun(f"slurms/{simulation}_{pdb}_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)

                    # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
                    jobIdList.append(jobId)
    if args.mode == 2:
        pdb_list = ["1bpi"]
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        simulation = "run_without_er_gamma_noCysCys_2"
        start_from ="extended"

        force_file = "forces_setup_gamma_noCysCys.py"
        subMode_list = []
        subMode_list += [0, 1, 2, 3, 4, 5]
        # subMode_list += [10, 11, 12, 13, 14, 15]
        # subMode_list += [20, 21, 22, 23, 24, 25]
        subMode_list += [30, 31, 32, 33, 34, 35]
        # subMode_list = [26]
        # subMode_list += [100, 101, 102, 103, 104, 105, 106]
        # subMode_list += [110, 111, 112, 113, 114, 115, 116]
        # subMode_list += [120, 121, 122, 123, 124, 125, 126]
        # subMode_list += [130, 131, 132, 133, 134, 135]
        for i in range(20):
            for pdb in pdb_list:
                for subMode in subMode_list:
                    print(pdb)
                    cmd = f"python mm_run.py setups/{pdb}/{start_from} --to {simulation}/{pdb}/{subMode}_{i} -f {force_file} --subMode {subMode} --platform CUDA --reportFrequency 4000 -s 1e7 --tempStart 800 --tempEnd 200"
                    # jobId = slurmRun(f"slurms/run_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                    jobId = slurmRun(f"slurms/{simulation}_{pdb}_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)

                    # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
                    jobIdList.append(jobId)

if args.day == "nov25":
    # native cys, 1fs3
    pdb_list = ["1fs3"]
    if args.mode == 1:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        # simulation = "run2_k10"
        # simulation = "run3"

        simulation = "run_with_ho_with_and_without_er"
        start_from ="extended"

        force_file = "forces_setup.py"
        subMode_list = []
        # subMode_list += [0, 1, 2, 3, 4, 5]
        # subMode_list += [10, 11, 12, 13, 14, 15]
        # subMode_list += [20, 21, 22, 23, 24, 25]
        subMode_list += [30, 31, 32, 33, 34, 35]
        # subMode_list = [26]
        # subMode_list += [100, 101, 102, 103, 104, 105, 106]
        # subMode_list += [110, 111, 112, 113, 114, 115, 116]
        # subMode_list += [120, 121, 122, 123, 124, 125, 126]
        subMode_list += [130, 131, 132, 133, 134, 135]
        for i in range(10):
            for pdb in pdb_list:
                for subMode in subMode_list:
                    print(pdb)
                    cmd = f"python mm_run.py setups/{pdb}/{start_from} --to {simulation}/{pdb}/{subMode}_{i} -f {force_file} --subMode {subMode} --platform CUDA --reportFrequency 4000 -s 1e7 --tempStart 800 --tempEnd 200"
                    # jobId = slurmRun(f"slurms/run_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                    jobId = slurmRun(f"slurms/{simulation}_{pdb}_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)

                    # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
                    jobIdList.append(jobId)
    if args.mode == 2:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        # simulation = "run2_k10"
        # simulation = "run3"

        simulation = "run_without_er_gamma_noCysCys"
        start_from ="extended"

        force_file = "forces_setup_gamma_noCysCys.py"
        subMode_list = []
        # subMode_list += [0, 1, 2, 3, 4, 5]
        # subMode_list += [10, 11, 12, 13, 14, 15]
        # subMode_list += [20, 21, 22, 23, 24, 25]
        subMode_list += [30, 31, 32, 33, 34, 35]
        # subMode_list = [26]
        # subMode_list += [100, 101, 102, 103, 104, 105, 106]
        # subMode_list += [110, 111, 112, 113, 114, 115, 116]
        # subMode_list += [120, 121, 122, 123, 124, 125, 126]
        # subMode_list += [130, 131, 132, 133, 134, 135]
        for i in range(10):
            for pdb in pdb_list:
                for subMode in subMode_list:
                    print(pdb)
                    cmd = f"python mm_run.py setups/{pdb}/{start_from} --to {simulation}/{pdb}/{subMode}_{i} -f {force_file} --subMode {subMode} --platform CUDA --reportFrequency 4000 -s 1e7 --tempStart 800 --tempEnd 200"
                    # jobId = slurmRun(f"slurms/run_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                    jobId = slurmRun(f"slurms/{simulation}_{pdb}_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)

                    # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
                    jobIdList.append(jobId)
if args.day == "nov23":
    # native cys, 1fs3
    pdb_list = ["1fs3"]
    if args.mode == 1:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        # simulation = "run2_k10"
        # simulation = "run3"
        # simulation = "native_energy"
        simulation = "constant_temp"
        # start_from ="extended"
        start_from = "1fs3"
        subMode_list = []
        # subMode_list += [0, 1, 2, 3, 4, 5]
        # subMode_list += [10, 11, 12, 13, 14, 15]
        subMode_list += [26]
        subMode_list += [20, 21, 22, 23, 24, 25]
        # subMode_list += [100, 101, 102, 103, 104, 105, 106]
        # subMode_list += [110, 111, 112, 113, 114, 115, 116]
        # subMode_list += [120, 121, 122, 123, 124, 125, 126]
        for i in range(5):
            for pdb in pdb_list:
                for subMode in subMode_list:
                    print(pdb)
                    cmd = f"python mm_run.py setups/{pdb}/{start_from} --to {simulation}/{pdb}/{subMode}_{i} -f forces_setup.py --subMode {subMode} --platform CUDA --reportFrequency 4000 -s 1e7 --tempStart 400 --tempEnd 200"
                    # jobId = slurmRun(f"slurms/run_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                    jobId = slurmRun(f"slurms/{simulation}_{pdb}_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)

                    # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
                    jobIdList.append(jobId)

if args.day == "aug11":
    if args.mode == 1:
        # pdb = "1pv6"
        do("mkdir -p forces_recompute")
        for pdb in dataset["hybrid"]:
            forceLocation = "forces_recompute"
            do(f"cp forces/forces_setup_{pdb}.py {forceLocation}/forces_setup_{pdb}.py")
            with fileinput.FileInput(f"{forceLocation}/forces_setup_{pdb}.py", inplace=True) as file:
                for line in file:
                    fromLine = 'MembranePart ='
                    toLine = 'MembranePart =[1,2]\nMembranePart ='
                    tmp = line.replace(fromLine, toLine)
                    print(tmp, end='')

if args.day == "nov20":
    # native cys, 1fs3
    pdb_list = ["1fs3"]
    if args.mode == 1:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        # simulation = "run2_k10"
        # simulation = "run3"

        # simulation = "run_with_er"
        # start_from ="extended"

        # subMode_list = []
        # subMode_list += [0, 1, 2, 3, 4, 5]
        # subMode_list += [10, 11, 12, 13, 14, 15]
        # subMode_list += [20, 21, 22, 23, 24, 25]
        # subMode_list = [26]
        # subMode_list += [100, 101, 102, 103, 104, 105, 106]
        # subMode_list += [110, 111, 112, 113, 114, 115, 116]
        # subMode_list += [120, 121, 122, 123, 124, 125, 126]

        simulation = "native_energy"
        start_from = "native"
        subMode_list = [27]

        for i in range(5):
            for pdb in pdb_list:
                for subMode in subMode_list:
                    print(pdb)
                    cmd = f"python mm_run.py setups/{pdb}/{start_from} --to {simulation}/{pdb}/{subMode}_{i} -f forces_setup.py --subMode {subMode} --platform CUDA --reportFrequency 4000 -s 1e7 --tempStart 800 --tempEnd 200"
                    # jobId = slurmRun(f"slurms/run_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                    jobId = slurmRun(f"slurms/{simulation}_{pdb}_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)

                    # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
                    jobIdList.append(jobId)
    if args.mode == 2:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        # simulation = "run2_k10"
        # simulation = "run3"

        simulation = "run_with_er_and_gamma_noCysCys"
        start_from ="extended"

        force_file = "forces_setup_gamma_noCysCys.py"
        subMode_list = []
        # subMode_list += [0, 1, 2, 3, 4, 5]
        # subMode_list += [10, 11, 12, 13, 14, 15]
        # subMode_list += [20, 21, 22, 23, 24, 25]
        # subMode_list = [26]
        subMode_list += [100, 101, 102, 103, 104, 105, 106]
        subMode_list += [110, 111, 112, 113, 114, 115, 116]
        subMode_list += [120, 121, 122, 123, 124, 125, 126]

        for i in range(5):
            for pdb in pdb_list:
                for subMode in subMode_list:
                    print(pdb)
                    cmd = f"python mm_run.py setups/{pdb}/{start_from} --to {simulation}/{pdb}/{subMode}_{i} -f {force_file} --subMode {subMode} --platform CUDA --reportFrequency 4000 -s 1e7 --tempStart 800 --tempEnd 200"
                    # jobId = slurmRun(f"slurms/run_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                    jobId = slurmRun(f"slurms/{simulation}_{pdb}_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)

                    # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
                    jobIdList.append(jobId)

if args.day == "nov15":
    pdb_list = ["1igd", "2sni", "1snb", "3il8", "1ubi", "1pht", "1poh", "1tig", "2acy", "1frd", "1opc", "1rds", "3chy", "5nul"]
    if args.mode == 1:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        # simulation = "run2_k10"
        # simulation = "new_beta_pap"
        # subMode = 0
        # simulation = "old_beta_pap"
        # subMode = 1
        simulation = "no_beta_pap"
        subMode = -1
        for i in range(5):
            for pdb in pdb_list:
                print(pdb)
                cmd = f"python mm_run.py setups/{pdb}/extended --to {simulation}/{pdb}/{i} -f forces_setup.py --subMode {subMode} --platform CUDA --reportFrequency 4000 -s 1e7 --tempStart 1000 --tempEnd 200"
                # jobId = slurmRun(f"slurms/run_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                jobId = slurmRun(f"slurms/{simulation}_{pdb}_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)

                # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
                jobIdList.append(jobId)

if args.day == "nov14":
    if args.mode ==1:
        pdb_list = ["1igd", "2sni", "1snb", "3il8", "1ubi", "1pht", "1poh", "1tig", "2acy", "1frd", "1opc", "1rds", "3chy", "5nul"]
        # pdb_list = ["6btc", "6g57", "6q64", "pex22", "sos", "unk2"]
        # pdb_list = ["fmt", "6n7n", "2is1"]
        # pdb_list = ["6jdq", "5uak"]
        for pdb in pdb_list:
            print(pdb)
            do(f"python mm_evaluate_native.py setups/{pdb}/{pdb} --to native/{pdb}_complete --platform CPU")
            pass

if args.day == "nov11":
    pdb_list = ["1fs3"]
    if args.mode == 1:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        # simulation = "run2_k10"
        # simulation = "run3"
        simulation = "run_weaker_beta_pap"
        # subMode_list = [0, 1, 2, 3, 4, 5]
        # subMode_list += [10, 11, 12, 13, 14, 15]
        subMode_list = [20, 21, 22, 23, 24, 25]
        for i in range(5):
            for pdb in pdb_list:
                for subMode in subMode_list:
                    print(pdb)
                    cmd = f"python mm_run.py setups/{pdb}/extended --to {simulation}/{pdb}/{subMode}_{i} -f forces_setup.py --subMode {subMode} --platform CUDA --reportFrequency 4000 -s 1e7 --tempStart 1000 --tempEnd 200"
                    # jobId = slurmRun(f"slurms/run_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                    jobId = slurmRun(f"slurms/{simulation}_{pdb}_{i}_{subMode}.slurm", cmd, template=gpu_commons_slurm)

                    # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
                    jobIdList.append(jobId)

if args.day == "nov10":
    pdb_list = ["1fs3"]
    if args.mode == 1:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        # simulation = "run2_k10"
        simulation = "run2"
        for i in range(5):
            for pdb in pdb_list:
                print(pdb)
                cmd = f"python mm_run.py setups/{pdb}/extended --to {simulation}/{pdb}/{i} -f forces_setup.py --platform CUDA --reportFrequency 4000 -s 1e7 --tempStart 1000 --tempEnd 200"
                # jobId = slurmRun(f"slurms/run_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                jobId = slurmRun(f"slurms/{simulation}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)

                # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
                jobIdList.append(jobId)
    if args.mode == 2:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        for i in range(5):
            for pdb in pdb_list:
                print(pdb)
                cmd = f"python mm_analysis.py setups/1fs3/1fs3 -t run1/1fs3/{i}/movie.dcd -f forces_setup.py --platform CUDA"
                # jobId = slurmRun(f"slurms/run_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                jobId = slurmRun(f"slurms/run1_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)

                # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
                jobIdList.append(jobId)
    if args.mode == 3:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        simulation = "run2_k10"
        # simulation = "run2"
        for i in range(5):
            for pdb in pdb_list:
                print(pdb)
                cmd = f"python mm_run.py setups/{pdb}/extended --to {simulation}/{pdb}/{i} -f forces_setup_k10.py --platform CUDA --reportFrequency 4000 -s 1e7 --tempStart 1000 --tempEnd 200"
                # jobId = slurmRun(f"slurms/run_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                jobId = slurmRun(f"slurms/{simulation}_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)

                # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
                jobIdList.append(jobId)
if args.day == "nov08":
    # pdb_list = ["1tt8", "6btc", "6g57", "6q64", "pex22", "sos", "unk2"]
    # pdb_list += ["fmt", "6n7n", "2is1"]
    pdb_list = ["6jdq", "5uak"]
    if args.mode == 1:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        for i in range(5):
            for pdb in pdb_list:
                print(pdb)
                if pdb in ["fmt", "6n7n", "2is1", "6jdq", "5uak"]:
                    subMode = 2
                elif pdb == "unk2":
                    subMode = 0
                else:
                    subMode = 1
                cmd = f"python mm_run.py setups/{pdb}/{pdb} --to run2/{pdb}/{i} --platform CUDA --subMode {subMode} --reportFrequency 4000 -s 1e6 --tempStart 800 --tempEnd 200"
                # jobId = slurmRun(f"slurms/run_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                jobId = slurmRun(f"slurms/run2_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)

                # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
                jobIdList.append(jobId)

if args.day == "nov06":
    pdb_list = ["fmt", "6n7n", "2is1"]
    if args.mode == 1:
        for pdb in pdb_list:
            do(f"cp ../frag_er_template/{pdb}_frag.mem setups/{pdb}/frag_HE.mem")

if args.day == "nov03":
    pdb_list = ["1tt8", "6btc", "6g57", "6q64", "pex22", "sos", "unk2"]
    pdb_list = ["fmt", "6n7n", "2is1"]
    if args.mode == 1:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        for i in range(3):
            for pdb in pdb_list:
                print(pdb)
                if pdb in ["fmt", "6n7n", "2is1"]:
                    subMode = 2
                elif pdb == "unk2":
                    subMode = 0
                else:
                    subMode = 1
                cmd = f"python mm_run.py setups/{pdb}/{pdb} --to run1/{pdb}/{i} --platform CUDA --subMode {subMode} --reportFrequency 4000 -s 1e7 --tempStart 1000 --tempEnd 200"
                # jobId = slurmRun(f"slurms/run_{pdb}_{i}.slurm", cmd, template=gpu_base_slurm)
                jobId = slurmRun(f"slurms/run_{pdb}_{i}.slurm", cmd, template=gpu_commons_slurm)

                # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
                jobIdList.append(jobId)
        # waitForJobs(jobIdList, sleepInterval=300)
if args.day == "oct31":
    pdb_list = ["1tt8", "6btc", "6g57", "6q64", "pex22", "sos", "unk2"]
    if args.mode == 1:
        # move frag, er, template file into right location.
        for pdb in pdb_list:
            do(f"cp ../frag_er_template/{pdb}_extended/run1/frag_HE.mem setups/{pdb}/")
            do(f"cp ../frag_er_template/{pdb}_extended/run1/go_rnativeC* setups/{pdb}/")
            do(f"cp ../frag_er_template/{pdb}_extended/run1/rnative.dat setups/{pdb}/")
            do(f"cp ../frag_er_template/{pdb}_extended/run1/ssweight setups/{pdb}/")
        pass
    if args.mode ==2:
        # pdb_list = ["6btc", "6g57", "6q64", "pex22", "sos", "unk2"]
        # pdb_list = ["fmt", "6n7n", "2is1"]
        pdb_list = ["6jdq", "5uak"]
        for pdb in pdb_list:
            print(pdb)
            if pdb in ["fmt", "6n7n", "2is1", "6jdq", "5uak"]:
                subMode = 2
            elif pdb == "unk2":
                subMode = 0
            else:
                subMode = 1
            do(f"python mm_evaluate_native.py setups/{pdb}/{pdb} --to native/{pdb}_complete --platform CPU --subMode {subMode}")
            pass

if args.day == "oct30":
    # from jun20
    # also compute rotation.
    # relative k computation, change is from change upper or lower to shift center.(which is the original idea, but instead of shift pdb, now shift the membrane.)
    # relative k now also compute the burial term and the membrane term.
    if args.mode == 1:
        do(f"mkdir -p proteins_name_list")
        # do(f"mkdir -p decoys/shifted")
        do(f"mkdir -p gammas")
        do(f"mkdir -p outs")
        do("mkdir -p ../phis")
        do("mkdir slurms")
        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(3):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i)
        n = i
        i = 0
        jobIdList = []
        # exit()
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d oct30 -m 2 -l {proteins}", template=base_slurm)
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
        waitForJobs(jobIdList, sleepInterval=300)
    if args.mode == 2:
        print(datetime.now())
        from pyCodeLib import *
        proteinlist = readList(args.label)
        print(proteinlist)
        for name in proteinlist:
            print(name)
            input_pdb_filename = f"../database/dompdb/{name}"
            structure = parse_pdb(input_pdb_filename)
            res_list = get_res_list(structure)
            neighbor_list = get_neighbor_list(structure)
            z_m = 0
            z_m_low = -15
            z_m_high = 15
            inside_or_not_table = []
            for res in res_list:
                z_loc = res["CA"].get_coord()[-1]
                if z_loc > z_m_low and z_loc < z_m_high:
                    inside_or_not_table.append(1)
                else:
                    inside_or_not_table.append(0)

            phi = phi_relative_k_with_membrane_well(res_list, neighbor_list, parameter_list="")
            np.savetxt(f"../phis/phi_relative_k_with_membrane_well_{name}_native_rotation", [phi], fmt='%1.4f')
            phis_list = []
            q_list = []
            # rotation_all = [30] * 12
            # z_shift_all = [-20] + [2]*20
            rotation_all = np.arange(0, 180, 15)
            # z_shift_all = np.arange(-20, 20, 4)
            z_shift_all = np.arange(-40, 40, 4)
            rotation_axis=Vector(1,0,0)
            for z_shift in z_shift_all:
                for degree in rotation_all:
                    structure = parse_pdb(input_pdb_filename)
                    # degree = 360
                    # degree = 0
                    radian=math.radians(degree)
                    # Iterate through all atoms and rotate by 90 degress
                    rotation_matrix = rotaxis2m(radian, rotation_axis)
                    translation=(0, 0, z_shift)
                    translation_matrix = np.array(translation, 'f')
                    # translation = np.array((0, 0, 100), 'f')
                    # print(translation_matrix)
                    for atom in structure.get_atoms():
                        # print(atom.get_coord())
                        atom.transform(rotation_matrix, translation_matrix)
                        # print(atom.get_coord())
                    count = 0
                    res_list = get_res_list(structure)
                    neighbor_list = get_neighbor_list(structure)
                    for i, res in enumerate(res_list):
                        z_loc = res["CA"].get_coord()[-1]
                        if z_loc > z_m_low and z_loc < z_m_high:
                            is_inside = 1
                        else:
                            is_inside = 0
                        if is_inside == inside_or_not_table[i]:
                            count += 1
                    # print(count, count/len(inside_or_not_table))
                    # io = PDBIO()
                    # io.set_structure(structure)
                    # io.save(f"p_{z_shift}_{degree}.pdb")
                    q = count/len(inside_or_not_table)

                    phi = phi_relative_k_with_membrane_well(res_list, neighbor_list, parameter_list="", z_m_high=z_m+15, z_m_low=z_m-15)
                    # q = interaction_well(z_m, -5, 5, 0.2)
                    # q = 0
                    print(phi, q)
                    phis_list.append(phi)
                    q_list.append(q)
            np.savetxt(f"../phis/phi_relative_k_with_membrane_well_{name}_decoys_rotation", phis_list, fmt='%1.4f')
            np.savetxt(f"../phis/phi_relative_k_with_membrane_well_{name}_decoysQ_rotation", q_list, fmt='%1.4f')
        print("done", args.label)
        print(datetime.datetime.now())
    if args.mode == 4:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list"
        with open("protein_list") as f:
            a = f.readlines()
        random.shuffle(a)
        with open("random_100", "w") as out:
            for line in a[:100]:
                out.write(line)
        complete_proteins = "top_500"
        complete_proteins = "last_500"
        complete_proteins = "tiny"
        complete_proteins = "random_100"
        print(datetime.datetime.now())
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='shifted',
                                        num_decoys=240, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, withBiased=True)
    if args.mode == 5:
        a = "rotation"
        b = "shifted_rotation"
        c = glob.glob("../back_phis/*")
        for line in c:
            if "native" in line:
                toLine = line.replace("back_phis", "phis")
            else:
                toLine = line.replace("back_phis", "phis").replace(a, b)
            # print(toLine)
            do(f"cp {line} {toLine}")
            # break
    # if args.mode == 2:
    #     from pyCodeLib import *
    #     proteinlist = readList(args.label)
    #     print(proteinlist)
    #     for name in proteinlist:
    #         print(name)
    #         input_pdb_filename = f"../database/dompdb/{name}"
    #         structure = parse_pdb(input_pdb_filename)
    #         res_list = get_res_list(structure)
    #         neighbor_list = get_neighbor_list(structure)
    #         z_m = 0
    #         z_m_low = -15
    #         z_m_high = 15
    #         inside_or_not_table = []
    #         for res in res_list:
    #             z_loc = res["CA"].get_coord()[-1]
    #             if z_loc > z_m_low and z_loc < z_m_high:
    #                 inside_or_not_table.append(1)
    #             else:
    #                 inside_or_not_table.append(0)

    #         phi = phi_relative_k_with_membrane_well(res_list, neighbor_list, parameter_list="")
    #         np.savetxt(f"../phis/phi_relative_k_with_membrane_well_{name}_native_rotation", [phi], fmt='%1.4f')
    #         phis_list = []
    #         q_list = []
    #         # rotation_all = [30] * 12
    #         # z_shift_all = [-20] + [2]*20
    #         rotation_all = [60] * 6
    #         z_shift_all = [-20] + [4]*10
    #         rotation_axis=Vector(1,0,0)
    #         z = 0
    #         rad = 0
    #         for z_shift in z_shift_all:
    #             z += z_shift
    #             translation=(0, 0, z_shift*10)
    #             translation_matrix = np.array(translation, 'f')
    #             rotation_matrix = np.identity(3)
    #             # print(rotation_matrix)
    #             for atom in structure.get_atoms():
    #                 # print(atom.get_coord())
    #                 atom.transform(rotation_matrix, translation_matrix)
    #             for degree in rotation_all:
    #                 # degree = 360
    #                 # degree = 0
    #                 rad += degree
    #                 radian=math.radians(degree)
    #                 # Iterate through all atoms and rotate by 90 degress
    #                 rotation_matrix = rotaxis2m(radian, rotation_axis)

    #                 translation_matrix = np.array((0, 0, 0), 'f')
    #                 # translation = np.array((0, 0, 100), 'f')
    #                 # print(translation_matrix)
    #                 for atom in structure.get_atoms():
    #                     # print(atom.get_coord())
    #                     atom.transform(rotation_matrix, translation_matrix)
    #                     # print(atom.get_coord())
    #                 count = 0
    #                 # res_list = get_res_list(structure)
    #                 for i, res in enumerate(res_list):
    #                     z_loc = res["CA"].get_coord()[-1]
    #                     if z_loc > z_m_low and z_loc < z_m_high:
    #                         is_inside = 1
    #                     else:
    #                         is_inside = 0
    #                     if is_inside == inside_or_not_table[i]:
    #                         count += 1
    #                 # print(count, count/len(inside_or_not_table))
    #                 io = PDBIO()
    #                 io.set_structure(structure)
    #                 io.save(f"p_{z}_{rad}.pdb")
    #                 q = count/len(inside_or_not_table)
    #                 res_list = get_res_list(structure)
    #                 neighbor_list = get_neighbor_list(structure)
    #                 phi = phi_relative_k_with_membrane_well(res_list, neighbor_list, parameter_list="", z_m_high=z_m+15, z_m_low=z_m-15)
    #                 # q = interaction_well(z_m, -5, 5, 0.2)
    #                 # q = 0
    #                 print(phi, q)
    #                 phis_list.append(phi)
    #                 q_list.append(q)
    #         np.savetxt(f"../phis/phi_relative_k_with_membrane_well_{name}_decoys_rotation", phis_list, fmt='%1.4f')
    #         np.savetxt(f"../phis/phi_relative_k_with_membrane_well_{name}_decoysQ_rotation", q_list, fmt='%1.4f')
    #     print("done", args.label)

if args.day == "oct22":
    # from jun20
    # relative k computation, change is from change upper or lower to shift center.(which is the original idea, but instead of shift pdb, now shift the membrane.)
    # relative k now also compute the burial term and the membrane term.
    if args.mode == 1:
        do(f"mkdir -p proteins_name_list")
        # do(f"mkdir -p decoys/shifted")
        do(f"mkdir -p gammas")
        do(f"mkdir -p outs")
        do("mkdir -p ../phis")
        do("mkdir slurms")
        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(3):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i)
        n = i
        i = 0
        jobIdList = []
        # exit()
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d oct22 -m 2 -l {proteins}")
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
        waitForJobs(jobIdList, sleepInterval=300)
    if args.mode == 2:
        from pyCodeLib import *
        proteinlist = readList(args.label)
        print(proteinlist)
        limit = 30
        for name in proteinlist:
            print(name)
            input_pdb_filename = f"../database/dompdb/{name}"
            structure = parse_pdb(input_pdb_filename)
            res_list = get_res_list(structure)
            neighbor_list = get_neighbor_list(structure)

            phi = phi_relative_k_with_membrane_well(res_list, neighbor_list, parameter_list="")
            np.savetxt(f"../phis/phi_relative_k_with_membrane_well_{name}_native_4.5_6.5_5.0_10", [phi], fmt='%1.4f')
            phis_list = []
            q_list = []
            z_m_all = np.arange(-limit, limit, 1)
            for z_m in z_m_all:
                phi = phi_relative_k_with_membrane_well(res_list, neighbor_list, parameter_list="", z_m_high=z_m+15, z_m_low=z_m-15)
                q = interaction_well(z_m, -5, 5, 0.2)
                # q = 0
                print(phi, q)
                phis_list.append(phi)
                q_list.append(q)
            np.savetxt(f"../phis/phi_relative_k_with_membrane_well_{name}_decoys_shifted_4.5_6.5_5.0_10", phis_list, fmt='%1.4f')
            np.savetxt(f"../phis/phi_relative_k_with_membrane_well_{name}_decoysQ_shifted_4.5_6.5_5.0_10", q_list, fmt='%1.4f')
        print("done", args.label)
    if args.mode == 4:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list"
        print(datetime.datetime.now())
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='shifted',
                                        num_decoys=240, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False)


if args.day == "oct20":
    # from jun20
    # relative k computation, change is from change upper or lower to shift center.(which is the original idea, but instead of shift pdb, now shift the membrane.)
    if args.mode == 1:
        do(f"mkdir -p proteins_name_list")
        # do(f"mkdir -p decoys/shifted")
        do(f"mkdir -p gammas")
        do(f"mkdir -p outs")
        do("mkdir -p ../phis")
        do("mkdir slurms")
        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(3):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i)
        n = i
        i = 0
        jobIdList = []
        # exit()
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d oct20 -m 2 -l {proteins}")
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
        waitForJobs(jobIdList, sleepInterval=300)
    if args.mode == 2:
        from pyCodeLib import *
        proteinlist = readList(args.label)
        print(proteinlist)
        limit = 30
        for name in proteinlist:
            input_pdb_filename = f"../database/dompdb/{name}"
            structure = parse_pdb(input_pdb_filename)
            res_list = get_res_list(structure)
            neighbor_list = get_neighbor_list(structure)

            phi = phi_relative_k_well(res_list, neighbor_list, parameter_list="")
            np.savetxt(f"../phis/phi_relative_k_well_{name}_native_4.5_6.5_5.0_10", [phi], fmt='%1.4f')
            phis_list = []
            q_list = []
            z_m_all = np.arange(5, limit, 1)
            for z_m in z_m_all:
                phi = phi_relative_k_well(res_list, neighbor_list, parameter_list="", z_m_high=z_m+15, z_m_low=z_m-15)
                q = interaction_well(z_m, -5, 5, 0.2)
                # q = 0
                print(phi, q)
                phis_list.append(phi)
                q_list.append(q)
            # np.savetxt(f"../phis/{name}_phis_high", phis_list, fmt='%1.4f')
            z_m_all = np.arange(-limit, -5, 1)
            for z_m in z_m_all:
                phi = phi_relative_k_well(res_list, neighbor_list, parameter_list="", z_m_high=z_m+15, z_m_low=z_m-15)
                q = interaction_well(z_m, -5, 5, 0.2)
                # q = 0
                print(phi, q)
                phis_list.append(phi)
                q_list.append(q)
            np.savetxt(f"../phis/phi_relative_k_well_{name}_decoys_shifted_4.5_6.5_5.0_10", phis_list, fmt='%1.4f')
            np.savetxt(f"../phis/phi_relative_k_well_{name}_decoysQ_shifted_4.5_6.5_5.0_10", q_list, fmt='%1.4f')
        print("done", args.label)
    if args.mode == 4:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list"
        print(datetime.datetime.now())
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='shifted',
                                        num_decoys=50, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False)


if args.day == "oct17":
    if args.mode == 1:
        # get DMP prediction
        # name = "beta_2_adrenergic_receptor"
        for name in dataset["membrane"]:
            folder = f"DMP/{name}"
            do(f"mkdir -p {folder}")
            cd(folder)
            do(f"cp ../../setups/{name}/{name}.fasta .")
            do(f"gg_server.py -d dmp -l {name}.fasta")
            cd("../..")

            # get TM prediction
            folder = f"TM_pred/{name}"
            do(f"mkdir -p {folder}")
            cd(folder)
            do(f"cp ../../setups/{name}/{name}.fasta .")
            do(f"/projects/pw8/wl45/topology_prediction/PureseqTM_Package/PureseqTM_proteome.sh -i {name}.fasta")
            cd("../..")

            # get secondary structure prediction
            folder = f"secondary/{name}"
            do(f"mkdir -p {folder}")
            cd(folder)
            do(f"cp ../../setups/{name}/{name}.fasta .")
            do(f"/projects/pw8/wl45/build/Predict_Property/Predict_Property.sh -i {name}.fasta")
            cd("../..")
    if args.mode == 2:
        # convert Porter5 format
        # pdb = "cannabinoid_receptor"
        # name = "serotonin_1A_receptor"
        # name = "cannabinoid_receptor"

        for name in dataset["membrane"]:
            toPre = f"setups/{name}"
            # from_secondary = f"secondary/{name}/{name}_PROP/{name}.ss3"
            # to_ssweight = f"{toPre}/ssweight"
            # print("convert ssweight")

            # data = pd.read_csv(from_secondary, comment="#", names=["i", "Res", "ss3", "Helix", "Sheet", "Coil"], sep="\s+")
            # # print(data)
            # with open(to_ssweight, "w") as out:
            #     for i, line in data.iterrows():
            #         if line["ss3"] == "H":
            #             out.write("1.0 0.0\n")
            #         if line["ss3"] == "E":
            #             out.write("0.0 1.0\n")
            #         if line["ss3"] == "C":
            #             out.write("0.0 0.0\n")

            # fasta_file = f"../original_fasta_files/{name}.fasta"
            # DMP_file = f"DMP/{name}/{name}.deepmetapsicov.con"
            # convertDMPToInput(name, DMP_file, fasta_file)
            # do(f"cp ~/opt/gremlin/protein/{name}/DMP/go_rnativeC* {name}/setup/")

            # name = "serotonin_1A_receptor"
            print("convert predictedZim")
            topo_name = f"TM_pred/{name}/{name}_topo"
            get_PredictedZim(topo_name, f"{toPre}/PredictedZim")
            get_PredictedZimSide(topo_name, f"{toPre}/PredictedZimSide")
    if args.mode ==3:
        cmd = "python mm_run_with_pulling_start.py setups/2bg9/2bg9 --to test1 --subMode 0 --platform CPU"

if args.day == "oct16":
    # from apr31
    # membrane contact optimization.
    from pyCodeLib import *
    import warnings
    warnings.filterwarnings('ignore')
    # n_decoys = 1000
    n_decoys = 2000
    # n_decoys = 4000
    # n_decoys = 8000
    separateDecoysNum = -1
    proteins = "protein_list"
    if args.mode == 1:
        do("mkdir proteins_name_list")
        do("mkdir slurms")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p outs")
        do("mkdir -p ../phis")
        proteins = f"protein_list"
        # generate_decoy_sequences(proteins, separateDecoysNum=separateDecoysNum, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
    if args.mode == 2:
        # time.sleep(16000)
        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(1):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i)
        n = i
        i = 0
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d oct16 -m 222 -l {proteins}")
            # jobId = slurmRun(f"slurms/compute_z_{i}.slurm", f"python3 compute_z_score.py {proteins}")
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
        # exit()
        waitForJobs(jobIdList, sleepInterval=300)
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d oct16 -m 4"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    # if args.mode == 22:
    #     protein = args.label
    #     do(f"python3 ~/opt/compute_phis.py -m 0 {protein}")
    #     # do("python3 ~/opt/compute_phis.py -m 0 test_protein")
    if args.mode == 222:
        proteins = args.label
        generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
        evaluate_phis_over_training_set_for_native_structures_Wei(proteins, "phi_list.txt",
                    decoy_method='shuffle', max_decoys=1e+10, tm_only=False, num_processors=1, separateDecoysNum=separateDecoysNum)
        print(proteins, "done")
    if args.mode == 3:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d oct16 -m 4"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 4:
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list"
        print(datetime.datetime.now())
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='shuffle',
                                        num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False)
    if args.mode == 5:
        # time.sleep(16000)
        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(1):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i)
        n = i
        i = 0
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d oct16 -m 222 -l {proteins}")
            jobId = slurmRun(f"slurms/compute_z_{i}.slurm", f"python3 compute_z_score.py {proteins}")
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            # do(f"cat {proteins} >> iter0_complete.txt")


if args.day == "oct09":
    if args.mode == 1:
        # get DMP prediction
        name = "beta_2_adrenergic_receptor"
        folder = f"DMP/{name}"
        do(f"mkdir -p {folder}")
        cd(folder)
        do(f"cp ../../../original_fasta_files/{name}.fasta .")
        do(f"gg_server.py -d dmp -l {name}.fasta")
        cd("../..")

        # get TM prediction
        folder = f"TM_pred/{name}"
        do(f"mkdir -p {folder}")
        cd(folder)
        do(f"cp ../../../original_fasta_files/{name}.fasta .")
        do(f"/projects/pw8/wl45/topology_prediction/PureseqTM_Package/PureseqTM_proteome.sh -i {name}.fasta")
        cd("../..")

        # get secondary structure prediction
        folder = f"secondary/{name}"
        do(f"mkdir -p {folder}")
        cd(folder)
        do(f"cp ../../../original_fasta_files/{name}.fasta .")
        do(f"/projects/pw8/wl45/build/Predict_Property/Predict_Property.sh -i {name}.fasta")
        cd("../..")

if args.day == "oct05":
    if args.mode == 1:
        # get DMP prediction
        name = "serotonin_1A_receptor"
        folder = f"DMP/{name}"
        do(f"mkdir -p {folder}")
        cd(folder)
        do(f"cp ../../../original_fasta_files/{name}.fasta .")
        do(f"gg_server.py -d dmp -l {name}.fasta")
        cd("../..")
    if args.mode == 2:
        # get TM prediction
        name = "serotonin_1A_receptor"
        folder = f"TM_pred/{name}"
        do(f"mkdir -p {folder}")
        cd(folder)
        do(f"cp ../../../original_fasta_files/{name}.fasta .")
        do(f"/projects/pw8/wl45/topology_prediction/PureseqTM_Package/PureseqTM_proteome.sh -i {name}.fasta")
        cd("../..")
    if args.mode == 3:
        # get secondary structure prediction
        # name = "serotonin_1A_receptor"
        name = "beta_2_adrenergic_receptor"
        # name = "cannabinoid_receptor"
        folder = f"secondary/{name}"
        do(f"mkdir -p {folder}")
        cd(folder)
        do(f"cp ../../../original_fasta_files/{name}.fasta .")
        do(f"/projects/pw8/wl45/build/Predict_Property/Predict_Property.sh -i {name}.fasta")
# if args.day == "oct02":
#     if args.mode == 1:
#         # get DMP prediction
#         do("gg_server.py -d dmp -l beta_2_adrenergic_receptor.fasta")
if args.day == "sep25":
    if args.mode == 1:
        do("gg_server.py -d sep02 -m 1 -l cannabinoid_receptor.fasta")
    if args.mode == 2:
        # topo_name = "cannabinoid_receptor_topo"
        # get_PredictedZim(topo_name, "PredictedZim")
        name = "serotonin_1A_receptor"
        topo_name = f"TM_pred/{name}/{name}_topo"
        get_PredictedZim(topo_name, f"{name}/setup/PredictedZim")
        get_PredictedZimSide(topo_name, f"{name}/setup/PredictedZimSide")
if args.day == "sep19":
    # shuffle iter0.
    # ligand optimization.
    from pyCodeLib import *
    import warnings
    warnings.filterwarnings('ignore')
    n_decoys = 100
    separateDecoysNum = -1

    if args.mode == 44:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d jun10 -m 4"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 1:
        do("mkdir -p proteins_name_list")
        do("mkdir -p slurms")
        # do("mkdir -p data")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p outs")
        do("mkdir -p ../phis")
        proteins = f"protein_list"
        # generate_decoy_sequences(proteins, separateDecoysNum=separateDecoysNum, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
    if args.mode == 2:
        # time.sleep(36000)
        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(1):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i)
        n = i
        i = 0
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d sep19 -m 222 -l {proteins}", template=base_slurm)
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
        exit()
        waitForJobs(jobIdList, sleepInterval=300)
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d sep19 -m 4"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 22:
        protein = args.label
        do(f"python3 ~/opt/compute_phis.py -m 0 {protein}")
        # do("python3 ~/opt/compute_phis.py -m 0 test_protein")
    if args.mode == 222:
        proteins = args.label
        evaluate_phis_over_training_set_for_native_structures_Wei(proteins, "phi_list.txt",
                    decoy_method='shuffle', max_decoys=1e+10, tm_only=False, num_processors=1, separateDecoysNum=separateDecoysNum)
    if args.mode == 3:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d sep19 -m 4"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 4:
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list"
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='shuffle',
                                        num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False)



if args.day == "aug20":
    if args.mode == 1:
        # pdb = args.label
        if args.label[-6:] == ".fasta":
            pdb = args.label[:-6]
        else:
            pdb = args.label
        do("mkdir -p TM_pred")
        cd("TM_pred")
        os.system(f"cp ../{pdb}.fasta .")
        # pdb = "4rws"
        # pdb = args.label
        # do(f"/projects/pw8/wl45/topology_prediction/PureseqTM_Package/PureseqTM.sh -i ../setup/{pdb}/{pdb}.fasta")
        # do(f"/projects/pw8/wl45/topology_prediction/PureseqTM_Package/PureseqTM_proteome.sh -i ../setup/{pdb}/{pdb}.fasta")
        do(f"/projects/pw8/wl45/topology_prediction/PureseqTM_Package/PureseqTM_proteome.sh -i {pdb}.fasta")
        cd("..")

if args.day == "aug15":
    if args.mode == 1:
        for pdb in pdb_list:
            folder = pdb
            os.mkdir(folder)
            os.chdir(folder)
            os.system(f"cp ../fastas/{pdb}.fasta .")
            cmd = f"bash /projects/pw8/wl45/DeepMetaPSICOV/run_DMP.sh -i {pdb}.fasta"
            out = slurmRun(f"{pdb}.slurm", cmd)
            print(out)
            os.chdir("..")
if args.day == "aug11":
    if args.mode == 1:
        # pdb = "1pv6"
        do("mkdir -p forces_recompute")
        for pdb in dataset["hybrid"]:
            forceLocation = "forces_recompute"
            do(f"cp forces/forces_setup_{pdb}.py {forceLocation}/forces_setup_{pdb}.py")
            with fileinput.FileInput(f"{forceLocation}/forces_setup_{pdb}.py", inplace=True) as file:
                for line in file:
                    fromLine = 'MembranePart ='
                    toLine = 'MembranePart =[1,2]\nMembranePart ='
                    tmp = line.replace(fromLine, toLine)
                    print(tmp, end='')

if args.day == "aug10":
    if args.mode == 1:
        for pdb in dataset["hybrid"]:
            do(f"cp -r frag_database/{pdb}_HA/HA_combined.mem setup/{pdb}/")
            do(f"cp -r frag_database/{pdb}_HA/HA_frags setup/{pdb}/")
    if args.mode == 2:
        # pdb = "1pv6"
        for pdb in dataset["hybrid"]:
            forceLocation = "forces_HA"
            do(f"cp forces/forces_setup_{pdb}.py {forceLocation}/forces_setup_{pdb}.py")
            with fileinput.FileInput(f"{forceLocation}/forces_setup_{pdb}.py", inplace=True) as file:
                for line in file:
                    fromLine = 'fragment_memory_term(oa, frag_file_list_file="./frags.mem", npy_frag_table="./frags.npy", UseSavedFragTable=True),'
                    toLine = 'fragment_memory_term(oa, frag_file_list_file="./HA_combined.mem", npy_frag_table="./HA_combined.npy", UseSavedFragTable=True),'
                    tmp = line.replace(fromLine, toLine)
                    print(tmp, end='')
if args.day == "jul18":
    if args.mode == 1:
        pdb = "5d91"
        for i in range(3):
            cmd = f"mm_run_with_pulling_start.py setup/{pdb}/{pdb} --to cpu1/{pdb}/{i} -s 1e4 --platform CPU -t 1 -f energy_forces_with_er/forces_setup_{pdb}.py --subMode 3"
            jobId = slurmRun(f"slurms/{pdb}_{i}.slurm", cmd, memory=5)
            print(jobId)


if args.day == "jun29":
    pdb_list = ["2jo1", "2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91"]
    if args.mode == 1:
        for pdb in pdb_list:
            if pdb == "6e67A":
                continue
            print(pdb)
            do(f"cp -r frag_database/{pdb}_HA/HA_frags setup/{pdb}/")
            do(f"cp frag_database/{pdb}_HA/HA_combined.mem setup/{pdb}/")
    if args.mode == 2:
        for pdb in pdb_list:
            # do(f"cp forces_setup_{pdb}.py energy_forces/")
            # replace(f"energy_forces/forces_setup_{pdb}.py", "single_frags", "HA_combined")
            for i in range(3):
                cmd = f"python mm_analysis.py setup/{pdb}/{pdb} -t second_small_batch/{pdb}/{i}/movie.dcd --subMode 3 -f energy_forces/forces_setup_{pdb}.py -o frag_ha_energy.dat"
                slurmRun(f"{pdb}_{i}.slurm", cmd, memory=5)
    if args.mode == 3:
        pdb = "2xov_complete"
        for i in range(3):
            cmd = f"mm_analysis.py setup/2xov_complete/2xov_complete -t second_small_batch/2xov_complete/{i}/movie.dcd --subMode 3 -f forces_setup_2xov_complete.py"
            jobId = slurmRun(f"slurms/{pdb}_{i}.slurm", cmd, memory=5)
            print(jobId)

if args.day == "jun27":
    pdb_list = ["2jo1", "2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91"]
    # pdb_list = ["5xpd", "3kp9", "4a2n", "5d91"]
    # pdb_list = ["4xt3", "5uiw", "6akg"]
    if args.mode == 1:
        for pdb in pdb_list:
            cd("TM_pred")
            do(f"/projects/pw8/wl45/topology_prediction/PureseqTM_Package/PureseqTM.sh -i ../setup/{pdb}/{pdb}.fasta")
            cd("..")
            do(f"gg_server.py -d jun24 -m 2 -l {pdb}")
    # if args.mode == 2:

if args.day == "jun24":
    if args.mode == 1:
        # do("/projects/pw8/wl45/topology_prediction/PureseqTM_Package/PureseqTM.sh -i crystal_structure.fasta")
        print("This script will predict TM topology, create zimPosition, set up force-step file.")
        pdb = args.label
        cd("TM_pred")
        do(f"/projects/pw8/wl45/topology_prediction/PureseqTM_Package/PureseqTM.sh -i ../setup/{pdb}/{pdb}.fasta")
        cd("..")
        do(f"gg_server.py -d jun24 -m 2 -l TM_pred/6e67B_PureTM/{pdb}.prob")
    if args.mode == 2:
        # predict TM topology, create zimPosition, set up force-step file.

        pdb = args.label
        # pdb = "6e67A"
        # pdb = "6e67B"
        do(f"cp forces_setup.py forces_setup_{pdb}.py")
        # probFile =args.label
        # print(probFile)
        probFile = f"TM_pred/{pdb}_PureTM/{pdb}.prob"
        with open(f"{probFile}") as f:
            a = f.readlines()
        res_list = []
        first = None
        count = 1
        previousEnd = 0
        # print("g_all = [")
        zimOut = open("zimPosition", "w")
        out = "[\n"
        for i, line in enumerate(a[3:]):
            prob = float(line.strip().split()[3])
            res = "0" if prob < 0.5 else "1"
            o = "2" if res == "1" else "1"
            zimOut.write(o+"\n")
            if res == "0":
                if len(res_list) > 0:
                    # print(f"g{count} =", res_list)
                    print(res_list, ", ")
                    out += f"    {res_list},\n"
                    count += 1
                    last = res_list[-1]
                    first = res_list[0] if first is None else first
                    span = res_list[0] - previousEnd
                    if span > 30:
                        print(f"{pdb} Globular", previousEnd, res_list[0])
                        globular = list(range(previousEnd+10, res_list[0]-10))
                    previousEnd = last
                res_list = []
            if res == "1":
                res_list.append(i)
        n = len(a[3:])
        print(f"{pdb}: size {n}")
        span = n - previousEnd
        if span > 30:
            print(f"{pdb} Globular", previousEnd, n)
            globular = list(range(previousEnd+10, n-10))

        out += "]\n"
        zimOut.close()

        membranePart = []
        for i in range(first-5, last+5):
            if i not in globular:
                membranePart.append(i)
        # print("]")
        # replace(, "GALL", out)
        # , backup='.bak'
        with fileinput.FileInput(f"forces_setup_{pdb}.py", inplace=True) as file:
            for line in file:
                tmp = line.replace("GALL", out).replace("FIRST", str(first)).replace("LAST", str(last))
                tmp = tmp.replace("RESMEMB", f"{membranePart}")
                tmp = tmp.replace("RESGLOBULAR", f"{globular}")
                print(tmp, end='')
        do(f"cp zimPosition setup/{pdb}/PredictedZim")



# def PredictedZimAndForceSetupFile(pdb):
#     do(f"cp forces_setup.py forces_setup_{pdb}.py")
#     # probFile =args.label
#     # print(probFile)
#     probFile = f"TM_pred/{pdb}_PureTM/{pdb}.prob"
#     with open(f"{probFile}") as f:
#         a = f.readlines()
#     res_list = []
#     first = None
#     count = 1
#     previousEnd = 0
#     # print("g_all = [")
#     zimOut = open("zimPosition", "w")
#     out = "[\n"
#     for i, line in enumerate(a[3:]):
#         prob = float(line.strip().split()[3])
#         res = "0" if prob < 0.5 else "1"
#         o = "2" if res == "1" else "1"
#         zimOut.write(o+"\n")
#         if res == "0":
#             if len(res_list) > 0:
#                 # print(f"g{count} =", res_list)
#                 print(res_list, ", ")
#                 out += f"    {res_list},\n"
#                 count += 1
#                 last = res_list[-1]
#                 first = res_list[0] if first is None else first
#                 span = res_list[0] - previousEnd
#                 if span > 30:
#                     print(f"{pdb} Globular", previousEnd, res_list[0])
#                     globular = list(range(previousEnd+10, res_list[0]-10))
#                 previousEnd = last
#             res_list = []
#         if res == "1":
#             res_list.append(i)
#     n = len(a[3:])
#     print(f"{pdb}: size {n}")
#     span = n - previousEnd
#     if span > 30:
#         print(f"{pdb} Globular", previousEnd, n)
#         globular = list(range(previousEnd+10, n-10))

#     out += "]\n"
#     zimOut.close()

#     membranePart = []
#     for i in range(first-5, last+5):
#         if i not in globular:
#             membranePart.append(i)
#     # print("]")
#     # replace(, "GALL", out)
#     # , backup='.bak'
#     with fileinput.FileInput(f"forces_setup_{pdb}.py", inplace=True) as file:
#         for line in file:
#             tmp = line.replace("GALL", out).replace("FIRST", str(first)).replace("LAST", str(last))
#             tmp = tmp.replace("RESMEMB", f"{membranePart}")
#             tmp = tmp.replace("RESGLOBULAR", f"{globular}")
#             print(tmp, end='')
#     do(f"cp zimPosition setup/{pdb}/PredictedZim")

# def PredictedZimAndForceSetupFile(pdb, forceLocation="."):
#     do(f"mkdir -p {forceLocation}")
#     do(f"cp forces_setup.py {forceLocation}/forces_setup_{pdb}.py")
#     # probFile =args.label
#     # print(probFile)
#     # probFile = f"TM_pred/{pdb}_PureTM/{pdb}.prob"
#     topFile = f"TM_pred/{pdb}_topo"
#     with open(f"{topFile}") as f:
#         a = f.readlines()
#     assert len(a) == 3
#     res_list = []
#     first = None
#     count = 1
#     previousEnd = 0
#     # print("g_all = [")
#     zimOut = open("zimPosition", "w")
#     out = "[\n"
#     topo = a[2].strip()

#     cutoff = 30 if pdb != "1py6" else 15
#     linkerSize = 10 if pdb != "1py6" else 5
#     for i, res in enumerate(topo):
#         o = "2" if res == "1" else "1"
#         zimOut.write(o+"\n")
#         if res == "0":
#             if len(res_list) > 0:
#                 # print(f"g{count} =", res_list)
#                 print(res_list, ", ")
#                 out += f"    {res_list},\n"
#                 count += 1
#                 last = res_list[-1]
#                 first = res_list[0] if first is None else first
#                 span = res_list[0] - previousEnd
#                 if span > cutoff:
#                     print(f"{pdb} Globular", previousEnd, res_list[0])
#                     globular = list(range(previousEnd+linkerSize, res_list[0]-linkerSize))
#                 previousEnd = last
#             res_list = []
#         if res == "1":
#             res_list.append(i)
#     n = len(topo)
#     print(f"{pdb}: size {n}")
#     span = n - previousEnd
#     if span > cutoff:
#         print(f"{pdb} Globular", previousEnd, n)
#         globular = list(range(previousEnd+linkerSize, n-linkerSize))

#     out += "]\n"
#     zimOut.close()

#     membranePart = []
#     for i in range(first-5, last+5):
#         if i not in globular:
#             membranePart.append(i)
#     # print("]")
#     # replace(, "GALL", out)
#     # , backup='.bak'
#     with fileinput.FileInput(f"{forceLocation}/forces_setup_{pdb}.py", inplace=True) as file:
#         for line in file:
#             tmp = line.replace("GALL", out).replace("FIRST", str(first)).replace("LAST", str(last))
#             tmp = tmp.replace("RESMEMB", f"{membranePart}")
#             tmp = tmp.replace("RESGLOBULAR", f"{globular}")
#             print(tmp, end='')
#     do(f"cp zimPosition setup/{pdb}/PredictedZim")

if args.day == "jun20":
    if args.mode == 1:
        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(3):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i)
        n = i
        i = 0
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d jun20 -m 2 -l {proteins}")
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
        waitForJobs(jobIdList, sleepInterval=300)
    if args.mode == 2:
        from pyCodeLib import *
        proteinlist = readList(args.label)
        limit = 60
        for name in proteinlist:
            input_pdb_filename = f"../database/dompdb/{name}"
            structure = parse_pdb(input_pdb_filename)
            res_list = get_res_list(structure)
            neighbor_list = get_neighbor_list(structure)

            phi = phi_relative_k_well(res_list, neighbor_list, parameter_list="")
            np.savetxt(f"../phis/phi_relative_k_well_{name}_native_4.5_6.5_5.0_10", [phi], fmt='%1.4f')
            phis_list = []
            q_list = []
            offset_all = np.arange(1, limit, 1)
            for offset in offset_all:
                phi = phi_relative_k_well(res_list, neighbor_list, parameter_list="", z_m_high=offset)
                q = interaction_well(offset, 12, 18, 1)
                print(phi)
                phis_list.append(phi)
                q_list.append(q)
            # np.savetxt(f"../phis/{name}_phis_high", phis_list, fmt='%1.4f')
            offset_all = np.arange(-limit, -1, 1)
            for offset in offset_all:
                phi = phi_relative_k_well(res_list, neighbor_list, parameter_list="", z_m_low=offset)
                q = interaction_well(offset, -18, -12, 1)
                print(phi)
                phis_list.append(phi)
                q_list.append(q)
            np.savetxt(f"../phis/phi_relative_k_well_{name}_decoys_shifted_4.5_6.5_5.0_10", phis_list, fmt='%1.4f')
            np.savetxt(f"../phis/phi_relative_k_well_{name}_decoysQ_shifted_4.5_6.5_5.0_10", q_list, fmt='%1.4f')

if args.day == "jun14":
    if args.mode == 1:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        pdb_list = dataset["membrane"]
        # folder = "with_restart"
        # folder = "jun02"
        # folder = "jun03"
        # folder = "original"
        # folder = "original_compare_speed"
        folder = "with_pulling"
        totalRuns = 20
        subMode = 1


        # folder = "jun02_original"
        for pdb in pdb_list:
            initial_structure = f"{pdb}"

            toPath = f"{folder}/{pdb}/native"
            cmd = ""
            cmd += f"python mm_run.py setup/{pdb}/{pdb} --to {toPath} -m 1 -s 4e2 -p CPU -t 1 --tempStart 800 --tempEnd 200 --subMode {subMode}\n"
            cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t {toPath}/movie.dcd --subMode {subMode}\n"
            cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t {toPath}/native.pdb -o native.dat --subMode {subMode}"
            jobId = slurmRun(f"slurms/{folder}_{pdb}_native.slurm", cmd, memory=3)
            for i in range(totalRuns):
                toPath = f"{folder}/{pdb}/{i}_0"
                cmd = ""
                cmd += f"python mm_run_with_pulling_start.py setup/{pdb}/{initial_structure} --to {toPath} -m 1 -s 1e5 -p CPU -t 1 --tempStart 800 --tempEnd 600 --subMode {subMode}\n"
                cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t {toPath}/movie.dcd --subMode {subMode}"
                jobId = slurmRun(f"slurms/{folder}_{pdb}_{i}.slurm", cmd, memory=3)
                jobIdList.append(jobId)
                # exit()
        waitForJobs(jobIdList, sleepInterval=100)
        for pdb in pdb_list:
            initial_structure = f"{pdb}"

            for i in range(totalRuns):
                toPath = f"{folder}/{pdb}/{i}_1"
                cmd = ""
                cmd += f"python mm_run.py setup/{pdb}/{initial_structure} --to {toPath} -m 1 -s 1e5 -p CPU -t 1 --tempStart 600 --tempEnd 200 --subMode {subMode} --fromCheckPoint {folder}/{pdb}/{i}_0/checkpnt.chk\n"
                cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t {toPath}/movie.dcd --subMode {subMode}"
                jobId = slurmRun(f"slurms/{folder}_{pdb}_{i}_2.slurm", cmd, memory=3)
                jobIdList.append(jobId)
                # exit()


if args.day == "jun12":
    if args.mode == 1:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        pdb_list = dataset["membrane"]
        # folder = "with_restart"
        # folder = "jun02"
        # folder = "jun03"
        folder = "original"
        folder = "original_compare_speed"
        totalRuns = 20
        subMode = 1


        # folder = "jun02_original"
        for pdb in pdb_list:
            initial_structure = f"{pdb}"

            toPath = f"{folder}/{pdb}/native"
            cmd = ""
            cmd += f"python mm_run.py setup/{pdb}/{pdb} --to {toPath} -m 1 -s 4e2 -p CPU -t 1 --tempStart 800 --tempEnd 200 --subMode {subMode}\n"
            cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t {toPath}/movie.dcd --subMode {subMode}\n"
            cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t {toPath}/native.pdb -o native.dat --subMode {subMode}"
            jobId = slurmRun(f"slurms/{folder}_{pdb}_native.slurm", cmd, memory=3)
            for i in range(totalRuns):
                toPath = f"{folder}/{pdb}/{i}_0"
                cmd = ""
                cmd += f"python mm_run.py setup/{pdb}/{initial_structure} --to {toPath} -m 1 -s 1e5 -p CPU -t 1 --tempStart 800 --tempEnd 600 --subMode {subMode}\n"
                cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t {toPath}/movie.dcd --subMode {subMode}"
                jobId = slurmRun(f"slurms/{folder}_{pdb}_{i}.slurm", cmd, memory=3)
                jobIdList.append(jobId)
                # exit()
        waitForJobs(jobIdList, sleepInterval=100)
        for pdb in pdb_list:
            initial_structure = f"{pdb}"

            for i in range(totalRuns):
                toPath = f"{folder}/{pdb}/{i}_1"
                cmd = ""
                cmd += f"python mm_run.py setup/{pdb}/{initial_structure} --to {toPath} -m 1 -s 1e5 -p CPU -t 1 --tempStart 600 --tempEnd 200 --subMode {subMode} --fromCheckPoint {folder}/{pdb}/{i}_0/checkpnt.chk\n"
                cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t {toPath}/movie.dcd --subMode {subMode}"
                jobId = slurmRun(f"slurms/{folder}_{pdb}_{i}_2.slurm", cmd, memory=3)
                jobIdList.append(jobId)
                # exit()


if args.day == "jun10":
    # multi chain iter 0 membrane.
    from pyCodeLib import *
    import warnings
    warnings.filterwarnings('ignore')
    n_decoys = 2000

    if args.mode == 44:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d jun10 -m 4"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 11:
        n = 10
        # n = 1000
        do("mkdir -p proteins_name_list")
        do("mkdir -p slurms")
        do("mkdir -p decoys/multiShuffle")
        do("mkdir -p gammas")
        do('mkdir -p outs')
        do("mkdir -p ../phis")
        with open("protein_list") as f:
            content = f.readlines()
        for line in content:
            name = line.strip()
            generate_multiShuffle(name, num_decoys=n, nameMode=1)
    if args.mode == 111:
        # filter out protein with small MSAs.
        with open("protein_list") as f:
            content = f.readlines()
        with open("protein_list_filtered", "w") as out:
            for line in content:
                name = line.strip()
                with open(f"alignments/{name}.seqs") as f:
                    a = f.readlines()
                if len(a) < 100:
                    pass
                else:
                    out.write(line)

    if args.mode == 1:
        # do("cp ../phi_list.txt .")
        # do("cp ~/opt/optimization/phi_list_contact.txt phi_list.txt")

        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(2):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i)
        n = i
        i = 0
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d jun10 -m 2 -l {proteins}")
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
        waitForJobs(jobIdList, sleepInterval=300)
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d apr31 -m 4"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")

    if args.mode == 2:
        proteins = args.label
        sampleK=2
        # sampleK=1000
        evaluate_phis_over_training_set_for_native_structures_Wei(proteins, "phi_list.txt", decoy_method='multiShuffle', max_decoys=1e+10, tm_only=False, num_processors=1, multi_seq=True, sampleK=sampleK)

    if args.mode == 22:
        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        jobIdList = []
        for i in range(n):
            sub_list_file = f"proteins_name_list/proteins_name_list_{i}.txt"
            jobId = slurmRun(f"slurms/run_AB_{i}.slurm", f"python3 ~/opt/gg_server.py -d jun10 -m 3 -l {sub_list_file}", memory=5)
            jobIdList.append(jobId)
        # for testing toymodel
    if args.mode == 3:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # n = len(simulation_location_list)
        # n = 2
        # n = len(new_simulation_list)

        # import time
        # time.sleep(4000)
        # with open(f"proteins_name_list.txt", "w") as out:
        #     for data in ["new", "old"]:
        #         pdb_list, steps = dataset[data]
        #         for p in pdb_list:
        #             name = p.lower()[:4]
        #             out.write(f"{name}\n")
        # complete_proteins = "proteins_name_list/proteins_name_list_tiny.txt"
        # complete_proteins = "proteins_name_list.txt"
        # complete_proteins = "tiny_list.txt"
        complete_proteins = args.label
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='multiShuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=True)

        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        out = calculate_A_B_and_gamma_wl45_parallel(complete_proteins, "phi_list.txt", decoy_method='multiShuffle',
                                        num_decoys=1000, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=True)
        print("Success", out)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_xl23(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False, mode=1, simulation_location_list=simulation_location_list)
    # if args.mode == 4:
    #     with open(f"slurms/run_on_scavenge.slurm", "w") as out:
    #         out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr13 -m 3"))
    #     replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
    #     do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 4:
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list"
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='multiShuffle',
                                        num_decoys=1000, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=True)
if args.day == "jun09":
    if args.mode == 1:
        def pdbToFasta(pdb, pdbLocation, fastaFile, chains="A"):
            import textwrap
            from Bio.PDB.PDBParser import PDBParser
            from Bio.PDB.Polypeptide import three_to_one
            # pdb = "1r69"
            # pdbLocation = "/Users/weilu/Research/server/may_2019/family_fold/1r69.pdb"
            # chains = "A"
            # fastaFile = "/Users/weilu/Research/server/may_2019/family_fold/1r69.fasta"
            s = PDBParser().get_structure("X", pdbLocation)
            # m = s[0]  # model 0
            seq = ""
            with open(fastaFile, "w") as out:
                chain = "A"
                out.write(f">{pdb.upper()}:{chain.upper()}|PDBID|CHAIN|SEQUENCE\n")
                seq = ""
                for res in s.get_residues():
                    resName = three_to_one(res.resname)
                    seq += resName
                out.write("\n".join(textwrap.wrap(seq, width=80))+"\n")
            return seq

        # a = glob.glob("cleaned_pdbs/*.pdb")
        jobIdList = []
        a = glob.glob("../membrane_only_contact_optimization/database/dompdb/*.pdb")
        do("mkdir -p fasta")
        for line in a:
            print(line)
            pdb = line.split("/")[-1].split(".")[0]
            pdbToFasta(pdb, line, f"fasta/{pdb}.fasta")
            cmd = f"/projects/pw8/wl45/hh-suite/build/bin/hhblits -i fasta/{pdb}.fasta -d ../uniclust30_2018_08/uniclust30_2018_08 -o {pdb} -oa3m {pdb}.a3m -n 2"
            # do(cmd)
            jobId = slurmRun(f"slurms/{pdb}.slurm", cmd, memory=10)
            jobIdList.append(jobId)
    if args.mode == 2:
        with open("protein_list_unfiltered") as f:
            content = f.readlines()
        for line1 in content:
            a = line1.strip()
            print(a)
            data = get_MSA_data(f"../../MSA/{a}.a3m")
            with open(f"alignments/{a}_filtered_0.05.seqs", "w") as out:
                for line in data:
                    out.write(line+"\n")

if args.day == "jun06":

    pdb_list = ["2bg9", "1j4n", "1py6_SD", "2bl2", "1rhz", "1iwg", "2ic8", "1pv6", "1occ", "1kpl", "2bs2", "1py6", "1u19"]
    folder = "first"
    steps = 30
    if args.mode == 1:
        print("hi")
        for pdb in pdb_list:
            do(f"cp -r all_simulations/{pdb} {folder}/")
            cd(f"{folder}/{pdb}/")
            do(f"run.py -n 10 {pdb} --commons 2 -s {steps}")
            cd("../..")
    if args.mode == 2:
        for pp in pdb_list:
            p = pp
            cd(f"{folder}/{pp}/")
            do(f"run.py {p} -n -1 -s 1 --start crystal --commons 2")
            cd("../../")
    if args.mode == 3:
        for pp in pdb_list:
            p = pp
            cd(f"{folder}/{p}/simulation/native")
            do(f"rerun.py {p} -t -m 1")
            cd("../../../../")


############################    -2019-    ############################
'''
if args.day == "may30":
    if args.mode == 1:
        # iteratively generate gammas
        # first get pdbs.
        folder = "with_restart"
        database = "../../database"
        pdb_list = dataset["may13"]
        for pdb in pdb_list:
            if pdb in ["4cpv", "2mhr", "1mba", "2fha"]:
                continue
            print(pdb)
            all_data = []
            for i in range(20):
                cd(f"{pdb}/{i}_1/")
                do("python /projects/pw8/wl45/openawsem/helperFunctions/convertOpenmmTrajectoryToStandardMovie.py movie.pdb")
                tmp = pd.read_csv("info.dat", sep="\s+")
                tmp = tmp.assign(Run=i, Protein=pdb)
                all_data.append(tmp)
                cd("../../")
                do(f"mkdir -p {database}/{folder}_{pdb}_{i}")
                do(f"cp {pdb}/{i}_1/* {database}/{folder}_{pdb}_{i}/")
            data = pd.concat(all_data)
            data.to_csv(f"{database}/Q_{folder}_{pdb}")
        # process Qw info.

    if args.mode == 444:
        # generate decoys
        database_location = "../../database"
        do("mkdir -p decoys/lammps")
        import io
        from Bio.PDB.PDBParser import PDBParser
        simulation_location, name = args.label.split("__")
        simulation_location_name = f"{simulation_location}_{name}"

        def getStructures(x, all_movies):
            index = int(x["index"])+1
            run = int(x["Run"])

            start = index * size
            end = (index + 1) * size
            f = io.StringIO("".join(all_movies[run][start:end]))
            parser = PDBParser()
            return parser.get_structure(f"{index}", f)

        a = pd.read_csv(f"{database_location}/Q_{simulation_location_name}", index_col=0).query(f"Rank < {decoy_n*3}")
        sampled = a.sample(decoy_n)
        all_movies = {}
        for i in sampled["Run"].unique():
            with open(f"{database_location}/{simulation_location_name}_{i}/movie.pdb") as f:
                movie = f.readlines()
            all_movies[i] = movie
        size = 0
        for line in movie:
            size += 1
            if line == "ENDMDL\n":
                break
        print(simulation_location_name, size)
        sampled["structure"] = sampled.apply(getStructures, all_movies=all_movies, axis=1)
        sampled["Qw"] = sampled[" Qw"].round(3)
        sampled.drop(" Qw", axis=1)
        sampled.to_pickle(f"decoys/lammps/{name}_{simulation_location}.pkl")

# if args.day == "may29":
#     if args.mode == 1:
#         folder = "may29"
#         pdb = "1mba"
#         with open("cuda_run.sh", "w") as out:
#             for i in range(20):
#                 out.write(f"python mm_run.py setup/{pdb}/extended --to {folder}/{pdb}/{i}_0 -m 1 -s 3e6 -p CUDA --tempStart 800 --tempEnd 200\n")

if args.day == "may28":
    if args.mode == 2:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        pdb_list = dataset["may13"]
        # folder = "with_restart"
        # folder = "jun02"
        # folder = "jun03"
        folder = "frags_original"
        totalRuns = 20
        subMode = 1
        # folder = "jun02_original"
        for pdb in pdb_list:
            if pdb == "1mba" or pdb == "2fha":
                continue
            toPath = f"{folder}/{pdb}/native"
            cmd = ""
            cmd += f"python mm_run.py setup/{pdb}/{pdb} --to {toPath} -m 1 -s 1e2 -p CPU -t 1 --tempStart 800 --tempEnd 200 --subMode {subMode}\n"
            cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t {toPath}/movie.dcd --subMode {subMode}\n"
            cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t {toPath}/native.pdb -o native.dat --subMode {subMode}"
            jobId = slurmRun(f"slurms/{folder}_{pdb}_native.slurm", cmd, memory=3)
            for i in range(totalRuns):
                toPath = f"{folder}/{pdb}/{i}_0"
                cmd = ""
                cmd += f"python mm_run.py setup/{pdb}/extended --to {toPath} -m 1 -s 1e5 -p CPU -t 1 --tempStart 800 --tempEnd 600 --subMode {subMode}\n"
                cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t {toPath}/movie.dcd --subMode {subMode}"
                jobId = slurmRun(f"slurms/{folder}_{pdb}_{i}.slurm", cmd, memory=3)
                jobIdList.append(jobId)
                # exit()
        waitForJobs(jobIdList, sleepInterval=100)
        for pdb in pdb_list:
            if pdb == "1mba" or pdb == "2fha":
                continue
            for i in range(totalRuns):
                toPath = f"{folder}/{pdb}/{i}_1"
                cmd = ""
                cmd += f"python mm_run.py setup/{pdb}/extended --to {toPath} -m 1 -s 1e5 -p CPU -t 1 --tempStart 600 --tempEnd 200 --subMode {subMode} --fromCheckPoint {folder}/{pdb}/{i}_0/checkpnt.chk\n"
                cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t {toPath}/movie.dcd --subMode {subMode}"
                jobId = slurmRun(f"slurms/{folder}_{pdb}_{i}_2.slurm", cmd, memory=3)
                jobIdList.append(jobId)
                # exit()

    # if args.mode == 2:
    #     do("mkdir -p slurms")
    #     do("mkdir -p outs")
    #     jobIdList = []
    #     pdb_list = dataset["may13"]
    #     # folder = "with_restart"
    #     # folder = "jun02"
    #     folder = "jun02_original"
    #     totalRuns = 20
    #     for pdb in pdb_list:
    #         if pdb == "1mba" or pdb == "2fha":
    #             continue
    #         toPath = f"{folder}/{pdb}/native"
    #         cmd = ""
    #         cmd += f"python mm_run_original.py setup/{pdb}/{pdb} --to {toPath} -m 1 -s 1e2 -p CPU -t 1 --tempStart 800 --tempEnd 200\n"
    #         cmd += f"python mm_analysis_original.py setup/{pdb}/{pdb} -t {toPath}/movie.dcd"
    #         cmd += f"python mm_analysis_original.py setup/{pdb}/{pdb} -t {toPath}/native.pdb -o native.dat"
    #         jobId = slurmRun(f"slurms/{pdb}_native.slurm", cmd, memory=3)
    #         for i in range(totalRuns):
    #             toPath = f"{folder}/{pdb}/{i}_0"
    #             cmd = ""
    #             cmd += f"python mm_run_original.py setup/{pdb}/extended --to {toPath} -m 1 -s 1e5 -p CPU -t 1 --tempStart 600 --tempEnd 400\n"
    #             cmd += f"python mm_analysis_original.py setup/{pdb}/{pdb} -t {toPath}/movie.dcd"
    #             jobId = slurmRun(f"slurms/{pdb}_{i}_{folder}.slurm", cmd, memory=3)
    #             jobIdList.append(jobId)
    #             # exit()
    #     waitForJobs(jobIdList, sleepInterval=100)
    #     for pdb in pdb_list:
    #         if pdb == "1mba" or pdb == "2fha":
    #             continue
    #         for i in range(totalRuns):
    #             toPath = f"{folder}/{pdb}/{i}_1"
    #             cmd = ""
    #             cmd += f"python mm_run_original.py setup/{pdb}/extended --to {toPath} -m 1 -s 1e5 -p CPU -t 1 --tempStart 400 --tempEnd 200 --fromCheckPoint {folder}/{pdb}/{i}_0/checkpnt.chk\n"
    #             cmd += f"python mm_analysis_original.py setup/{pdb}/{pdb} -t {toPath}/movie.dcd"
    #             jobId = slurmRun(f"slurms/{pdb}_{i}_{folder}_2.slurm", cmd, memory=3)
    #             jobIdList.append(jobId)
    #             # exit()
    if args.mode == 1:
        do("mkdir -p slurms")
        do("mkdir -p outs")
        jobIdList = []
        pdb_list = dataset["may13"]
        # folder = "with_restart"
        # folder = "jun02"
        # folder = "jun03"
        folder = "frags"
        totalRuns = 20
        subMode = 0
        # folder = "jun02_original"
        for pdb in pdb_list:
            if pdb == "1mba" or pdb == "2fha":
                continue
            toPath = f"{folder}/{pdb}/native"
            cmd = ""
            cmd += f"python mm_run.py setup/{pdb}/{pdb} --to {toPath} -m 1 -s 1e2 -p CPU -t 1 --tempStart 800 --tempEnd 200 --subMode {subMode}\n"
            cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t {toPath}/movie.dcd --subMode {subMode}\n"
            cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t {toPath}/native.pdb -o native.dat --subMode {subMode}"
            jobId = slurmRun(f"slurms/{pdb}_native.slurm", cmd, memory=3)
            for i in range(totalRuns):
                toPath = f"{folder}/{pdb}/{i}_0"
                cmd = ""
                cmd += f"python mm_run.py setup/{pdb}/extended --to {toPath} -m 1 -s 1e5 -p CPU -t 1 --tempStart 800 --tempEnd 600 --subMode {subMode}\n"
                cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t {toPath}/movie.dcd --subMode {subMode}"
                jobId = slurmRun(f"slurms/{pdb}_{i}.slurm", cmd, memory=3)
                jobIdList.append(jobId)
                # exit()
        waitForJobs(jobIdList, sleepInterval=100)
        for pdb in pdb_list:
            if pdb == "1mba" or pdb == "2fha":
                continue
            for i in range(totalRuns):
                toPath = f"{folder}/{pdb}/{i}_1"
                cmd = ""
                cmd += f"python mm_run.py setup/{pdb}/extended --to {toPath} -m 1 -s 1e5 -p CPU -t 1 --tempStart 600 --tempEnd 200 --subMode {subMode} --fromCheckPoint {folder}/{pdb}/{i}_0/checkpnt.chk\n"
                cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t {toPath}/movie.dcd --subMode {subMode}"
                jobId = slurmRun(f"slurms/{pdb}_{i}_2.slurm", cmd, memory=3)
                jobIdList.append(jobId)
                # exit()

if args.day == "may26":
    if args.mode == 1:
        do("mkdir -p slurms")
        jobIdList = []
        pdb_list = dataset["may13"]
        for pdb in pdb_list:
            for i in range(20):
                cmd = f"python mm_run.py setup/{pdb}/extended --to fast_10f_timeStep/{pdb}/{i} -m 1 -s 1e5 -p CPU -t 1\n"
                cmd += f"python mm_analysis.py setup/{pdb}/{pdb} -t fast_10f_timeStep/{pdb}/{i}/movie.dcd"
                jobId = slurmRun(f"slurms/{pdb}_{i}.slurm", cmd, memory=3)
                jobIdList.append(jobId)
                # exit()
if args.day == "may21":
    if args.mode == 1:
        # make sure all dump files exists.
        a = glob.glob("/scratch/wl45/may_2019/database/*")
        for one in a:
            if "Q_" in one:
                continue
            if os.path.exists(os.path.join(one, "dump.lammpstrj")):
                print(one)
    if args.mode == 2:
        splitPDB(".", "movie.pdb")

if args.day == "may19":
    if args.mode == 1:
        newFolder = "optimization_4_correct"
        do("gg_server.py -d may13 -m 1")
        do("gg_server.py -d may13 -m 3")
        cd("../iter1_optimization_decoys_2000")
        do(f"mkdir -p {newFolder}")
        cd(newFolder)
        do("gg_server.py -d may13 -m 4")
        do("gg_server.py -d may13 -m 5")
        do("gg_server.py -d may13 -m 6")

if args.day == "may13":
    # pdb_list = dataset["combined"]
    formatName = False
    pdb_list = dataset["may13"]
    Run = 30
    # decoy_n = 1000
    decoy_n = 2000
    folder_list = []
    folder_list_1 = ["original", "multi_iter0"]
    folder_list_2 = ["multi_constant_tc_frag", "multi_constant_tc", "original_fragMemory", "multi_iter0_fragMemory"]
    folder_list_3 = ["iter1_30", "iter1_90", "iter1_30_frag", "iter1_90_frag"]
    folder_list_4 = ["iter2_30", "iter2_90", "iter2_30_frag", "iter2_90_frag"]
    folder_list_5 = ["iter3_30", "iter3_90", "iter3_30_frag", "iter3_90_frag"]
    folder_list_6 = ["iter1_30_correct", "iter1_30_correct_frag", "iter1_80_correct", "iter1_80_correct_frag"]
    folder_list_7 = ["iter2_30_correct", "iter2_30_correct_frag", "iter2_80_correct", "iter2_80_correct_frag"]
    folder_list_8 = ["iter3_30_correct", "iter3_30_correct_frag", "iter3_80_correct", "iter3_80_correct_frag"]
    folder_list = folder_list_1 + folder_list_2
    # folder_list = folder_list_2
    # folder_list += folder_list_3
    # folder_list += folder_list_4
    # folder_list += folder_list_5
    folder_list += folder_list_6
    folder_list += folder_list_7
    folder_list += folder_list_8
    # folder_list = folder_list_1
    # folder_list = ["original"]

    if args.mode == 1:
        # Convert Dump File to Pdbs
        jobIdList = []
        for simulation_location in folder_list:
            cd(simulation_location)
            for pp in pdb_list:
                if formatName:
                    p = pp.lower()[:4]
                else:
                    p = pp
                cd(f"{p}/simulation")
                proteins = p
                jobId = slurmRun(f"convertDumpToPdbs.slurms", f"python3 ~/opt/gg_server.py -d may13 -m 11 -l {proteins}")
                jobIdList.append(jobId)
                cd("../..")
            cd("..")
        waitForJobs(jobIdList, sleepInterval=100)
    if args.mode == 2:
        # Transport Pdbs to database folder
        jobIdList = []
        for simulation_location in folder_list:
            cd(simulation_location)
            for pp in pdb_list:
                if formatName:
                    p = pp.lower()[:4]
                else:
                    p = pp
                cd(f"{p}/simulation")
                proteins = f"/scratch/wl45/may_2019/database/{simulation_location}_{p}"
                jobId = slurmRun(f"transportPDBs.slurm", f"python3 ~/opt/gg_server.py -d may13 -m 22 -l {proteins}", memory=3)
                jobIdList.append(jobId)
                cd("../..")
            cd("..")
        waitForJobs(jobIdList, sleepInterval=200)
    if args.mode == 11:
        name = args.label
        for i in range(30):
            cd(f"{i}/0/")
            do(f"cp ../{name}.seq .")
            do(f"python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie {name}.seq")
            cd("../..")
    if args.mode == 22:
        # python3 ~/opt/gg_server.py -d may06 -m 22 -l /scratch/wl45/may_2019/database/original_1r69
        cmd = args.label
        for i in range(30):
            pre = f"{cmd}_{i}/"
            do(f"mkdir -p {pre}")
            do(f"cp -r {i}/0/* {pre}/")
            fileName = "movie.pdb"
            # splitPDB(pre, fileName)
            # do(f"rm {pre}movie.*")
    if args.mode == 3:
        # process Qw info.
        print("Processing Qw files")
        # print(len(simulation_location_list))
        for simulation_location in folder_list:
            print(simulation_location)
            for p in pdb_list:
                # name = p.lower()[:4]
                name = p
                complete_Q = []
                for i in range(Run):
                    if not os.path.exists(f"../database/{simulation_location}_{name}_{i}/dump.lammpstrj"):
                        print(f"File not exist {simulation_location}_{name}_{i}/dump.lammpstrj")
                        continue
                    if not os.path.exists(f"../database/{simulation_location}_{name}_{i}/movie.pdb"):
                        print(f"File not exist {simulation_location}_{name}_{i}/movie.pdb")
                        continue
                    try:
                        Q = pd.read_csv(f"../database/{simulation_location}_{name}_{i}/wham.dat")[" Qw"]
                    except:
                        print(f"File not exist {simulation_location}_{name}_{i}/wham.dat")
                        continue
                    n = len(Q)
                    if n != 2000 and n != 1000 and n != 750:
                        print(f"database/{simulation_location}_{name}_{i},  {n}")
                    if n == 0:
                        # print("!! Using previous i.")
                        print(f"An error at {i}")
                        continue

                    complete_Q.append(pd.DataFrame(Q).assign(Run=i))
                data = pd.concat(complete_Q).reset_index(drop=False)
                data["Rank"] = data["index"].rank(ascending=False)
                data.to_csv(f"../database/Q_{simulation_location}_{name}")


    if args.mode == 4:
        do("mkdir -p slurms")
        jobIdList = []
        for simulation_location in folder_list:
            for p in pdb_list:
                # name = p.lower()[:4]
                name = p
                cmd = f"python3 ~/opt/gg_server.py -d may13 -m 444 -l {simulation_location}__{name}"
                jobId = slurmRun(f"slurms/decoys_{name}_{simulation_location}.slurm", cmd, memory=10)
                jobIdList.append(jobId)
        waitForJobs(jobIdList, sleepInterval=200)

    if args.mode == 444:
        # generate decoys
        database_location = "../../database"
        do("mkdir -p decoys/lammps")
        import io
        from Bio.PDB.PDBParser import PDBParser
        simulation_location, name = args.label.split("__")
        simulation_location_name = f"{simulation_location}_{name}"

        def getStructures(x, all_movies):
            index = int(x["index"])+1
            run = int(x["Run"])

            start = index * size
            end = (index + 1) * size
            f = io.StringIO("".join(all_movies[run][start:end]))
            parser = PDBParser()
            return parser.get_structure(f"{index}", f)

        a = pd.read_csv(f"{database_location}/Q_{simulation_location_name}", index_col=0).query(f"Rank < {decoy_n*3}")
        sampled = a.sample(decoy_n)
        all_movies = {}
        for i in sampled["Run"].unique():
            with open(f"{database_location}/{simulation_location_name}_{i}/movie.pdb") as f:
                movie = f.readlines()
            all_movies[i] = movie
        size = 0
        for line in movie:
            size += 1
            if line == "END\n":
                break
        print(simulation_location_name, size)
        sampled["structure"] = sampled.apply(getStructures, all_movies=all_movies, axis=1)
        sampled["Qw"] = sampled[" Qw"].round(3)
        sampled.drop(" Qw", axis=1)
        sampled.to_pickle(f"decoys/lammps/{name}_{simulation_location}.pkl")

    # if args.mode == 44:
    #     # generate decoys
    #     database_location = "../../database"
    #     do("mkdir -p decoys/lammps")
    #     import io
    #     from Bio.PDB.PDBParser import PDBParser

    #     def getStructures(x, all_movies):
    #         index = int(x["index"])+1
    #         run = int(x["Run"])

    #         start = index * size
    #         end = (index + 1) * size
    #         f = io.StringIO("".join(all_movies[run][start:end]))
    #         parser = PDBParser()
    #         return parser.get_structure(f"{index}", f)

    #     for simulation_location in folder_list:
    #         for p in pdb_list:
    #             # name = p.lower()[:4]
    #             name = p
    #             a = pd.read_csv(f"{database_location}/Q_{simulation_location}_{name}", index_col=0).query(f"Rank < {decoy_n*3}")
    #             sampled = a.sample(decoy_n)
    #             all_movies = {}
    #             for i in sampled["Run"].unique():
    #                 with open(f"{database_location}/{simulation_location}_{name}_{i}/movie.pdb") as f:
    #                     movie = f.readlines()
    #                 all_movies[i] = movie
    #             size = 0
    #             for line in movie:
    #                 size += 1
    #                 if line == "END\n":
    #                     break
    #             print(simulation_location, p, size)
    #             sampled["structure"] = sampled.apply(getStructures, all_movies=all_movies, axis=1)
    #             sampled["Qw"] = sampled[" Qw"].round(3)
    #             sampled.drop(" Qw", axis=1)
    #             sampled.to_pickle(f"decoys/lammps/{name}_{simulation_location}.pkl")
    # if args.mode == 4:
    #     # generate decoys
    #     database_location = "../../database"
    #     do("mkdir -p decoys/lammps")
    #     for simulation_location in folder_list:
    #         for p in pdb_list:
    #             # name = p.lower()[:4]
    #             name = p
    #             a = pd.read_csv(f"{database_location}/Q_{simulation_location}_{name}", index_col=0).query(f"Rank < {decoy_n*3}")
    #             sampled = a.sample(decoy_n)
    #             with open(f"decoys/lammps/{name}_{simulation_location}.decoys", "w") as out:
    #                 for i, item in sampled.iterrows():
    #                     location = f"{database_location}/{simulation_location}_{name}_{int(item['Run'])}/frame{int(item['index']+1)}"
    #                     qw = np.round(item[' Qw'], 3)
    #                     if not os.path.exists(location+".pdb"):
    #                         print(location, "not exists")
    #                         location = pre
    #                         qw = preQw
    #                     out.write(f"{location} {qw}\n")
    #                     pre = location
    #                     preQw = qw
    if args.mode == 5:
        # Compute Phis
        print("Computing Phis")
        do("mkdir -p proteins_name_list")
        do("mkdir -p slurms")
        do("mkdir -p ../phis")
        do("cp ../../gammas/for_simulation/iteration_multi_constraint_tc_constant_burial_gamma.dat burial_gamma.dat")
        do("cp ../../gammas/for_simulation/iteration_multi_constraint_tc_constant_gamma.dat gamma.dat")
        do("cp ~/opt/optimization/phi_list_contact.txt .")
        do("cp ~/opt/optimization/group_normalization.txt .")
        do("rm phi_list.txt")
        do("cat phi_list_contact.txt >> phi_list.txt")
        do("cat group_normalization.txt >> phi_list.txt")
        jobIdList = []
        for simulation_location in folder_list:
            for p in pdb_list:
                # name = p.lower()[:4]
                name = p
                with open(f"proteins_name_list/proteins_name_list_{name}_{simulation_location}.txt", "w") as out:
                        out.write(f"{name}_{simulation_location}\n")
                cmd = f"python3 ~/opt/compute_phis.py -m 3 proteins_name_list/proteins_name_list_{name}_{simulation_location}.txt"
                jobId = slurmRun(f"slurms/run_{name}_{simulation_location}.slurm", cmd, memory=10)
                jobIdList.append(jobId)
        waitForJobs(jobIdList, sleepInterval=200)
    if args.mode == 6:
        # time.sleep(4000)
        do("mkdir -p gammas")
        runFile = f"slurms/run_on_scavenge_2.slurm"
        with open(runFile, "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d may13 -m 66 -l group_normalization.txt"))
        replace(runFile, "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch {runFile}")

        runFile = f"slurms/run_on_scavenge_1.slurm"
        with open(runFile, "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d may13 -m 66 -l phi_list_contact.txt"))
        replace(runFile, "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch {runFile}")
    if args.mode == 66:
        # for testing toymodel
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # n = len(simulation_location_list)
        # n = 2
        folder_list = folder_list_1 + folder_list_2
        # folder_list += folder_list_3
        # folder_list += folder_list_4
        # folder_list += folder_list_5
        folder_list += folder_list_6
        folder_list += folder_list_7
        folder_list += folder_list_8
        n = len(folder_list)
        print(n, folder_list)
        phi_file = args.label
        # import time
        # time.sleep(4000)
        with open(f"proteins_name_list.txt", "w") as out:
            for p in pdb_list:
                # name = p.lower()[:4]
                name = p
                out.write(f"{name}\n")
        complete_proteins = "proteins_name_list.txt"
        simulation_location_list_dic = defaultdict(list)
        for p in pdb_list:
            # name = p.lower()[:4]
            name = p
            simulation_location_list_dic[name] += folder_list

        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, phi_file, decoy_method='lammps',
                                        num_decoys=n*decoy_n, noise_filtering=True, jackhmmer=False, read=False, mode=2, withBiased=True, simulation_location_list_dic=simulation_location_list_dic)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_xl23(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False, mode=1, simulation_location_list=simulation_location_list)
    if args.mode == 7:
        do("mkdir -p database/dompdb")
        do("mkdir -p optimization")
        for p in pdb_list:
            # name = p.lower()[:4]
            name = p
            do(f"cp ../database/dompdb/{name}.pdb database/dompdb/")



if args.day == "may11":
    if args.mode == 1:
        formatName = True
        pdb_list = dataset["combined"]
        folder_list = ["original", "multi_iter0"]
        # folder_list = ["multi_constant_tc_frag", "multi_constant_tc", "original_fragMemory", "multi_iter0_fragMemory"]
        jobIdList = []
        for simulation_location in folder_list:
            for pp in pdb_list:
                print(pp)
                if formatName:
                    p = pp.lower()[:4]
                else:
                    p = pp
                for i in range(30):
                    cd(f"{simulation_location}_{p}_{i}")
                    proteins = p
                    jobId = slurmRun(f"convertDumpToPdbs.slurms", f"python3 ~/opt/gg_server.py -d may11 -m 2 -l {proteins}")
                    jobIdList.append(jobId)
                    cd("..")
        waitForJobs(jobIdList, sleepInterval=300)
    if args.mode == 2:
        protein = args.label
        a = glob.glob("frame*.pdb")
        a.sort(key=natural_keys)
        print(len(a))
        with open("ff_energy_smaller.dat", "w") as out:
            for frame in a[::20]:
                cmd = f"python3 ~/opt/compute_energy.py {frame} -l /scratch/wl45/may_2019/family_fold/ff_contact/{protein}/"
                line = getFromTerminal(cmd)
                out.write(line)
    if args.mode == 3:
        pdb_list = dataset["combined"]
        folder_list = ["original", "multi_iter0"]
        for location in folder_list:
            cd(location)
            do("pwd")
            for pp in pdb_list:
                p = pp.lower()[:4]
                name = p
                cd(f"{p}/simulation/native/rerun")
                # do(f"cp ../{name}.seq .")
                # do(f"python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie {name}.seq")
                # fileName = "movie.pdb"
                # splitPDB("./", fileName)
                cmd = f"python3 ~/opt/compute_energy.py frame0.pdb -l /scratch/wl45/may_2019/family_fold/ff_contact/{name}/ > newContact.dat"
                line = do(cmd)
                cd("../../../..")
            cd("..")
if args.day == "may10":
    if args.mode == 1:
        a = glob.glob("*.out")
        for line in a:
            a = do(f"grep 'mm_run' {line}", get=True, show=False)
            try:
                b = a.split(",")[5]
                print(b)
                do(f"grep 'hours' {line}")
            except:
                pass

if args.day == "may09":
    if args.mode == 1:
        for cpu in [1, 2, 4, 8, 16]:
            for thread in [1, 2, 4, 8, 16]:
                name = f"without_beta_test_cpu_{cpu}_t_{thread}"
                cmd = f"python3 mm_run.py 1r69 -m 1 --to {name} --platform CPU -t {thread}"
                slurmRun(f"{name}.slurm", cmd, memory=2)
                replace(f"{name}.slurm", "#SBATCH --cpus-per-task=1", f"#SBATCH --cpus-per-task={cpu}")
if args.day == "may07":
    if args.mode == 1:
        splitPDB(".", "movie.pdb")


if args.day == "may05":
    folder_list = ["original", "multi_iter0"]
    folder_list = ["noContact", "multi_iter0_fragMemory", "original_fragMemory"]
    if args.mode == 1:
        pdb_list = dataset["combined"]

        for location in folder_list:
            cd(location)
            do("pwd")
            for pp in pdb_list:
                p = pp.lower()[:4]
                # do(f"run.py {p} -n -1 -s 1 --start crystal --commons 2")
                # cd(f"{p}/simulation/native")
                for i in range(30):
                    cd(f"{p}/simulation/{i}")
                    do(f"rerun.py {p} -t -m 1")
                    cd("../../..")
            cd("..")
    if args.mode == 2:
        pdb_list = dataset["combined"]
        for location in folder_list:
            cd(location)
            do("pwd")
            for pp in pdb_list:
                p = pp.lower()[:4]
                cd(p)
                do(f"run.py {p} -n -1 -s 1 --start crystal --commons 2")
                cd("..")
            cd("..")
    if args.mode == 3:
        pdb_list = dataset["combined"]
        for location in folder_list:
            cd(location)
            do("pwd")
            for pp in pdb_list:
                p = pp.lower()[:4]
                # do(f"run.py {p} -n -1 -s 1 --start crystal --commons 2")
                cd(f"{p}/simulation/native")
                do(f"rerun.py {p} -t -m 1")
                cd("../../..")
            cd("..")
if args.day == "may04":
    # multi chain iter 0
    # data = pd.read_csv("chosen.csv", index_col=0)
    data = pd.read_csv("data_info_3.csv", index_col=0).query("Problematic != 4")
    if args.mode == 1:
        # do("cp ../phi_list.txt .")
        # do("cp ~/opt/optimization/phi_list_contact.txt phi_list.txt")
        do("mkdir proteins_name_list")
        do("mkdir slurms")
        do("mkdir -p decoys/multiShuffle")
        do("mkdir -p gammas")
        do("mkdir -p ../phis")
        complete_data_filtered = data.query("Problematic != 4").reset_index(drop=True)
        to_location = "./"
        n = len(complete_data_filtered)
        number_of_runs = int(np.ceil(n/5))
        perRun = 5
        count = 0
        for i in range(number_of_runs):
            with open(to_location+f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                cc = 0
                while count < n and cc < perRun:
                    fullName = complete_data_filtered.iloc[count]["FullName"]
                    out.write(fullName+"\n")
                    cc += 1
                    count += 1
        print(number_of_runs)
    if args.mode == 2:
        do("mkdir -p decoys/multiShuffle")
        # for i, line in data.iterrows():
        #     fullName = line["FullName"]
        #     generate_multiShuffle(fullName, num_decoys=1000)
        jobIdList = []
        for i in range(362):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 4 {proteins}")
            jobIdList.append(jobId)
        # waitForJobs(jobIdList, sleepInterval=300)
        # with open(f"slurms/run_on_scavenge.slurm", "w") as out:
        #     out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d may04 -m 3"))
        # replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        # do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 22:
        jobIdList = []
        for i in range(362):
            sub_list_file = f"proteins_name_list/proteins_name_list_{i}.txt"
            jobId = slurmRun(f"slurms/run_AB_{i}.slurm", f"python3 ~/opt/gg_server.py -d may04 -m 3 -l {sub_list_file}", memory=5)
            jobIdList.append(jobId)
        # for testing toymodel
    if args.mode == 3:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # n = len(simulation_location_list)
        # n = 2
        # n = len(new_simulation_list)

        # import time
        # time.sleep(4000)
        # with open(f"proteins_name_list.txt", "w") as out:
        #     for data in ["new", "old"]:
        #         pdb_list, steps = dataset[data]
        #         for p in pdb_list:
        #             name = p.lower()[:4]
        #             out.write(f"{name}\n")
        # complete_proteins = "proteins_name_list/proteins_name_list_tiny.txt"
        # complete_proteins = "proteins_name_list.txt"
        # complete_proteins = "tiny_list.txt"
        complete_proteins = args.label
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='multiShuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=True)

        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        out = calculate_A_B_and_gamma_wl45_parallel(complete_proteins, "phi_list.txt", decoy_method='multiShuffle',
                                        num_decoys=1000, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=True)
        print("Success", out)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_xl23(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False, mode=1, simulation_location_list=simulation_location_list)
    # if args.mode == 4:
    #     with open(f"slurms/run_on_scavenge.slurm", "w") as out:
    #         out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr13 -m 3"))
    #     replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
    #     do(f"sbatch slurms/run_on_scavenge.slurm")

if args.day == "may02":
    # relative k
    offset_all = np.arange(-50, 50, 1)
    n_decoys = len(offset_all)
    if args.mode == 1:
        do(f"mkdir -p proteins_name_list")
        do(f"mkdir -p decoys/shifted")
        do(f"mkdir -p gammas")
        do(f"mkdir -p outs")
        do("mkdir -p ../phis")
        do("mkdir slurms")
    if args.mode == 2:
        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(1):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i)
        n = i
        i = 0
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d may02 -m 22 -l {proteins}")
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
        waitForJobs(jobIdList, sleepInterval=300)
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d may02 -m 3"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 22:
        proteinlist = readList(args.label)
        print(proteinlist)
        # generate decoys file and decoys.
        for name in proteinlist:
            print(name)
            os.system(f"mkdir -p decoyData/{name}")
            for offset in offset_all:
                duplicate_pdb(f"../database/dompdb/{name}.pdb", f"decoyData/{name}/{name}_{offset}.pdb", offset_z=offset, new_chain=-1)
            original_inside_or_not_table = get_inside_or_not_table(f"../database/dompdb/{name}.pdb")
            with open(f"decoys/shifted/{name}.decoys", "w") as out:
                for offset in offset_all:
                    inside_or_not_table = get_inside_or_not_table(f"decoyData/{name}/{name}_{offset}.pdb")
                    ratio = np.sum(np.array(inside_or_not_table) == np.array(original_inside_or_not_table))/len(original_inside_or_not_table)
                    q = ratio
                    out.write(f"decoyData/{name}/{name}_{offset} {q}\n")
            # cleanPdb([f"{name}"], chain="A", formatName=False)
        do(f"python3 ~/opt/compute_phis.py -m 6 {args.label}")
    if args.mode == 3:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list"
        n_decoys = 118
        # complete_proteins = "top70"
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='shifted',
                                        num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, withBiased=True)
if args.day == "apr31":
    # four body helix
    from pyCodeLib import *
    import warnings
    warnings.filterwarnings('ignore')
    n_decoys = 1000
    n_decoys = 2000
    # n_decoys = 4000
    # n_decoys = 8000
    separateDecoysNum = 2
    # do("cp -r transmembrane/* ~/Research/server/may_2019/four_body_helix/database/dompdb/")
    if args.mode == 1:
        do("mkdir proteins_name_list")
        do("mkdir slurms")
        do("mkdir data")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p outs")
        do("mkdir -p ../phis")
        proteins = f"protein_list"
        generate_decoy_sequences(proteins, separateDecoysNum=separateDecoysNum, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
    if args.mode == 2:
        # time.sleep(16000)
        with open("protein_list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(2):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i)
        n = i
        i = 0
        jobIdList = []
        for i in range(n):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d apr31 -m 222 -l {proteins}")
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
        waitForJobs(jobIdList, sleepInterval=300)
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d apr31 -m 4"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 22:
        protein = args.label
        do(f"python3 ~/opt/compute_phis.py -m 0 {protein}")
        # do("python3 ~/opt/compute_phis.py -m 0 test_protein")
    if args.mode == 222:
        proteins = args.label
        evaluate_phis_over_training_set_for_native_structures_Wei(proteins, "phi_list.txt",
                    decoy_method='shuffle', max_decoys=1e+10, tm_only=False, num_processors=1, separateDecoysNum=separateDecoysNum)
    if args.mode == 3:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr31 -m 4"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 4:
        # complete_proteins = "iter0.txt"
        complete_proteins = "protein_list"
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='shuffle',
                                        num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False)

if args.day == "apr30":
    from pyCodeLib import *
    import warnings
    warnings.filterwarnings('ignore')
    n_decoys = 1000
    # iteration 0
    if args.mode == 11:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr30 -m 5"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")

    if args.mode == 1:
        do("mkdir proteins_name_list")
        do("mkdir slurms")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p ../phis")
        do("cp ~/opt/optimization/phi* .")
        do("cp phi_")
        with open("../database/cath-dataset-nonredundant-S20Clean.list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(20):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i-1)

    if args.mode == 2:
        # do("cp ~/opt/optimization/phi_list_* .")
        # do("cp phi_list_contact.txt phi_list.txt")
        from pyCodeLib import *
        i = 0
        jobIdList = []
        do("rm iter0.txt")
        for i in range(712):
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d apr30 -m 5 -l {proteins}")
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
        waitForJobs(jobIdList, sleepInterval=300)

        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr30 -m 6"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        replace(f"slurms/run_on_scavenge.slurm", "scavenge", "commons")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 22:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr30 -m 6"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        replace(f"slurms/run_on_scavenge.slurm", "scavenge", "commons")
        replace(f"slurms/run_on_scavenge.slurm", "--time=04:00:00", "--time=23:50:00")

        do(f"sbatch slurms/run_on_scavenge.slurm")

    if args.mode == 3:
        do("python ~/opt/compute_phis.py proteins_name_list/proteins_name_list_0.txt -m 0")
    if args.mode == 5:
        proteins = args.label
        generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
        do(f"python3 ~/opt/compute_phis.py -m 0 {proteins}")
    if args.mode == 6:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        complete_proteins = "iter0.txt"
        # complete_proteins = "iter_partial.txt"
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='shuffle',
                                        num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False)
    if args.mode == 7:
        complete_proteins = "../database/cath-dataset-nonredundant-S20Clean.list"
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
                                num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
    if args.mode == 77:
        complete_proteins = "../database/cath-dataset-nonredundant-S20Clean.list"
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
                                num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=1)
    if args.mode == 8:
        name = args.label
        subset = read_column_from_file(name, 1)
        print(subset)
        complete_proteins = "../database/cath-dataset-nonredundant-S20Clean.list"
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
                                        num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=subset, subset_index=args.label.split('/')[-1].split('.')[0], read=0)
    if args.mode == 9:
        jobIdList = []
        for i in range(712):
            # proteins = f"../database/cath-dataset-nonredundant-S20Clean.list"
            proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
            # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
            jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/gg_server.py -d apr30 -m 8 -l {proteins}", memory=5)
            # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
            jobIdList.append(jobId)
            do(f"cat {proteins} >> iter0_complete.txt")
        waitForJobs(jobIdList, sleepInterval=300)
    if args.mode == 10:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr30 -m 7"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        # replace(f"slurms/run_on_scavenge.slurm", "scavenge", "commons")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    # if args.mode == 6:
    #     # existing check
    #     with open("../database/cath-dataset-nonredundant-S20Clean.list", "r") as f:
    #         for line in f:
    #             # print("!", line.strip(), "!")
    #             a = line.strip()
    #             if os.path.exists(f"../database/S20_seq/{a}.seq") and os.path.exists(f"../database/dompdb/{a}.pdb"):
    #                 # print("yes")
    #                 pass
    #             else:
    #                 print(f"../database/S20_seq/{a}.seq")
    #                 print(f"../database/dompdb/{a}.pdb")
    #                 print("no")
    #             # break


if args.day == "apr24":
    if args.mode == 1:
            fileName = "movie.pdb"
            splitPDB("./", fileName)


if args.day == "apr20":
    if args.mode == 1:
        do("cp ../../weighted/clean/gamma.dat .")
        do("cp ../../weighted/clean/membrane_gamma_rescaled.dat .")
if args.day == "apr18":
    if args.mode == 1:
        a_list = glob.glob("*.out")
        do("rm proteins_name_list.txt")
        for a in a_list:
            line = getFromTerminal(f"grep 'phi_relative_k_well' {a} | wc -l")
            if line.strip() == "21":
                first_line = getFromTerminal(f"grep '/home/wl45/opt/compute_phis.py' {a}")
                # splited = (first_line.strip().split(",")[-1]).strip()[1:-2]
                splited = (first_line.strip().split(" ")[-1])
                print("!"+splited+"!")
                # break
                do(f"cat {splited} >> proteins_name_list.txt")
    if args.mode == 2:
        # generate gammas
        num_decoys = len(np.arange(-50, 50, 2))
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "proteins_name_list/proteins_name_list.txt"
        complete_proteins = "top1000"
        # complete_proteins = "top10"
        # complete_proteins = "middle1000"
        # complete_proteins = "chosen"
        # complete_proteins = "last1000"
        # complete_proteins = "proteins_name_list.txt"
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='shifted',
                                num_decoys=num_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, withBiased=True)



if args.day == "apr19":
    if args.mode == 1:
        d = pd.read_csv("seq_info.csv", index_col=0)
        pdb_list = d.query("length < 150 and index % 2 == 0")["protein"].tolist()
        # cd("original_with_minimization/")
        cd("multi_iter0_with_minimization")
        do("pwd")
        for p in pdb_list:
            cd(p)
            do(f"run.py {p} -n -1 -s 1 --start crystal --commons 2")
            cd("..")
        cd("..")
    if args.mode == 2:
        d = pd.read_csv("seq_info.csv", index_col=0)
        pdb_list = d.query("length < 150 and index % 2 == 0")["protein"].tolist()
        # cd("original_with_minimization/")
        for location in ["multi_iter0_with_minimization", "original_with_minimization"]:
            cd(location)
            do("pwd")
            for p in pdb_list:
                # do(f"run.py {p} -n -1 -s 1 --start crystal --commons 2")
                cd(f"{p}/simulation/native")
                do(f"rerun.py {p} -t -m 1")
                cd("../../..")
            cd("..")
    if args.mode == 3:
        d = pd.read_csv("seq_info.csv", index_col=0)
        pdb_list = d.query("length < 150 and index % 2 == 0")["protein"].tolist()
        # cd("original_with_minimization/")
        for location in ["multi_iter0_with_minimization", "original_with_minimization"]:
            cd(location)
            do("pwd")
            for p in pdb_list:
                # do(f"run.py {p} -n -1 -s 1 --start crystal --commons 2")
                # cd(f"{p}/simulation/native")
                for i in range(10):
                    cd(f"{p}/simulation/{i}")
                    do(f"rerun.py {p} -t -m 1")
                    cd("../../..")
            cd("..")



if args.day == "apr16":
    # pre = "/scratch/wl45/april_second_2019/rough_approximation/noClean/"
    # database = f"{pre}/../database/dompdb/"
    if args.mode == 1:
        # clean pdbs
        a = glob.glob(f"{pre}/../database/original_pdbs/*.pdb")
        print(len(a))
        random.seed(1)
        subset = a
        # subset = random.sample(a, 5000)
        content =[one.split("/")[-1][:4] for one in subset]
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"{pre}/proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(20):
                    if pos < n:
                        out.write(content[pos]+"\n")
                    pos += 1
                i += 1
        print(i-1)
        runs = i
        # Cleanup the pdbs
        jobIdList = []
        do("mkdir -p slurms")
        for i in range(runs):
            with open(f"slurms/run_{i}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr16 -m 4 -l proteins_name_list/proteins_name_list_{i}.txt"))
            # replace(f"run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=5G")
            # replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=40G")
            # do(f"sbatch run_on_scavenge.slurm")
            a = getFromTerminal(f"sbatch {pre}slurms/run_{i}.slurm")
            jobId = a.split(" ")[-1].strip()
            # print("!"+jobId+"!")
            jobIdList.append(jobId)
        waitForJobs(jobIdList, sleepInterval=100)

    if args.mode == 2:
        do(f"mkdir -p proteins_name_list")
        do(f"mkdir -p decoys/shifted")
        do(f"mkdir -p gammas")
        do(f"mkdir -p ../phis")
        # generate phis.
        a = glob.glob(f"../database/dompdb/*.pdb")
        print(len(a))
        random.seed(1)
        subset = random.sample(a, 4000)
        content =[one.split("/")[-1][:4] for one in subset]
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(20):
                    if pos < n:
                        out.write(content[pos]+"\n")
                    pos += 1
                i += 1
        print(i-1)
        runs = i
        jobIdList = []
        do("mkdir -p slurms")
        for i in range(runs):
            with open(f"slurms/get_phi_{i}.slurm", "w") as out:
                # out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr16 -m 5 -l proteins_name_list/proteins_name_list_{i}.txt"))
                out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr16 -m 6 -l proteins_name_list/proteins_name_list_{i}.txt"))
            # replace(f"run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=5G")
            # replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=40G")
            # do(f"sbatch run_on_scavenge.slurm")
            a = getFromTerminal(f"sbatch slurms/get_phi_{i}.slurm")
            jobId = a.split(" ")[-1].strip()
            # print("!"+jobId+"!")
            jobIdList.append(jobId)
        waitForJobs(jobIdList, sleepInterval=100)

    if args.mode == 3:
        # generate gammas
        # a = glob.glob("/scratch/wl45/april_second_2019/compare/database/cleaned_pdbs/*.pdb")
        # a = glob.glob(f"../database/dompdb/*.pdb")
        # random.seed(1)
        # subset = random.sample(a, 2000)
        # # subset = subset[500:]
        # # subset = subset[:500]
        # subset = random.sample(subset, 1000)
        num_decoys = len(np.arange(-50, 50, 2))

        # 6, 45, 52, 74
        for i in range(50,100):
            if i not in [6, 45, 52, 74]:
                do(f"cat proteins_name_list/proteins_name_list_{i}.txt >> proteins_name_list.txt")

        # with open(f"proteins_name_list/proteins_name_list.txt", "w") as out:
        #     for one in subset:
        #         name = one.split("/")[-1][:4]
        #         out.write(f"{name}\n")
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # complete_proteins = "proteins_name_list/proteins_name_list.txt"
        complete_proteins = "proteins_name_list.txt"
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='shifted',
                                num_decoys=num_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False)

    if args.mode == 4:
        proteinlist = readList(args.label)
        cleanPdb1(proteinlist, chain="first", formatName=False, toFolder=f"{pre}/../database/cleaned_pdbs/", source=f"{pre}/../database/original_pdbs/")

    if args.mode == 5:
        proteinlist = readList(args.label)
        offset_all = np.arange(-50, 50, 2)
        print(proteinlist)
        # generate decoys file and decoys.
        for name in proteinlist:
            print(name)
            os.system(f"mkdir -p decoyData/{name}")
            for offset in offset_all:
                duplicate_pdb(f"../database/dompdb/{name}.pdb", f"decoyData/{name}/{name}_{offset}.pdb", offset_z=offset, new_chain="A")

            with open(f"decoys/shifted/{name}.decoys", "w") as out:
                for offset in offset_all:
                    out.write(f"decoyData/{name}/{name}_{offset}\n")
            # cleanPdb([f"{name}"], chain="A", formatName=False)
        do(f"python3 ~/opt/compute_phis.py -m 5 {args.label}")

    if args.mode == 6:
        proteinlist = readList(args.label)
        offset_all = np.arange(-50, 50, 2)
        print(proteinlist)
        # generate decoys file and decoys.
        for name in proteinlist:
            print(name)
            os.system(f"mkdir -p decoyData/{name}")
            for offset in offset_all:
                duplicate_pdb(f"../database/dompdb/{name}.pdb", f"decoyData/{name}/{name}_{offset}.pdb", offset_z=offset, new_chain="A")
            original_inside_or_not_table = get_inside_or_not_table(f"../database/dompdb/{name}.pdb")
            with open(f"decoys/shifted/{name}.decoys", "w") as out:
                for offset in offset_all:
                    inside_or_not_table = get_inside_or_not_table(f"decoyData/{name}/{name}_{offset}.pdb")
                    ratio = np.sum(np.array(inside_or_not_table) == np.array(original_inside_or_not_table))/len(original_inside_or_not_table)
                    q = ratio
                    out.write(f"decoyData/{name}/{name}_{offset} {q}\n")
            # cleanPdb([f"{name}"], chain="A", formatName=False)
        do(f"python3 ~/opt/compute_phis.py -m 6 {args.label}")


if args.day == "apr16_2":
    pre = "/scratch/wl45/april_second_2019/compare/noClean/"
    database = "/scratch/wl45/april_second_2019/compare/database/dompdb/"

    if args.mode == 1:
        import glob
        a = glob.glob("original_pdbs/*")
        print(len(a))
        for i, one in enumerate(a):
            if i % 100 == 0:
                print(i, one)
            name = one.split("/")[-1][:4]
            do(f"mv {one} original_pdbs/{name}.pdb")
    if args.mode == 4:
        import glob
        a = glob.glob("original_pdbs_with_HETATM/*.pdb")
        print(len(a))
        for i, one in enumerate(a):
            if i % 100 == 0:
                print(i, one)
            name = one.split("/")[-1][:4]
            do(f"grep -v HETATM {one} > original_pdbs/{name}.pdb")

# if args.day == "apr15":
#     pre = "/scratch/wl45/april_second_2019/compare/noClean/"
#     if args.mode == 1:
#         import glob
#         # a = glob.glob("/scratch/wl45/april_second_2019/optimization_k_relative/database/original_pdbs/*.pdb")
#         # a = glob.glob("/scratch/wl45/april_second_2019/compare/database/old/*.pdb")
#         a = glob.glob("/scratch/wl45/april_second_2019/compare/database/cleaned_pdbs/*.pdb")
#         do(f"mkdir -p {pre}/proteins_name_list")
#         do(f"mkdir -p {pre}/decoys/shifted")
#         do(f"mkdir -p {pre}/gammas")
#         print(len(a))
#         random.seed(1)
#         # subset = random.sample(a, 100)
#         # subset = a
#         subset = a[:10]
#         offset_all = np.arange(-40, 40, 5)
#         print(subset)
#         # copy pdbs
#         for one in subset:
#             name = one.split("/")[-1][:4]
#             # source = "/scratch/wl45/april_second_2019/optimization_k_relative/database/original_pdbs/"
#             source = "/scratch/wl45/april_second_2019/compare/database/cleaned_pdbs/"
#             to = "/scratch/wl45/april_second_2019/compare/database/dompdb/"
#             do(f"cp {source}/{name}.pdb {to}/ ")
#         # create protein name list
#         with open(f"{pre}/proteins_name_list/proteins_name_list.txt", "w") as out:
#             for one in subset:
#                 name = one.split("/")[-1][:4]
#                 out.write(f"{name}\n")
#         # generate decoys file and decoys.
#         for one in subset:
#             name = one.split("/")[-1][:4]
#             print(name)
#             os.system(f"mkdir -p {pre}/decoyData/{name}")
#             for offset in offset_all:
#                 duplicate_pdb(f"{to}/{name}.pdb", f"{pre}/decoyData/{name}/{name}_{offset}.pdb", offset_z=offset, new_chain="A")

#             with open(f"{pre}/decoys/shifted/{name}.decoys", "w") as out:
#                 for offset in offset_all:
#                     out.write(f"decoyData/{name}/{name}_{offset}\n")

#             # cleanPdb([f"{name}"], chain="A", formatName=False)
#         do("python3 ~/opt/compute_phis.py -m 5 proteins_name_list/proteins_name_list.txt ")
#     if args.mode == 2:
#         from pyCodeLib import *
#         import warnings
#         warnings.filterwarnings('ignore')
#         complete_proteins = "proteins_name_list/proteins_name_list.txt"
#         A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='shifted',
#                                 num_decoys=20, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False)

#     if args.mode == 3:
#         name = "1kn9"
#         name = "1e6c"
#         # cleanPdb([f"{name}"], chain="A", formatName=False)
#         cleanPdb([f"{name}"], chain="first", formatName=False, toFolder=f"{pre}/../database/cleaned_pdbs/", source=f"{pre}/../database/original_pdbs/")


# if args.day == "apr13":
#     # multi chain iter 0
#     # data = pd.read_csv("chosen.csv", index_col=0)
#     data = pd.read_csv("data_info_3.csv", index_col=0).query("Problematic != 4")
#     if args.mode == 1:
#         do("cp ../phi_list.txt .")
#         do("mkdir proteins_name_list")
#         do("mkdir slurms")
#         do("mkdir -p decoys")
#         do("mkdir -p gammas")
#         do("mkdir -p phis")
#     if args.mode == 2:
#         do("mkdir -p decoys/multiShuffle")
#         for i, line in data.iterrows():
#             fullName = line["FullName"]
#             generate_multiShuffle(fullName, num_decoys=1000)
#         jobIdList = []
#         for i in range(362):
#             with open(f"slurms/run_{i}.slurm", "w") as out:
#                 out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 4 proteins_name_list/proteins_name_list_{i}.txt"))
#                 # out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 1 proteins_name_list/proteins_name_list_{name}.txt"))
#             # replace(f"slurms/run_{i}.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=20G")
#             a = getFromTerminal(f"sbatch slurms/run_{i}.slurm")
#             jobId = a.split(" ")[-1].strip()
#             print("!"+jobId+"!")
#             jobIdList.append(jobId)
#         waitForJobs(jobIdList, sleepInterval=300)
#         with open(f"slurms/run_on_scavenge.slurm", "w") as out:
#             out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr13 -m 3"))
#         replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
#         do(f"sbatch slurms/run_on_scavenge.slurm")
#         # for testing toymodel
#     if args.mode == 3:
#         from pyCodeLib import *
#         import warnings
#         warnings.filterwarnings('ignore')
#         # n = len(simulation_location_list)
#         # n = 2
#         # n = len(new_simulation_list)

#         # import time
#         # time.sleep(4000)
#         # with open(f"proteins_name_list.txt", "w") as out:
#         #     for data in ["new", "old"]:
#         #         pdb_list, steps = dataset[data]
#         #         for p in pdb_list:
#         #             name = p.lower()[:4]
#         #             out.write(f"{name}\n")
#         # complete_proteins = "proteins_name_list/proteins_name_list_tiny.txt"
#         complete_proteins = "proteins_name_list/proteins_name_list.txt"
#         # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
#         #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
#         A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='multiShuffle',
#                                         num_decoys=1000, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=True)
#         # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_xl23(complete_proteins, "phi_list.txt", decoy_method='lammps',
#         #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False)
#         # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='lammps',
#         #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False, mode=1, simulation_location_list=simulation_location_list)
#     if args.mode == 4:
#         with open(f"slurms/run_on_scavenge.slurm", "w") as out:
#             out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr13 -m 3"))
#         replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
#         do(f"sbatch slurms/run_on_scavenge.slurm")

if args.day == "apr12":
    pre = "/scratch/wl45/april_2019/database/"
    cluster_folder = "/scratch/wl45/april_2019/cluster_results/"
    simulation_location_list = ["iter3_with_rg", "iter2_with_rg", "iter1_with_rg"]
    pdb_list = dataset["combined"]
    if args.mode == 1:
        jobIdList = []
        for simulation_location in simulation_location_list:
            print(simulation_location)
            for p in pdb_list:
                name = p.lower()[:4]
                print(name)
                group = f"{simulation_location}_{name}"
                with open(f"run_{group}.slurm", "w") as out:
                    out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr12 -m 2 -l {simulation_location}-{name}"))
                # replace(f"run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=5G")
                # replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=40G")
                # do(f"sbatch run_on_scavenge.slurm")
                a = getFromTerminal(f"sbatch run_{group}.slurm")
                jobId = a.split(" ")[-1].strip()
                # print("!"+jobId+"!")
                jobIdList.append(jobId)
        waitForJobs(jobIdList, sleepInterval=100)
    if args.mode == 2:
        print(args.label.split("-"))
        all_data = []
        simulation_location, name = args.label.split("-")
        t1 = f"{cluster_folder}{simulation_location}_{name}"
        data = pd.read_csv(f"{pre}/Q_{simulation_location}_{name}", index_col=0)
        a = data.groupby("Run").apply(lambda grp: grp.nlargest(1, "index")).reset_index(drop=True)
        for i, item in a.iterrows():
            for simulation_location_2 in simulation_location_list:
                t2 = f"{cluster_folder}{simulation_location_2}_{name}"
                data2 = pd.read_csv(f"{pre}/Q_{simulation_location_2}_{name}", index_col=0)
                b = data2.groupby("Run").apply(lambda grp: grp.nlargest(1, "index")).reset_index(drop=True)
                for j, item2 in b.iterrows():
                    p1 = f"{t1}/run{i}.pdb"
                    p2 = f"{t2}/run{j}.pdb"
                    aligned_length, rmsd, tmscore, seqid = get_aligned_info(p1, p2)
                    all_data.append([name, simulation_location, i, simulation_location_2, j, aligned_length, rmsd, tmscore, seqid])
                    # break
        data = pd.DataFrame(all_data, columns=["protein", "sim1", "i", "sim2", "j", "aligned_length", "rmsd", "tmscore", "seqid"])
        data.to_csv(f"{simulation_location}_{name}.csv")
if args.day == "apr07":
    pdb_list = ["2lep"]

    gammaSource = "../complete_gammas/for_simulation"
    # if args.label == "1":
    #     simulation_location = "native"
    #     iteration = "0_normalized"
    #     gamma = f"iteration_{iteration}_gamma.dat"
    #     burial = f"iteration_{iteration}_burial_gamma.dat"
    if args.label == "1":
        # simulation_location = "start_from_unfolded_2"
        simulation_location = "strengthen_beta"
        gamma = f"original_gamma.dat"
        burial = f"original_burial_gamma.dat"
        startFrom = ""
    if args.label == "2":
        simulation_location = "start_from_native_2"
        gamma = f"original_gamma.dat"
        burial = f"original_burial_gamma.dat"
        startFrom = "--start crystal"
    if args.mode == 1:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            # do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            do(f"cp -r all_simulations/{name}/ {simulation_location}/")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            # replace(f"{pre}/{name}_multi.in", "minimize", "#minimize")

            do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
            do(f"cp {gammaSource}/{gamma} {pre}/gamma.dat")
            # print(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
            do(f"cp {pre}/burial_gamma.dat {pre}/original_burial_gamma.dat")
            do(f"cp {gammaSource}/{burial} {pre}/burial_gamma.dat")
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            # steps = 40
            # steps = 2
            steps = 80
            cd(name)
            # do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            # do(f"run.py -n 20 {name} --commons 2 -s {steps} --runs 1")
            do(f"run.py -n 20 {name} --commons 2 -s {steps} --runs 1 {startFrom}")
            # do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")
        cd("..")


if args.day == "apr02":
    pdb_list = dataset["combined"]

    gammaSource = "../complete_gammas/for_simulation"
    if args.label == "1":
        simulation_location = "native"
        iteration = "0_normalized"
        gamma = f"iteration_{iteration}_gamma.dat"
        burial = f"iteration_{iteration}_burial_gamma.dat"
    if args.label == "2":
        simulation_location = "compare_native_with_original_gamma"
        gamma = f"original_gamma.dat"
        burial = f"original_burial_gamma.dat"
    if args.mode == 1:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            # do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            do(f"cp -r all_simulations/{name}/ {simulation_location}/")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            # replace(f"{pre}/{name}_multi.in", "minimize", "#minimize")

            do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
            do(f"cp {gammaSource}/{gamma} {pre}/gamma.dat")
            # print(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
            do(f"cp {pre}/burial_gamma.dat {pre}/original_burial_gamma.dat")
            do(f"cp {gammaSource}/{burial} {pre}/burial_gamma.dat")
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            # steps = 40
            steps = 2
            # steps = 80
            cd(name)
            # do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            do(f"run.py -n 1 {name} --commons 2 -s {steps} --runs 1 --start crystal")
            # do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")
        cd("..")
    if args.mode == 2:
        # cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            cd(f"{name}/simulation/0/")
            do("cp energy.log back_energy.log")
            do(f"rerun.py {name} -t")
            cd("rerun")
            do(f"cp ../{name}.seq .")
            do(f"python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie {name}.seq")
            cd("..")
            cd(f"../../../")
    if args.mode == 3:
        splitPDB("./", "movie.pdb")

if args.day == "apr01":
    print(args.day)
    pdb_list_dic = {"../iterative_optimization_old_set":"old",
                    "../iterative_optimization_new_temp_range":"new",
                    "../iterative_optimization_biased_sampling":"new"}
    pdb_list_dic_rev = {"old":"iterative_optimization_old_set",
                            "new":"iterative_optimization_new_temp_range"}

    iteration_source_dic = {"bias_2":"../iterative_optimization_biased_sampling",
                            "bias_old_gamma":"../iterative_optimization_biased_sampling",
                            "iter1_with_bias_96percent":"../iterative_optimization_new_temp_range",
                            "iter1_with_bias_98percent":"../iterative_optimization_new_temp_range",
                            "new_iter1_0":"../iterative_optimization_new_temp_range",
                            "new_iter1_90":"../iterative_optimization_new_temp_range",
                            "new_iter1_96":"../iterative_optimization_new_temp_range",
                            "new_iter1_98":"../iterative_optimization_new_temp_range",
                            "new_iter1_combined_on_B":"../iterative_optimization_new_temp_range",
                            "new_iter2_8":"../iterative_optimization_new_temp_range",
                            "new_iter2_10":"../iterative_optimization_new_temp_range",
                            "old_new_iter2_8":"../iterative_optimization_new_temp_range",
                            "new_iter3_10":"../iterative_optimization_old_set",
                            "single":"../iterative_optimization_old_set",
                            "iter4_30":"../iterative_optimization_old_set",
                            "iter4_6":"../iterative_optimization_old_set",
                            "iter4_13":"../iterative_optimization_old_set",
                            "iter5_30":"../iterative_optimization_old_set",
                            "iter6_30":"../iterative_optimization_old_set",
                            "noFrag":"../iterative_optimization_old_set",
                            "iter0_normalized_noFrag":"../iterative_optimization_combined_train_set",
                            "iter1_normalized_noFrag":"../iterative_optimization_combined_train_set",
                            "iter2_normalized_noFrag":"../iterative_optimization_combined_train_set",
                            "iter3_normalized_noFrag":"../iterative_optimization_combined_train_set",
                            "iter4_normalized_noFrag":"../iterative_optimization_combined_train_set",
                            "iter3_normalized_noFrag_90":"../iterative_optimization_combined_train_set",
                            "iter5_normalized_noFrag":"../iterative_optimization_combined_train_set_with_frag",
                            "original":"../iterative_optimization_combined_train_set_with_frag",
                            "iter6_normalized_noFrag":"../iterative_optimization_combined_train_set_with_frag",
                            "iter0":"../iterative_optimization_combined_train_set_with_frag",
                            "without_contact":"../iterative_optimization_combined_train_set_with_frag",
                            "original_with_rg":"../iterative_optimization_combined_train_set_with_frag",
                            "iter1_with_rg":"../iterative_optimization_combined_train_set_with_frag",
                            "iter6_with_rg":"../iterative_optimization_combined_train_set_with_frag",
                            "iter2_with_rg":"../iterative_optimization_combined_train_set_with_frag",
                            "iter3_with_rg":"../iterative_optimization_combined_train_set_with_frag",
                            "iter2_with_rg_90":"../iterative_optimization_combined_train_set_with_frag",
                            "iter3_with_rg_less_frag":"../iterative_optimization_combined_train_set_with_frag",
                            "multi_iter0":"../iterative_optimization_combined_train_set_with_frag",
                            "multi_iter1":"../iterative_optimization_combined_train_set_with_frag",
                            "multi_groupedNorm":"../iterative_optimization_combined_train_set_with_frag",
                            "multi_iter2":"../iterative_optimization_combined_train_set_with_frag",
                            "multi_groupedNorm_check":"../iterative_optimization_combined_train_set_with_frag",
                            "multi_iter1_correct":"../iterative_optimization_combined_train_set_with_frag",
                            "multi_iter1_correct_30":"../iterative_optimization_combined_train_set_with_frag",
                            "multi_group_iter1":"../iterative_optimization_combined_train_set_with_frag",


                            }
    pdb_list_dic = {"../iterative_optimization_old_set":"old",
                    "../iterative_optimization_new_temp_range":"new",
                    "../iterative_optimization_biased_sampling":"new",
                    "../iterative_optimization_combined_train_set":"combined",
                    "../iterative_optimization_combined_train_set_with_frag":"combined"}
    # new_simulation_list = ["iter1_with_bias_96percent", "new_iter2_10"]
    # old_protein_simulation_list = ["single", "new_iter3_10"]

    # new_simulation_list = ["bias_2","bias_old_gamma", "iter1_with_bias_96percent", "iter1_with_bias_98percent", "new_iter2_10", "new_iter1_90", "new_iter2_8", "old_new_iter2_8"]
    # old_protein_simulation_list = ["noFrag", "iter6_30", "iter5_30", "single", "new_iter3_10", "iter4_30", "iter4_6", "iter4_13"]
    # combined_simulation_list = ["iter5_normalized_noFrag", "original", "iter0_normalized_noFrag", "iter1_normalized_noFrag", "iter2_normalized_noFrag", "iter3_normalized_noFrag", "iter4_normalized_noFrag", "iter3_normalized_noFrag_90"]
    # new_data = ["iter5_normalized_noFrag", "original"]


    new_simulation_list = []
    old_protein_simulation_list = []
    combined_simulation_list = ["multi_iter1_correct", "multi_iter1_correct_30", "multi_group_iter1", "multi_iter2", "multi_groupedNorm_check", "multi_iter1", "iter3_with_rg_less_frag", "multi_iter0", "iter3_with_rg", "iter2_with_rg_90", "iter2_with_rg", "iter0", "without_contact", "original_with_rg", "iter1_with_rg", "iter6_with_rg"]
    combined_simulation_list = ["multi_groupedNorm_check", "multi_iter1", "iter3_with_rg_less_frag", "multi_iter0", "iter3_with_rg", "iter2_with_rg_90", "iter2_with_rg", "iter0", "without_contact", "original_with_rg", "iter1_with_rg", "iter6_with_rg"]
    # combined_simulation_list = ["multi_iter2"]
    combined_simulation_list = ["multi_iter0"]
    # multi_iter2 will be special

    # new_data = ["multi_iter1"]
    # combined_simulation_list = ["multi_groupedNorm", "multi_iter2", "multi_groupedNorm_check"]
    # combined_simulation_list = ["multi_groupedNorm"]
    new_data = ["multi_iter1_correct", "multi_iter1_correct_30", "multi_group_iter1"]

    simulation_location_list_dic = defaultdict(list)
    for p in dataset["new"]:
        name = p.lower()[:4]
        simulation_location_list_dic[name] += new_simulation_list
    for p in dataset["old"]:
        name = p.lower()[:4]
        simulation_location_list_dic[name] += old_protein_simulation_list
    for p in dataset["combined"]:
        name = p.lower()[:4]
        simulation_location_list_dic[name] += combined_simulation_list

    simulation_location_list = new_data
    simulation_location_list = []
    cwd = os.getcwd()
    print(cwd)
    Run = 30
    decoy_n = 1000
    if args.mode == 1:
        do("mkdir -p proteins_name_list")
        do("mkdir -p slurms")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p phis")
        if not os.path.exists("phi_list.txt"):
            do("cp ../phi_list.txt .")
        # Convert Dump File to Pdbs
        jobIdList = []
        for simulation_location in simulation_location_list:
            print(simulation_location)
            # cd(iteration_source_dic[simulation_location]+"/"+simulation_location)
            iteration_source = iteration_source_dic[simulation_location]
            pdb_list = dataset[pdb_list_dic[iteration_source]]
            cd(iteration_source+"/"+simulation_location)
            for p in pdb_list:
                name = p.lower()[:4]
                print(name)
                cd(f"{name}/simulation/")
                with open(f"run_on_scavenge.slurm", "w") as out:
                    out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d mar30 -m 33 -l {name}"))
                # replace(f"run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=5G")
                # replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=40G")
                # do(f"sbatch run_on_scavenge.slurm")
                a = getFromTerminal(f"sbatch run_on_scavenge.slurm")
                jobId = a.split(" ")[-1].strip()
                # print("!"+jobId+"!")
                jobIdList.append(jobId)
                cd("../../")
            cd("..")
        waitForJobs(jobIdList, sleepInterval=100)


        # simulation_location_list = combined_simulation_list

        # Transport Pdbs to database folder
        print("Transporting Pdbs to database")
        jobIdList = []
        for simulation_location in simulation_location_list:
            print(simulation_location)
            # cd(iteration_source_dic[simulation_location]+"/"+simulation_location)
            iteration_source = iteration_source_dic[simulation_location]
            pdb_list = dataset[pdb_list_dic[iteration_source]]
            cd(iteration_source+"/"+simulation_location)
            for p in pdb_list:
                name = p.lower()[:4]
                print(name)
                cd(f"{name}/simulation/")
                with open(f"run_on_scavenge.slurm", "w") as out:
                    out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d mar30 -m 44 -l {simulation_location}_{name}"))
                replace(f"run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=5G")
                # replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=40G")
                # do(f"sbatch run_on_scavenge.slurm")
                a = getFromTerminal(f"sbatch run_on_scavenge.slurm")
                jobId = a.split(" ")[-1].strip()
                # print("!"+jobId+"!")
                jobIdList.append(jobId)
                cd("../../")
            cd("..")
        waitForJobs(jobIdList, sleepInterval=100)

        # process Qw info.
        print("Processing Qw files")
        # print(len(simulation_location_list))
        cd(cwd)
        os.system("pwd")
        for simulation_location in simulation_location_list:
            print(simulation_location)
            iteration_source = iteration_source_dic[simulation_location]
            pdb_list = dataset[pdb_list_dic[iteration_source]]
            for p in pdb_list:
                name = p.lower()[:4]
                complete_Q = []
                for i in range(Run):
                    if not os.path.exists(f"../database/{simulation_location}_{name}_{i}/dump.lammpstrj"):
                        print(f"File not exist {simulation_location}_{name}_{i}/dump.lammpstrj")
                        continue
                    try:
                        Q = pd.read_csv(f"../database/{simulation_location}_{name}_{i}/wham.dat")[" Qw"]
                    except:
                        print(f"File not exist {simulation_location}_{name}_{i}/wham.dat")
                        continue
                    n = len(Q)
                    if n != 2000 and n != 1000 and n != 750:
                        print(f"database/{simulation_location}_{name}_{i},  {n}")
                    if n == 0:
                        # print("!! Using previous i.")
                        print(f"An error at {i}")
                        continue

                    complete_Q.append(pd.DataFrame(Q).assign(Run=i))
                data = pd.concat(complete_Q).reset_index(drop=False)
                data["Rank"] = data["index"].rank(ascending=False)
                data.to_csv(f"../database/Q_{simulation_location}_{name}")

        print("generating decoys file")
        cd(cwd)
        os.system("pwd")
        do("mkdir -p decoys/lammps")

        # generate decoy files
        simulation_location_list = combined_simulation_list
        # simulation_location_list = new_data

        print(len(simulation_location_list))
        for simulation_location in simulation_location_list:
            print(simulation_location)
            iteration_source = iteration_source_dic[simulation_location]
            pdb_list = dataset[pdb_list_dic[iteration_source]]
            for p in pdb_list:
                name = p.lower()[:4]
                a = pd.read_csv(f"../database/Q_{simulation_location}_{name}", index_col=0).query(f"Rank < {decoy_n*3}")
                sampled = a.sample(decoy_n)
                with open(f"decoys/lammps/{name}_{simulation_location}.decoys", "w") as out:
                    for i, item in sampled.iterrows():
                        out.write(f"../database/{simulation_location}_{name}_{int(item['Run'])}/frame{int(item['index'])} {np.round(item[' Qw'], 3)}\n")


        # Compute Phis
        print("Computing Phis")
        cd(cwd)
        os.system("pwd")
        jobIdList = []
        for simulation_location in simulation_location_list:
            iteration_source = iteration_source_dic[simulation_location]
            pdb_list = dataset[pdb_list_dic[iteration_source]]
            for p in pdb_list:
                name = p.lower()[:4]
                with open(f"proteins_name_list/proteins_name_list_{name}_{simulation_location}.txt", "w") as out:
                        out.write(f"{name}_{simulation_location}\n")
                with open(f"slurms/run_{name}_{simulation_location}.slurm", "w") as out:
                    out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 3 proteins_name_list/proteins_name_list_{name}_{simulation_location}.txt"))
                    # out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 1 proteins_name_list/proteins_name_list_{name}.txt"))
                replace(f"slurms/run_{name}_{simulation_location}.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=10G")
                # do(f"sbatch slurms/run_{name}_{simulation_location}.slurm")
                a = getFromTerminal(f"sbatch slurms/run_{name}_{simulation_location}.slurm")
                jobId = a.split(" ")[-1].strip()
                jobIdList.append(jobId)
        waitForJobs(jobIdList, sleepInterval=200)

        # exit()
        # # do("cp phis/* ../phis/")
        # with open(f"slurms/run_on_scavenge.slurm", "w") as out:
        #     out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr01 -m 2"))
        # replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        # do(f"sbatch slurms/run_on_scavenge.slurm")
        runFile = f"slurms/run_on_scavenge_2.slurm"
        with open(runFile, "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr01 -m 5 -l phi_list_relative_k.txt"))
        replace(runFile, "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch {runFile}")

        runFile = f"slurms/run_on_scavenge_1.slurm"
        with open(runFile, "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr01 -m 5 -l phi_list_contact.txt"))
        replace(runFile, "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch {runFile}")

    if args.mode == 2:
        # for testing toymodel
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # n = len(simulation_location_list)
        # n = 2
        n = len(new_simulation_list) + len(combined_simulation_list)
        # import time
        # time.sleep(4000)
        with open(f"proteins_name_list.txt", "w") as out:
            for data in ["new", "old"]:
                pdb_list = dataset[data]
                for p in pdb_list:
                    name = p.lower()[:4]
                    out.write(f"{name}\n")
        complete_proteins = "proteins_name_list.txt"

        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='lammps',
                                        num_decoys=n*decoy_n, noise_filtering=True, jackhmmer=False, read=False, mode=2, withBiased=False, simulation_location_list_dic=simulation_location_list_dic)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_xl23(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False, mode=1, simulation_location_list=simulation_location_list)
    if args.mode == 3:
        do("mkdir -p proteins_name_list")
        do("mkdir -p slurms")
        do("mkdir -p decoys/lammps")
        do("mkdir -p gammas")
        do("mkdir -p phis")
        do("cp ../phi_list.txt .")
        simulation_location_list = combined_simulation_list

        print(len(simulation_location_list))
        for simulation_location in simulation_location_list:
            print(simulation_location)
            iteration_source = iteration_source_dic[simulation_location]
            pdb_list = dataset[pdb_list_dic[iteration_source]]
            for p in pdb_list:
                name = p.lower()[:4]
                a = pd.read_csv(f"../database/Q_{simulation_location}_{name}", index_col=0).query(f"Rank < {decoy_n*3}")
                sampled = a.sample(decoy_n)
                with open(f"decoys/lammps/{name}_{simulation_location}.decoys", "w") as out:
                    for i, item in sampled.iterrows():
                        out.write(f"../database/{simulation_location}_{name}_{int(item['Run'])}/frame{int(item['index'])} {np.round(item[' Qw'], 3)}\n")


        # Compute Phis
        print("Computing Phis")
        cd(cwd)
        os.system("pwd")
        jobIdList = []
        for simulation_location in simulation_location_list:
            iteration_source = iteration_source_dic[simulation_location]
            pdb_list = dataset[pdb_list_dic[iteration_source]]
            for p in pdb_list:
                name = p.lower()[:4]
                with open(f"proteins_name_list/proteins_name_list_{name}_{simulation_location}.txt", "w") as out:
                        out.write(f"{name}_{simulation_location}\n")
                with open(f"slurms/run_{name}_{simulation_location}.slurm", "w") as out:
                    out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 3 proteins_name_list/proteins_name_list_{name}_{simulation_location}.txt"))
                    # out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 1 proteins_name_list/proteins_name_list_{name}.txt"))
                replace(f"slurms/run_{name}_{simulation_location}.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=10G")
                # do(f"sbatch slurms/run_{name}_{simulation_location}.slurm")
                a = getFromTerminal(f"sbatch slurms/run_{name}_{simulation_location}.slurm")
                jobId = a.split(" ")[-1].strip()
                jobIdList.append(jobId)
        waitForJobs(jobIdList, sleepInterval=200)

        # # do("cp phis/* ../phis/")
        # with open(f"slurms/run_on_scavenge.slurm", "w") as out:
        #     out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr01 -m 2"))
        # replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        # do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 4:
        runFile = f"slurms/run_on_scavenge_2.slurm"
        with open(runFile, "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr01 -m 5 -l phi_list_relative_k.txt"))
        replace(runFile, "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch {runFile}")

        runFile = f"slurms/run_on_scavenge_1.slurm"
        with open(runFile, "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d apr01 -m 5 -l phi_list_contact.txt"))
        replace(runFile, "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch {runFile}")
    if args.mode == 5:
        # for testing toymodel
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # n = len(simulation_location_list)
        # n = 2
        n = len(new_simulation_list) + len(combined_simulation_list)
        phi_file = args.label
        # import time
        # time.sleep(4000)
        with open(f"proteins_name_list.txt", "w") as out:
            for data in ["new", "old"]:
                pdb_list = dataset[data]
                for p in pdb_list:
                    name = p.lower()[:4]
                    out.write(f"{name}\n")
        complete_proteins = "proteins_name_list.txt"

        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, phi_file, decoy_method='lammps',
                                        num_decoys=n*decoy_n, noise_filtering=True, jackhmmer=False, read=False, mode=2, withBiased=False, simulation_location_list_dic=simulation_location_list_dic)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_xl23(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False, mode=1, simulation_location_list=simulation_location_list)
if args.day == "mar31":
    # pdb_list, steps = dataset["old"]
    pdb_list = dataset["combined"]

    # simulation_location = "inverseBurial"
    gammaSource = "../complete_gammas/for_simulation"
    # gammaSource = "../optimization_iter1"

    if args.label == "1":
        simulation_location = "iter0_normalized_noFrag"
        iteration = "0_normalized"
        gamma = f"iteration_{iteration}_gamma.dat"
        burial = f"iteration_{iteration}_burial_gamma.dat"

    if args.label == "2":
        simulation_location = "iter1_normalized_noFrag"
        iteration = "iter_1"
        percent = 30
        gamma = f"iteration_{iteration}_gamma_{percent}.dat"
        burial = f"iteration_{iteration}_burial_gamma_{percent}.dat"

    if args.label == "3":
        simulation_location = "original"
        gamma = f"original_gamma.dat"
        burial = f"original_burial_gamma.dat"

    if args.label == "4":
        i = 2
        simulation_location = f"iter{i}_normalized_noFrag"
        iteration = f"iter_{i}"
        percent = 30
        gamma = f"iteration_{iteration}_gamma_{percent}.dat"
        burial = f"iteration_{iteration}_burial_gamma_{percent}.dat"

    if args.label == "5":
        i = 3
        simulation_location = f"iter{i}_normalized_noFrag"
        iteration = f"iter_{i}"
        percent = 30
        gamma = f"iteration_{iteration}_gamma_{percent}.dat"
        burial = f"iteration_{iteration}_burial_gamma_{percent}.dat"

    if args.label == "6":
        i = 3
        iteration = f"iter_{i}"
        percent = 90
        simulation_location = f"iter{i}_normalized_noFrag_{percent}"
        gamma = f"iteration_{iteration}_gamma_{percent}.dat"
        burial = f"iteration_{iteration}_burial_gamma_{percent}.dat"

    if args.label == "7":
        i = 4
        simulation_location = f"iter{i}_normalized_noFrag"
        iteration = f"iter_{i}"
        percent = 30
        gamma = f"iteration_{iteration}_gamma_{percent}.dat"
        burial = f"iteration_{iteration}_burial_gamma_{percent}.dat"

    if args.label == "8":
        i = 5
        simulation_location = f"iter{i}_normalized_noFrag"
        iteration = f"iter_{i}"
        percent = 30
        gamma = f"iteration_{iteration}_gamma_{percent}.dat"
        burial = f"iteration_{iteration}_burial_gamma_{percent}.dat"

    if args.label == "9":
        i = 6
        simulation_location = f"iter{i}_normalized_noFrag"
        iteration = f"iter_{i}"
        percent = 30
        gamma = f"iteration_{iteration}_gamma_{percent}.dat"
        burial = f"iteration_{iteration}_burial_gamma_{percent}.dat"

    if args.label == "10":
        i = "7"
        simulation_location = f"iter{i}_normalized"
        iteration = f"iter_{i}"
        percent = 30
        gamma = f"iteration_{iteration}_gamma_{percent}.dat"
        burial = f"iteration_{iteration}_burial_gamma_{percent}.dat"

    if args.label == "11":
        i = "7"
        simulation_location = f"iter{i}_normalized_90"
        iteration = f"iter_{i}"
        percent = 90
        gamma = f"iteration_{iteration}_gamma_{percent}.dat"
        burial = f"iteration_{iteration}_burial_gamma_{percent}.dat"

    if args.mode == 1:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            # do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            do(f"cp -r all_simulations/{name}/ {simulation_location}/")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            # if simulation_location == "newContactWithBurial":

            do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
            do(f"cp {gammaSource}/{gamma} {pre}/gamma.dat")
            # print(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
            do(f"cp {pre}/burial_gamma.dat {pre}/original_burial_gamma.dat")
            do(f"cp {gammaSource}/{burial} {pre}/burial_gamma.dat")

    if args.mode == 2:
        cd(simulation_location)
        for p in pdb_list:
            steps = returnSteps(p)
            name = p.lower()[:4]
            # steps = 40
            cd(name)
            do(f"run.py -n 10 {name} --commons 2 -s {steps} --runs 2")
            cd("..")

if args.day == "mar30":
    print(args.day)
    # dataset = {"old":("1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "), 40),
    #             "new":("1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "), 80),
    #             "test":(['t089', 't120', 't251', 'top7', '1ubq', 't0766', 't0778', 't0782',
    #                     't0792', 't0803', 't0815', 't0833', 't0842', 't0844'], 40)}

    # , 't0846
    pdb_list_dic = {"../iterative_optimization_old_set":"old",
                    "../iterative_optimization_new_temp_range":"new",
                    "../iterative_optimization_biased_sampling":"new"}
    pdb_list_dic_rev = {"old":"iterative_optimization_old_set",
                            "new":"iterative_optimization_new_temp_range"}

    iteration_source_dic = {"bias_2":"../iterative_optimization_biased_sampling",
                            "bias_old_gamma":"../iterative_optimization_biased_sampling",
                            "iter1_with_bias_96percent":"../iterative_optimization_new_temp_range",
                            "iter1_with_bias_98percent":"../iterative_optimization_new_temp_range",
                            "new_iter1_0":"../iterative_optimization_new_temp_range",
                            "new_iter1_90":"../iterative_optimization_new_temp_range",
                            "new_iter1_96":"../iterative_optimization_new_temp_range",
                            "new_iter1_98":"../iterative_optimization_new_temp_range",
                            "new_iter1_combined_on_B":"../iterative_optimization_new_temp_range",
                            "new_iter2_8":"../iterative_optimization_new_temp_range",
                            "new_iter2_10":"../iterative_optimization_new_temp_range",
                            "old_new_iter2_8":"../iterative_optimization_new_temp_range",
                            "new_iter3_10":"../iterative_optimization_old_set",
                            "single":"../iterative_optimization_old_set",
                            "iter4_30":"../iterative_optimization_old_set",
                            "iter4_6":"../iterative_optimization_old_set",
                            "iter4_13":"../iterative_optimization_old_set",
                            "iter5_30":"../iterative_optimization_old_set",
                            "iter6_30":"../iterative_optimization_old_set",
                            "noFrag":"../iterative_optimization_old_set",
                            "iter0_normalized_noFrag":"../iterative_optimization_combined_train_set",
                            "iter1_normalized_noFrag":"../iterative_optimization_combined_train_set",
                            "iter1_normalized_noFrag":"../iterative_optimization_combined_train_set",
                            }
    pdb_list_dic = {"../iterative_optimization_old_set":"old",
                    "../iterative_optimization_new_temp_range":"new",
                    "../iterative_optimization_biased_sampling":"new",
                    "../iterative_optimization_combined_train_set":"combined"}
    # new_simulation_list = ["iter1_with_bias_96percent", "new_iter2_10"]
    # old_protein_simulation_list = ["single", "new_iter3_10"]
    new_simulation_list = ["bias_2","bias_old_gamma", "iter1_with_bias_96percent", "iter1_with_bias_98percent", "new_iter2_10", "new_iter1_90", "new_iter2_8", "old_new_iter2_8"]
    old_protein_simulation_list = ["noFrag", "iter6_30", "iter5_30", "single", "new_iter3_10", "iter4_30", "iter4_6", "iter4_13"]
    combined_simulation_list = ["iter0_normalized_noFrag", "iter1_normalized_noFrag"]


    simulation_location_list_dic = defaultdict(list)
    for p in dataset["new"]:
        name = p.lower()[:4]
        simulation_location_list_dic[name] += new_simulation_list
    for p in dataset["old"]:
        name = p.lower()[:4]
        simulation_location_list_dic[name] += old_protein_simulation_list
    for p in dataset["combined"]:
        name = p.lower()[:4]
        simulation_location_list_dic[name] += combined_simulation_list


    if args.mode == 1:
        for d_set in ["old", "new"]:
            pdb_list, steps = dataset[d_set]
            native_source = pdb_list_dic_rev[d_set]
            # native_source = "iterative_optimization_old_set"
            do("mkdir -p database/dompdb")
            for p in pdb_list:
                name = p.lower()[:4]
                do(f"mkdir -p database/{name}")
                do(f"cp -r {native_source}/native/{name}/simulation/0/rerun/* database/{name}/")
                pre = f"database/{name}/"
                fileName = "movie.pdb"
                splitPDB(pre, fileName)
                do(f"cp database/{name}/frame0.pdb database/dompdb/{name}.pdb")
    if args.mode == 2:
        do("mkdir proteins_name_list")
        do("mkdir slurms")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p phis")
        do("cp ../phi_list.txt .")
    if args.mode == 3:
        # simulation_location_list = ["single", "new_iter3_10"]
        # simulation_location_list = ["iter4_30", "iter4_13", "iter4_6"]
        simulation_location_list = ["iter5_30"]
        # simulation_location_list = ["iter6_30"]
        simulation_location_list = ["noFrag"]
        simulation_location_list = ["iter0_normalized_noFrag", "iter1_normalized_noFrag"]
        # Go To list of simulation folder
        for simulation_location in simulation_location_list:
            print(simulation_location)
            # cd(iteration_source_dic[simulation_location]+"/"+simulation_location)
            iteration_source = iteration_source_dic[simulation_location]
            pdb_list = dataset[pdb_list_dic[iteration_source]]
            cd(iteration_source+"/"+simulation_location)
            for p in pdb_list:
                name = p.lower()[:4]
                print(name)
                cd(f"{name}/simulation/")
                with open(f"run_on_scavenge.slurm", "w") as out:
                    out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d mar30 -m 333 -l {name}"))
                # replace(f"run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=5G")
                # replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=40G")
                do(f"sbatch run_on_scavenge.slurm")
                cd("../../")
            cd("..")
        # for simulation_location in simulation_location_list:
        #     print(simulation_location)
        #     # cd(iteration_source_dic[simulation_location]+"/"+simulation_location)
        #     iteration_source = iteration_source_dic[simulation_location]
        #     pdb_list, steps = dataset[pdb_list_dic[iteration_source]]
        #     cd(iteration_source+"/"+simulation_location)
        #     for p in pdb_list:
        #         name = p.lower()[:4]
        #         print(name)
        #         cd(f"{name}/simulation/")
        #         with open(f"run_on_scavenge.slurm", "w") as out:
        #             out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d mar30 -m 33 -l {name}"))
        #         # replace(f"run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=5G")
        #         # replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=40G")
        #         do(f"sbatch run_on_scavenge.slurm")
        #         cd("../../")
        #     cd("..")
    if args.mode == 33:
        name = args.label
        for i in range(30):
            cd(f"{i}/0/")
            do(f"cp ../{name}.seq .")
            do(f"python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie {name}.seq")
            cd("../..")
    if args.mode == 333:
        name = args.label
        rerun = 1
        for i in range(10):
            cd(f"{i}/{rerun}/")
            do(f"cp ../{name}.seq .")
            do(f"python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie {name}.seq")
            cd("../..")
    if args.mode == 4:
        simulation_location_list = new_simulation_list + old_protein_simulation_list
        simulation_location_list = list(iteration_source_dic.keys())
        simulation_location_list = ["noFrag"]
        simulation_location_list = ["iter0_normalized_noFrag", "iter1_normalized_noFrag"]
        for simulation_location in simulation_location_list:
            print(simulation_location)
            # cd(iteration_source_dic[simulation_location]+"/"+simulation_location)
            iteration_source = iteration_source_dic[simulation_location]
            pdb_list = dataset[pdb_list_dic[iteration_source]]
            cd(iteration_source+"/"+simulation_location)
            for p in pdb_list:
                name = p.lower()[:4]
                print(name)
                cd(f"{name}/simulation/")
                with open(f"run_on_scavenge.slurm", "w") as out:
                    out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d mar30 -m 444 -l {simulation_location}_{name}"))
                replace(f"run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=5G")
                # replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=40G")
                do(f"sbatch run_on_scavenge.slurm")
                cd("../../")
            cd("..")
        # for simulation_location in simulation_location_list:
        #     print(simulation_location)
        #     # cd(iteration_source_dic[simulation_location]+"/"+simulation_location)
        #     iteration_source = iteration_source_dic[simulation_location]
        #     pdb_list, steps = dataset[pdb_list_dic[iteration_source]]
        #     cd(iteration_source+"/"+simulation_location)
        #     for p in pdb_list:
        #         name = p.lower()[:4]
        #         print(name)
        #         cd(f"{name}/simulation/")
        #         with open(f"run_on_scavenge.slurm", "w") as out:
        #             out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d mar30 -m 44 -l {simulation_location}_{name}"))
        #         replace(f"run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=5G")
        #         # replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=40G")
        #         do(f"sbatch run_on_scavenge.slurm")
        #         cd("../../")
        #     cd("..")
    if args.mode == 44:
        cmd = args.label
        for i in range(30):
            pre = f"/scratch/wl45/april_2019/database/{cmd}_{i}/"
            do(f"mkdir -p {pre}")
            do(f"cp -r {i}/0/* {pre}/")
            fileName = "movie.pdb"
            splitPDB(pre, fileName)
            do(f"rm {pre}movie.*")
    if args.mode == 444:
        cmd = args.label
        rerun = 1
        for i in range(10):
            pre = f"/scratch/wl45/april_2019/database/{cmd}_{i}/"
            do(f"cp -r {i}/{rerun}/ {pre}")
            fileName = "movie.pdb"
            splitPDB(pre, fileName)
            do(f"rm {pre}movie.*")
    if args.mode == 5:
        # generate decoy file
        simulation_location_list = new_simulation_list + old_protein_simulation_list
        simulation_location_list = ["iter0_normalized_noFrag", "iter1_normalized_noFrag"]
        # simulation_location_list = ["noFrag"]
        print("generating decoys file")
        do("mkdir -p decoys/lammps")
        print(len(simulation_location_list))
        Run = 10
        for simulation_location in simulation_location_list:
            print(simulation_location)
            iteration_source = iteration_source_dic[simulation_location]
            pdb_list = dataset[pdb_list_dic[iteration_source]]
            for p in pdb_list:
                name = p.lower()[:4]
                complete_Q = []
                for i in range(Run):
                    Q = pd.read_csv(f"../database/{simulation_location}_{name}_{i}/wham.dat")[" Qw"]
                    n = len(Q)
                    if n != 2000 and n != 1000 and n != 750:
                        print(f"database/{simulation_location}_{name}_{i},  {n}")
                    if n == 0:
                        # print("!! Using previous i.")
                        print(f"An error at {i}")
                        continue

                    complete_Q.append(pd.DataFrame(Q).assign(Run=i))
                data = pd.concat(complete_Q).reset_index(drop=False)
                data["Rank"] = data["index"].rank(ascending=False)
                data.to_csv(f"../database/Q_{simulation_location}_{name}")
    if args.mode == 6:
        simulation_location_list = new_simulation_list + old_protein_simulation_list
        simulation_location_list += ["iter0_normalized_noFrag", "iter1_normalized_noFrag"]
        print("generating decoys file")
        do("mkdir -p decoys/lammps")
        decoy_n = 500
        print(len(simulation_location_list))
        for simulation_location in simulation_location_list:
            print(simulation_location)
            iteration_source = iteration_source_dic[simulation_location]
            pdb_list = dataset[pdb_list_dic[iteration_source]]
            for p in pdb_list:
                name = p.lower()[:4]
                a = pd.read_csv(f"../database/Q_{simulation_location}_{name}", index_col=0).query(f"Rank < {decoy_n*3}")
                sampled = a.sample(decoy_n)
                with open(f"decoys/lammps/{name}_{simulation_location}.decoys", "w") as out:
                    for i, item in sampled.iterrows():
                        out.write(f"../database/{simulation_location}_{name}_{int(item['Run'])}/frame{int(item['index'])} {np.round(item[' Qw'], 3)}\n")
        # for simulation_location in simulation_location_list:
        #     print(simulation_location)
        #     iteration_source = iteration_source_dic[simulation_location]
        #     pdb_list, steps = dataset[pdb_list_dic[iteration_source]]
        #     for p in pdb_list:
        #         name = p.lower()[:4]
        #         count = 0
        #         a = pd.read_csv(f"../database/Q_{simulation_location}_{name}", index_col=0).query(f"Rank < {decoy_n*3}")
        #         sampled = a.sample(decoy_n)
        #         with open(f"decoys/lammps/{name}_{simulation_location}.decoys", "w") as out:
        #             for i, item in sampled.iterrows():
        #                 out.write(f"../database/{simulation_location}_{name}_{int(item['Run'])}/frame{int(item['index'])} {np.round(item[' Qw'], 3)}\n")
    if args.mode == 7:
        # simulation_location_list = ["single", "new_iter3_10"]
        simulation_location_list = new_simulation_list + old_protein_simulation_list
        simulation_location_list += ["iter0_normalized_noFrag", "iter1_normalized_noFrag"]
        do("mkdir -p slurms")
        do("mkdir -p proteins_name_list")
        for simulation_location in simulation_location_list:
            iteration_source = iteration_source_dic[simulation_location]
            pdb_list = dataset[pdb_list_dic[iteration_source]]
            for p in pdb_list:
                name = p.lower()[:4]
                with open(f"proteins_name_list/proteins_name_list_{name}_{simulation_location}.txt", "w") as out:
                        out.write(f"{name}_{simulation_location}\n")
                with open(f"slurms/run_{name}_{simulation_location}.slurm", "w") as out:
                    out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 3 proteins_name_list/proteins_name_list_{name}_{simulation_location}.txt"))
                    # out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 1 proteins_name_list/proteins_name_list_{name}.txt"))
                replace(f"slurms/run_{name}_{simulation_location}.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=10G")
                do(f"sbatch slurms/run_{name}_{simulation_location}.slurm")
    if args.mode == 8:
        # for testing toymodel
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # n = len(simulation_location_list)
        # n = 2
        n = len(new_simulation_list) + len(combined_simulation_list)
        decoy_n = 500
        # import time
        # time.sleep(4000)
        with open(f"proteins_name_list.txt", "w") as out:
            for data in ["new", "old"]:
                pdb_list = dataset[data]
                for p in pdb_list:
                    name = p.lower()[:4]
                    out.write(f"{name}\n")
        complete_proteins = "proteins_name_list.txt"

        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='lammps',
                                        num_decoys=n*decoy_n, noise_filtering=True, jackhmmer=False, read=False, mode=2, simulation_location_list_dic=simulation_location_list_dic)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_xl23(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False, mode=1, simulation_location_list=simulation_location_list)
    if args.mode == 9:
        # import time
        # time.sleep(7200)
        do("mkdir -p gammas")
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d mar30 -m 8"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")

if args.day == "mar29":
    print(args.day)
    dataset = {"old":("1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "), 40),
                "new":("1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "), 80),
                "test":(['t089', 't120', 't251', 'top7', '1ubq', 't0766', 't0778', 't0782',
                        't0792', 't0803', 't0815', 't0833', 't0842', 't0844', 't0846'], 40)}    # pdb_list, steps = dataset["old"]
    # pdb_list, steps = dataset["test"]
    # pdb_list, steps = dataset["test"]
    pdb_list, steps = dataset["new"]
    # pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    # pdb_list = "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", ")

    simulation_location = f"single"

    # percent = 30
    # iteration = 7
    # simulation_location = f"iter7_2"
    # gammaSource = f"../fix_gamma_mix_error"
    # gamma = f"iteration_test_7_2_gamma_30.dat"
    # burial = f"iteration_test_7_2_burial_gamma_30.dat"
    if args.mode == 1:
        for p in pdb_list:
            name = p.lower()[:4]
            # name = p
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            # do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            do(f"cp -r all_simulations/{name}/ {simulation_location}/")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            # if simulation_location == "newContactWithBurial":
            if simulation_location != "single":
                do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
                do(f"cp {gammaSource}/{gamma} {pre}/gamma.dat")
                # print(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
                do(f"cp {pre}/burial_gamma.dat {pre}/original_burial_gamma.dat")
                do(f"cp {gammaSource}/{burial} {pre}/burial_gamma.dat")

    if args.mode == 2:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            # name = p
            # steps = 40
            cd(name)
            do(f"run.py -n 10 {name} --commons 2 -s {steps} --runs 1")
            cd("..")


if args.day == "mar28":
    if args.mode == 1:
        for i in range(2, 5):
            do(f"rm -r optimization_weighted_by_q_iter{i}/database/")
            do(f"cp -r optimization_weighted_by_q_iter{i}/gammas/ ../april_2019/old_gammas/optimization_weighted_by_q_iter{i}_gamma")
if args.day == "mar27":
    print(args.day)
    dataset = {"old":("1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "), 40),
                "new":("1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "), 80),
                "test":(['1A2J', '1A3H', '1A5Y', '1A8Q', '1AGY', '1AKZ', '1AUZ', '1B1A', '1B31', '1B8X', '1BCO', '1BN6', '1BOH', '1BOI'], 40)}
    pdb_list, steps = dataset["old"]
    # pdb_list, steps = dataset["test"]
    # pdb_list, steps = dataset["new"]

    percent = 30
    # simulation_location = f"withoutContact"
    # simulation_location = f"single"


    # simulation_location = f"iter7_2"
    # gammaSource = f"../fix_gamma_mix_error"
    # gamma = f"iteration_test_7_2_gamma_30.dat"
    # burial = f"iteration_test_7_2_burial_gamma_30.dat"


    simulation_location = f"noFrag"
    # gammaSource = f"../fix_gamma_mix_error"
    # gamma = f"iteration_multiSeq_gamma_30.dat"
    # burial = f"iteration_multiSeq_burial_gamma_30.dat"
    # gamma = f"iteration_gamma_combined_on_B.dat"
    # burial = f"iteration_burial_gamma_combined_on_B.dat"
    if args.mode == 1:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Water\]", "\[Water\]-")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Burial\]", "\[Burial\]-")
            # if simulation_location == "newContactWithBurial":

            # if simulation_location != "single":
            #     do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
            #     do(f"cp {gammaSource}/{gamma} {pre}/gamma.dat")
            #     # print(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
            #     do(f"cp {pre}/burial_gamma.dat {pre}/original_burial_gamma.dat")
            #     do(f"cp {gammaSource}/{burial} {pre}/burial_gamma.dat")

    if args.mode == 2:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            # steps = 40
            cd(name)
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")


if args.day == "mar26":
    # data = pd.read_csv("chosen.csv", index_col=0)
    data = pd.read_csv("data_info_3.csv", index_col=0).query("Problematic != 4")
    if args.mode == 1:
        do("cp ../phi_list.txt .")
        do("mkdir -p database/dompdb")
        do("mkdir -p database/S20_seq")
        do("mkdir proteins_name_list")
        do("mkdir slurms")
        do("mkdir -p subsetgammas")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p phis")
    if args.mode == 2:
        do("mkdir -p decoys/multiShuffle")
        for i, line in data.iterrows():
            fullName = line["FullName"]
            generate_multiShuffle(fullName, num_decoys=1000)
        jobIdList = []
        for i in range(181):
            with open(f"slurms/run_{i}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 4 proteins_name_list/proteins_name_list_{i}.txt"))
                # out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 1 proteins_name_list/proteins_name_list_{name}.txt"))
            # replace(f"slurms/run_{i}.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=20G")
            a = getFromTerminal(f"sbatch slurms/run_{i}.slurm")
            jobId = a.split(" ")[-1].strip()
            print("!"+jobId+"!")
            jobIdList.append(jobId)
        waitForJobs(jobIdList, sleepInterval=300)
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d mar26 -m 3"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
        # for testing toymodel
    if args.mode == 3:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # n = len(simulation_location_list)
        # n = 2
        # n = len(new_simulation_list)

        # import time
        # time.sleep(4000)
        # with open(f"proteins_name_list.txt", "w") as out:
        #     for data in ["new", "old"]:
        #         pdb_list, steps = dataset[data]
        #         for p in pdb_list:
        #             name = p.lower()[:4]
        #             out.write(f"{name}\n")
        # complete_proteins = "proteins_name_list/proteins_name_list_tiny.txt"
        complete_proteins = "proteins_name_list/proteins_name_list.txt"
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
        #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='multiShuffle',
                                        num_decoys=1000, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=True)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_xl23(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False, mode=1, simulation_location_list=simulation_location_list)
    if args.mode == 12:
        # import time
        # time.sleep(7200)
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d mar26 -m 3"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
if args.day == "mar25":
    if args.mode == 1:
        do("mkdir -p decoys/multiShuffle")
        data = pd.read_csv("chosen.csv", index_col=0)
        for i, line in data.iterrows():
            generate_decoy_sequences(line["FullName"], num_decoys=5000)
    if args.mode == 2:
        time.sleep(60)
    if args.mode == 3:
        # a = getFromTerminal("sbatch run_on_scavenge.slurm")
        a = getFromTerminal("gg_server.py -d mar25 -m 4")
        jobId = a.split(" ")[-1].strip()
        print("!"+jobId+"!")
        jobIdList = [jobId]
        previousJobNotFinished = True
        while previousJobNotFinished:
            print("Waiting for previous jobs", datetime.now())
            time.sleep(10)
            previousJobNotFinished = False
            a = getFromTerminal("squeue -u wl45")
            for jobId in jobIdList:
                if jobId in a:
                    previousJobNotFinished = True
        print("Continue Next Script")
    if args.mode == 4:
        do("sbatch run_on_scavenge.slurm")
if args.day == "mar11":
    print(args.day)
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    # simulation_location = "noFrag"
    # simulation_location = "newFrag"
    # simulation_location = "top20"
    simulation_location = "top1"
    # simulation_location = "top3"
    # simulation_location = "fragment"

    if args.mode == 1:
        fragLibrary = f"fragment_memory_{simulation_location}"
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r ../fragment_memory/fraglib {pre}/")
            do(f"cp {fragLibrary}/{name}.mem {pre}/frags.mem")

            # replace(f"{simulation_location}/{name}/{name}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Water\]", "\[Water\]-")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Burial\]", "\[Burial\]-")
            # do(f"rm all_simulations/{name}/{name}/gamma.dat")
    if args.mode == 2:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            steps = 80
            cd(name)
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")

if args.day == "mar09":
    print(args.day)
    dataset = {"old":("1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "), 40),
                "new":("1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "), 80),
                "test":(['1A2J', '1A3H', '1A5Y', '1A8Q', '1AGY', '1AKZ', '1AUZ', '1B1A', '1B31', '1B8X', '1BCO', '1BN6', '1BOH', '1BOI'], 40)}
    # pdb_list, steps = dataset["old"]
    # pdb_list, steps = dataset["test"]
    pdb_list, steps = dataset["new"]
    # pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    # pdb_list = "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", ")

    simulation_location = f"frag"
    if args.mode == 1:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            # if simulation_location == "newContactWithBurial":

    if args.mode == 2:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            # steps = 40
            cd(name)
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")


if args.day == "mar06":
    print(args.day)
    dataset = {"old":("1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "), 40),
                "new":("1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "), 80),
                "test":(['1A2J', '1A3H', '1A5Y', '1A8Q', '1AGY', '1AKZ', '1AUZ', '1B1A', '1B31', '1B8X', '1BCO', '1BN6', '1BOH', '1BOI'], 40)}
    pdb_list, steps = dataset["old"]
    # pdb_list, steps = dataset["test"]
    # pdb_list, steps = dataset["new"]
    submode = 0
    submode = 1
    submode = 2
    submode = 3
    submode = 4
    submode = 5
    submode = 6
    # pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    # pdb_list = "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", ")
    # simulation_location = "noFrag"
    # simulation_location = "newFrag"
    # simulation_location = "newContactWithBurial"
    # simulation_location = "iter1_with_bias_96percent"
    # percent = 90
    # percent = 8
    # percent = 10
    # percent = 96
    # simulation_location = f"new_iter2_{percent}"
    # simulation_location = f"new_iter3_{percent}"
    # simulation_location = f"iter0"
    # simulation_location = f"new_iter1_{percent}"
    # simulation_location = f"new_iter1_combined_on_B"
    # simulation_location = "single"
    # simulation_location = "inverseBurial"
    # gammaSource = "../optimization"
    # gammaSource = "../optimization_iter1"
    # gammaSource = "../optimization_with_biased_iter1"
    # gammaSource = "../optimization_weighted_by_q_iter1"
    # gammaSource = "../optimization_weighted_by_q_iter2"
    # gammaSource = "../optimization_weighted_by_q_iter3"
    # gamma = f"iteration_gamma_{percent}.dat"
    # burial = f"iteration_burial_gamma_{percent}.dat"
    if submode == 0:
        simulation_location = f"iter0"
        gammaSource = "../optimization"
        gamma = f"iteration_gamma.dat"
        burial = f"iteration_burial_gamma.dat"
    elif submode == 1:
        # percent = 13
        # percent = 22
        # percent = 5
        # percent = 6
        percent = 30

        # iteration = 5
        # iteration = 6
        iteration = 7
        simulation_location = f"iter{iteration}_{percent}"
        gammaSource = f"../optimization_weighted_by_q_iter{iteration}"
        gamma = f"iteration_{iteration}_gamma_{percent}.dat"
        burial = f"iteration_{iteration}_burial_gamma_{percent}.dat"
    elif submode == 2:
        percent = 30
        iteration = 7
        simulation_location = f"iter7_2"
        gammaSource = f"../fix_gamma_mix_error"
        gamma = f"iteration_test_7_2_gamma_30.dat"
        burial = f"iteration_test_7_2_burial_gamma_30.dat"
    elif submode == 3:
        percent = 30
        iteration = 7
        simulation_location = f"iter7_compare"
        gammaSource = f"../fix_gamma_mix_error"
        gamma = f"iteration_test_7_gamma_30.dat"
        burial = f"iteration_test_7_burial_gamma_30.dat"
    elif submode == 4:
        percent = 30
        iteration = 7
        simulation_location = f"iter7_fast"
        gammaSource = f"../fix_gamma_mix_error"
        gamma = f"iteration_test_7_3_gamma_30.dat"
        burial = f"iteration_test_7_3_burial_gamma_30.dat"
    elif submode == 5:
        percent = 30
        iteration = 7
        simulation_location = f"iter7_normalized"
        gammaSource = f"../fix_gamma_mix_error"
        gamma = f"iteration_7_normalized_gamma_30.dat"
        burial = f"iteration_7_normalized_burial_gamma_30.dat"
    elif submode == 6:
        percent = 30
        simulation_location = f"multiSeq"
        gammaSource = f"../fix_gamma_mix_error"
        gamma = f"iteration_multiSeq_gamma_30.dat"
        burial = f"iteration_multiSeq_burial_gamma_30.dat"
    # gamma = f"iteration_gamma_combined_on_B.dat"
    # burial = f"iteration_burial_gamma_combined_on_B.dat"
    if args.mode == 1:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            # if simulation_location == "newContactWithBurial":
            if simulation_location != "single":
                do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
                do(f"cp {gammaSource}/{gamma} {pre}/gamma.dat")
                # print(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
                do(f"cp {pre}/burial_gamma.dat {pre}/original_burial_gamma.dat")
                do(f"cp {gammaSource}/{burial} {pre}/burial_gamma.dat")

    if args.mode == 2:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            # steps = 40
            cd(name)
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")


if args.day == "mar05":
    # pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    dataset = {"old":("1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "), 40),
                "new":("1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "), 80),
                "test":(['1A2J', '1A3H', '1A5Y', '1A8Q', '1AGY', '1AKZ', '1AUZ', '1B1A', '1B31', '1B8X', '1BCO', '1BN6', '1BOH', '1BOI'], 40)}
    # optimization_folder = "optimization_weighted_by_q_iter1"
    # native_source = "../iterative_optimization_new_temp_range"
    # gammaSource = "../optimization"
    # simulation_location_list = ["bias_2", "bias_old_gamma", "iter1_with_bias_96percent", "iter1_with_bias_98percent"]
    # iteration_source_dic = {"bias_2":"../iterative_optimization_biased_sampling",
    #                         "bias_old_gamma":"../iterative_optimization_biased_sampling",
    #                         "iter1_with_bias_96percent":"../iterative_optimization_new_temp_range",
    #                         "iter1_with_bias_98percent":"../iterative_optimization_new_temp_range"}
    # optimization_folder = "optimization_weighted_by_q_iter2_improved"
    # optimization_folder = "optimization_weighted_by_q_iter5"
    optimization_folder = "optimization_weighted_by_q_iter6"
    optimization_folder = "optimization_weighted_by_q_iter7"
    gammaSource = "../optimization"
    # simulation_location_list = ["bias_2", "bias_old_gamma", "iter1_with_bias_96percent", "iter1_with_bias_98percent",
    #             "new_iter1_0", "new_iter1_90", "new_iter1_96", "new_iter1_98", "new_iter1_combined_on_B",
    #             "new_iter2_8", "new_iter2_10", "old_new_iter2_8", "single", "new_iter3_10"]
    # simulation_location_list = [
    #             "new_iter2_8", "new_iter2_10", "old_new_iter2_8"]
    iteration_source_dic = {"bias_2":"../iterative_optimization_biased_sampling",
                            "bias_old_gamma":"../iterative_optimization_biased_sampling",
                            "iter1_with_bias_96percent":"../iterative_optimization_new_temp_range",
                            "iter1_with_bias_98percent":"../iterative_optimization_new_temp_range",
                            "new_iter1_0":"../iterative_optimization_new_temp_range",
                            "new_iter1_90":"../iterative_optimization_new_temp_range",
                            "new_iter1_96":"../iterative_optimization_new_temp_range",
                            "new_iter1_98":"../iterative_optimization_new_temp_range",
                            "new_iter1_combined_on_B":"../iterative_optimization_new_temp_range",
                            "new_iter2_8":"../iterative_optimization_new_temp_range",
                            "new_iter2_10":"../iterative_optimization_new_temp_range",
                            "old_new_iter2_8":"../iterative_optimization_new_temp_range",
                            "new_iter3_10":"../iterative_optimization_old_set",
                            "single":"../iterative_optimization_old_set",
                            "iter4_30":"../iterative_optimization_old_set",
                            "iter4_6":"../iterative_optimization_old_set",
                            "iter4_13":"../iterative_optimization_old_set",
                            "iter5_30":"../iterative_optimization_old_set",
                            "iter6_30":"../iterative_optimization_old_set",
                            }
    pdb_list_dic = {"../iterative_optimization_old_set":"old",
                    "../iterative_optimization_new_temp_range":"new",
                    "../iterative_optimization_biased_sampling":"new"}
    # new_simulation_list = ["iter1_with_bias_96percent", "new_iter2_10"]
    # old_protein_simulation_list = ["single", "new_iter3_10"]
    new_simulation_list = ["bias_2","bias_old_gamma", "iter1_with_bias_96percent", "iter1_with_bias_98percent", "new_iter2_10", "new_iter1_90", "new_iter2_8", "old_new_iter2_8"]
    old_protein_simulation_list = ["noFrag", "iter6_30", "iter5_30", "single", "new_iter3_10", "iter4_30", "iter4_6", "iter4_13"]
    simulation_location_list_dic = {}
    for p in dataset["new"][0]:
        name = p.lower()[:4]
        simulation_location_list_dic[name] = new_simulation_list
    for p in dataset["old"][0]:
        name = p.lower()[:4]
        simulation_location_list_dic[name] = old_protein_simulation_list

    if args.mode == 1:
        do(f"mkdir {optimization_folder}")
        cd(optimization_folder)
        do("cp ../phi_list.txt .")
        do("mkdir -p database/dompdb")
        do("mkdir proteins_name_list")
        do("mkdir slurms")
        do("mkdir -p subsetgammas")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p phis")
    if args.mode == 2:
        pdb_list, steps = dataset["old"]
        simulation_location = "native"
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            replace(f"{pre}/{name}_multi.in", "minimize", "#minimize")

            do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
            do(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
            # print(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
            do(f"cp {pre}/burial_gamma.dat {pre}/original_burial_gamma.dat")
            do(f"cp {gammaSource}/iteration_burial_gamma.dat {pre}/burial_gamma.dat")
    if args.mode == 3:
        pdb_list, steps = dataset["old"]
        simulation_location = "native"
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            # steps = 40
            steps = 2
            # steps = 80
            cd(name)
            # do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            do(f"run.py -n 1 {name} --commons 2 -s {steps} --runs 1 --start crystal")
            # do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")
    if args.mode == 4:
        pdb_list, steps = dataset["old"]
        cd("native")
        for p in pdb_list:
            name = p.lower()[:4]
            cd(f"{name}/simulation/0/")
            do("cp energy.log back_energy.log")
            do(f"rerun.py {name} -t")
            cd(f"../../../")
    if args.mode == 5:
        pdb_list, steps = dataset["old"]
        simulation_location = "native"
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            cd(f"{name}/simulation/0/rerun/")
            do(f"cp ../{name}.seq .")
            do(f"python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie {name}.seq")
            cd("../../../../")
    if args.mode == 6:
        pdb_list, steps = dataset["old"]
        # native_source = "../iterative_optimization_new_temp_range"
        native_source = "../iterative_optimization_old_set"

        do("mkdir -p database/dompdb")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p database/{name}")
            do(f"cp -r {native_source}/native/{name}/simulation/0/rerun/* database/{name}/")
            pre = f"database/{name}/"
            fileName = "movie.pdb"
            splitPDB(pre, fileName)
            do(f"cp database/{name}/frame0.pdb database/dompdb/{name}.pdb")

        pdb_list, steps = dataset["new"]
        native_source = "../iterative_optimization_new_temp_range"
        do("mkdir -p database/dompdb")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p database/{name}")
            do(f"cp -r {native_source}/native/{name}/simulation/0/rerun/* database/{name}/")
            pre = f"database/{name}/"
            fileName = "movie.pdb"
            splitPDB(pre, fileName)
            do(f"cp database/{name}/frame0.pdb database/dompdb/{name}.pdb")
    if args.mode == 7:
        # simulation_location_list = ["single", "new_iter3_10"]
        # simulation_location_list = ["iter4_30", "iter4_13", "iter4_6"]
        simulation_location_list = ["iter5_30"]
        # simulation_location_list = ["iter6_30"]
        simulation_location_list = ["noFrag"]
        # Go To list of simulation folder
        for simulation_location in simulation_location_list:
            print(simulation_location)
            # cd(iteration_source_dic[simulation_location]+"/"+simulation_location)
            iteration_source = iteration_source_dic[simulation_location]
            pdb_list, steps = dataset[pdb_list_dic[iteration_source]]
            cd(iteration_source+"/"+simulation_location)
            for p in pdb_list:
                name = p.lower()[:4]
                print(name)
                cd(f"{name}/simulation/")
                with open(f"run_on_scavenge.slurm", "w") as out:
                    out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d mar05 -m 77 -l {name}"))
                # replace(f"run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=5G")
                # replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=40G")
                do(f"sbatch run_on_scavenge.slurm")
                cd("../../")
            cd("..")
        # for simulation_location in simulation_location_list:
        #     print(simulation_location)
        #     # cd(iteration_source_dic[simulation_location]+"/"+simulation_location)
        #     iteration_source = iteration_source_dic[simulation_location]
        #     cd(iteration_source+"/"+simulation_location)
        #     for p in pdb_list:
        #         name = p.lower()[:4]
        #         print(name)
        #         for i in range(30):
        #             cd(f"{name}/simulation/{i}/0/")
        #             with open(f"run_on_scavenge.slurm", "w") as out:
        #                 out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d mar05 -m 77 -l {name}"))
        #             # replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=40G")
        #             do(f"sbatch run_on_scavenge.slurm")
        #             cd("../../../../")
        #     cd("..")
        # for simulation_location in simulation_location_list:
        #     print(simulation_location)
        #     cd(iteration_source_dic[simulation_location]+"/"+simulation_location)
        #     for p in pdb_list:
        #         name = p.lower()[:4]
        #         print(name)
        #         for i in range(30):
        #             cd(f"{name}/simulation/{i}/0/")
        #             do(f"cp ../{name}.seq .")
        #             do(f"python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie {name}.seq")
        #             cd("../../../../")
        #     cd("..")
    if args.mode == 77:
        name = args.label
        for i in range(30):
            cd(f"{i}/0/")
            do(f"cp ../{name}.seq .")
            do(f"python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie {name}.seq")
            cd("../..")
        # name = args.label
        # do(f"cp ../{name}.seq .")
        # do(f"python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie {name}.seq")
    if args.mode == 8:
        simulation_location_list = new_simulation_list + old_protein_simulation_list
        simulation_location_list = list(iteration_source_dic.keys())
        simulation_location_list = ["bias_2","bias_old_gamma"]
        for simulation_location in simulation_location_list:
            print(simulation_location)
            # cd(iteration_source_dic[simulation_location]+"/"+simulation_location)
            iteration_source = iteration_source_dic[simulation_location]
            pdb_list, steps = dataset[pdb_list_dic[iteration_source]]
            cd(iteration_source+"/"+simulation_location)
            for p in pdb_list:
                name = p.lower()[:4]
                print(name)
                cd(f"{name}/simulation/")
                with open(f"run_on_scavenge.slurm", "w") as out:
                    out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d mar05 -m 88 -l {simulation_location}_{name}"))
                replace(f"run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=5G")
                # replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=40G")
                do(f"sbatch run_on_scavenge.slurm")
                cd("../../")
            cd("..")
            # time.sleep(60)
            # for p in pdb_list:
            #     name = p.lower()[:4]
            #     print(name)
            #     for i in range(30):
            #         cd(f"{name}/simulation/{i}/")
            #         with open(f"run_on_scavenge.slurm", "w") as out:
            #             out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d mar05 -m 88 -l {simulation_location}_{name}_{i}"))
            #         # replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=40G")
            #         do(f"sbatch run_on_scavenge.slurm")
            #         cd("../../../")
            # cd("..")
            # time.sleep(60)
        # for simulation_location in simulation_location_list:
        #     print(simulation_location)
        #     iteration_source = iteration_source_dic[simulation_location]
        #     for p in pdb_list:
        #         print(p)
        #         name = p.lower()[:4]
        #         for i in range(30):
        #             do(f"cp -r {iteration_source}/{simulation_location}/{name}/simulation/{i}/0/ database/{simulation_location}_{name}_{i}")
        #             pre = f"database/{simulation_location}_{name}_{i}/"
        #             fileName = "movie.pdb"
        #             splitPDB(pre, fileName)
    if args.mode == 88:
        cmd = args.label
        for i in range(30):
            pre = f"/scratch/wl45/april_2019/database/{cmd}_{i}/"
            do(f"cp -r {i}/0/ {pre}")
            fileName = "movie.pdb"
            splitPDB(pre, fileName)
            do(f"rm {pre}movie.*")
        # cmd = args.label
        # for i in range(30):
        #     do(f"cp -r {i}/0/ ../../../../{optimization_folder}/database/{cmd}_{i}")
        #     pre = f"../../../../{optimization_folder}/database/{cmd}_{i}/"
        #     fileName = "movie.pdb"
        #     splitPDB(pre, fileName)
        #     do(f"rm {pre}movie.*")
    if args.mode == 9:
        # generate decoy file
        simulation_location_list = new_simulation_list + old_protein_simulation_list
        print("generating decoys file")
        do("mkdir -p decoys/lammps")
        decoy_n = 200
        for simulation_location in simulation_location_list:
            iteration_source = iteration_source_dic[simulation_location]
            pdb_list, steps = dataset[pdb_list_dic[iteration_source]]
            for p in pdb_list:
                name = p.lower()[:4]
                count = 0
                with open(f"decoys/lammps/{name}_{simulation_location}.decoys", "w") as out:
                    for i in range(30):
                        Q = pd.read_csv(f"database/{simulation_location}_{name}_{i}/wham.dat")[" Qw"].values
                        n = len(Q)
                        if n != 2000 and n != 1000:
                            print(f"database/{simulation_location}_{name}_{i},  {n}")
                        if n == 0:
                            # print("!! Using previous i.")
                            print(f"An error at {i}")
                            continue
                        for j in range(decoy_n):
                            # out.write(f"database/{name}_{i}/frame{j+1801}\n")
                            try:
                                out.write(f"database/{simulation_location}_{name}_{i}/frame{j+n-decoy_n+1} {Q[j+n-decoy_n]}\n")
                            except:
                                print(f"!!!database/{simulation_location}_{name}_{i}\n")
                        preQ = Q
                        prei = i
                        count += 1
                    if count != 30:
                        print(f"!! Using previous i. for {30-count}")
                    n = len(preQ)
                    Q = preQ
                    i = prei
                    while count < 30:
                        for j in range(decoy_n):
                            # out.write(f"database/{name}_{i}/frame{j+1801}\n")
                            try:
                                out.write(f"database/{simulation_location}_{name}_{i}/frame{j+n-decoy_n+1} {Q[j+n-decoy_n]}\n")
                            except:
                                print(f"!!!database/{simulation_location}_{name}_{i}\n")
                        count += 1
        # for p in pdb_list:
        #     name = p.lower()[:4]
        #     for simulation_location in simulation_location_list:
        #         with open(f"decoys/lammps/{name}_{simulation_location}.decoys", "w") as out:
        #             for i in range(30):
        #                 Q = pd.read_csv(f"database/{simulation_location}_{name}_{i}/wham.dat")[" Qw"].values
        #                 n = len(Q)
        #                 if n != 2000:
        #                     print(f"database/{simulation_location}_{name}_{i},  {n}")
        #                 for j in range(decoy_n):
        #                     # out.write(f"database/{name}_{i}/frame{j+1801}\n")
        #                     try:
        #                         out.write(f"database/{simulation_location}_{name}_{i}/frame{j+n-decoy_n+1} {Q[j+n-decoy_n]}\n")
        #                     except:
        #                         print(f"!!!database/{simulation_location}_{name}_{i}\n")
        # for p in pdb_list:
        #     name = p.lower()[:4]
        #     with open(f"decoys/lammps/{name}.decoys", "w") as out:
        #         for simulation_location in simulation_location_list:
        #             for i in range(30):
        #                 Q = pd.read_csv(f"database/{simulation_location}_{name}_{i}/wham.dat")[" Qw"].values
        #                 n = len(Q)
        #                 if n != 2000:
        #                     print(f"database/{simulation_location}_{name}_{i},  {n}")
        #                 for j in range(decoy_n):
        #                     # out.write(f"database/{name}_{i}/frame{j+1801}\n")
        #                     try:
        #                         out.write(f"database/{simulation_location}_{name}_{i}/frame{j+n-decoy_n+1} {Q[j+n-decoy_n]}\n")
        #                     except:
        #                         print(f"!!!database/{simulation_location}_{name}_{i}\n")
    if args.mode == 99:
        # generate decoy file
        print("generating decoys file")
        do("mkdir -p decoys/lammps")
        decoy_n = 100
        for p in pdb_list:
            name = p.lower()[:4]
            with open(f"decoys/lammps/{name}.decoys", "w") as out:
                for simulation_location in simulation_location_list:
                    for i in range(30):
                        Q = pd.read_csv(f"database/{simulation_location}_{name}_{i}/wham.dat")[" Qw"].values
                        for j in range(decoy_n):
                            # out.write(f"database/{name}_{i}/frame{j+1801}\n")
                            out.write(f"database/{simulation_location}_{name}_{i}/frame{j*2+1801} {Q[j*2+1800]}\n")
    if args.mode == 10:
        # simulation_location_list = ["single", "new_iter3_10"]
        simulation_location_list = new_simulation_list + old_protein_simulation_list

        for simulation_location in simulation_location_list:
            iteration_source = iteration_source_dic[simulation_location]
            pdb_list, steps = dataset[pdb_list_dic[iteration_source]]
            for p in pdb_list:
                name = p.lower()[:4]
                with open(f"proteins_name_list/proteins_name_list_{name}_{simulation_location}.txt", "w") as out:
                        out.write(f"{name}_{simulation_location}\n")
                with open(f"slurms/run_{name}_{simulation_location}.slurm", "w") as out:
                    out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 3 proteins_name_list/proteins_name_list_{name}_{simulation_location}.txt"))
                    # out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 1 proteins_name_list/proteins_name_list_{name}.txt"))
                replace(f"slurms/run_{name}_{simulation_location}.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=20G")
                do(f"sbatch slurms/run_{name}_{simulation_location}.slurm")
            # with open(f"proteins_name_list/proteins_name_list_{name}.txt", "w") as out:
            #     out.write(f"{name}\n")
            # with open(f"slurms/run_{name}.slurm", "w") as out:
            #     out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 2 proteins_name_list/proteins_name_list_{name}.txt"))
            #     # out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 1 proteins_name_list/proteins_name_list_{name}.txt"))
            # replace(f"slurms/run_{name}.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=50G")
            # do(f"sbatch slurms/run_{name}.slurm")
    if args.mode == 11:
        # for testing toymodel
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # n = len(simulation_location_list)
        # n = 2
        n = len(new_simulation_list)
        # import time
        # time.sleep(4000)
        with open(f"proteins_name_list.txt", "w") as out:
            for data in ["new", "old"]:
                pdb_list, steps = dataset[data]
                for p in pdb_list:
                    name = p.lower()[:4]
                    out.write(f"{name}\n")
        complete_proteins = "proteins_name_list.txt"

        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='lammps',
                                        num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False, mode=2, simulation_location_list_dic=simulation_location_list_dic)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_xl23(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False)
        # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='lammps',
        #                                 num_decoys=n*6000, noise_filtering=True, jackhmmer=False, read=False, mode=1, simulation_location_list=simulation_location_list)
    if args.mode == 12:
        # import time
        # time.sleep(7200)
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d mar05 -m 11"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
        do(f"sbatch slurms/run_on_scavenge.slurm")

if args.day == "feb27":
    print(args.day)
    dataset = {"old":("1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "), 40),
                "new":("1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "), 80),
                "test":(['1A2J', '1A3H', '1A5Y', '1A8Q', '1AGY', '1AKZ', '1AUZ', '1B1A', '1B31', '1B8X', '1BCO', '1BN6', '1BOH', '1BOI'], 80)}
    # pdb_list, steps = dataset["old"]
    pdb_list, steps = dataset["new"]
    # pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    # pdb_list = "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", ")
    # simulation_location = "noFrag"
    # simulation_location = "newFrag"
    # simulation_location = "newContactWithBurial"
    # simulation_location = "single"
    # simulation_location = "iter1"
    simulation_location = "bias_2"
    # simulation_location = "inverseBurial"
    gammaSource = "../optimization"
    # gammaSource = "../optimization_iter1"

    if args.mode == 1:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            if simulation_location != "single":
                do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
                do(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
                # print(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
                do(f"cp {pre}/burial_gamma.dat {pre}/original_burial_gamma.dat")
                do(f"cp {gammaSource}/iteration_burial_gamma.dat {pre}/burial_gamma.dat")

    if args.mode == 2:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            # steps = 40
            cd(name)
            # do(f"run.py -n 3 {name} --commons 1 -s {steps} --runs 1")
            fileName = f"{name}/{name}_multi.in"
            with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
                for line in file:
                    tmp = line.replace("K_RG", "1")  # remove in future.
                    tmp = tmp.replace("TARGET_RG", str(computeRg(f"{name}/crystal_structure.pdb")))  # change temp, remove in future
                    print(tmp, end='')
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1 --bias 1")
            cd("..")


if args.day == "feb25_2":
    print(args.day)
    dataset = {"old":("1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "), 40),
                "new":("1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "), 80),
                "test":(['1A2J', '1A3H', '1A5Y', '1A8Q', '1AGY', '1AKZ', '1AUZ', '1B1A', '1B31', '1B8X', '1BCO', '1BN6', '1BOH', '1BOI'], 40)}
    # pdb_list, steps = dataset["old"]
    # pdb_list, steps = dataset["test"]
    pdb_list, steps = dataset["new"]
    # pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    # pdb_list = "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", ")
    # simulation_location = "noFrag"
    # simulation_location = "newFrag"
    # simulation_location = "newContactWithBurial"
    # simulation_location = "iter1_with_bias_96percent"
    simulation_location = "iter1_with_bias_98percent"
    # simulation_location = "single"
    # simulation_location = "inverseBurial"
    # gammaSource = "../optimization"
    # gammaSource = "../optimization_iter1"
    gammaSource = "../optimization_with_biased_iter1"

    if args.mode == 1:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            # if simulation_location == "newContactWithBurial":
            if simulation_location != "single":
                do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
                do(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
                # print(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
                do(f"cp {pre}/burial_gamma.dat {pre}/original_burial_gamma.dat")
                do(f"cp {gammaSource}/iteration_burial_gamma.dat {pre}/burial_gamma.dat")

    if args.mode == 2:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            # steps = 40
            cd(name)
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")


if args.day == "feb25":

    # get the native phi by starting from the crystal structure.
    print(args.day)
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    # pdb_list = "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", ")
    # pdb_list = ['1A2J', '1A3H', '1A5Y', '1A8Q', '1AGY', '1AKZ', '1AUZ', '1B1A', '1B31', '1B8X', '1BCO', '1BN6', '1BOH', '1BOI']
    # simulation_location = "noFrag"
    # simulation_location = "newFrag"
    simulation_location = "newContactWithBurial"
    simulation_location = "bias_2"
    # simulation_location = "single"
    # simulation_location = "native"
    # iteration_source = "../iterative_optimization_new_temp_range"
    iteration_source = "../iterative_optimization_biased_sampling"
    # simulation_location = "inverseBurial"
    gammaSource = "../optimization"
    # individual_gammas_randomized_decoy=read_all_gammas("phi_list.txt", complete_proteins, training_decoy_method="shuffle", noise_filtering=True)
    if args.mode == 1:
        simulation_location = "native"
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")

            do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
            do(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
            # print(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
            do(f"cp {pre}/burial_gamma.dat {pre}/original_burial_gamma.dat")
            do(f"cp {gammaSource}/iteration_burial_gamma.dat {pre}/burial_gamma.dat")
    if args.mode == 2:
        simulation_location = "native"
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            # steps = 40
            steps = 5
            # steps = 80
            cd(name)
            # do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            do(f"run.py -n 1 {name} --commons 2 -s {steps} --runs 1 --start crystal")
            # do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")
    if args.mode == 3:
        cd("native")
        for p in pdb_list:
            name = p.lower()[:4]
            cd(f"{name}/simulation/0/")
            do("cp energy.log back_energy.log")
            do(f"rerun.py {name} -t")
            cd(f"../../../")
    if args.mode == 4:
        simulation_location = "native"
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            cd(f"{name}/simulation/0/rerun/")
            do(f"cp ../{name}.seq .")
            do(f"python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie {name}.seq")
            cd("../../../../")
    if args.mode == 5:
        do("mkdir -p database/dompdb")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p database/{name}")
            do(f"cp -r {iteration_source}/native/{name}/simulation/0/rerun/* database/{name}/")
            pre = f"database/{name}/"
            fileName = "movie.pdb"
            splitPDB(pre, fileName)
            do(f"cp database/{name}/frame0.pdb database/dompdb/{name}.pdb")
    if args.mode == 6:
        # Go To simulation folder
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            print(name)
            for i in range(30):
                cd(f"{name}/simulation/{i}/0/")
                do(f"cp ../{name}.seq .")
                do(f"python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie {name}.seq")
                cd("../../../../")
    if args.mode == 7:
        for p in pdb_list:
            name = p.lower()[:4]
            for i in range(30):
                do(f"cp -r {iteration_source}/{simulation_location}/{name}/simulation/{i}/0/ database/{name}_{i}")
                pre = f"database/{name}_{i}/"
                fileName = "movie.pdb"
                splitPDB(pre, fileName)
    if args.mode == 8:
        # generate decoy file
        print("generating decoys file")
        do("mkdir -p decoys/lammps")
        decoy_n = 200
        for p in pdb_list:
            name = p.lower()[:4]
            with open(f"decoys/lammps/{name}.decoys", "w") as out:
                for i in range(30):
                    Q = pd.read_csv(f"database/{name}_{i}/wham.dat")[" Qw"].values
                    for j in range(decoy_n):
                        # out.write(f"database/{name}_{i}/frame{j+1801}\n")
                        out.write(f"database/{name}_{i}/frame{j+1801} {Q[j+1800]}\n")
    if args.mode == 9:
        do("mkdir proteins_name_list")
        do("mkdir slurms")
        do("mkdir -p subsetgammas")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p phis")
    if args.mode == 10:
        for p in pdb_list:
            name = p.lower()[:4]
            with open(f"proteins_name_list/proteins_name_list_{name}.txt", "w") as out:
                out.write(f"{name}\n")
            with open(f"slurms/run_{name}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 2 proteins_name_list/proteins_name_list_{name}.txt"))
                # out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 1 proteins_name_list/proteins_name_list_{name}.txt"))
            do(f"sbatch slurms/run_{name}.slurm")
    if args.mode == 11:
        # for testing toymodel
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # import time
        # time.sleep(4000)
        with open(f"proteins_name_list.txt", "w") as out:
            for p in pdb_list:
                name = p.lower()[:4]
                out.write(f"{name}\n")
        complete_proteins = "proteins_name_list.txt"
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_xl23(complete_proteins, "phi_list.txt", decoy_method='lammps',
                                        num_decoys=6000, noise_filtering=True, jackhmmer=False, read=False)
    if args.mode == 12:
        with open(f"slurms/run_on_scavenge.slurm", "w") as out:
            out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d feb25 -m 11"))
        replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=40G")
        do(f"sbatch slurms/run_on_scavenge.slurm")
    if args.mode == 13:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"cp -r ../optimization_iter1/database/{name} database/")
    # if args.mode == 10:
    #     for p in pdb_list:
    #         name = p.lower()[:4]
    #         for i in range(30):
    #             with open(f"proteins_name_list/proteins_name_list_{name}_{i}.txt", "w") as out:
    #                 out.write(f"{name}_{i}\n")
    #             with open(f"slurms/run_{name}_{i}.slurm", "w") as out:
    #                 out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 1 proteins_name_list/proteins_name_list_{name}_{i}.txt"))
    #             do(f"sbatch slurms/run_{name}_{i}.slurm")
    # if args.mode == 8:
    #     # generate decoy file
    #     print("generating decoys file")
    #     do("mkdir -p decoys/lammps")
    #     decoy_n = 200
    #     for p in pdb_list:
    #         name = p.lower()[:4]
    #         for i in range(30):
    #             with open(f"decoys/lammps/{name}_{i}.decoys", "w") as out:
    #                 for j in range(decoy_n):
    #                     out.write(f"database/{name}_{i}/frame{j+1801}\n")
    # if args.mode == 10:
    #     with open("proteins_name_list.txt", "w") as out:
    #         for p in pdb_list:
    #             name = p.lower()[:4]
    #             out.write(f"{name}\n")
    #     do("python3 ~/opt/compute_phis.py ./proteins_name_list.txt -m 1")
    # if args.mode == 4:
    #     cd(simulation_location)
    #     for p in pdb_list:
    #         name = p.lower()[:4]
    #         print(name)
    #         for i in range(1):
    #             cd(f"{name}/simulation/{i}/0/")
    #             do(f"cp ../{name}.seq .")
    #             do(f"python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie {name}.seq")
    #             cd("../../../../")
    # if args.mode == 4:
    #     # clean up
    #     # do("rm database/S20_seq/*")
    #     # do("rm database/dompdb/*")
    #     do("mkdir -p decoys")
    #     do("mkdir -p gammas")
    #     do("mkdir -p phis")
    #     do("mkdir -p database")
    #     # do("cp -r ../optimization_iter1/database/dompdb database/")
    #     # do("cp ../optimization_iter1/proteins_name_list.txt .")
    #     # do("cp ../optimization_iter1/phi_list.txt .")
    # if args.mode == 5:
    #     for p in pdb_list:
    #         name = p.lower()[:4]
    #         for i in range(30):
    #             do(f"cp -r ../iterative_optimization_2/{simulation_location}/{name}/simulation/{i}/0/ database/{name}_{i}")
    #             pre = f"database/{name}_{i}/"
    #             fileName = "movie.pdb"
    #             splitPDB(pre, fileName)
    # if args.mode == 6:
    #     # generate decoy file
    #     print("generating decoys file")
    #     for p in pdb_list:
    #         name = p.lower()[:4]
    #         with open(f"decoys/lammps/{name}.decoys", "w") as out:
    #             for i in range(30):
    #                 for j in range(decoy_n):
    #                     out.write(f"database/{name}_{i}/frame{j+1}\n")
    # if args.mode == 7:
    #     with open("proteins_name_list.txt", "w") as out:
    #         for p in pdb_list:
    #             name = p.lower()[:4]
    #             out.write(f"{name}\n")
    #     do("python3 ~/opt/compute_phis.py proteins_name_list.txt")


if args.day == "feb23":
    print(args.day)
    dataset = {"old":("1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "), 40),
                "new":("1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "), 80)}
    pdb_list, steps = dataset["old"]
    # pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    pdb_list = "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", ")
    # simulation_location = "noFrag"
    # simulation_location = "newFrag"
    simulation_location = "newContactWithBurial"
    # simulation_location = "single"
    # simulation_location = "inverseBurial"
    gammaSource = "../optimization"

    if args.mode == 1:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")

            do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
            do(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
            # print(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
            do(f"cp {pre}/burial_gamma.dat {pre}/original_burial_gamma.dat")
            do(f"cp {gammaSource}/iteration_burial_gamma.dat {pre}/burial_gamma.dat")

    if args.mode == 2:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            steps = 40
            cd(name)
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")


if args.day == "feb22":
    print(args.day)
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    # simulation_location = "noFrag"
    # simulation_location = "newFrag"
    # simulation_location = "newContactWithBurial"
    simulation_location = "inverseBurial"
    gammaSource = "../optimization"
    # simulation_location = "single"
    if args.mode == 1:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")

            do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
            do(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")

            do(f"cp {pre}/burial_gamma.dat {pre}/original_burial_gamma.dat")
            do(f"cp {gammaSource}/iteration_burial_gamma.dat {pre}/burial_gamma.dat")

    if args.mode == 2:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            steps = 80
            cd(name)
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")

if args.day == "feb21":
    if args.mode == 1:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # import time
        # time.sleep(4000)
        complete_proteins = "database/cath-dataset-nonredundant-S20Clean.list"
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
                                        num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
if args.day == "feb20":
    # for first round of optimization
    if args.mode == 1:
        with open("database/cath-dataset-nonredundant-S20Clean.list") as f:
            content = f.readlines()
        pos = 0
        i = 0
        n = len(content)
        # n = 100  # for testing
        while pos < n:
            with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
                for ii in range(20):
                    if pos < n:
                        out.write(content[pos])
                    pos += 1
                i += 1
        print(i-1)
        do("mkdir proteins_name_list")
        do("mkdir slurms")
        do("mkdir -p subsetgammas")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p phis")
    if args.mode == 2:
        n = 712
        for i in range(n):
            with open(f"slurms/run_{i}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py proteins_name_list/proteins_name_list_{i}.txt"))
            do(f"sbatch slurms/run_{i}.slurm")
    if args.mode == 3:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # import time
        # time.sleep(4000)
        complete_proteins = "database/cath-dataset-nonredundant-S20Clean.list"
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
                                        num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
        # then do gg_server.py -d feb18 -m 2
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

    if args.mode == 5:
        # for testing toymodel
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # import time
        # time.sleep(4000)
        complete_proteins = "proteins_name_list_tiny.txt"
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_xl23(complete_proteins, "phi_list.txt", decoy_method='shuffle',
                                        num_decoys=1000, noise_filtering=True, jackhmmer=False)
        # individual_gammas_randomized_decoy=read_all_gammas("phi_list.txt", complete_proteins, training_decoy_method="shuffle", noise_filtering=True)

if args.day == "feb19":
    print(args.day)
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    # simulation_location = "noFrag"
    # simulation_location = "newFrag"
    simulation_location = "newContact"
    simulation_location = "single"
    if args.mode == 1:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            do(f"cp ../jan_optimization/iteration_gamma.dat {pre}/")
            do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
            do(f"cp {pre}/iteration_gamma.dat {pre}/gamma.dat")
    if args.mode == 2:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            steps = 80
            cd(name)
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")

if args.day == "feb18_2":
    print(args.day)
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    # simulation_location = "noFrag"
    # simulation_location = "newFrag"
    # simulation_location = "top20"
    # simulation_location = "top1"
    simulation_location = "top3"
    # simulation_location = "fragment"
    if args.mode == 1:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            steps = 80
            cd(name)
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")
    if args.mode == 2:
        fragLibrary = f"fragment_memory_{simulation_location}"
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r ../fragment_memory/fraglib {pre}/")
            do(f"cp {fragLibrary}/{name}.mem {pre}/frags.mem")

            # replace(f"{simulation_location}/{name}/{name}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Water\]", "\[Water\]-")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Burial\]", "\[Burial\]-")
            # do(f"rm all_simulations/{name}/{name}/gamma.dat")

if args.day == "feb18":
    if args.mode == 1:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        name = f"proteins_name_list/proteins_name_list_{args.label}.txt"
        subset = read_column_from_file(name, 1)
        print(subset)
        # import time
        # time.sleep(4000)
        complete_proteins = "database/cath-dataset-nonredundant-S20Clean.list"
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
                                        num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=subset, subset_index=args.label, read=0)
    if args.mode == 2:
        n = 712
        for i in range(n):
            with open(f"slurms/run_on_scavenge_{i}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python3 ~/opt/gg_server.py -d feb18 -m 1 -l {i}"))
            replace(f"slurms/run_on_scavenge_{i}.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=20G")
            do(f"sbatch slurms/run_on_scavenge_{i}.slurm")
    if args.mode == 3:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # import time
        # time.sleep(4000)
        complete_proteins = "database/cath-dataset-nonredundant-S20Clean.list"
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
                                        num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=1)

if args.day == "feb13":
    print(args.day)
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    decoy_n = 20
    # simulation_location = "iter1"
    simulation_location = "iter2"
    if args.mode == 1:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            print(name)
            for i in range(30):
                cd(f"{name}/simulation/{i}/0/")
                do(f"cp ../{name}.seq .")
                do(f"python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie {name}.seq")
                cd("../../../../")
    if args.mode == 2:
        # clean up
        # do("rm database/S20_seq/*")
        # do("rm database/dompdb/*")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p phis")
        do("mkdir -p database")
        # do("cp -r ../optimization_iter1/database/dompdb database/")
        # do("cp ../optimization_iter1/proteins_name_list.txt .")
        # do("cp ../optimization_iter1/phi_list.txt .")
    if args.mode == 3:
        for p in pdb_list:
            name = p.lower()[:4]
            for i in range(30):
                do(f"cp -r ../iterative_optimization_2/{simulation_location}/{name}/simulation/{i}/0/ database/{name}_{i}")
                pre = f"database/{name}_{i}/"
                fileName = "movie.pdb"
                splitPDB(pre, fileName)
    if args.mode == 4:
        # generate decoy file
        print("generating decoys file")
        for p in pdb_list:
            name = p.lower()[:4]
            with open(f"decoys/lammps/{name}.decoys", "w") as out:
                for i in range(30):
                    for j in range(decoy_n):
                        out.write(f"database/{name}_{i}/frame{j+1}\n")
    if args.mode == 5:
        with open("proteins_name_list.txt", "w") as out:
            for p in pdb_list:
                name = p.lower()[:4]
                out.write(f"{name}\n")
        do("python3 ~/opt/compute_phis.py proteins_name_list.txt")
    if args.mode == 6:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # import time
        # time.sleep(4000)
        # complete_proteins = "database/cath-dataset-nonredundant-S20Clean.list"
        complete_proteins = "proteins_name_list.txt"
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_xl23(complete_proteins, "phi_list.txt", decoy_method='lammps', num_decoys=decoy_n, noise_filtering=True, jackhmmer=False)
        # individual_gammas_randomized_decoy=read_all_gammas("phi_list.txt", complete_proteins, training_decoy_method="shuffle", noise_filtering=True)
    if args.mode == 7:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            do(f"cp ../optimization_iter2/iteration_gamma.dat {pre}/")
            do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
            do(f"cp {pre}/iteration_gamma.dat {pre}/gamma.dat")
    if args.mode == 8:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            # steps = 40
            # steps = 1
            steps = 80
            cd(name)
            # do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            # do(f"run.py -n 10 {name} --commons 2 -s {steps} --runs 1 --start crystal")
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")



if args.day == "feb07":
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    decoy_n = 20
    if args.mode == 1:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"cp -r ../iterative_optimization/native/{name}/simulation/0/rerun database/{name}")
            pre = f"database/{name}/"
            fileName = "movie.pdb"
            splitPDB(pre, fileName)
            do(f"cp database/{name}/frame0.pdb database/dompdb/{name}.pdb")
    if args.mode == 2:
        for p in pdb_list:
            name = p.lower()[:4]
            pre = f"database/{name}/"
            fileName = "movie.pdb"
            splitPDB(pre, fileName)
            do(f"cp database/{name}/frame0.pdb database/dompdb/{name}.pdb")
    if args.mode == 3:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"cp database/{name}/frame0.pdb database/dompdb/{name}.pdb")
    if args.mode == 4:
        # generate decoy file
        print("generating decoys file")
        for p in pdb_list:
            name = p.lower()[:4]
            with open(f"decoys/lammps/{name}.decoys", "w") as out:
                for i in range(decoy_n):
                    out.write(f"database/{name}/frame{i+1}\n")
    if args.mode == 5:
        # clean up
        # do("rm database/S20_seq/*")
        # do("rm database/dompdb/*")
        do("mkdir -p decoys")
        do("mkdir -p gammas")
        do("mkdir -p phis")
    if args.mode == 6:
        with open("proteins_name_list.txt", "w") as out:
            for p in pdb_list:
                name = p.lower()[:4]
                out.write(f"{name}\n")
        do("python3 ~/opt/compute_phis.py proteins_name_list.txt")
    if args.mode == 7:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # import time
        # time.sleep(4000)
        # complete_proteins = "database/cath-dataset-nonredundant-S20Clean.list"
        complete_proteins = "proteins_name_list.txt"
        A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_xl23(complete_proteins, "phi_list.txt", decoy_method='lammps', num_decoys=decoy_n, noise_filtering=True, jackhmmer=False)
        # individual_gammas_randomized_decoy=read_all_gammas("phi_list.txt", complete_proteins, training_decoy_method="shuffle", noise_filtering=True)
    if args.mode == 8:
        for i in range(475):
            with open(f"run_{i}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py proteins_name_list_{i}.txt"))
            do(f"sbatch run_{i}.slurm")
    if args.mode == 9:
        simulation_location = "native"
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            cd(f"{name}/simulation/0/rerun/")
            do(f"cp ../{name}.seq .")
            do(f"python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie {name}.seq")
            cd("../../../../")

if args.day == "feb06":
    print(args.day)
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    # simulation_location = "noFrag"
    # simulation_location = "newContact_noFrag"
    # simulation_location = "native"
    # simulation_location = "iter1"
    simulation_location = "filtered_gamma_iter1"
    if args.mode == 1:
        print('hi')
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            do(f"cp ../optimization_iter1/iteration_gamma.dat {pre}/")
            do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
            do(f"cp {pre}/iteration_gamma.dat {pre}/gamma.dat")
    if args.mode == 2:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            # steps = 40
            # steps = 1
            steps = 40
            cd(name)
            # do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            # do(f"run.py -n 10 {name} --commons 2 -s {steps} --runs 1 --start crystal")
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")
    if args.mode == 3:
        cd("native")
        for p in pdb_list:
            name = p.lower()[:4]
            cd(f"{name}/simulation/0/")
            do("cp energy.log back_energy.log")
            do(f"rerun.py {name} -t")
            cd(f"../../../")

if args.day == "jan29":
    print(args.day)
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    # simulation_location = "noFrag"
    # simulation_location = "newFrag"
    simulation_location = "top5_noeven"
    if args.mode == 1:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            steps = 40
            cd(name)
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")
    if args.mode == 2:
        fragLibrary = "fragment_memory_top5_noeven"
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r fragment_memory/fraglib {pre}/")
            do(f"cp {fragLibrary}/{name}.mem {pre}/frags.mem")

            # replace(f"{simulation_location}/{name}/{name}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Water\]", "\[Water\]-")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Burial\]", "\[Burial\]-")
            # do(f"rm all_simulations/{name}/{name}/gamma.dat")

if args.day == "jan23":
    print(args.day)
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    # simulation_location = "noFrag"
    # simulation_location = "newFrag"
    simulation_location = "all_simulations"
    if args.mode == 1:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            steps = 40
            cd(name)
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")

if args.day == "jan22":
    print(args.day)
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    # simulation_location = "noFrag"
    # simulation_location = "newFrag"
    simulation_location = "top5"
    if args.mode == 1:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            steps = 40
            cd(name)
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")
    if args.mode == 2:
        fragLibrary = "fragment_memory_top5"
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r fragment_memory/fraglib {pre}/")
            do(f"cp {fragLibrary}/{name}.mem {pre}/frags.mem")

            # replace(f"{simulation_location}/{name}/{name}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Water\]", "\[Water\]-")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Burial\]", "\[Burial\]-")
            # do(f"rm all_simulations/{name}/{name}/gamma.dat")

if args.day == "jan20":
    print(args.day)
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    # simulation_location = "noFrag"
    simulation_location = "FragOnly"
    if args.mode == 1:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            steps = 40
            cd(name)
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")
    if args.mode == 2:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            pre = f"{simulation_location}/{name}/{name}"
            # replace(f"{simulation_location}/{name}/{name}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            replace(f"{pre}/fix_backbone_coeff.data", "\[Water\]", "\[Water\]-")
            replace(f"{pre}/fix_backbone_coeff.data", "\[Burial\]", "\[Burial\]-")
            # do(f"rm all_simulations/{name}/{name}/gamma.dat")

if args.day == "jan17":
    print(args.day)
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    # simulation_location = "noFrag"
    # simulation_location = "newContact_noFrag"
    simulation_location = "newContact_singleFrag"
    if args.mode == 1:
        cd(simulation_location)
        for p in pdb_list:
            name = p.lower()[:4]
            steps = 40
            cd(name)
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")
    if args.mode == 2:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            replace(f"{simulation_location}/{name}/{name}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            # replace(f"{simulation_location}/{name}/{name}/fix_backbone_coeff.data", "\[Water\]", "\[Water\]-")
            # replace(f"{simulation_location}/{name}/{name}/fix_backbone_coeff.data", "\[Burial\]", "\[Burial\]-")
            # do(f"rm all_simulations/{name}/{name}/gamma.dat")
    if args.mode == 3:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {simulation_location}/{name}")
            pre = f"{simulation_location}/{name}/{name}"
            do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
            # replace(f"{pre}/fix_backbone_coeff.data", "\[Fragment_Memory_Table\]", "\[Fragment_Memory_Table\]-")
            do(f"cp ../optimization/iteration_gamma.dat {pre}/")
            do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
            do(f"cp {pre}/iteration_gamma.dat {pre}/gamma.dat")

if args.day == "jan16":
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    if args.mode == 1:
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            if name == "2fha" or name == "1mba":
                steps = 20
            else:
                steps = 40
            cd(name)
            do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
            cd("..")


if args.day == "jan10":
    pdb_list = "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", ")
    if args.mode == 1:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p all_simulations/{name}")
            do(f"cp -r ../iterative_optimization_4/all_simulations/{name}/{name} all_simulations/{name}/")
            replace(f"all_simulations/{name}/{name}/fix_backbone_coeff.data", "\[Water\]", "\[Water\]-")
            replace(f"all_simulations/{name}/{name}/fix_backbone_coeff.data", "\[Burial\]", "\[Burial\]-")
            do(f"rm all_simulations/{name}/{name}/gamma.dat")

if args.day == "jan06":
    if args.mode == 1:
        cd("all_simulations")
        name = "2fha"
        if name == "2fha" or name == "1mba":
            steps = 20
        else:
            steps = 40
        cd(name)
        do(f"run.py -n 20 {name} --commons 2 -s {steps} --runs 3")
        cd("..")

'''
'''
#####----------------------------2018--------------------------------------------
if args.day == "dec30":
    pdb_list = "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", ")
    if args.mode == 1:
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p all_simulations/{name}")
            do(f"cp -r ../iterative_optimization_4/all_simulations/{name}/{name} all_simulations/{name}/")
            do(f"cp ../optimization/iteration_gamma.dat all_simulations/{name}/{name}/")
            do(f"cp all_simulations/{name}/{name}/gamma.dat all_simulations/{name}/{name}/original_gamma.dat")
            do(f"cp all_simulations/{name}/{name}/iteration_gamma.dat all_simulations/{name}/{name}/gamma.dat")
    if args.mode == 2:
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

if args.day == "dec11":
    if args.mode == 1:
        all_results = pd.read_csv("filtered.csv", index_col=0)
        for i, line in all_results.iterrows():
            print(i, line["name"], line["folder"])

            name = line["name"]
            folder = line["folder"]
            do(f"mkdir -p {name}")
            from_file = f"/scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/{line['name']}/refinement/post-processing/1PCbias/sequentialPC/towardsPC1/lowTstructure/lowTstructure{line['folder']}.pdb"
            os.system(f"cp {from_file} {name}/{name}_{folder}_{i%5}.pdb")

if args.day == "dec09":
    if args.mode == 1:
        for i in range(10):
            for percore in [1]:
                for cpu in [1]:
                    for thread in [-1]:
                        do(f"python3 mm_sbatch.py -c {cpu} -t {thread} --percore {percore} --memory 32 -i {i}")
if args.day == "dec06":
    if args.mode == 1:
        for i in range(10):
            for percore in [1]:
                for cpu in [2]:
                    for thread in [-1]:
                        do(f"python3 mm_sbatch.py -c {cpu} -t {thread} --percore {percore} --memory 16 -i {i}")
    if args.mode == 2:
        all_results = pd.read_csv("all_results.csv", index_col=0)
        for i, line in all_results.query("result == 'picked_top5'").sort_values(["name", "prob"], ascending=False).reset_index(drop=True).iterrows():
            print(i, line["name"], line["folder"])
            name = line["name"]
            folder = line["folder"]
            from_file = f"/scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/{line['name']}/refinement/post-processing/1PCbias/sequentialPC/towardsPC1/lowTstructure/lowTstructure{line['folder']}.pdb"
            os.system(f"cp {from_file} {name}_{folder}_{i%5}.pdb")
if args.day == "dec05":
    if args.mode == 1:
        for cpu in [8]:
            for thread in [-1, 1, 2, 4, 8, 16, 32, 48]:
                for percore in [2]:
                    do(f"python3 mm_sbatch.py -c {cpu} -t {thread} --percore {percore}")
    if args.mode == 2:
        for percore in [1, 2]:
            for cpu in [4]:
                for thread in [-1, 1, 2, 4, 8, 16, 32, 48]:
                    do(f"python3 mm_sbatch.py -c {cpu} -t {thread} --percore {percore} --memory 8")
    if args.mode == 3:
        name_list = ["tr884-halfDIHE", "tr872-halfDIHE", "tr948-halfDIHE", "tr894", "tr882", "tr594", "tr869", "tr898", "tr862", "tr877", "tr872", "tr885", "tr866", "tr868", "tr884", "tr895", "tr896", "tr870", "tr921", "tr922", "tr891", "tr948", "tr947"]
        for name in name_list:
            do(f"mkdir {name}")
            # do(f"cp /scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/selection/{name}/PCselection/* {name}/")
            do(f"cp /scratch/xl23/home/xl23/notsJob/gromacs/all-atom/aawsem/pca/results/selection/PCselection/{name}/* {name}/")

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
    if args.mode == 5:
        # do("mkdir -p decoys")
        # do("mkdir -p gammas")
        # do("mkdir -p phis")
        for i in range(400, 475):
            with open(f"run_{i}.slurm", "w") as out:
                out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py proteins_name_list_{i}.txt"))
            do(f"sbatch run_{i}.slurm")
if args.day == "oct26":
    if args.mode == 1:
        from pyCodeLib import *
        import warnings
        warnings.filterwarnings('ignore')
        # import time
        # time.sleep(4000)
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
