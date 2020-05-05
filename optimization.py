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
# from small_script.variable_test import variable_test
# from small_script.variable_test2 import variable_test2
import subprocess
# from small_script.myFunctions import compute_theta_for_each_helix
# from small_script.myFunctions import *
from collections import defaultdict
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
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
with open('cmd_optimization.txt', 'a') as f:
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

base_casp14_slurm = '''\
#!/bin/bash
#SBATCH --reservation=casp14
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
##SBATCH --gres=gpu:1
##SBATCH --gres=gpu:1
#SBATCH --time=1:00:00
##SBATCH --mem-per-cpu=10G
#SBATCH --export=ALL
#SBATCH -o outs/slurm-%j.out
#SBATCH --mail-user=luwei0917@gmail.com
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
def getFromTerminal(CMD):
    return subprocess.Popen(CMD,stdout=subprocess.PIPE,shell=True).communicate()[0].decode()

def split_proteins_name_list(pdb_per_txt=1):
    with open("protein_list") as f:
        content = f.readlines()
    pos = 0
    i = 0
    n = len(content)
    # n = 100  # for testing
    while pos < n:
        with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
            for ii in range(pdb_per_txt):
                if pos < n:
                    out.write(content[pos])
                pos += 1
            i += 1
    print(i)
    n = i
    return n

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

def waitForJobs(jobIdList, sleepInterval=30, userName="wl45"):
    from datetime import datetime as dt
    if len(jobIdList) == 0:
        return
    previousJobNotFinished = True
    while previousJobNotFinished:
        print(f"Waiting for previous jobs {jobIdList}", dt.now())
        previousJobNotFinished = False
        a = getFromTerminal(f"squeue -u {userName}")
        for jobId in jobIdList:
            if jobId in a:
                previousJobNotFinished = True
        if previousJobNotFinished:
            time.sleep(sleepInterval)
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
# data = pd.read_csv("~/Research/database/training_set_apr11_2020.csv")
# specific_decoys = data.query("Length < 150 and Length > 70").reset_index(drop=True)
# pdb_list = specific_decoys["Protein"].to_list()
# pdb_list = [a.lower() for a in pdb_list]
# skip_pdb_list = ["1puc", "1skz"]
# skip_pdb_list += ["1msc", "1fmb", "1gvp", "2tgi", "1whi", "1baj", "1rmd", "1div"]  #dimer.
# skip_pdb_list += ["1aqe"]  # lots of ligand
# skip_pdb_list += ['1by2', '1rcb', '1hyp', '3lzt', '3cyr', '7rsa', '1rzl', '1b6e', '1poc', '1bea', '1poa']   # has at least 6 CYS.
# filtered_pdb_list = [x for x in pdb_list if x not in skip_pdb_list]
dataset["apr11_2020"] = ['1hoe', '1tif', '1vcc', '1by9', '1bdo', '451c', '1cc5', '1bb9', '1pht', '1opd', '1a32', '1ptf', '1cyo', '1tig', '1ctj', '1fna', '1who', '2cbp', '2acy', '1plc', '1bm8', '1opc', '3vub', '1tul', '1kte', '1erv', '1btn', '1a1x', '1bkf', '1ycc', '1sfp', '1kpf', '2mcm', '2pii', '1a6f', '1tmy', '2a0b', '1mai', '1neu', '1dun', '2sak', '1dhn', '1cxc', '1bgf', '1bqk', '3pyp', '1bfg', '1opy', '1rlw', '1rie', '3chy', '1cpq', '1pdo', '1hmt', '1htp', '1c52', '1kuh', '1crb', '1aqt', '2end', '5nul', '1pne', '1lcl', '2sns', '1flp', '1tfe', '1ax8', '1pkp', '1rss', '1jon', '1vls', '1lba', '1aly', '1mba', '2hbg', '1akr', '1osa']



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

def slurmRun(slurmFileName, cmd, template=scavenge_slurm, memory=1, thread=1, runOnServer=True):
    os.system("mkdir -p outs")
    if not runOnServer:
        do(cmd)
        return "runOnLocalMachine"
    with open(slurmFileName, "w") as out:
        out.write(template.format(cmd))
        # out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 1 proteins_name_list/proteins_name_list_{name}.txt"))
    replace(slurmFileName, "#SBATCH --mem-per-cpu=1G", f"#SBATCH --mem-per-cpu={memory}G")
    replace(slurmFileName, "#SBATCH --cpus-per-task=1", f"#SBATCH --cpus-per-task={thread}")

    a = getFromTerminal(f"sbatch {slurmFileName}")
    jobId = a.split(" ")[-1].strip()
    return jobId


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
        filtered_gamma_ = np.loadtxt(name)
        if len(filtered_gamma_) == 840:
            filtered_gamma = filtered_gamma_
            n = 4
            ax_title_list=["Direct_H", "Direct_P", "High density(protein)", "Low density(water)"]
        else:
            n = 3
            filtered_gamma = np.zeros(630)
            for i in range(min(len(filtered_gamma_), 630)):
                filtered_gamma[i] = filtered_gamma_[i]
            ax_title_list=["Direct", "High density(protein)", "Low density(water)"]
        figureName = f"{save_gamma_pre}/figures/{trial_name}_cutoff{cutoff_i}_equal_legend"
        title = f"{trial_name}_cutoff{cutoff_i}"
        print(cutoff_i, name, figureName)
        show_together_v2(filtered_gamma, figureName, title=title, inferBound=2, n=n, ax_title_list=ax_title_list)

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


def get_z_scores(trial_name):
    from pyCodeLib import validate_hamiltonian_wei
    import warnings
    warnings.filterwarnings('ignore')
    # print("mar08")
    cutoff = 400
    # cutoff = 600
    pre = "."
    gamma_pre = f"{pre}/saved_gammas"
    # trial_name = args.label
    do(f"cp -r gammas back_gammas_{trial_name}")
    gamma_file_name = f"{gamma_pre}/{trial_name}_cutoff{cutoff}_impose_Aprime_constraint"
    # gamma_file_name = f"{pre}/iter_2_30"
    data = validate_hamiltonian_wei("phi_list.txt", "protein_list_complete", gamma_file_name, "openMM", 100, mode=0)
    dataFile = f"optimization_2020.csv"
    data.to_csv(dataFile)

def get_weighted_Q(preDecoyBiasName, decoyBiasName):
    from pyCodeLib import read_phi_list, read_column_from_file, get_parameters_string
    phi_list_file_name = "phi_list.txt"
    training_set_file = "protein_list_complete"
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
            fromFile = f"../phis/{phi}_{protein}_{preDecoyBiasName}_{decoy_method}_{parameters_string}"
            a = np.loadtxt(fromFile)
            # print(a)
            # b = (1-a) * normalized_Z_weight
            # toFile = f"../phis/{phi}_{protein}_{decoyBiasName}_{decoy_method}_{parameters_string}"
            # np.savetxt(toFile, b, fmt='%.4f')
            if preDecoyBiasName == "decoysQ":
                b = (1-a) * normalized_Z_weight
            elif preDecoyBiasName == "deocyQZBias":
                b = a * normalized_Z_weight
            else:
                b = a * normalized_Z_weight
            # b = (1-a) * normalized_Z_weight_sq
            # b = (1-a) * normalized_Z_weight_half
            toFile = f"../phis/{phi}_{protein}_{decoyBiasName}_{decoy_method}_{parameters_string}"
            np.savetxt(toFile, b, fmt='%.4f')

def getSeqFromPDB(fileLocation):
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

def optimization_setupDatabase(pdb_list, fromFolder, toFolder):
    do(f"mkdir -p {toFolder}/database/dompdb")
    do(f"mkdir -p {toFolder}/database/S20_seq")
    for pdb in pdb_list:
        do(f"cp {fromFolder}/{pdb}.pdb {toFolder}/database/dompdb/")

    for pdb in pdb_list:
        fileLocation = f"{toFolder}/database/dompdb/{pdb}.pdb"
        seq = getSeqFromPDB(fileLocation)
        opt_pdbName = pdb
        fileLocation = f"{toFolder}/database/S20_seq/{opt_pdbName}.seq"
        with open(fileLocation, "w") as out:
            out.write(seq+"\n")

def optimization_setupFolder(phi_list, pdb_list):
    do(f"cp {phi_list} phi_list.txt")
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

def optimization_setupDecoys_shuffle(n_decoys, pdb_per_txt=2, runOnServer=True, template=base_slurm):
    n = split_proteins_name_list(pdb_per_txt=pdb_per_txt)
    jobIdList = []
    for i in range(n):
        proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
        # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
        # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 -d genDecoy -m 1 {proteins}", template=base_slurm, memory=10)
        jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/generate_decoys.py -n {n_decoys} -m 1 {proteins}", template=template, memory=10, runOnServer=runOnServer)
        # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 {proteins}", template=base_slurm, memory=10, runOnServer=runOnServer)
        # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
        jobIdList.append(jobId)
        do(f"cat {proteins} >> iter0_complete.txt")
    waitForJobs(jobIdList, sleepInterval=60)


def optimization_convertPDB(base_path, simulationType, folder_list, pdb_list, run_n):
    cd(f"{base_path}/{simulationType}")
    jobIdList = []
    for folder in folder_list:
        cd(folder)
        for pdb in pdb_list:
            for i in range(run_n):
                cd(f"{pdb}/{i}")
                # do()
                cmd = "python /projects/pw8/wl45/openawsem/helperFunctions/convertOpenmmTrajectoryToStandardMovie.py movie.pdb"
                jobId = slurmRun(f"convert.slurm", cmd, template=scavenge_slurm)
                jobIdList.append(jobId)
                cd("../..")
        cd("..")
    waitForJobs(jobIdList, sleepInterval=30)

def optimization_setupDecoys(optimizationBasePath, optimizationFolder, folder_list, pdb_list, simulationType, base_path, run_n, lastN_frame,
                    phi_list="../phi_list_environment_complete.txt"):
    do(f"mkdir -p {optimizationBasePath}/{optimizationFolder}")
    cd(f"{optimizationBasePath}/{optimizationFolder}")
    do(f"cp {phi_list} phi_list.txt")
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
# if args.mode == 3:
    outFile = read_simulation_info(folder_list, pdb_list, simulationType, base_path, run_n)
    # today = datetime.today().strftime('%m-%d')
    outFile = f"{simulationType}.csv"
# if args.mode == 13:
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

def optimization_compute_phi_and_gamma(complete_folder_list, pdb_list, runOnServer=True, phi_cmd="python3 ~/opt/compute_phis.py -m 7 {}", gamma_cmd="python3 ~/opt/gg_server.py -d apr11 -m 44", template=base_slurm):
    n = split_proteins_name_list(pdb_per_txt=2)
    jobIdList = []
    for i in range(n):
        proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
        # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
        jobId = slurmRun(f"slurms/run_{i}.slurm", phi_cmd.format(proteins), template=template, memory=10, runOnServer=runOnServer)
        # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
        jobIdList.append(jobId)
        do(f"cat {proteins} >> iter0_complete.txt")
    waitForJobs(jobIdList, sleepInterval=60)

    with open("protein_list_complete", "w") as out:
        for folder in complete_folder_list:
            for pdb in pdb_list:
                out.write(f"{pdb}_{folder}\n")
    jobIdList = []
    jobId = slurmRun(f"slurms/run_on_scavenge.slurm", f"{gamma_cmd}", template=template, memory=60, runOnServer=runOnServer)
    jobIdList.append(jobId)
    waitForJobs(jobIdList, sleepInterval=60)


def optimization_z_bias(iteration, c=-50):
    # mode = 1 for this environment.
    do(f"optimization_analyze.py {iteration} --proteinList protein_list_complete -c {c} -m 1")
    trial_name = iteration
    get_z_scores(trial_name)
    get_weighted_Q("decoysQ", "deocyQZBias")
# if args.mode == 7:
    decoyBiasName = "deocyQZBias"
    jobIdList = []
    jobId = slurmRun(f"slurms/get_gamma_{decoyBiasName}.slurm", f"python3 ~/opt/gg_server.py -d apr11 -m 77 -l {decoyBiasName}", template=base_slurm, memory=60)
    # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
    jobIdList.append(jobId)
    waitForJobs(jobIdList, sleepInterval=60)
# if args.mode == 8:
    # i = "iter4_shift_center_z_weighted"
    trial_name = f"{iteration}_zBias_1"
    # i = args.label
    do(f"optimization_analyze.py {trial_name} --proteinList protein_list_complete -c {c} -m 1")

    get_z_scores(trial_name)
    get_weighted_Q("deocyQZBias", "decoyQZBias_2")

    decoyBiasName = "decoyQZBias_2"
    jobIdList = []
    jobId = slurmRun(f"slurms/get_gamma_{decoyBiasName}.slurm", f"python3 ~/opt/gg_server.py -d apr11 -m 77 -l {decoyBiasName}", template=base_slurm, memory=60)
    # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
    jobIdList.append(jobId)
    waitForJobs(jobIdList, sleepInterval=60)

    trial_name = f"{iteration}_zBias_2"
    print(trial_name)
    # i = args.label
    do(f"optimization_analyze.py {trial_name} --proteinList protein_list_complete -c {c} -m 1")


def gamma_to_plot(name, figureName, title):
    filtered_gamma_ = np.loadtxt(name)
    if len(filtered_gamma_) == 840:
        filtered_gamma = filtered_gamma_
        n = 4
        ax_title_list=["Direct_H", "Direct_P", "High density(protein)", "Low density(water)"]
    else:
        n = 3
        filtered_gamma = np.zeros(630)
        for i in range(min(len(filtered_gamma_), 630)):
            filtered_gamma[i] = filtered_gamma_[i]
        ax_title_list=["Direct", "High density(protein)", "Low density(water)"]
    # print(cutoff_i, name, figureName)
    show_together_v2(filtered_gamma, figureName, title=title, inferBound=2, n=n, ax_title_list=ax_title_list)



if args.day == "example":

    phi_list = "/home/wl45/opt/optimization/phi_list_contact.txt"
    # phi_list = "phi_list.txt"
    # pdb_list = dataset["apr11_2020"]
    # data = pd.read_csv("/home/wl45/dataset/has_structures_small_dataset_cleaned.csv")
    # pdb_list = data["Protein"].to_list()
    pdb_list = ['3s5rA00', '1ba5A00', '4fdyA01', '1b43A02', '1zzhA02']

    base_path = "/scratch/wl45/may_week1_2020"
    # simulationType = "mass_iterative_run_traditional"
    optimizationBasePath = f"{base_path}/example_optimization"

    if args.mode == 1:
        optimizationFolder = "iter0_example"
        runOnServer = False
        # lastN_frame = 50
        # run_n = 2
        n_decoys = 20
        trial_name = optimizationFolder
        opt_folder = "~/opt"
        phi_cmd=f"python3 {opt_folder}/compute_phis.py -m 0 " + "{}"
        # phi_cmd="python3 ~/opt/compute_phis.py -m 7 {proteins}"
        gamma_cmd=f"python3 {opt_folder}/generate_gamma.py -m 1 -n {n_decoys}"
        complete_folder_list = "protein_list"

        # start from the pdbs in source folder, move them to database folder.
        fromFolder = f"{optimizationBasePath}/source"
        optimization_setupDatabase(pdb_list, fromFolder, optimizationBasePath)
        # setup the optimization folder.
        do(f"mkdir -p {optimizationBasePath}/{optimizationFolder}")
        cd(f"{optimizationBasePath}/{optimizationFolder}")
        optimization_setupFolder(phi_list, pdb_list)
        optimization_setupDecoys_shuffle(n_decoys, runOnServer=runOnServer)
        optimization_compute_phi_and_gamma(complete_folder_list, pdb_list,
                runOnServer=runOnServer, phi_cmd=phi_cmd, gamma_cmd=gamma_cmd)
        do(f"optimization_analyze.py {trial_name}")
        do(f"gg_server.py -d plot -l {trial_name}")

    if args.mode == 3:
        optimizationFolder = "iter0_test_server"
        runOnServer = True
        # lastN_frame = 50
        # run_n = 2
        n_decoys = 100
        trial_name = optimizationFolder
        phi_cmd="python3 ~/opt/compute_phis.py -m 0 {}"
        # phi_cmd="python3 ~/opt/compute_phis.py -m 7 {proteins}"
        gamma_cmd=f"python3 ~/opt/generate_gamma.py -m 1 -n {n_decoys}"
        complete_folder_list = "protein_list"

        # start from the pdbs in source folder, move them to database folder.
        fromFolder = f"{optimizationBasePath}/source"
        optimization_setupDatabase(pdb_list, fromFolder, optimizationBasePath)
        # setup the optimization folder.
        do(f"mkdir -p {optimizationBasePath}/{optimizationFolder}")
        cd(f"{optimizationBasePath}/{optimizationFolder}")
        optimization_setupFolder(phi_list, pdb_list)
        optimization_setupDecoys_shuffle(n_decoys, runOnServer=runOnServer, template=base_casp14_slurm)
        optimization_compute_phi_and_gamma(complete_folder_list, pdb_list,
                runOnServer=runOnServer, phi_cmd=phi_cmd, gamma_cmd=gamma_cmd, template=base_casp14_slurm)
        do(f"optimization_analyze.py {trial_name}")
        do(f"gg_server.py -d plot -l {trial_name}")
    # if args.mode == 1:
    #     random.seed(5)
    #     selected = random.sample(pdb_list, 5)
    #     print(selected)
    #     for pdb in selected:
    #         do(f"cp ../cath_dataset_shuffle_optimization/database/dompdb/{pdb}.pdb source/")

