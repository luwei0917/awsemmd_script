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
parser.add_argument("-l", "--label", type=str, default="label")
parser.add_argument("-t", "--test", action="store_true", default=False)
args = parser.parse_args()

if args.test:
    do = print
else:
    do = os.system
cd = os.chdir

def replace(TARGET, FROM, TO):
    do("sed -i.bak 's/{}/{}/g' {}".format(FROM,TO,TARGET))

def getFromTerminal(CMD):
    return subprocess.Popen(CMD,stdout=subprocess.PIPE,shell=True).communicate()[0].decode()

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

# name = "T0953S2"
# name = "T0954"
# name = "T0955"
# name = "T0956"
# name = "T0957S1"
# name = "T0957S2"
# name = "T0958"
# name = "T0960"
name = "T0959"
if args.mode == 1:
    # name = "T0953S1"
    do(f"scp -r wl45@davinci.rice.edu:/scratch/xl23/casp13/human/result/{name}/post-processing {name}")

    # copy from one location
    my_from = f"/scratch/mc70/CASP13/{name}/run1/"
    my_to = f"{name}/awsem"
    cmd = "rsync -a --exclude='dump.lammpstrj' --exclude='slurm-*' --exclude='run.pdb' {} {}".format(my_from, my_to)
    print(cmd)
    os.system(cmd)
    # another possible location
    my_from = f"wl45@davinci.rice.edu:/scratch/mc70/CASP13/{name}/run1/"
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
    # queue = "interactive"
    queue = "ctbp"
    # queue = "commons"
    model_list = glob.glob(f"model.*")
    for model in model_list:
        print(model)
        i = int(model.split(".")[-1])
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
if args.mode == 4:
    do(f"python3 /home/wl45/opt/native_structure_finder.py {name}")
