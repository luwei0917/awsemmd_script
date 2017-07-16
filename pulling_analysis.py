#!/usr/bin/env python3
# @Author: Wei Lu <weilu>
# @Date:   25-Jan-2017
# @Email:  wl45@rice.edu
# @Last modified by:   weilu
# @Last modified time: 25-Jan-2017
# @Copyright: Free


import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
import imp
from myPersonalFunctions import *
import glob
import numpy as np
import datetime
import pickle
# Useful codes
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# tr " " "\n"
# sed 1d
# sort -u -k 3
# sed -e 's/+T//'
mypath = os.environ["PATH"]
os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
my_env = os.environ.copy()

parser = argparse.ArgumentParser(description="Analysis code, need run multiple times")

# parser.add_argument("protein", help="the name of protein")
# parser.add_argument("template", help="the name of template file")
parser.add_argument("-t", "--test", help="test ", action="store_true", default=False)
parser.add_argument("--pulling", action="store_true", default=False)
parser.add_argument("--pulling2", action="store_true", default=False)
parser.add_argument("--qnqc", action="store_true", default=False)
parser.add_argument("--mutation", action="store_true", default=False)
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("--protein", default="2xov")
parser.add_argument("--dimension", type=int, default=2)
parser.add_argument("-f", "--freeEnergy", action="store_true", default=False)
parser.add_argument("-m", "--mode", type=int, default=1)
parser.add_argument("--rate", action="store_true", default=False)
parser.add_argument("-s", "--save", action="store_true", default=False)
parser.add_argument("-r", "--reproduce", default=None)
parser.add_argument("--force", type=float, default=1.0)
args = parser.parse_args()

if(args.reproduce):
    print("Reproducing!")
    with open(args.reproduce, "rb") as f:
        args = pickle.load(f)
        print(args)


if(args.save):
    print(os.getcwd())
    print("Saving")
    args.save = False
    with open("args"+datetime.datetime.now().strftime("%m%d-%H%M"), "wb") as f:
        pickle.dump(args, f)

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir



if(args.test):
    print("hello world")

if(args.rate):
    file_name = "test"
    with open("data") as f:
        names = next(f)


def move_data_to_wham(temp_list):
    for temp in temp_list:
        do("cp ../data/{}/data data".format(temp))
        do("awk '{print $1}' data > qn_t%i" % (temp))
        do("awk '{print $2}' data > qc_t%i" % (temp))
        do("awk '{print $3}' data > q_t%i" % (temp))
        do("awk '{print $4}' data > energy_t%i" % (temp))


def write_simulation_list(temp_list):
    with open("T_list", "w") as f:
        for temp in temp_list:
            f.write(str(temp)+"\n")
    with open("sim_list", "w") as f:
        for temp in temp_list:
            f.write("t"+str(temp)+"\n")


def get_total_x(temp_list):
    x_list = ["q", "qn", "qc", "energy"]
    for x in x_list:
        for temp in temp_list:
            do("cat {0}_t{1} >> {0}_total".format(x, temp))



freeEnergy = """\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=23:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun python2 ~/opt/pulling_compute-pmf.py {}
"""

if(args.mode == 1):
    nsample = 600
<<<<<<< HEAD
    force_list = np.linspace(0.5, 2, 0.1)
=======
<<<<<<< HEAD
    force_list = np.arange(0.5, 2, 0.1)
=======
    force_list = np.linespace(0.5, 2, 0.1)
>>>>>>> ddfb20ed2b095da89efb6979e848751d49f87505
>>>>>>> 1db191c4c33e171cb199ab2f9d1bced1a8590c7b
    for force in force_list:
        # force = 1
        temp_arg = "-f {} -nsamples {}".format(force, nsample)
        folder_name = "force_{}".format(force)
        do("mkdir "+folder_name)
        cd(folder_name)
        do("make_metadata.py -m 1")
        arg = "-b 2 -e 1 -d 1 -v1 2 -v1n 30 " + temp_arg
        with open("freeEnergy.slurm", "w") as f:
            f.write(freeEnergy.format(arg))
        do("sbatch freeEnergy.slurm")
