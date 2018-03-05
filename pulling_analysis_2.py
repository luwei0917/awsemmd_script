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
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("--protein", default="2xov")
parser.add_argument("--dimension", type=int, default=2)
parser.add_argument("-m", "--mode", type=int, default=1)
parser.add_argument("-s", "--save", action="store_true", default=False)
parser.add_argument("-r", "--reproduce", default=None)
parser.add_argument("--force", type=float, default=1.0)
parser.add_argument("-p", "--patch", type=int, default=1)
parser.add_argument("--commons", type=int, default=0)
parser.add_argument("--nsample", type=int, default=2500)
parser.add_argument("--submode", type=int, default=-1)
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
#SBATCH --mem-per-cpu=30G
#SBATCH --time=23:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun python2 ~/opt/pulling_compute-pmf.py {}
"""

if args.commons:
    freeEnergy = freeEnergy.replace("ctbp-common", "commons")
    freeEnergy = freeEnergy.replace("--time=23:00:00", "--time=08:00:00")


if args.mode >= 9:
    print("two d")
    nsample = args.nsample
    # force_list = [0.0, 0.1, 0.2]
    force_list = [0.0]
    # force_list = [0.0, 0.1, 0.2]
    for force in force_list:
        # force = 1
        temp_arg = "-f {} -nsamples {}".format(force, nsample)
        folder_name = "force_{}".format(force)
        do("mkdir -p "+folder_name)
        cd(folder_name)
        # do("make_metadata.py -m 1")
        do("cp ../metadatafile .")
        if args.submode == 5 and args.mode == 9:
            arg = "-b 3 -e 1 -d 1 " + temp_arg
            arg += " -v1 3 -v1n 50"

        if args.submode == 5 and args.mode == 10:
            arg = "-b 3 -e 1 -d 1 " + temp_arg
            arg += " -v1 2 -v1n 50"
        if args.submode == 5 and args.mode == 11:
            arg = "-b 3 -e 1 -d 2 " + temp_arg
            arg += " -v1 3 -v1n 30 -v2 2 -v2n 30"

        if args.submode == 5 and args.mode == 12:
            arg = "-b 3 -e 1 -d 1 " + temp_arg
            arg += " -v1 4 -v1n 50"

        if args.submode == 5 and args.mode == 13:
            arg = "-b 3 -e 1 -d 2 " + temp_arg
            arg += " -v1 4 -v1n 30 -v2 2 -v2n 30"
        if args.submode == 5 and args.mode == 14:
            arg = "-b 3 -e 1 -d 2 " + temp_arg
            arg += " -v1 4 -v1n 30 -v2 3 -v2n 30"
        arg += " -ti 20 -st 380 -et 540 -ev 5-185 -pb y -ss y"
        with open("freeEnergy.slurm", "w") as f:
            f.write(freeEnergy.format(arg))
        do("sbatch freeEnergy.slurm")
        cd("..")
