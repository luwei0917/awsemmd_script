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
parser.add_argument("--force", type=int, default=0)
parser.add_argument("-p", "--patch", type=int, default=1)
parser.add_argument("--commons", type=int, default=0)
parser.add_argument("--nsample", type=int, default=2500)
parser.add_argument("--submode", type=int, default=-1)
parser.add_argument("--subsubmode", type=int, default=-1)
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
#SBATCH --mem-per-cpu=20G
#SBATCH --time=23:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun python2 ~/opt/pulling_compute-pmf.py {}
"""

if args.commons:
    freeEnergy = freeEnergy.replace("ctbp-common", "commons")
    # freeEnergy = freeEnergy.replace("--time=23:00:00", "--time=08:00:00")

print("two d")
nsample = args.nsample
if args.force == 0:
    force_list = [0.0]
elif args.force == 1:
    force_list = [0.05, 0.02, 0.1, 0.2, 0.0]
elif args.force == 2:
    force_list = [0.2]
elif args.force == 2:
    force_list = [0.2, 0.3, 0.4]
elif args.force == 3:
    force_list = [0.13, 0.14, 0.15, 0.16]
elif args.force == 4:
    force_list = [0.0, 0.05, 0.1]
elif args.force == 5:
    force_list = [0.08, 0.09, 0.1]
elif args.force == 6:
    force_list = [0.1, 0.13, 0.07]
elif args.force == 7:
    force_list = [0.1, 0.15, 0.2, 0.25]
elif args.force == 8:
    force_list = [0.1, 0.2, 0.0]
elif args.force == 9:
    force_list = [0.1, 0.0]
elif args.force == 10:
    force_list = [0.0, 0.05, 0.1, 0.2, 0.3]
# force_list = [0.0, 0.1, 0.2]

# force_list = [0.0, 0.1, 0.2]
for force in force_list:
    # force = 1
    temp_arg = "-f {} -nsamples {}".format(force, nsample)
    folder_name = "force_{}".format(force)
    do("mkdir -p "+folder_name)
    cd(folder_name)
    # do("make_metadata.py -m 1")
    do("cp ../metadatafile .")
    if args.mode ==1:
        arg = "-b 2 -e 1 -d 1 " + temp_arg
        arg += " -v1 2 -v1n 50"

    if args.mode == 9:
        arg = "-b 3 -e 1 -d 1 " + temp_arg
        arg += " -v1 3 -v1n 50"
    if args.mode == 91:
        arg = "-b 3 -e 1 -d 1 " + temp_arg
        arg += " -v1 6 -v1n 50"
    if args.mode == 92:
        arg = "-b 3 -e 1 -d 1 " + temp_arg
        arg += " -v1 4 -v1n 50"
    if args.mode == 93:
        arg = "-b 3 -e 1 -d 1 " + temp_arg
        arg += " -v1 5 -v1n 50"
    if args.mode == 94:
        arg = "-b 3 -e 1 -d 1 " + temp_arg
        arg += " -v1 8 -v1n 50"
    if args.mode == 10:
        arg = "-b 3 -e 1 -d 1 " + temp_arg
        arg += " -v1 2 -v1n 50"
    if args.mode == 11:
        arg = "-b 3 -e 1 -d 2 " + temp_arg
        arg += " -v1 3 -v1n 30 -v2 2 -v2n 30"

    if args.mode == 12:
        arg = "-b 3 -e 1 -d 1 " + temp_arg
        arg += " -v1 4 -v1n 50"

    if args.mode == 13:
        arg = "-b 3 -e 1 -d 2 " + temp_arg
        arg += " -v1 4 -v1n 30 -v2 2 -v2n 30"

    if args.mode == 14:
        arg = "-b 3 -e 1 -d 2 " + temp_arg
        arg += " -v1 4 -v1n 30 -v2 3 -v2n 30"
    if args.mode == 15:
        arg = "-b 3 -e 1 -d 2 " + temp_arg
        arg += " -v1 4 -v1n 30 -v2 5 -v2n 30"

    if args.mode == 16:  # energy at 6
        arg = "-b 3 -e 6 -d 2 " + temp_arg
        arg += " -v1 4 -v1n 30 -v2 2 -v2n 30"
    if args.mode == 17:
        arg = "-b 3 -e 1 -d 2 " + temp_arg
        arg += " -v1 3 -v1n 40 -v2 5 -v2n 40"
    if args.mode == 18:
        arg = "-b 3 -e 1 -d 2 " + temp_arg
        arg += " -v1 4 -v1n 40 -v2 7 -v2n 40"
    if args.mode == 20:
        arg = "-b 3 -e 1 -d 2 " + temp_arg
        # arg += " -v1 3 -v1n 40 -v2 8 -v2n 40"
        arg += " -v1 8 -v1n 40 -v2 3 -v2n 40"
    if args.mode == 21:
        arg = "-b 3 -e 1 -d 2 " + temp_arg
        arg += " -v1 3 -v1n 40 -v2 5 -v2n 40"
    if args.mode == 22:
        arg = "-b 3 -e 1 -d 2 " + temp_arg
        arg += " -v1 3 -v1n 40 -v2 11 -v2n 40"
    if args.mode == 23:
        arg = "-b 3 -e 1 -d 2 " + temp_arg
        arg += " -v1 10 -v1n 40 -v2 5 -v2n 40"
    if args.mode == 24:
        arg = "-b 3 -e 1 -d 2 " + temp_arg
        arg += " -v1 8 -v1n 40 -v2 4 -v2n 40"
    if args.mode == 25:
        arg = "-b 3 -e 1 -d 2 " + temp_arg
        arg += " -v1 9 -v1n 40 -v2 4 -v2n 40"
    if args.mode == 26:
        arg = "-b 3 -e 1 -d 2 " + temp_arg
        arg += " -v1 10 -v1n 40 -v2 4 -v2n 40"
    if args.submode == 1:
        arg += " -ti 10 -st 350 -et 600 -p 5 -p 6 -p 7 -p 8 -p 9 -p 10 -p 11 -p 12 -pb y -ss y"
    if args.submode == 2:
        arg += " -ti 10 -st 350 -et 600 -ev 3 -pb y -ss y"
    if args.submode == 21 and not args.mode == 17:
        arg += " -ti 10 -st 360 -et 390 -ev 3 -pb y -ev 17-20 -p 6 -p 7 -ss y"
    if args.submode == 21 and args.mode == 17:
        arg += " -ti 10 -st 360 -et 390 -p 6 -p 7 -p 9 -p 10 -p 11 -p 12 -p 13 -p 14 -p 15 -p 16 -ss y"
    if args.submode == 22:
        arg += " -ti 10 -st 290 -et 350 -ev 3 -pb y"
    if args.submode == 23:
        # arg += " -ti 10 -st 260 -et 280 -ev 5-65 -pb y"
        if args.subsubmode != -1:
            start = 5 + args.subsubmode * 10
            end = start + 9
            if args.subsubmode == 17:
                end += 1
            arg += f" -ti 10 -st 260 -et 280 -ev {start}-{end} -pb y"
    if args.submode == 24:
        # arg += " -ti 10 -st 260 -et 280 -ev 5-65 -pb y"
        if args.subsubmode != -1:
            start = 5 + args.subsubmode * 5
            end = start + 4
            if args.subsubmode == 35:
                end += 1
            arg += f" -ti 10 -st 260 -et 280 -ev {start}-{end} -pb y"
            #  -p 186 -p 187 -p 188 -p 189 -p 190 -p 191 -p 192 -p 193
    # if args.submode == 24:
    #     arg += " -ti 10 -st 260 -et 280 -ev 66-116 -pb y"
    if args.submode == 25:
        arg += " -ti 10 -st 260 -et 280 -ev 117-185 -pb y"
    if args.submode == 26:
        arg += " -ti 10 -st 250 -et 350 -ev 5-8 -pb y"
    if args.submode == 27:
        arg += " -ti 10 -st 360 -et 390 -p 7 -ss y"
    if args.submode == 28:
        arg += " -ti 10 -st 360 -et 390 -p 9 -p 10 -ss y"
    if args.submode == 29:
        arg += " -ti 10 -st 330 -et 420 -p 10 -ss y -pb y -ev 20-23"
    if args.submode == 30:
        arg += " -ti 10 -st 330 -et 420 -p 13 -ss y -pb y -ev 23-26"
    if args.submode == 31:
        arg += " -ti 1 -st 373 -et 374 -p 6 -p 7 -p 8 -p 9 -ss y"
    if args.submode == 32:
        arg += " -ti 3 -st 370 -et 373 -p 6 -p 7 -p 8 -p 9 -ss y"
    if args.submode == 3:
        arg += " -ti 10 -st 400 -et 600 -ev 5 -pb y -ss y"
    if args.submode == 4:
        arg += " -ti 10 -st 250 -et 500"
    if args.submode == 5:
        arg += " -ti 20 -st 380 -et 540 -ev 5-185 -pb y -ss y"
    if args.submode == 6:
        arg += " -ti 2 -st 256 -et 264 -p 5 -p 6 -p 7 -p 8 -p 9 -p 10 -p 11 -p 12 -pb y"
    with open("freeEnergy.slurm", "w") as f:
        f.write(freeEnergy.format(arg))
    do("sbatch freeEnergy.slurm")
    cd("..")
