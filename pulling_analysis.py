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
import numpy
# Useful codes
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# tr " " "\n"
# sed 1d
# sort -u -k 3
# sed -e 's/+T//'
mypath = os.environ["PATH"]
os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
my_env = os.environ.copy()

parser = argparse.ArgumentParser(description="This is my playground for current project")

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
args = parser.parse_args()


if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir



if(args.test):
    folder_list = glob.glob("T_*")
    with open("complete_folder_list","w") as f:
        for i in folder_list:
            f.write(i+"\n")

if(args.mutation):
    print("Pulling Free energy batch compute")
    force_list = numpy.arange(1,2.5,0.1)
    dimension = 1
    mut_list = ["a206g", "l155a"]
    do("mkdir -p one_d_mutation")
    cd("one_d_mutation")
    for mut in mut_list:
        for i in range(1,3):
            for force in force_list:
                folder = str(i) + "_force_" + str(force) + "_" + mut
                do("mkdir -p "+folder)
                cd(folder)
                do("cp ../../folder_list .")
                cmd = "make_metadata.py --pulling --server -m {} --protein {}".format(i, mut)
                do(cmd)
                do("cp ~/opt/pulling/freeEnergy.slurm .")
                do(
                    "sed -i.bak 's/FORCE/" +
                    str(force) +
                    "/g' freeEnergy.slurm")
                do(
                    "sed -i.bak 's/DIMENSION/" +
                    str(dimension) +
                    "/g' freeEnergy.slurm")
                do("sbatch freeEnergy.slurm")
                cd("..")
    cd("..")


if(args.qnqc):
    array = []
    cwd = os.getcwd()
    print(cwd)
    with open('folder_list_jan16', 'r') as ins:
        for line in ins:
            target = line.strip('\n')
            t1 = target + "/simulation/0"
            array.append(t1)
    for i in array:
        os.chdir(i)
        os.system("pwd")
        # os.system("cp ~/opt/pulling/qnqc.slurm .")
        # os.system("sbatch qnqc.slurm")
        os.system("tail -n+3 energy.log > energy")
        os.system("head -n 6000 energy > energy_all")
        os.system("awk '{print $17}' energy_all > etotal_all")
        os.system("tail -n 4000 energy_all > energy_half")
        os.system("awk '{print $17}' energy_half > etotal")
        os.system("sed '/^#/ d' x.colvars.traj > test")
        os.system("awk 'NR % 10 == 1 ' test > x")
        os.system("head -n 6000 x > x_all")
        os.system("awk '{print $2}' x_all > myx_all")
        os.system("tail -n 4000 x_all > x_half")
        os.system("awk '{print $2}' x_half > myx_half")
        os.system("paste etotal myx_half > halfdata")
        os.system("head -n 6000 qn > qn_all")
        os.system("head -n 6000 qc > qc_all")
        os.system("head -n 6000 qc2 > qc2_all")
        os.system("tail -n 4000 qn_all > qn_half")
        os.system("tail -n 4000 qc_all > qc_half")
        os.system("tail -n 4000 qc2_all > qc2_half")
        os.system("paste qn_all qc_all qc2_all etotal_all myx_all > all_data")
        os.system("paste qn_half qc_half qc2_half etotal myx_half > halfdata")
        # uniq energy.log > test
        os.system("mv halfdata halfdata_back")
        os.system("awk '{print $0, $1-$2}' halfdata_back > halfdata")
        os.chdir(cwd)

if(args.pulling):
    print("Pulling Free energy batch compute")
    force_list = numpy.arange(1.5,2.5,0.1)
    # force_list = numpy.arange(1,1.2,0.1)
    dimension = args.dimension

    side_list = ["n_term", "c_term"]
    for side in side_list:
        do("mkdir -p "+side)
        cd(side)
        for i in range(1,2):
            for force in force_list:
                folder = str(i) + "_force_" + str(force)
                do("mkdir -p "+folder)
                cd(folder)
                do("cp ../../folder_list folder_list")
                cmd = "make_metadata.py --pulling3 --server -m {}".format(i)
                do(cmd)
                do("cp ~/opt/pulling/{0}.slurm freeEnergy.slurm".format(side))
                do(
                    "sed -i.bak 's/FORCE/" +
                    str(force) +
                    "/g' freeEnergy.slurm")
                do(
                    "sed -i.bak 's/DIMENSION/" +
                    str(dimension) +
                    "/g' freeEnergy.slurm")
                do("sbatch freeEnergy.slurm")
                cd("..")
                # cmd = "python2 ~/opt/pulling_compute-pmf.py -f {}".format(force)
        cd("..")
# if(args.pulling):
#     print("Pulling Free energy batch compute")
#     force_list = numpy.arange(1,2.5,0.1)
#     dimension = args.dimension
#     do("mkdir -p n_term")
#     cd("one_d")
#     for i in range(1,3):
#         for force in force_list:
#             folder = str(dimension) + "d_" + str(i) + "_force_" + str(force)
#             do("mkdir -p "+folder)
#             cd(folder)
#             do("cp ../../folder_list .")
#             cmd = "make_metadata.py --pulling3 --server -m {}".format(i)
#             do(cmd)
#             do("cp ~/opt/pulling/freeEnergy2.slurm freeEnergy.slurm")
#             do(
#                 "sed -i.bak 's/FORCE/" +
#                 str(force) +
#                 "/g' freeEnergy.slurm")
#             do(
#                 "sed -i.bak 's/DIMENSION/" +
#                 str(dimension) +
#                 "/g' freeEnergy.slurm")
#             do("sbatch freeEnergy.slurm")
#             cd("..")
#             # cmd = "python2 ~/opt/pulling_compute-pmf.py -f {}".format(force)
#     cd("..")

if(args.pulling2):
    print("Pulling Free energy batch compute")
    force_list = numpy.arange(1,2.5,0.1)
    dimension = args.dimension
    do("mkdir -p one_d")
    cd("one_d")
    for i in range(1,2):
        for force in force_list:
            folder = str(dimension) + "d_" + str(i) + "_force_" + str(force)
            do("mkdir -p "+folder)
            cd(folder)
            do("cp ../../folder_list .")
            cmd = "make_metadata.py --pulling3 --server -m {}".format(i)
            do(cmd)
            do("cp ~/opt/pulling/freeEnergy2.slurm freeEnergy.slurm")
            do(
                "sed -i.bak 's/FORCE/" +
                str(force) +
                "/g' freeEnergy.slurm")
            do(
                "sed -i.bak 's/DIMENSION/" +
                str(dimension) +
                "/g' freeEnergy.slurm")
            do("sbatch freeEnergy.slurm")
            cd("..")
            # cmd = "python2 ~/opt/pulling_compute-pmf.py -f {}".format(force)
    cd("..")
