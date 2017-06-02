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
#SBATCH --mem-per-cpu=5G
#SBATCH --time=23:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun python2 ~/opt/pulling_compute-pmf.py {}
"""

if(args.freeEnergy):
    arg = ""
    nsample = 2000
    force = 0
    temp_arg = "-f {} -nsamples {}".format(force, nsample)
    if(args.mode == 9):
        folder_name = "multi_temp_2"
        do("mkdir "+folder_name)
        cd(folder_name)
        temp_list = [135, 160, 185, 210]
        move_data_to_wham(temp_list)
        write_simulation_list(temp_list)
        get_total_x(temp_list)
        sim_list = 't135 t160 t185 t210'
        temp_list = '135 160 185 210'
        do("mult_calc_cv.sc . '{}' 20 '{}' 150 350 10 30 200 0 0.95 2xov q".format(sim_list, temp_list))
        cd("..")
    if(args.mode == 8):
        arg = "-b 3 -e 4 -d 2 -v1 1 -v1n 20 -v2 2 -v2n 20 -f 0 -nsamples 10000"
    if(args.mode == 7):
        arg = "-b 3 -e 4 -d 1 -v1 3 -v1n 30 -f 0 -nsamples 10000"
    if(args.mode == 6):
        arg = "-b 3 -e 4 -d 1 -v1 3 -v1n 30 -f 0 -nsamples 4000"
    if(args.mode == 1):
        arg = "-b 3 -e 4 -d 2 -v1 1 -v1n 30 -v2 2 -v2n 30 -f 1.7 -nsamples 2000"
    if(args.mode == 2):
        arg = "-b 3 -e 4 -d 2 -v1 1 -v1n 20 -v2 2 -v2n 20 -f 1.7 -nsamples 4000"
    if(args.mode == 3):
        arg = "-b 3 -e 4 -d 2 -v1 1 -v1n 20 -v2 2 -v2n 20 -f 0 -nsamples 4000"
    if(args.mode == 4):
        folder_name = "wham"
        cd(folder_name)
        do("rm all_wham.dat")
        for i in range(40):
            do("cat ../{}/halfdata >> all_wham.dat".format(i))
        os.system("awk '{print $3}' all_wham.dat > Qw_total")
        os.system("awk '{print $3}' all_wham.dat > p_total")
        os.system("awk '{print $1}' all_wham.dat > qn_total")
        os.system("awk '{print $2}' all_wham.dat > qc_total")
        os.system("awk '{print $4}' all_wham.dat > e_total")
        os.system("cp ~/opt/wham_analysis/*.m .")
        cd("..")
        # os.system("~/opt/script/wham/fused_calc_cv.sc {} top7 50 400 350 450 5 50 100 0 0.98".format(folder_name))
    if(args.mode == 5):
        arg = "-b 3 -e 4 -d 1 -v1 3 -v1n 40 -f {} -nsamples 4000".format(args.force)
    with open("freeEnergy.slurm", "w") as f:
        f.write(freeEnergy.format(arg))

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


#
# if(args.qnqc):
#     array = []
#     cwd = os.getcwd()
#     print(cwd)
#     with open('complete_folder_list', 'r') as ins:
#         for line in ins:
#             target = line.strip('\n')
#             t1 = "simulation/" + target + "/simulation/0"
#             array.append(t1)
#             t2 = "simulation/" + target + "/simulation/1"
#             array.append(t2)
#     for i in array:
#         os.chdir(i)
#         os.system("pwd")
#         if(args.mode == 2):         # cal qn, qc first.
#             os.system("cp ~/opt/pulling/qnqc.slurm .")
#             os.system("sbatch qnqc.slurm")
#         if(args.mode == 1):         # default mode. assemble halfdata
#             # do("tail -n+3 energy.log | awk '{print $NF}' > energy")
#             # do("sed '/^#/ d' x.colvars.traj | awk 'NR % 10 == 1'  | awk '{print $2}' > x")
#             # do("paste qn qc x energy | tail -n 2000 > halfdata")
#             do("paste qn qc x energy -d ',' > test_data")
#             # do("sed -i '1iqn,qc,x,energy' test_data")
#             # os.system("tail -n+3 energy.log > energy")
#             # os.system("head -n 6000 energy > energy_all")
#             # os.system("awk '{print $17}' energy_all > etotal_all")
#             # os.system("tail -n 4000 energy_all > energy_half")
#             # os.system("awk '{print $17}' energy_half > etotal")
#             # os.system("sed '/^#/ d' x.colvars.traj > test")
#             # os.system("awk 'NR % 10 == 1 ' test > x")
#             # os.system("head -n 6000 x > x_all")
#             # os.system("awk '{print $2}' x_all > myx_all")
#             # os.system("tail -n 4000 x_all > x_half")
#             # os.system("awk '{print $2}' x_half > myx_half")
#             # os.system("paste etotal myx_half > halfdata")
#             # os.system("head -n 6000 qn > qn_all")
#             # os.system("head -n 6000 qc > qc_all")
#             # os.system("head -n 6000 qc2 > qc2_all")
#             # os.system("tail -n 4000 qn_all > qn_half")
#             # os.system("tail -n 4000 qc_all > qc_half")
#             # os.system("tail -n 4000 qc2_all > qc2_half")
#             # os.system("paste qn_all qc_all qc2_all etotal_all myx_all > all_data")
#             # os.system("paste qn_half qc_half qc2_half etotal myx_half > halfdata")
#             # # uniq energy.log > test
#             # os.system("mv halfdata halfdata_back")
#             # os.system("awk '{print $0, $1-$2}' halfdata_back > halfdata")
#             os.chdir(cwd)
