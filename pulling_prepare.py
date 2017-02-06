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
import numpy

parser = argparse.ArgumentParser(
    description="Prepare the data for run and analysis. \
                Codes here only need run once")

# parser.add_argument("protein", help="the name of protein")
# parser.add_argument("template", help="the name of template file")
parser.add_argument("-t", "--test", help="test ", action="store_true", default=False)
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("--distance", action="store_true", default=False)
parser.add_argument("--replace", action="store_true", default=False)
parser.add_argument("--summary", action="store_true", default=False)
parser.add_argument("--make_metadata", action="store_true", default=False)
parser.add_argument("--qnqc", action="store_true", default=False)
parser.add_argument("--data", action="store_true", default=False)
parser.add_argument("-m", "--mode", type=int, default=1)
parser.add_argument("-s", "--switch", type=int, default=1)
args = parser.parse_args()


if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

# compute distance by "read dump file"
# if(args.test):
#     for i in range(40):
#         print(i)
#         cd(str(i))
#         do("read_dump_file.py")
#         cd("..")


def replace(TARGET, FROM, TO):
    do("sed -i.bak 's/{}/{}/g' {}".format(FROM,TO,TARGET))


def extract_data():
    # do("tail -n+3 energy.log | awk '{print $NF}' > etotal")
    do("awk '{print $17}' energy.dat  > etotal")
    do("head -n 6000 etotal | tail -n 2000 > etotal_half")
    do("head -n 6000 qn | tail -n 2000 > qn_half")
    do("head -n 6000 qc | tail -n 2000 > qc_half")
    do("head -n 6000 qo | tail -n 2000 > qo_half")
    do("paste qn qc etotal qo | tail -n 4000 > data")
    do("paste qn_half qc_half etotal_half qo_half > halfdata")


server_run = """\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=01:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun {}
"""

if(args.test):
    for i in range(40):
        print(i)
        cd(str(i))
        # do("cp ~/opt/pulling/qo.slurm .")
        # do("sbatch qo.slurm")
        extract_data()
        cd("..")

if(args.summary):
    if(args.mode == 1):
        with open("data", "w") as out:
            out.write("step, qw, run, energy\n")
            for i in range(40):
                print(i)
                with open(str(i)+"/halfdata") as f:
                    step = 0
                    for line in f:
                        step += 1
                        qn, qc, qw, energy = line.split()
                        out.write("{}, {}, run_{}, {}\n".format(step, qw, i, energy))
                        # out.write(str(n)+", "+qw+", run_"+str(i)+", "+energy+"\n"
    if(args.mode == 2):
        with open("data", "w") as out:
            out.write("step, qw, run, energy\n")
            for i in range(40):
                print(i)
                with open(str(i)+"/halfdata") as f:
                    step = 0
                    for line in f:
                        step += 1
                        qn, qc, energy, qw = line.split()
                        out.write("{}, {}, run_{}, {}\n".format(step, qw, i, energy))
                        # out.write(str(n)+", "+qw+", run_"+str(i)+", "+energy+"\n"
if(args.qnqc):
    if(args.mode == 1):
        n = 40
        temp = 350
        cwd = os.getcwd()
        for i in range(n):
            cd("simulation/{}/{}".format(temp, i))
            do("cp ~/opt/pulling/qnqc.slurm .")
            do("sbatch qnqc.slurm")
            # with open("server_run.slurm", "w") as f:
            #     f.write(server_run.format("read_dump_file.py"))
            # do("sbatch server_run.slurm")
            cd(cwd)
    if(args.mode == 2):
        array = []
        cwd = os.getcwd()
        print(cwd)
        with open('folder_list', 'r') as ins:
            for line in ins:
                target = line.strip('\n')
                t1 = "simulation/" + target + "/simulation/0"
                array.append(t1)
                # t2 = "simulation/" + target + "/simulation/1"
                # array.append(t2)
        for i in array:
            os.chdir(i)
            os.system("pwd")
            os.system("cp ~/opt/pulling/qnqc.slurm .")
            os.system("sbatch qnqc.slurm")
            os.chdir(cwd)

if(args.data):
    if(args.mode == 4):
        n = 40
        cwd = os.getcwd()
        for i in range(n):
            print(str(i))
            cd("simulation/300/{}".format(i))
            do("awk '{print $2}' wham.dat > qw")
            do("awk '{print $6}' wham.dat > energy")
            do("paste qn qc qw energy | tail -n 4000 > halfdata")
            cd(cwd)
    if(args.mode == 3):
        n = 40
        temp = 350
        cwd = os.getcwd()
        for i in range(n):
            print(str(i))
            cd("simulation/{}/{}".format(temp, i))
            do("head -n 5800 wham.dat | tail -n 4000 | awk '{print $2}' > qw")
            do("head -n 5800 wham.dat | tail -n 4000 | awk '{print $5}' > e")
            do("head -n 5800 qn | tail -n 4000 > qn_half")
            do("head -n 5800 qc | tail -n 4000 > qc_half")
            # do("paste qn qc | head -n 5800 | tail -n 4000 > qnqc")

            do("paste qn_half qc_half qw e > halfdata")
            cd(cwd)
    if(args.mode == 1):
        n = 40
        cwd = os.getcwd()
        for i in range(n):
            cd("simulation/300/{}".format(i))
            do("tail -n+3 energy.log | awk '{print $NF}' > energy")
            do("paste qn qc distance energy | tail 4000 > halfdata")
            cd(cwd)
    if(args.mode == 2):
        array = []
        cwd = os.getcwd()
        print(cwd)
        with open('folder_list', 'r') as ins:
            for line in ins:
                target = line.strip('\n')
                t1 = "simulation/" + target + "/simulation/0"
                array.append(t1)
                # t2 = "simulation/" + target + "/simulation/1"
                # array.append(t2)
        for i in array:
            os.chdir(i)
            os.system("pwd")
            # do("tail -n+3 energy.log | awk '{print $NF}' > energy")
            do("sed '/^#/ d' x.colvars.traj | awk 'NR % 10 == 1'  | awk '{print $2}' > x")
            do("paste qn qc x energy | tail -n 4000 > halfdata")
            do("paste qn qc x energy -d ',' | tail -n 4000 > test_data")
            # do("sed -i '1iqn,qc,x,energy' test_data")

            os.chdir(cwd)



if(args.make_metadata):
    if(args.mode == 4):
        kconstant = 400   # double the k constant
        temp = 350
        q0 = 0.12
        metadata = open("metadatafile", "w")
        for i in range(40):
            q = q0 + i*0.02
            target = "../simulation/300/" + str(i) + "/halfdata {} {} {:.2f}\n".format(temp, kconstant, q)
            # target = "../simulation/350/" + str(i) + "/data {} {} {:.2f}\n".format(temp, kconstant, q)
            metadata.write(target)
        metadata.close()
    if(args.mode == 1):
        kconstant = 2000   # double the k constant
        temp = 350
        q0 = 0.12
        metadata = open("metadatafile", "w")
        for i in range(40):
            q = q0 + i*0.02
            target = "../simulation/350/" + str(i) + "/halfdata {} {} {:.2f}\n".format(temp, kconstant, q)
            # target = "../simulation/350/" + str(i) + "/data {} {} {:.2f}\n".format(temp, kconstant, q)
            metadata.write(target)
        metadata.close()
    elif(args.mode == 2):
            kconstant = 0.02   # double the k constant
            metadata = open("metadatafile", "w")
            with open('folder_list', 'r') as ins:
                for line in ins:
                    target = line.strip(' \n')
                    temp = target.split("_")[1]
                    x = target.split("_")[3]
                    # print(temp)
                    t1 = "/scratch/wl45/project/freeEnergy_2xov/pullingDistance/simulation/" + target + "/simulation/0/halfdata {} {} {}\n".format(temp, kconstant, x)
                    metadata.write(t1)
                    # elif(args.mode == 2):
                    # t2 = "/scratch/wl45/freeEnergy_2xov/pullingDistance/simulation/" + target + "/simulation/1/halfdata {} {} {}\n".format(temp, kconstant, x)
                    # metadata.write(t2)
            metadata.close()
    elif(args.mode == 3):
            kconstant = 0.04   # double the k constant
            metadata = open("metadatafile", "w")
            with open('folder_list', 'r') as ins:
                for line in ins:
                    target = line.strip(' \n')
                    temp = target.split("_")[1]
                    x = target.split("_")[3]
                    # print(temp)
                    t1 = "/scratch/wl45/project/freeEnergy_2xov/pullingDistance_v3/simulation/" + target + "/simulation/0/halfdata {} {} {}\n".format(temp, kconstant, x)
                    metadata.write(t1)
                    # elif(args.mode == 2):
                    # t2 = "/scratch/wl45/freeEnergy_2xov/pullingDistance/simulation/" + target + "/simulation/1/halfdata {} {} {}\n".format(temp, kconstant, x)
                    # metadata.write(t2)
            metadata.close()

if(args.replace):
    target = "2xov.in"
    replace(target, "TEMPERATURE", "300")
    replace(target, "RANDOM", str(randint(1, 10**6)))

if(args.distance):
    do("read_dump_file.py")
