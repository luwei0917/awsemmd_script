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
parser.add_argument("--continue_run", action="store_true", default=False)
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
#SBATCH --time=02:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun {}
"""

if(args.test):
    cd("simulation")
    a = glob.glob("*")
    print(a)
    # for i in range(40):
    #     print(i)
    #     cd(str(i))
    #     do("paste qw energy | tail -n 4000 > q_e")
    #     # do("cp ~/opt/pulling/qo.slurm .")
    #     # do("sbatch qo.slurm")
    #     # extract_data()
    #     cd("..")

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
    if(args.mode == 3):
        cd("simulation")
        with open("data", "w") as out:
            out.write("step, qw, run\n")
            for i in range(20):
                print(i)
                with open(str(i)+"/wham.dat") as f:
                    next(f)
                    for line in f:
                        step, qw, *rest = line.split()
                        out.write("{}, {}, run_{}\n".format(step, qw, i))
                        # out.write(str(n)+", "+qw+", run_"+str(i)+", "+energy+"\n"
    if(args.mode == 4):
        n = 40
        with open("data", "w") as out:
            out.write("step, qn, qc, dis, qw, run, energy\n")
            for i in range(n):
                print(i)
                cd(str(i))
                do("awk '{print $2}' addforce.dat > dis")
                do("paste qn qc dis wham.dat| tail -n+2 | head -n 6000 > data")
                # do("paste qn qc dis wham.dat| tail -n+2 > data")
                cd("..")
                with open(str(i)+"/data") as f:
                    for line in f:
                        qn, qc, dis, step, qw, *rest, energy = line.split()
                        out.write("{}, {}, {}, {}, {}, run_{}, {}\n".format(step, qn, qc, dis, qw, i, energy))
                        # out.write(str(n)+", "+qw+", run_"+str(i)+", "+energy+"\n"
    if(args.mode == 5):
        with open("data", "w") as out:
            out.write("i, step, qw, target, run, energy\n")
            q = 0.1

            for i in range(40):
                print(i)
                q += 0.02
                count = i*6000
                with open(str(i)+"/wham.dat") as f:
                    next(f)
                    for line in f:
                        count += 1
                        step, qw, *rest, energy = line.split()
                        out.write("{}, {}, {}, {}, {}, {}\n".format(count, step, qw, q, i, energy))
                        # out.write(str(n)+", "+qw+", run_"+str(i)+", "+energy+"\n"
if(args.qnqc):
    if(args.mode == 4):
        n = 80
        # temp_list = [300,350]
        # temp_list = [250,275, 325]
        # temp_list = [200]
        # temp_list = [0, 1, 2]
        # temp_list = [3]
        # temp_list = [4, 5]
        # temp_list = [6]
        # temp_list = [7]
        # temp_list = [8, 9]
        temp_list = [2]
        temp_list = [0]
        temp_list = [1]
        run_list = [10, 11, 12, 13]
        # temp_list = ['300', '200', '250']
        cwd = os.getcwd()
        for temp in run_list:
            for i in range(n):
                cd("simulation/{}/{}".format(i, temp))
                do("cp ../2xov.pdb .")
                do("cp ~/opt/pulling/qnqc.slurm .")
                do("sbatch qnqc.slurm")
                # with open("server_run.slurm", "w") as f:
                #     f.write(server_run.format("read_dump_file.py"))
                # do("sbatch server_run.slurm")
                cd(cwd)
    if(args.mode == 1):
        n = 20
        # temp_list = [300,350]
        # temp_list = [250,275, 325]
        # temp_list = [200]
        temp_list = [135, 160, 185, 210]
        # temp_list = ['300', '200', '250']
        cwd = os.getcwd()
        for temp in temp_list:
            for i in range(n):
                cd("simulation/{}/{}".format(temp, i))
                do("cp 2xov_translocon.pdb 2xov.pdb")
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
                t1 = "simulation/" + target + "/"
                array.append(t1)
                # t2 = "simulation/" + target + "/simulation/1"
                # array.append(t2)
        for i in array:
            os.chdir(i)
            os.system("pwd")
            os.system("cp ~/opt/pulling/qnqc.slurm .")
            os.system("sbatch qnqc.slurm")
            os.chdir(cwd)
    if(args.mode == 3):
        n = 40
        cwd = os.getcwd()
        for i in range(n):
            # cd("simulation/{}".format(i))
            cd("{}".format(i))
            do("cp ~/opt/pulling/qnqc.slurm .")
            do("sbatch qnqc.slurm")
            # with open("server_run.slurm", "w") as f:
            #     f.write(server_run.format("read_dump_file.py"))
            # do("sbatch server_run.slurm")
            cd(cwd)
if(args.data):
    if(args.mode == 8):
        n = 80
        # temp_list = [300,350]
        # temp_list = [250,275, 325]
        # temp_list = [200]
        # run_list = [0, 1, 2]
        # run_list = [3]
        # run_list = [4, 5]
        # run_list = [6]
        run_list = [7, 8, 9]
        run_list = [0, 1, 2]
        run_list = [0]
        run_list = [0, 1]
        run_list = [10, 11, 12, 13]
        # temp_list = ['300', '200', '250']
        cwd = os.getcwd()
        for run in run_list:
            for i in range(n):
                cd("simulation/{}/{}".format(i, run))
                do("awk '{print $1}' addforce.dat > steps")
                do("awk '{print $2}' addforce.dat > distance")
                do("awk '{print $17}' energy.dat > energy")
                do("paste -d, steps distance qn qc energy > data")
                # do("paste -d, steps distance distance distance energy > data")
                cd(cwd)
                # with open("server_run.slurm", "w") as f:
                #     f.write(server_run.format("read_dump_file.py"))
                # do("sbatch server_run.slurm")
    if(args.mode == 7):
        n = 20
        temp_list = ['300', '200', '250']
        cwd = os.getcwd()
        for temp in temp_list:
            for i in range(n):
                print(str(i))
                cd("simulation/{}/{}".format(temp, i))
                do("awk '{print $2}' wham.dat > qw")
                do("awk '{print $6}' wham.dat > energy")
                do("paste qn qc qw energy | tail -n 2000 > data")
                # do("tail -n 2000 data > small_data")
                # do("head -n 5800 wham.dat | tail -n 4000 | awk '{print $2}' > qw")
                # do("head -n 5800 wham.dat | tail -n 4000 | awk '{print $5}' > e")
                # do("head -n 5800 qn | tail -n 4000 > qn_half")
                # do("head -n 5800 qc | tail -n 4000 > qc_half")
                # do("paste qn qc | head -n 5800 | tail -n 4000 > qnqc")
                cd(cwd)
    if(args.mode == 6):
        n = 20
        temp_list = [135, 160, 185, 210]
        cwd = os.getcwd()
        for temp in temp_list:
            for i in range(n):
                print(str(i))
                cd("simulation/{}/{}".format(temp, i))
                do("awk '{print $2}' wham.dat > qw")
                do("awk '{print $6}' wham.dat > energy")
                do("paste qn qc qw energy | tail -n 10000 > data")
                do("tail -n 2000 data > small_data")
                # do("head -n 5800 wham.dat | tail -n 4000 | awk '{print $2}' > qw")
                # do("head -n 5800 wham.dat | tail -n 4000 | awk '{print $5}' > e")
                # do("head -n 5800 qn | tail -n 4000 > qn_half")
                # do("head -n 5800 qc | tail -n 4000 > qc_half")
                # do("paste qn qc | head -n 5800 | tail -n 4000 > qnqc")
                cd(cwd)
    if(args.mode == 5):
        target = "all_halfdata"
        do("awk '{print $1}' %s > qn" % (target))
        do("awk '{print $2}' %s > qc" % (target))
        do("awk '{print $3}' %s > p_total" % (target))
        do("awk '{print $4}' %s > e_total" % (target))
    if(args.mode == 4):
        n = 40
        temp_list = [250, 275, 300, 325, 350]
        cwd = os.getcwd()
        target = "multi_temp"
        for temp in temp_list:
            cd("simulation/{}".format(temp))
            do("ls */halfdata | sort -g | xargs cat > data")
            do("awk '{print $1}' data > ../../%s/qn_t%i" % (target, temp))
            do("awk '{print $2}' data > ../../%s/qc_t%i" % (target, temp))
            do("awk '{print $3}' data > ../../%s/q_t%i" % (target, temp))
            do("awk '{print $4}' data > ../../%s/energy_t%i" % (target, temp))
            cd(cwd)
        # cd(target)
        # with open("sim_list", "w") as f:
        #     for temp in temp_list:
        #         f.write("t{}\n".format(temp))
        # with open("T_list", "w") as f:
        #     for temp in temp_list:
        #         f.write("{}\n".format(temp))
        # cd(cwd)
    if(args.mode == 3):
        n = 40
        # temp_list = [250, 275, 300, 325, 350]
        temp_list = [200]
        temp_list = [135, 160, 185, 210]
        cwd = os.getcwd()
        for temp in temp_list:
            for i in range(n):
                print(str(i))
                cd("simulation/{}/{}".format(temp, i))
                do("awk '{print $2}' wham.dat > qw")
                do("awk '{print $6}' wham.dat > energy")
                do("paste qn qc qw energy | tail -n 4000 > halfdata")
                # do("head -n 5800 wham.dat | tail -n 4000 | awk '{print $2}' > qw")
                # do("head -n 5800 wham.dat | tail -n 4000 | awk '{print $5}' > e")
                # do("head -n 5800 qn | tail -n 4000 > qn_half")
                # do("head -n 5800 qc | tail -n 4000 > qc_half")
                # do("paste qn qc | head -n 5800 | tail -n 4000 > qnqc")
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
                t1 = "simulation/" + target + "/"
                array.append(t1)
                # t2 = "simulation/" + target + "/simulation/1"
                # array.append(t2)
        for i in array:
            os.chdir(i)
            os.system("pwd")
            do("tail -n+3 energy.log | awk '{print $NF}' > energy")
            do("sed '/^#/ d' x.colvars.traj | awk '{print $2}' > x")
            do("awk '{print $2}' wham.dat > qw")
            # do("sed '/^#/ d' x.colvars.traj | awk 'NR % 10 == 1'  | awk '{print $2}' > x")
            do("paste qn qc x energy qw| tail -n 4000 > halfdata")
            do("paste qn qc x energy qw -d ',' | tail -n 4000 > test_data")
            # do("sed -i '1iqn,qc,x,energy' test_data")

            os.chdir(cwd)



if(args.make_metadata):
    if(args.mode == 5):
        kconstant = 300   # double the k constant

        q0 = 0.0
        metadata = open("metadatafile", "w")
        for i in range(20):
            q = q0 + i*0.05
            # temp_list = [135, 160, 185, 210]
            temp_list = [160]
            for temp in temp_list:
                target = "../simulation/{}/".format(temp) + str(i) + "/small_data {} {} {:.2f}\n".format(temp, kconstant, q)
                metadata.write(target)
        metadata.close()
    if(args.mode == 4):
        kconstant = 800   # double the k constant

        q0 = 0.12
        metadata = open("metadatafile", "w")
        for i in range(40):
            q = q0 + i*0.02
            temp_list = [250, 275, 300, 325, 350]
            for temp in temp_list:
                target = "../simulation/300/" + str(i) + "/halfdata {} {} {:.2f}\n".format(temp, kconstant, q)
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
                    t1 = "/scratch/wl45/project/freeEnergy_2xov/pullingDistance/complete0.02/simulation/" + target + "/halfdata {} {} {}\n".format(temp, kconstant, x)
                    # t1 = "/scratch/wl45/project/freeEnergy_2xov/pullingDistance/simulation/" + target + "/simulation/0/halfdata {} {} {}\n".format(temp, kconstant, x)
                    metadata.write(t1)
                    # elif(args.mode == 2):
                    # t2 = "/scratch/wl45/freeEnergy_2xov/pullingDistance/simulation/" + target + "/simulation/1/halfdata {} {} {}\n".format(temp, kconstant, x)
                    # metadata.write(t2)
            metadata.close()
    elif(args.mode == 3):
            kconstant = 0.04   # double the k constant
            metadata = open("metadatafile", "w")
            with open('../folder_list', 'r') as ins:
                for line in ins:
                    target = line.strip(' \n')
                    temp = target.split("_")[1]
                    x = target.split("_")[3]

                    cwd = os.getcwd()
                    t1 = cwd + "/../simulation/" + target + "/halfdata {} {} {}\n".format(temp, kconstant, x)
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


if(args.continue_run):
    folder_name = "continue_simulation_2"
    folder_name = "continue_simulation"
    do("mkdir {}".format(folder_name))
    # do("mkdir continue_simulation")
    cd(folder_name)
    n = 40
    simulation_steps = 6*1000*1000
    protein_name = "2xov"
    for i in range(n):
        do("mkdir {}".format(i))
        cd(str(i))
        do("cp -r ../../2xov/* .")
        # do("cp ../../continue_simulation/{0}/restart.12000000 .".format(i))
        do("cp ../../simulation/{0}/restart.6000000 .".format(i))
        do("cp ~/opt/pulling/2xov_continue_run.in 2xov.in")
        # do(
        #     "sed -i.bak 's/START_FROM/'" +
        #     "12000000" +
        #     "'/g' "+protein_name+".in")
        do(
            "sed -i.bak 's/START_FROM/'" +
            "6000000" +
            "'/g' "+protein_name+".in")
        seed(datetime.now())
        do(  # replace RANDOM with a radnom number
            "sed -i.bak 's/RANDOM/'" +
            str(randint(1, 10**6)) +
            "'/g' "+protein_name+".in")
        do(  # replace SIMULATION_STEPS with specific steps
            "sed -i.bak 's/SIMULATION_STEPS/'" +
            str(simulation_steps) +
            "'/g' "+protein_name+".in")
        cd("..")
