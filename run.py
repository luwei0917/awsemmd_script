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
from time import sleep
import fileinput

# from run_parameter import *
parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    automatic copy the template file, \
    run simulation and analysis")

parser.add_argument("protein", help="The name of the protein")
parser.add_argument("--name", default="simulation", help="Name of the simulation")
parser.add_argument("-n", "--number", type=int, default=1,
                    help="# of simulation run, default: 1")
parser.add_argument("--restart", type=int, default=0,
                    help="start from? default: 0")
parser.add_argument("--runs", type=int, default=1,
                    help="then do how many runs?, default: 1")
parser.add_argument("-s", "--steps", type=int, default=50,
                    help="How many steps in unit of hundred thousand(not millions),\
                    per run, default: 50")
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("-m", "--mode", type=int, default=2)
parser.add_argument("-i", "--inplace", action="store_true", default=False)
parser.add_argument("-f", "--force", type=float, default=1.0)
parser.add_argument("--start", default="native")
parser.add_argument("--commons", type=int, default=0)
args = parser.parse_args()

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

if args.mode == 3:
    run_slurm = '''\
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
echo "My job ran on:"
echo $SLURM_NODELIST
srun /home/wl45/build/lammps-16Mar18/src/lmp_serial -in {}_{}.in
    '''

if args.mode == 2:
    run_slurm = '''\
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
echo "My job ran on:"
echo $SLURM_NODELIST
srun /home/wl45/build/awsem_lipid_fluctuations/src/lmp_serial -in {}_{}.in
    '''
if args.mode == 1:
    run_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --constraint=opath
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun /home/wl45/build/awsem_lipid_fluctuations/src/lmp_serial -in {}_{}.in
    '''

if args.commons == 1:
    run_slurm = run_slurm.replace("ctbp-common", "commons")
proteinName = args.protein.strip("/.")
def set_up():
    seed(datetime.now())
    with open("my_loop_submit.bash", "w") as f:
        steps = args.steps*1e5
        runs = args.runs
        restart = args.restart
        for ii in range(runs):
            i = ii + restart
            with open("run_{}.slurm".format(i), "w") as r:
                r.write(run_slurm.format(proteinName, i))
            if(i != restart):
                dependency = "--dependency=afterany:$jobid"
            else:
                dependency = ""
            f.write("jobid=`sbatch "+dependency+" run_{}.slurm".format(i) + " | tail -n 1 | awk '{print $4}'`\necho $jobid\n")
            if i == 0:
                start_from = "read_data data.{}".format(proteinName)
                if args.start == "extended":
                    start_from = "read_restart restart.extended"
                if args.start == "crystal":
                    start_from = "read_data data.crystal"
            else:
                start_from = "read_restart restart." + str(int(steps*i))
            do("cp {0}_multi.in {0}_{1}.in".format(proteinName, i))
            fileName = "{0}_{1}.in".format(proteinName, i)
            # backbone_file = "fix_backbone_coeff_{}.data".format(args.name)
            with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
                for line in file:
                    tmp = line.replace("read_data data.crystal", "START_FROM")  # remove in future.
                    tmp = tmp.replace("langevin 800 800", "langevin 300 300")  # change temp, remove in future
                    tmp = tmp.replace("langevin 800 300", "langevin 300 300")  # change temp, remove in future
                    tmp = tmp.replace("START_FROM", start_from)
                    tmp = tmp.replace("MY_FORCE", str(args.force))
                    # tmp = tmp.replace("fix_backbone_coeff_er.data", backbone_file)
                    if i == 0:
                        tmp = tmp.replace("RESET_TIME_OR_NOT", "reset_timestep	0")
                    print(tmp, end='')

            os.system(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/NUMBER/'" +
                str(int(i)) +
                "'/g' {}_{}.in".format(proteinName, i))
            os.system("mkdir -p {}".format(i))
            os.system(  # replace RANDOM with a radnom number
                "sed -i.bak 's/RANDOM/'" +
                str(randint(1, 10**6)) +
                "'/g' *.in")



def batch_run():
    if(platform.system() == 'Darwin'):
        os.system("/Users/weilu/bin/lmp_serial < "+protein_name+".in")
        # os.system("/Users/weilu/Research/Build/lammps-9Oct12_modified/src/lmp_serial \
        # < "+protein_name+".in")
    elif(platform.system() == 'Linux'):
        os.system("bash my_loop_submit.bash")
        sleep(0.2)  # Time in seconds.
    else:
        print("system unkown")

if(args.inplace):
    print("inplace")
    set_up()
    # batch_run()
    do("~/build/brian/z_dependence/lmp_serial -in {}_multi.in".format(proteinName))
else:
    n = args.number
    cwd = os.getcwd()
    for i in range(n):
        if args.restart == 0:
            do("mkdir -p " + args.name)
            do("cp -r {} {}/{}".format(proteinName, args.name, i))
        cd(args.name + "/"+str(i))
        set_up()
        batch_run()
        cd(cwd)

# print("hello world")
