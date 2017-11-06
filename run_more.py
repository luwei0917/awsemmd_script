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

parser.add_argument("-n", "--number", type=int, default=1,
                    help="# of simulation run, default: 1")
parser.add_argument("--restart", type=int, default=0,
                    help="start from? default: 0")
parser.add_argument("--runs", type=int, default=1,
                    help="then do how many runs?, default: 1")
parser.add_argument("-s", "--steps", type=int, default=5,
                    help="How many steps in unit of million,\
                    per run, default: 5")
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("-m", "--mode", type=int, default=1)
parser.add_argument("-i", "--inplace", type=int, default=0)
args = parser.parse_args()

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir
# TODO:
# add clean command.
# test analysis, and fullfill missing anaylsis.

if args.mode == 0:
    run_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun ~/build/brian/src/lmp_serial -in 2xov_{}.in
    '''
elif args.mode == 1:
    run_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun ~/build/brian/bug_fix_jun04/lmp_serial -in 2xov_{}.in
    '''

def set_up():
    seed(datetime.now())
    with open("my_loop_submit.bash", "w") as f:
        steps = args.steps*1e6
        runs = args.runs
        restart = args.restart
        for ii in range(runs):
            i = ii + restart
            with open("run_{}.slurm".format(i), "w") as r:
                r.write(run_slurm.format(i))
            if(i != restart):
                dependency = "--dependency=afterany:$jobid"
            else:
                dependency = ""
            f.write("jobid=`sbatch "+dependency+" run_{}.slurm".format(i) + " | tail -n 1 | awk '{print $4}'`\necho $jobid\n")
            if i == 0:
                start_from = "read_data data.2xov"
            else:
                start_from = "read_restart restart." + str(int(steps*i))
            with fileinput.FileInput("2xov_multi.in", inplace=True, backup='.bak') as file:
                for line in file:
                    print(line.replace("reset_timestep", "# reset_timestep"), end='')
            do("cp 2xov_multi.in 2xov_{}.in".format(i))
            fileName = "2xov_{}.in".format(i)

            with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
                for line in file:
                    print(line.replace("START_FROM", start_from), end='')
            os.system(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/NUMBER/'" +
                str(int(i)) +
                "'/g' 2xov_{}.in".format(i))
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
    batch_run()
else:
    n = args.number
    cwd = os.getcwd()
    for i in range(n):
        if args.restart == 0:
            do("mkdir -p simulation")
            do("cp -r 2xov simulation/{}".format(i))
        cd("simulation/"+str(i))
        set_up()
        batch_run()
        cd(cwd)

# print("hello world")
