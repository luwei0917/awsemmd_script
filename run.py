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

# from run_parameter import *
parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    automatic copy the template file, \
    run simulation and analysis")

parser.add_argument("template", help="the name of template file")
parser.add_argument("-n", "--number", type=int, default=20,
                    help="Number of simulation run")
parser.add_argument("-s", "--steps", type=int, default=8,
                    help="Simulation steps in unit of million,\
                    default is 8 million, -1 means test run")
parser.add_argument("-r", "--read", help="Read from config file",
                    action="store_true")
parser.add_argument("-ws", "--warmSteps", type=int, default=20,
                    help="Warmup Simulation steps in unit of hundred thousand,\
                    default is 2 million")
parser.add_argument("-t", "--test", help="test mode",
                    action="store_true")
parser.add_argument("-c", "--copy",
                    help="copy the restart file before run",
                    action="store_true")
parser.add_argument("-o", "--offAuto", help="turn off from Read from \
                    config file", action="store_true", default=False)
parser.add_argument("-i", "--inplace", help="change in this folder",
                    action="store_true", default=False)
parser.add_argument("-m", "--mode",
                    help="mode 2 is dependence run",
                    type=int, default=1)
args = parser.parse_args()
# TODO:
# add clean command.
# test analysis, and fullfill missing anaylsis.

# protein_name = args.template.split('_', 1)[-1].strip('/')
n = args.number
protein_name = args.template.strip('/')
if args.steps == -1:  # smallest run for debug.
    simulation_steps = 10**4
    warm_up_steps = 10**4
    n = 1  # also set
elif args.test:  # test run
    simulation_steps = 50 * 10**3
    warm_up_steps = 50 * 10**3
else:
    simulation_steps = args.steps * 10**6
    warm_up_steps = args.warmSteps * 10**5

config = open('config.py', 'w')
config.write("protein_name = '%s'\nnumber_of_run = %d\nsimulation_steps = %d\n\
warm_up_steps = %d\n" % (protein_name, n, simulation_steps, warm_up_steps))
config.close()
if(not args.offAuto):
    exec(open("variables.dat").read())
    print(TSTART, TEND)


def set_up():
    seed(datetime.now())
    os.system(  # replace SIMULATION_STEPS with specific steps
        "sed -i.bak 's/WARM_UP_STEPS/'" +
        str(warm_up_steps) +
        "'/g' "+protein_name+".in")
    os.system(  # replace RANDOM with a radnom number
        "sed -i.bak 's/RANDOM/'" +
        str(randint(1, 10**6)) +
        "'/g' "+protein_name+".in")
    os.system(  # replace SIMULATION_STEPS with specific steps
        "sed -i.bak 's/SIMULATION_STEPS/'" +
        str(simulation_steps) +
        "'/g' "+protein_name+".in")
    if args.steps == -1:
        os.system(  # replace TEMPERATURE with specific steps
            "sed -i.bak 's/Q0/'" +
            str(0.5) +
            "'/g' fix_qbias_coeff.data")
        os.system(  # replace TEMPERATURE with specific steps
            "sed -i.bak 's/TEMPERATURE/'" +
            str(350) +
            "'/g' "+protein_name+".in")
    if(not args.offAuto):
            os.system(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/TSTART/'" +
                str(TSTART) +
                "'/g' "+protein_name+".in")
            os.system(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/TEND/'" +
                str(TEND) +
                "'/g' "+protein_name+".in")



def set_up2():
    seed(datetime.now())
    os.system(  # replace SIMULATION_STEPS with specific steps
        "sed -i.bak 's/WARM_UP_STEPS/'" +
        str(warm_up_steps) +
        "'/g' *.in")
    os.system(  # replace RANDOM with a radnom number
        "sed -i.bak 's/RANDOM/'" +
        str(randint(1, 10**6)) +
        "'/g' *.in")
    os.system(  # replace SIMULATION_STEPS with specific steps
        "sed -i.bak 's/SIMULATION_STEPS/'" +
        str(simulation_steps) +
        "'/g' *.in")
    if args.steps == -1:
        os.system(  # replace TEMPERATURE with specific steps
            "sed -i.bak 's/Q0/'" +
            str(0.5) +
            "'/g' fix_qbias_coeff.data")
        os.system(  # replace TEMPERATURE with specific steps
            "sed -i.bak 's/TEMPERATURE/'" +
            str(350) +
            "'/g' "+protein_name+".in")
    if(not args.offAuto):
            os.system(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/TSTART/'" +
                str(TSTART) +
                "'/g' "+protein_name+".in")
            os.system(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/TEND/'" +
                str(TEND) +
                "'/g' "+protein_name+".in")


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
srun ~/bin/lmp_serial -in 2xov_{}.in
'''

def set_up3():
    seed(datetime.now())
    with open("my_loop_submit.bash", "w") as f:

        for i in range(2,4):
            with open("run_{}.slurm".format(i), "w") as r:
                r.write(run_slurm.format(i))
            if(i != 2):
                dependency = "--dependency=afterany:$jobid"
            else:
                dependency = ""
            f.write("jobid=`sbatch "+dependency+" run_{}.slurm".format(i) + " | tail -n 1 | awk '{print $4}'`\necho $jobid\n")
            os.system("cp 2xov_multi.in 2xov_{}.in".format(i))
            os.system(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/START_FROM_N/'" +
                str(int(5e6*(i))) +
                "'/g' 2xov_{}.in".format(i))
            os.system(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/NUMBER/'" +
                str(int(i)) +
                "'/g' 2xov_{}.in".format(i))
            os.system("mkdir -p {}".format(i))
    os.system(  # replace RANDOM with a radnom number
        "sed -i.bak 's/RANDOM/'" +
        str(randint(1, 10**6)) +
        "'/g' *.in")
    if(not args.offAuto):
            os.system(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/TSTART/'" +
                str(TSTART) +
                "'/g' "+protein_name+".in")
            os.system(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/TEND/'" +
                str(TEND) +
                "'/g' "+protein_name+".in")


def set_up3_real():
    seed(datetime.now())
    with open("my_loop_submit.bash", "w") as f:
        with open("run_0.slurm", "w") as r:
            r.write(run_slurm.format(0))
        f.write("jobid=`sbatch run_0.slurm | tail -n 1 | awk '{print $4}'`\necho $jobid\n")
        os.system("mkdir 0")
        for i in range(1,4):
            with open("run_{}.slurm".format(i), "w") as r:
                r.write(run_slurm.format(i))
            f.write("jobid=`sbatch --dependency=afterany:$jobid run_{}.slurm".format(i) + " | tail -n 1 | awk '{print $4}'`\necho $jobid\n")
            os.system("cp 2xov_multi.in 2xov_{}.in".format(i))
            os.system(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/START_FROM_N/'" +
                str(int(5e6*(i))) +
                "'/g' 2xov_{}.in".format(i))
            os.system(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/NUMBER/'" +
                str(int(i)) +
                "'/g' 2xov_{}.in".format(i))
            os.system("mkdir {}".format(i))
    os.system(  # replace RANDOM with a radnom number
        "sed -i.bak 's/RANDOM/'" +
        str(randint(1, 10**6)) +
        "'/g' *.in")
    if(not args.offAuto):
            os.system(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/TSTART/'" +
                str(TSTART) +
                "'/g' "+protein_name+".in")
            os.system(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/TEND/'" +
                str(TEND) +
                "'/g' "+protein_name+".in")


def batch_run():
    if(platform.system() == 'Darwin'):
        os.system("/Users/weilu/bin/lmp_serial < "+protein_name+".in")
        # os.system("/Users/weilu/Research/Build/lammps-9Oct12_modified/src/lmp_serial \
        # < "+protein_name+".in")
    elif(platform.system() == 'Linux'):
        os.system("cp ~/opt/run.slurm run.slurm")
        os.system(  # replace PROTEIN with pdb name
            "sed -i.bak 's/PROTEIN/'" +
            protein_name +
            "'/g' run.slurm")
        os.system("sbatch run.slurm")
        sleep(0.2)  # Time in seconds.
    else:
        print("system unkown")



def batch_run2():
    if(platform.system() == 'Darwin'):
        os.system("/Users/weilu/bin/lmp_serial < "+protein_name+".in")
        # os.system("/Users/weilu/Research/Build/lammps-9Oct12_modified/src/lmp_serial \
        # < "+protein_name+".in")
    elif(platform.system() == 'Linux'):
        os.system("bash loopsubmit.bash")
        sleep(0.2)  # Time in seconds.
    else:
        print("system unkown")


def batch_run3():
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
    set_up()
    batch_run()
else:
    for i in range(n):
        # simulation set up
        os.system("mkdir -p simulation/"+str(i))
        os.system("cp -r "+protein_name+"/* simulation/"+str(i))
        os.chdir("simulation/"+str(i))
        if(args.mode == 3):
            set_up3()
            batch_run3()
            os.chdir("../..")
        elif(args.mode == 2):
            set_up2()
            batch_run2()
            os.chdir("../..")
        else:
            set_up()
            batch_run()
            os.chdir("../..")

# print("hello world")
