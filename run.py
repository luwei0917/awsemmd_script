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
# from run_parameter import *
parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        automatic copy the template file, \
        run simulation and analysis")

parser.add_argument("template", help="the name of template file")
parser.add_argument("-n", "--number", type=int, default=20,
                    help="Number of simulation run")
parser.add_argument("-s", "--steps", type=int, default=4,
                    help="Simulation steps in unit of million,\
                    default is 4 million, -1 means test run")
parser.add_argument("-r", "--read", help="Read from config file",
                    action="store_true")
parser.add_argument("-ws", "--warmSteps", type=int, default=1,
                    help="Warmup Simulation steps in unit of hundred thousand,\
                    default is 1 hundred thousand")
args = parser.parse_args()
# TODO:
# add clean command.
# test analysis, and fullfill missing anaylsis.

# protein_name = args.template.split('_', 1)[-1].strip('/')
n = args.number
protein_name = args.template.strip('/')
if args.steps != -1:
    simulation_steps = args.steps * 10**6
    warm_up_steps = args.warmSteps * 10**5
else:  # -1 means a test run with 10000 steps
    simulation_steps = 10**4
    warm_up_steps = 10**3
    n = 1  # also set n to be 1

# imp.load_source('run_paramter.py', '')
if(args.read):
    exec (open("myconfig.py").read())
    print(n, x, y, type(y))
    sys.exit(0)
for i in range(n):
    seed(datetime.now())
# simulation set up
    os.system("mkdir -p simulation/"+str(i))
    os.system("cp -r "+args.template+"* simulation/"+str(i))
    os.chdir("simulation/"+str(i))
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
# if(platform.system() == 'Darwin'):
#     os.system("/Users/weilu/Documents/lammps-9Oct12_modified/src/lmp_serial \
#     < "+protein_name+".in")
    if(platform.system() == 'Darwin'):
        os.system("lmp_serial \
        < "+protein_name+".in")
    elif(platform.system() == 'Linux'):
        os.system("cp ~/opt/run.slurm .")
        os.system(  # replace PROTEIN with pdb name
                "sed -i.bak 's/PROTEIN/'" +
                protein_name +
                "'/g' run.slurm")
        os.system("sbatch run.slurm")
    else:
        print("system unkown")
    os.chdir("../..")

# print("hello world")
