#!/Users/weilu/anaconda/envs/3.5/bin/python3
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime

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
args = parser.parse_args()

# TODO:
# add clean command.
# test analysis, and fullfill missing anaylsis.


# protein_name = args.template.split('_', 1)[-1].strip('/')
n = args.number
protein_name = args.template.strip('/')
if args.steps != -1:
    simulation_steps = args.steps * 10**6
    warm_up_steps = 10**5
else:  # -1 means a test run with 10000 steps
    simulation_steps = 10**4
    warm_up_steps = 10**3
    n = 1  # also set n to be 1
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
    if(platform.system() == 'Darwin'):
        os.system("lmp_serial < "+protein_name+".in")
    elif(platform.system() == 'Linux'):
        os.system(  # replace RANDOM with a radnom number
            "sed 's/NUMBER/'" +
            str(i) +
            "'/g' " "~/opt/run.slurm > run.slurm_"+str(i))
        os.system("sbatch run.slurm_"+str(i))
    else:
        print("system unkown")
    os.chdir("../..")

# print("hello world")
