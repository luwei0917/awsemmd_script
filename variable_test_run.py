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
        do see the difference variable make \
        run simulation")

parser.add_argument("template", help="the name of template file")
args = parser.parse_args()

# protein_name = args.template.split('_', 1)[-1].strip('/')
protein_name = args.template.strip('/')
simulation_steps = 4 * 10**6
warm_up_steps = 10 * 10**5

seed(datetime.now())
folder_name = ""
rg_cylindrical_spring_constants = [1, 0.1, 0.01, 0.001]


for SpringConstant in rg_cylindrical_spring_constants:
    # simulation set up
    folder_name += "SpringConstant"+str(SpringConstant)+"/"
    os.system("mkdir -p "+folder_name)
    os.system("cp -r "+args.template+"* "+folder_name)
    os.chdir(folder_name)
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
    os.system(  # replace SpringConstant with specific steps
            "sed -i.bak 's/SpringConstant/'" +
            str(SpringConstant) +
            "'/g' "+protein_name+".in")
# if(platform.system() == 'Darwin'):
#     os.system("/Users/weilu/Documents/lammps-9Oct12_modified/src/lmp_serial \
#     < "+protein_name+".in")
    if(platform.system() == 'Darwin'):
        os.system("/Users/weilu/Documents/lammps-9Oct12_modified/src/lmp_serial \
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
    os.chdir("..")
    folder_name = ""

# print("hello world")
