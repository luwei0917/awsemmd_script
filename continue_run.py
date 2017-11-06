#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import platform
import subprocess

mypath = os.environ["PATH"]
os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
my_env = os.environ.copy()

parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    continues previous simulation")

parser.add_argument("template", help="the name of template file")
parser.add_argument("-n", "--number", type=int, default=20,
                    help="Number of simulation run")
parser.add_argument("-m", "--movie", help="generate the movie",
                    action="store_true")
parser.add_argument("-p", "--plotOnly", help="only generate the plot",
                    action="store_true")
args = parser.parse_args()

list_of_max_q = []

n = args.number
protein_name = args.template.strip('/')

for i in range(n):
    os.chdir("continue_simulation/"+str(i))

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
    os.chdir("../..")
