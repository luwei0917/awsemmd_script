#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import glob

mypath = os.environ["PATH"]
os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
my_env = os.environ.copy()

parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        automatically analysis the simulation")

parser.add_argument("template", help="the name of template file")
# parser.add_argument("-n", "--number", type=int, default=20,
                    #  help="Number of simulation run")
parser.add_argument("-o", "--offAuto", help="turn off from Read from \
                    config file", action="store_true")
args = parser.parse_args()

# n = args.number
protein_name = args.template.strip('/')

# folder_list = glob.glob("*")
# print(folder_list)
# sys.exit()
if(not args.offAuto):
    exec (open("config.py").read())
    n = number_of_run
    steps = simulation_steps
os.system("mkdir -p "+protein_name+"/lowest_energy")
for i in range(n):
    # move
    os.chdir("analysis/"+str(i))
    # os.system("cp chosen.pdb ../../../weilu/"+folder+"/best_q/"+str(i)+".pdb")
    # os.system("cp ~/opt/plot_scripts/print_chosen.pml .")
    # os.system("/usr/local/bin/pymol -qc -r print_chosen.pml")
    # os.system("cp chosen.png ../../results/chosen_"+str(i)+".png")
    # os.system("cp final.png ../../results/final_"+str(i)+".png")
    # os.system("cp final.pdb ../../results/final_"+str(i)+".pdb")
    # os.system("cp final.txt ../../results/final_"+str(i)+".txt")
    os.system("cp lowest_energy.pdb \
        ../../"+protein_name+"/lowest_energy/lowest_energy_" + str(i)+".pdb")
    os.chdir("../..")
