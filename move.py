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

# parser.add_argument("template", help="the name of template file")
# parser.add_argument("-n", "--number", type=int, default=20,
                    #  help="Number of simulation run")
parser.add_argument("-o", "--offAuto", help="turn off from Read from \
                    config file", action="store_true")
args = parser.parse_args()

# n = args.number
# protein_name = args.template.strip('/')

# folder_list = glob.glob("*")
# print(folder_list)
# sys.exit()
folder_list = ['T089_ha', 'T089_he', 'T120_ha', 'T120_he', 'T251_ha', 'T251_he', 'top7_ha', 'top7_he', '1UBQ_ha', '1UBQ_he']
os.system("mkdir -p weilu")
for folder in folder_list:
    os.chdir(folder)
    if(not args.offAuto):
        exec (open("config.py").read())
        n = number_of_run
        steps = simulation_steps
    os.system("mkdir -p ../weilu/"+folder+"/lowest_energy")
    os.system("mkdir -p ../weilu/"+folder+"/best_q")
    os.system("sort analysis/list_of_max_q > ../weilu/q_"+folder+".dat")
    for i in range(n):
        # move
        os.chdir("analysis/"+str(i))
        os.system("cp chosen.pdb ../../../weilu/"+folder+"/best_q/"+str(i)+".pdb")
        # os.system("cp ~/opt/plot_scripts/print_chosen.pml .")
        # os.system("/usr/local/bin/pymol -qc -r print_chosen.pml")
        # os.system("cp chosen.png ../../results/chosen_"+str(i)+".png")
        # os.system("cp final.png ../../results/final_"+str(i)+".png")
        # os.system("cp final.pdb ../../results/final_"+str(i)+".pdb")
        # os.system("cp final.txt ../../results/final_"+str(i)+".txt")
        # os.system("cp lowest_energy.pdb \
        #     ../../results/lowest_energy/lowest_energy_" + str(i)+".pdb")
        os.chdir("../..")
    os.chdir("..")
