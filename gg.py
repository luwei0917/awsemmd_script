#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import imp

mypath = os.environ["PATH"]
os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
my_env = os.environ.copy()

parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        automatically analysis the simulation")

# parser.add_argument("template", help="the name of template file")
args = parser.parse_args()

result_folder = "WeiLu_Aug_07"
os.system("mkdir -p "+result_folder)
protein_list = ['T089', 'T120', 'T251', 'top7', '1UBQ']
# sublist = ['_ha', '_he']
sublist = ['_lp', '_he_lp']
folder_list = []
for protein in protein_list:
    for sub in sublist:
        folder_list += [protein+sub]
print(folder_list)
# exit(1)

for folder in folder_list:
    print(folder)
    os.chdir(folder)
    exec (open("config.py").read())
    n = number_of_run
    steps = simulation_steps
    os.system("mkdir -p ../{}/".format(result_folder)+folder+"/best_q")
    os.system("sort analysis/list_of_max_q > ../{}/q_".format(result_folder)+folder+".dat")
    for i in range(n):
        # move
        os.chdir("analysis/"+str(i))
        os.system("cp chosen.pdb ../../../{}/".format(result_folder) + folder+"/best_q/"+str(i)+".pdb")
        os.chdir("../..")
    os.chdir("..")
# Useful codes
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# tr " " "\n"
# sed 1d
