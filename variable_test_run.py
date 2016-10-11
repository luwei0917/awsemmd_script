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
protein_name = args.template.strip('/')
# protein_name = args.template.split('_', 1)[-1].strip('/')

folder_list = open('folder_list', 'w')
add_force_strengths = [-10, -5, -1]
for ForceStrength in add_force_strengths:
    folder_name = "ForceStrength"+str(ForceStrength)
    folder_list.write(folder_name+"\n")
    os.system("mkdir -p " + folder_name)
    os.system("cp -r "+protein_name+" "+folder_name+"/")
    os.chdir(folder_name)
    os.chdir(protein_name)
    os.system(  # replace ForceStrength with specific steps
        "sed -i.bak 's/ForceStrength/'" +
        str(ForceStrength) +
        "'/g' "+protein_name+".in")
    os.chdir("..")
    os.system("run.py "+protein_name + "-n 1")
    os.chdir("..")
# n = 5
# membrane_k = [1, 2, 3]
# rg_cylindrical_spring_constants = [1, 0.1]
# for i in range(n):
#     for MemK in membrane_k:
#         add_force_strengths = range(-3-MemK, -1-MemK)
#         pre_pre_folder_name = "MemK"+str(MemK)
#         for ForceStrength in add_force_strengths:
#             pre_folder_name = pre_pre_folder_name+"ForceStrength"+str(ForceStrength)
#             for SpringConstant in rg_cylindrical_spring_constants:
#                 # simulation set up
#                 folder_name = pre_folder_name+"SpringConstant"+str(SpringConstant)+"_"+str(i)+"/"
#                 os.system("mkdir -p "+folder_name)
#                 folder_list.write(folder_name+"\n")
#                 os.system("cp -r "+args.template+"* "+folder_name)
#                 os.chdir(folder_name)
#                 os.system(  # replace SIMULATION_STEPS with specific steps
#                     "sed -i.bak 's/WARM_UP_STEPS/'" +
#                     str(warm_up_steps) +
#                     "'/g' "+protein_name+".in")
#                 os.system(  # replace RANDOM with a radnom number
#                         "sed -i.bak 's/RANDOM/'" +
#                         str(randint(1, 10**6)) +
#                         "'/g' "+protein_name+".in")
#                 os.system(  # replace SIMULATION_STEPS with specific steps
#                         "sed -i.bak 's/SIMULATION_STEPS/'" +
#                         str(simulation_steps) +
#                         "'/g' "+protein_name+".in")
#                 os.system(  # replace SpringConstant with specific steps
#                         "sed -i.bak 's/SpringConstant/'" +
#                         str(SpringConstant) +
#                         "'/g' "+protein_name+".in")
#                 os.system(  # replace ForceStrength with specific steps
#                         "sed -i.bak 's/ForceStrength/'" +
#                         str(ForceStrength) +
#                         "'/g' "+protein_name+".in")
#                 os.system(  # replace ForceStrength with specific steps
#                         "sed -i.bak 's/MemK/'" +
#                         str(MemK) +
#                         "'/g' fix_backbone_coeff.data")
#             # if(platform.system() == 'Darwin'):
#             #     os.system("/Users/weilu/Documents/lammps-9Oct12_modified/src/lmp_serial \
#             #     < "+protein_name+".in")
#                 if(platform.system() == 'Darwin'):
#                     os.system("/Users/weilu/Documents/lammps-9Oct12_modified/src/lmp_serial \
#                     < "+protein_name+".in")
#                 elif(platform.system() == 'Linux'):
#                     os.system("cp ~/opt/run.slurm .")
#                     os.system(  # replace PROTEIN with pdb name
#                             "sed -i.bak 's/PROTEIN/'" +
#                             protein_name +
#                             "'/g' run.slurm")
#                     os.system("sbatch run.slurm")
#                 else:
#                     print("system unkown")
#                 os.chdir("..")

folder_list.close()
# print("hello world")
