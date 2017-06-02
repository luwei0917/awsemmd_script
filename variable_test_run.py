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
import numpy as np
import fileinput
# from run_parameter import *
parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    do see the difference variable make \
    run simulation")

parser.add_argument("template", help="the name of template file")
parser.add_argument("-d", "--debug", action="store_true", default=False)
args = parser.parse_args()
protein_name = args.template.strip('/')

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

# protein_name = args.template.split('_', 1)[-1].strip('/')
# os.system("cp ~/opt/variable_test_run.py .")

# run_slurm = '''\
# #!/bin/bash
# #SBATCH --job-name=CTBP_WL
# #SBATCH --account=ctbp-common
# #SBATCH --partition=ctbp-common
# #SBATCH --ntasks=1
# #SBATCH --mem-per-cpu=1G
# #SBATCH --time=1-00:00:00
# #SBATCH --mail-user=luwei0917@gmail.com
# #SBATCH --mail-type=FAIL
# echo "My job ran on:"
# echo $SLURM_NODELIST
# srun ~/build/brian/adjustable_z_dependence/lmp_serial -in 2xov_{}.in
# '''

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
srun ~/build/brian/adjustable_z_dependence/lmp_serial -in 2xov_{}.in
'''

fileName = "2xov_multi.in"
start_from = "read_data data.2xov"
# rg_list = [0, 1, 5, 10]
# force_list = [2.0]
# memb_k_list = [0, 1, 5, 10]

# rg_list = [0, 1, 2, 5]
# force_list = [0.0, 1.0, 2.0, 3.0]
# memb_k_list = [0, 1, 2, 5]

# rg_list = [0, 1, 2, 5]
# force_list = [0.0, 3.0]
# memb_k_list = [0, 1, 2, 5]

# rg_list = [0, 0.1, 1, 5, 10]
# force_list = [0.0, 3.0]
# memb_k_list = [0, 0.1, 1, 5, 10]

# rg_list = [0, 0.1, 1]
# force_list = ["ramp"]
# memb_k_list = [0, 0.1, 1]

# rg_list = [0, 0.1, 0.5, 1, 2]
# rg_list = [3, 4]
# force_list = ["ramp"]
# memb_k_list = [0, 0.1, 1, 2, 5, 10]

rg_list = [1.5, 2, 2.5, 3, 3.5]
force_list = ["ramp"]
memb_k_list = [1.5, 2, 2.5, 3, 3.5, 4]

# rg_list = [0.01]
for memb_k in memb_k_list:
    for force in force_list:
        for rg in rg_list:
            i = 0
            folder_name = "memb_{}_force_{}_rg_{}".format(memb_k, force, rg)
            do("cp -r 2xov " + folder_name)
            cd(folder_name)
            # fixFile = "fix_backbone_coeff_go.data"
            fixFile = "fix_backbone_coeff_single.data"
            with fileinput.FileInput(fixFile, inplace=True, backup='.bak') as file:
                for line in file:
                    print(line.replace("MY_MEMB_K", str(memb_k)), end='')
            with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
                for line in file:
                    print(line.replace("MY_FORCE", str(force)), end='')
            with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
                for line in file:
                    print(line.replace("MY_RG", str(rg)), end='')
            with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
                for line in file:
                    print(line.replace("START_FROM", start_from), end='')
            do("cp 2xov_multi.in 2xov_{}.in".format(i))
            do(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/NUMBER/'" +
                str(int(i)) +
                "'/g' 2xov_{}.in".format(i))
            do("mkdir -p {}".format(i))
            do(  # replace RANDOM with a radnom number
                "sed -i.bak 's/RANDOM/'" +
                str(randint(1, 10**6)) +
                "'/g' *.in")
            with open("run_{}.slurm".format(i), "w") as r:
                r.write(run_slurm.format(i))
            do("sbatch run_0.slurm")
            cd("..")
# folder_list = open('folder_list', 'w')
# distance_list = np.arange(20, 350, 5)
# # temp_list = np.arange(250, 400, 50)
# temp_list = [300]
# folder_name = ""
# cwd = os.getcwd()
# os.system("mkdir -p simulation")
# for temp in temp_list:
#     pre_folder_name = "T_"+str(temp)
#     for distance in distance_list:
#         folder_name = pre_folder_name + "_D_"+str(distance)
#         folder_list.write(folder_name+"\n")
#         os.system("mkdir -p simulation/" + folder_name)
#         os.system("cp -r "+protein_name+"/* simulation/"+folder_name+"/")
#         os.chdir("simulation")
#         os.chdir(folder_name)
#         os.system(
#             "sed -i.bak 's/TEMPERATURE/'" +
#             str(temp) +
#             "'/g' "+protein_name+".in")
#         os.system(
#             "sed -i.bak 's/DISTANCE/'" +
#             str(distance) +
#             "'/g' colvars.x")
#         # os.system(
#         #     "sed -i.bak 's/MGamma/'" +
#         #     str(MGamma) +
#         #     "'/g' fix_backbone_coeff.data")
#         do("run.py " + protein_name + " -o -s 5 -i")
#         # os.system("run.py " + protein_name + "/ -o -n 1 -s 6 -i")
#         cd(cwd)

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

# folder_list.close()
# print("hello world")
