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
parser.add_argument("--rerun",
                    type=int, default=1)
parser.add_argument("-m", "--mode", type=int, default=2)
parser.add_argument("--model", default="single")
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

# if args.mode == 1:
#     run_slurm = '''\
# #!/bin/bash
# #SBATCH --job-name=CTBP_WL
# #SBATCH --account=ctbp-common
# #SBATCH --partition=ctbp-common
# #SBATCH --ntasks=1
# #SBATCH --threads-per-core=1
# #SBATCH --mem-per-cpu=1G
# #SBATCH --time=1-00:00:00
# #SBATCH --mail-user=luwei0917@gmail.com
# #SBATCH --mail-type=FAIL
# echo "My job ran on:"
# echo $SLURM_NODELIST
# srun ~/build/brian/z_dependence/lmp_serial -in 2xov_{}.in
# '''
# if args.mode == 2:
#     run_slurm = '''\
# #!/bin/bash
# #SBATCH --job-name=CTBP_WL
# #SBATCH --account=ctbp-common
# #SBATCH --partition=ctbp-common
# #SBATCH --ntasks=1
# #SBATCH --threads-per-core=1
# #SBATCH --mem-per-cpu=1G
# #SBATCH --time=1-00:00:00
# #SBATCH --mail-user=luwei0917@gmail.com
# #SBATCH --mail-type=FAIL
# echo "My job ran on:"
# echo $SLURM_NODELIST
# srun /home/wl45/build/awsem_new_membrane/src/lmp_serial -in 2xov_{}.in
# '''

fileName = "2xov_multi.in"
# if args.rerun == 0:
#     start_from = "read_data data.2xov"
# if args.rerun == 1:
#     start_from = "read_restart restart.extended"
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

rg_list = [0, 0.1, 1, 2, 4, 8]
# rg_list = [0, 0.01, 0.04, 0.08, 0.1, 0.2, 0.5, 1]
force_list = ["ramp"]
memb_k_list = [0, 1, 2, 4, 8, 16]
i = args.rerun

simulation_steps = 5e7
# Defaults
# rg_list = [0]
# memb_k_list = [0]
# force_list = ["ramp"]
# force_ramp_rate_list = [1]
# repeat = 80


# rg_list = [0.08]
# force_list = [0.0, 0.02]
# force_list = [0.04, 0.06, 0.08]
# force_list = [0.045, 0.05, 0.055]
# force_list = [0.03, 0.07]
# force_list = [0.1]
# memb_k_list = [1]
# force_ramp_rate_list = [1]



# rg_list = [0, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.6, 1, 5]
# memb_k_list = [0, 1, 2, 5, 10, 20]
# force_list = ["ramp"]
# force_ramp_rate_list = [1]

# rg_list = [0, 0.02, 0.08, 0.16, 0.6, 1, 2, 5]
# memb_k_list = [0, 1, 2, 5, 10]
# force_list = ["ramp"]
# force_ramp_rate_list = [1]
# repeat = 2

# rg_list = [0, 0.04, 0.08]
# memb_k_list = [0, 1, 2, 5, 10]
# force_list = ["ramp"]
# force_ramp_rate_list = [1]
# repeat = 20

# rg_list = [0.08]
# memb_k_list = [1]
# force_list = [0.02, 0.03, 0.04]
# force_ramp_rate_list = [1]
# repeat = 100

# rg_list = [0.08]
# memb_k_list = [1]
# # force_list = [0.02, 0.03, 0.04, 0.05, 0.06]
# force_list = [0.025, 0.035, 0.045]

# force_list = [0.02, 0.03]
force_ramp_rate_list = [1]
repeat = 2
# repeat = 100


# rg_list = [0, 0.1, 0.5, 1, 1.5, 2]
# memb_k_list = [0, 1, 2, 3, 4]
# force_list = ["force_ramp"]
# force_ramp_rate_list = [1]
# repeat = 2


# force_ramp_rate_list = [1, 5, 10, 20]
# force_ramp_rate_list = [30, 40, 50, 60, 70, 80, 90, 100, 500, 1000]
# force_ramp_rate_list = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 500, 1000]
# force_ramp_rate_list = [1]
for force_ramp_rate in force_ramp_rate_list:
    for memb_k in memb_k_list:
        for force in force_list:
            for rg in rg_list:
                folder_name = "memb_{}_rg_{}".format(memb_k, rg)
                # folder_name = "memb_{}_force_{}_rg_{}".format(memb_k, force, rg)
                # folder_name = "rate_{}".format(force_ramp_rate)
                # folder_name = "force_{}".format(force)
                # if memb_k == 0 and rg == 0:
                #     continue
                print(folder_name)

                do("mkdir "+folder_name)
                do("cp -r 2xov " + folder_name + "/")
                cd(folder_name + "/2xov")
                if args.model == "go":
                    fixFile = "fix_backbone_coeff_go.data"
                if args.model == "single":
                    fixFile = "fix_backbone_coeff_single.data"
                with fileinput.FileInput(fixFile, inplace=True, backup='.bak') as file:
                    for line in file:
                        print(line.replace("MY_MEMB_K", str(memb_k)), end='')
                with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
                    for line in file:
                        tmp = line
                        tmp = tmp.replace("MY_FORCE", str(force))
                        tmp = tmp.replace("MY_RG", str(rg))
                        tmp = tmp.replace("RATE", str(force_ramp_rate))
                        tmp = tmp.replace("SIMULATION_STEPS", str(int(simulation_steps/force_ramp_rate)))
                        print(tmp, end='')
                cd("..")
                do("run.py -m 2 2xov -n {}".format(repeat))
                # do("run.py -m 2 2xov --start extended -n {}".format(repeat))
                cd("..")


                # do(  # replace SIMULATION_STEPS with specific steps
                #     "sed -i.bak 's/NUMBER/'" +
                #     str(int(i)) +
                #     "'/g' 2xov_{}.in".format(i))
                # do("mkdir -p {}".format(i))
                # do(  # replace RANDOM with a radnom number
                #     "sed -i.bak 's/RANDOM/'" +
                #     str(randint(1, 10**6)) +
                #     "'/g' *.in")
                # with open("run_{}.slurm".format(i), "w") as r:
                #     r.write(run_slurm.format(i))
                # do("sbatch " + "run_{}.slurm".format(i))
                # cd("..")




# force_list = ["ramp"]
# memb_k_list = [0, 0.5, 1, 2, 3, 5]
#
# for memb_k in memb_k_list:
#     for force in force_list:
#         rg = memb_k
#         i = 0
#         folder_name = "memb_{}_force_{}_rg_{}".format(memb_k, force, rg)
#         do("cp -r 2xov " + folder_name)
#         cd(folder_name)
#         # fixFile = "fix_backbone_coeff_go.data"
#         fixFile = "fix_backbone_coeff_single.data"
#         with fileinput.FileInput(fixFile, inplace=True, backup='.bak') as file:
#             for line in file:
#                 print(line.replace("MY_MEMB_K", str(memb_k)), end='')
#         with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
#             for line in file:
#                 print(line.replace("MY_FORCE", str(force)), end='')
#         with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
#             for line in file:
#                 print(line.replace("MY_RG", str(rg)), end='')
#
#         do("cp 2xov_multi.in 2xov_{}.in".format(i))
#         with fileinput.FileInput("2xov_{}.in".format(i), inplace=True, backup='.bak') as file:
#             for line in file:
#                 print(line.replace("START_FROM", start_from), end='')
#         do(  # replace SIMULATION_STEPS with specific steps
#             "sed -i.bak 's/NUMBER/'" +
#             str(int(i)) +
#             "'/g' 2xov_{}.in".format(i))
#         do("mkdir -p {}".format(i))
#         do(  # replace RANDOM with a radnom number
#             "sed -i.bak 's/RANDOM/'" +
#             str(randint(1, 10**6)) +
#             "'/g' *.in")
#         with open("run_{}.slurm".format(i), "w") as r:
#             r.write(run_slurm.format(i))
#         do("sbatch run_0.slurm")
#         cd("..")


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
