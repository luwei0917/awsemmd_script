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

args = parser.parse_args()
# TODO:
# add clean command.
# test analysis, and fullfill missing anaylsis.

# protein_name = args.template.split('_', 1)[-1].strip('/')
n = args.number
protein_name = args.template.strip('/')


warm_up_steps = 2*10**6
simulation_steps = 6*10**6

# temp_list = [400, 500]
# temp_list = [250, 275, 325]
# temp_list = [300]
temp_list = [200]
n = 40

if(platform.system() == 'Darwin'):
    print("Are you sure?")
    exit()
elif(platform.system() == 'Linux'):
    pass
else:
    print("system unkown")

config = open('config.py', 'w')
config.write("protein_name = '%s'\nnumber_of_run = %d\nsimulation_steps = %d\n\
warm_up_steps = %d\nn = %d\n" % (protein_name, n, simulation_steps, warm_up_steps, n))
config.close()


# temp_list = [300, 350, 400]
# print(temp)

# simulation set up
for temp in temp_list:
    seed(datetime.now())
    os.system("mkdir -p simulation/"+str(temp))
    os.chdir("simulation/"+str(temp))
    q_bias_step = 0.02
    q0 = 0.1
    for i in range(n):
        q0 += q_bias_step
        os.system("mkdir -p "+str(i))
        os.system("cp -r ../../"+args.template+"/* "+str(i))
        os.chdir(str(i))
        os.system(  # replace TEMPERATURE with specific steps
            "sed -i.bak 's/Q0/'" +
            str(q0) +
            "'/g' fix_qbias_coeff.data")
        # os.system(  # replace TEMPERATURE with specific steps
        #     "sed -i.bak 's/Q0/'" +
        #     str(q0) +
        #     "'/g' fix_qbias_coeff1.data")
        os.system(  # replace TEMPERATURE with specific steps
            "sed -i.bak 's/TEMPERATURE/'" +
            str(temp) +
            "'/g' "+protein_name+".in")
        os.system(  # replace TEMPERATURE with specific steps
            "sed -i.bak 's/TSTART/'" +
            str(temp) +
            "'/g' "+protein_name+".in")
        os.system(  # replace TEMPERATURE with specific steps
            "sed -i.bak 's/TEND/'" +
            str(temp) +
            "'/g' "+protein_name+".in")
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
    os.chdir("../..")

# for temp in range(200, 300, 50):
#     # print(temp)
#     seed(datetime.now())
# # simulation set up
#     os.system("mkdir -p simulation/"+str(temp))
#
#     os.chdir("simulation/"+str(temp))
#
#     for q0_percent in range(10, 95, 5):
#         q0 = q0_percent/100.0
#         os.system("mkdir -p "+str(q0))
#         os.system("cp -r ../../"+args.template+"* "+str(q0))
#         os.chdir(str(q0))
#         os.system(  # replace TEMPERATURE with specific steps
#             "sed -i.bak 's/Q0/'" +
#             str(q0) +
#             "'/g' fix_qbias_coeff.data")
#         os.system(  # replace TEMPERATURE with specific steps
#             "sed -i.bak 's/TEMPERATURE/'" +
#             str(temp) +
#             "'/g' "+protein_name+".in")
#         os.system(  # replace SIMULATION_STEPS with specific steps
#             "sed -i.bak 's/WARM_UP_STEPS/'" +
#             str(warm_up_steps) +
#             "'/g' "+protein_name+".in")
#         os.system(  # replace RANDOM with a radnom number
#                 "sed -i.bak 's/RANDOM/'" +
#                 str(randint(1, 10**6)) +
#                 "'/g' "+protein_name+".in")
#         os.system(  # replace SIMULATION_STEPS with specific steps
#                 "sed -i.bak 's/SIMULATION_STEPS/'" +
#                 str(simulation_steps) +
#                 "'/g' "+protein_name+".in")
#         if(platform.system() == 'Darwin'):
#             os.system("/Users/weilu/Documents/lammps-9Oct12_modified/src/lmp_serial \
#             < "+protein_name+".in")
#         elif(platform.system() == 'Linux'):
#             os.system("cp ~/opt/run.slurm .")
#             os.system(  # replace PROTEIN with pdb name
#                     "sed -i.bak 's/PROTEIN/'" +
#                     protein_name +
#                     "'/g' run.slurm")
#             os.system("sbatch run.slurm")
#         else:
#             print("system unkown")
#         os.chdir("..")
#     os.chdir("../..")

sys.exit(0)


# print("hello world")
