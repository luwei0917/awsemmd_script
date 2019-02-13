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
from time import sleep

# from run_parameter import *
parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    automatic copy the template file, \
    run simulation and analysis")
parser.add_argument("protein", help="The name of the protein")
parser.add_argument("-r", "--run", help="test mode",
                    action="store_true")
parser.add_argument("-s", "--see", help="test mode",
                    action="store_true")
# parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("-m", "--mode", type=int, default=0)
parser.add_argument("-d", "--day", type=str, default="someday")
parser.add_argument("-t", "--test", action="store_true", default=False)
args = parser.parse_args()

do = os.system
cd = os.chdir

def replace(TARGET, FROM, TO):
    do("sed -i.bak 's/{}/{}/g' {}".format(FROM,TO,TARGET))
def extra(fileName, offset=0):
    replace(fileName, "1 1 30 5", "1 1 30 {}".format(offset))
def rerun(proteinName, extra=extra, offset=0):
    do(f"cp {proteinName}_0.in rerun_{offset}.in")
    fileName = "rerun_{}.in".format(offset)
    replace(fileName, "fix               1 all nve", "")
    replace(fileName, "fix               2 all langevin", "#")
    replace(fileName, "run", "#")
    # replace(fileName, "0\/", "recompute_offset_{}\/".format(offset))
    replace(fileName, "0\/", "rerun\/")
    do("mkdir rerun")
    do("cp 0/dump.lammpstrj rerun/")
    replace(fileName, "minimize", "# minimize")
    replace(fileName, "dump		1 all atom 4000", "#")
    replace(fileName, "dump_modify	1 sort id", "")
    replace(fileName, "restart         100000 restart", "rerun rerun\/dump.lammpstrj dump x y z")

    # replace(fileName, "1 1 30 5", "1 1 30 0")
    # extra(fileName, offset)
    slurm = "rerun_{}.slurm".format(offset)
    do("cp run_0.slurm " + slurm)
    replace(slurm, "#SBATCH --time=04:00:00", "#SBATCH --time=00:10:00")
    replace(slurm, f"{proteinName}_0.in", "rerun_{}.in".format(offset))
    if args.test:
        do("/home/wl45/build/sep03/src/lmp_serial -in rerun_0.in")
    else:
        do("sbatch " + slurm)
    # do("cp ~/opt/2xov_eval/run.slurm " + slurm)
    # replace(slurm, "2xov_eval", "rerun_{}".format(offset))
    # replace(slurm, "commons", "ctbp-common")
    # do("mkdir recompute_offset_{}".format(offset))
    # do("sbatch " + slurm)

proteinName = args.protein.strip("/.").lower()
rerun(proteinName)
# parser.add_argument("template", help="the name of template file")
# parser.add_argument("-n", "--number", type=int, default=20,
#                     help="Number of simulation run")
# parser.add_argument("-s", "--steps", type=int, default=8,
#                     help="Simulation steps in unit of million,\
#                     default is 8 million, -1 means test run")
# parser.add_argument("-r", "--read", help="Read from config file",
#                     action="store_true")
# parser.add_argument("-ws", "--warmSteps", type=int, default=20,
#                     help="Warmup Simulation steps in unit of hundred thousand,\
#                     default is 2 million")
# parser.add_argument("-t", "--test", help="test mode",
#                     action="store_true")
# parser.add_argument("-c", "--copy",
#                     help="copy the restart file before run",
#                     action="store_true")
# parser.add_argument("-o", "--offAuto", help="turn off from Read from \
#                     config file", action="store_true", default=False)
# args = parser.parse_args()
# # TODO:
# # add clean command.
# # test analysis, and fullfill missing anaylsis.
#
# # protein_name = args.template.split('_', 1)[-1].strip('/')
# n = args.number
# protein_name = args.template.strip('/')
# if args.steps == -1:  # smallest run for debug.
#     simulation_steps = 10**5
#     warm_up_steps = 10**4
#     n = 1  # also set
# elif args.test:  # test run
#     simulation_steps = 50 * 10**3
#     warm_up_steps = 50 * 10**3
# else:
#     simulation_steps = args.steps * 10**6
#     warm_up_steps = args.warmSteps * 10**5
#
# config = open('config.py', 'w')
# config.write("protein_name = '%s'\nnumber_of_run = %d\nsimulation_steps = %d\n\
# warm_up_steps = %d\n" % (protein_name, n, simulation_steps, warm_up_steps))
# config.close()
# if(not args.offAuto):
#     exec(open("variables.dat").read())
#     print(TSTART, TEND)
# for i in range(n):
#     seed(datetime.now())
# # simulation set up
#     if(args.copy):
#         os.system("cp restart/{}/melt.4000000 1qjp/".format(i))
#     name = "rerun_" + str(i)
#     os.system("mkdir -p simulation/"+name)
#     os.system("cp -r "+args.template+"* simulation/"+name)
#     os.system("cp simulation/"+str(i)+"/dump.lammpstrj simulation/"+name)
#     os.chdir("simulation/"+name)
#     os.system("cp ~/opt/pulling/2xov_rerun.in 2xov.in")
#     os.system(  # replace SIMULATION_STEPS with specific steps
#         "sed -i.bak 's/WARM_UP_STEPS/'" +
#         str(warm_up_steps) +
#         "'/g' "+protein_name+".in")
#     os.system(  # replace RANDOM with a radnom number
#         "sed -i.bak 's/RANDOM/'" +
#         str(randint(1, 10**6)) +
#         "'/g' "+protein_name+".in")
#     os.system(  # replace SIMULATION_STEPS with specific steps
#         "sed -i.bak 's/SIMULATION_STEPS/'" +
#         str(simulation_steps) +
#         "'/g' "+protein_name+".in")
#     if args.steps == -1:
#         os.system(  # replace TEMPERATURE with specific steps
#             "sed -i.bak 's/Q0/'" +
#             str(0.5) +
#             "'/g' fix_qbias_coeff.data")
#         os.system(  # replace TEMPERATURE with specific steps
#             "sed -i.bak 's/TEMPERATURE/'" +
#             str(350) +
#             "'/g' "+protein_name+".in")
#     if(not args.offAuto):
#             os.system(  # replace SIMULATION_STEPS with specific steps
#                 "sed -i.bak 's/TSTART/'" +
#                 str(TSTART) +
#                 "'/g' "+protein_name+".in")
#             os.system(  # replace SIMULATION_STEPS with specific steps
#                 "sed -i.bak 's/TEND/'" +
#                 str(TEND) +
#                 "'/g' "+protein_name+".in")
# # if(platform.system() == 'Darwin'):
# #     os.system("/Users/weilu/Documents/lammps-9Oct12_modified/src/lmp_serial \
# #     < "+protein_name+".in")
#     if(platform.system() == 'Darwin'):
#         os.system("/Users/weilu/Research/bin/lmp_serial < "+protein_name+".in")
#         # os.system("/Users/weilu/Research/Build/lammps-9Oct12_modified/src/lmp_serial \
#         # < "+protein_name+".in")
#     elif(platform.system() == 'Linux'):
#         os.system("cp ~/opt/run.slurm run.slurm")
#         os.system(  # replace PROTEIN with pdb name
#             "sed -i.bak 's/PROTEIN/'" +
#             protein_name +
#             "'/g' run.slurm")
#         os.system("sbatch run.slurm")
#         sleep(0.2)  # Time in seconds.
#     else:
#         print("system unkown")
#     os.chdir("../..")

# print("hello world")
