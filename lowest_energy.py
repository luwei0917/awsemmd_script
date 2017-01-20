#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import imp

parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        find the lowest energy frame and cooresponding pdb")


parser.add_argument("protein", help="The name of the protein")
# parser.add_argument("template", help="the name of template file")
# parser.add_argument("-n", "--number", type=int, default=20,
#                     help="Number of simulation run")
# parser.add_argument("-m", "--movie", type=int, default=-1,
#                     help="generate the movie,defalut is all")
# parser.add_argument("-p", "--plotOnly", help="only generate the plot",
#                     action="store_true")
# parser.add_argument("-s", "--steps", type=int, default=4,
#                     help="Simulation steps in unit of million,\
#                     default is 4 million, -1 means test run")
parser.add_argument("-o", "--offAuto", help="turn off from Read from \
                    config file", action="store_true")
args = parser.parse_args()

list_of_lowest_potential_energy = []
protein_name = args.protein.split('.')[0]
# n = args.number
# steps = args.steps*1000*1000
# if args.steps == -1:
#     n = 1  # also set n to be 1 ,this is for debug
#     steps = 10*1000
# # imp.load_source('run_paramter.py', '')
if(not args.offAuto):
    exec (open("config.py").read())
    # print(n, x, y, type(y))
    n = number_of_run
    steps = simulation_steps
    # print(n, steps)
    # sys.exit(0)

for i in range(n):
    os.chdir("simulation/"+str(i))
    sys.stdout = open("lowest_energy.txt", "w")
    record_time = 0
    with open('energy.log') as input_data:
        # Skips text before the beginning of the interesting block:
        next(input_data)
        record_time = 0
        lowest_potential_energy = 0
        last_potential_energy = 0
        next(input_data)
        for line in input_data:
            # time, q = line.strip().split()

            time = int(line.strip().split()[0])
            potenial_energy = float(line.strip().split()[-1])
            # print(time, potenial_energy)
            if(potenial_energy < lowest_potential_energy):
                record_time = time
                lowest_potential_energy = potenial_energy
            last_potential_energy = potenial_energy
        list_of_lowest_potential_energy += [(lowest_potential_energy, record_time, last_potential_energy)]
    time_step = record_time

    print('ITEM: TIMESTEP')
    with open('dump.lammpstrj') as input_data:
        # Skips text before the beginning of the interesting block:
        for line in input_data:
            if line.strip() == str(time_step):
                print(line.strip())  # Or whatever test is needed
                break
        # Reads text until the end of the block:
        for line in input_data:  # This keeps reading the file
            if line.strip() == 'ITEM: TIMESTEP':
                break
            print(line.strip())
    sys.stdout.close()
    os.system(
        "python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py \
        lowest_energy.txt lowest_energy "+protein_name+".seq")
    os.chdir("../..")
sys.stdout = open("list_of_lowest_potential_energy", "w")
# store the globally loweest energy frame
global_min_idx = -1
energy_temp = 0
for idx, q in enumerate(list_of_lowest_potential_energy):
    print(q[0], q[1], q[2], idx)  # max q, timestep of max q, last q
    if(q[0] < energy_temp):
        energy_temp = q[0]
        global_min_idx = idx
sys.stdout.close()
os.system("cp simulation/{}/lowest_energy.pdb global_lowest_energy.pdb".format(global_min_idx))
