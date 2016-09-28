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
parser.add_argument("-m", "--movie", type=int, default=-2,
                    help="generate the movie,defalut is none")
args = parser.parse_args()

list_of_max_q = []

exec(open("config.py").read())
# print(n, x, y, type(y))
n = number_of_run
steps = simulation_steps
# print(n, steps)
# sys.exit(0)

# protein_name = args.template.strip('/')

os.system("mkdir -p analysis")
for i in range(n):
    # analysis
    os.system("mkdir -p analysis/"+str(i))
    os.chdir("analysis/"+str(i))
    # move necessary file into analysis folder
    sys.stdout = open("chosen.txt", "w")
    os.system("mv ../../simulation/"+str(i)+"/dump.lammpstrj .")
    os.system("mv ../../simulation/"+str(i)+"/wham.dat .")
    os.system("mv ../../simulation/"+str(i)+"/energy.dat .")
    record_time = 0
    with open('wham.dat') as input_data:
        # Skips text before the beginning of the interesting block:
        record_time = 0
        max_q = 0
        last_q = 0
        next(input_data)
        for line in input_data:
            # time, q = line.strip().split()

            time = int(line.strip().split()[0])
            q = float(line.strip().split()[1])
            if(q > max_q):
                record_time = time
                max_q = q
            last_q = q
        list_of_max_q += [(max_q, record_time, last_q)]
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
    if(args.movie == -1 or args.movie == i):
        os.system(
            "python2 ~/opt/script/BuildAllAtomsFromLammps.py \
            dump.lammpstrj movie")
    os.chdir("../..")
sys.stdout = open("analysis/list_of_max_q", "w")
for idx, q in enumerate(list_of_max_q):
    print(q[0], q[1], q[2], idx)  # max q, timestep of max q, last q
sys.stdout.close()
