#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess

mypath = os.environ["PATH"]
os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
my_env = os.environ.copy()

parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        automatically analysis the simulation")

parser.add_argument("template", help="the name of template file")
parser.add_argument("-n", "--number", type=int, default=20,
                    help="Number of simulation run")
parser.add_argument("-s", "--steps", type=int, default=4,
                    help="Simulation steps in unit of million,\
                    default is 4 million, -1 means test run")
args = parser.parse_args()

list_of_max_q = []

n = args.number
protein_name = args.template.strip('/')

os.system("mkdir -p results")
for i in range(n):
    # analysis
    os.system("mkdir -p analysis/"+str(i))
    os.chdir("analysis/"+str(i))

    sys.stdout = open("final.txt", "w")
    print('ITEM: TIMESTEP')
    time_step = args.steps*1000*1000
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

    os.chdir("../..")

os.system("> results/cross_q")
for i in range(n):
    for j in range(0, i):
        os.system("echo '"+str(i)+" "+str(j)+" \c' >> results/cross_q")
        os.system("python2 ~/opt/CalcQValue.py analysis/"+str(i)+"/final \
        analysis/"+str(j)+"/final.txt >> results/cross_q ")
