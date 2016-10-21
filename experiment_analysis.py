#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
from experiment_functions import *

mypath = os.environ["PATH"]
os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
my_env = os.environ.copy()

parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        automatically analysis the simulation")

args = parser.parse_args()

exec(open("config.py").read())
n = number_of_run
steps = simulation_steps
# protein_name = protein_name

clean()
os.chdir("analysis")

for i in range(n):
    os.chdir(str(i))
    gagb()
    os.chdir("..")

os.chdir("../results")
os.system("paste ga.dat gb.dat > gagb.dat")
os.system("cp ~/opt/gagb/q_ga-gb.gp .")
os.system("gnuplot q_ga-gb.gp")
os.system("open q_ga-gb.pdf")
os.chdir("..")
# os.system("mkdir -p results")
# for i in range(n):
#     # analysis
#     os.system("mkdir -p analysis/"+str(i))
#     os.chdir("analysis/"+str(i))
#
#     sys.stdout = open("final.txt", "w")
#     print('ITEM: TIMESTEP')
#     time_step = simulation_steps
#     with open('dump.lammpstrj') as input_data:
#         # Skips text before the beginning of the interesting block:
#         for line in input_data:
#             if line.strip() == str(time_step):
#                 print(line.strip())  # Or whatever test is needed
#                 break
#         # Reads text until the end of the block:
#         for line in input_data:  # This keeps reading the file
#             if line.strip() == 'ITEM: TIMESTEP':
#                 break
#             print(line.strip())
#     sys.stdout.close()
#
#     os.chdir("../..")


# os.system("> results/cross_q")
# for i in range(n):
#     for j in range(0, i):
#         os.system("echo '"+str(i)+" "+str(j)+" \c' >> results/cross_q")
#         os.system("python2 ~/opt/CalcQValue.py analysis/"+str(i)+"/final \
#         analysis/"+str(j)+"/final.txt >> results/cross_q ")
