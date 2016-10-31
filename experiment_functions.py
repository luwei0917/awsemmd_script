#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess


def clean():
    os.system("rm results/ga.dat")
    os.system("rm results/gb.dat")
    os.system("rm results/gaa.dat")
    os.system("rm results/ga_highest.dat")
    os.system("rm results/gb_highest.dat")

def gagb():
    os.system("cp ~/opt/gagb/*.pdb .")
    os.system("awk '{if ((!((($1>0 && $1<25) || ($1>159 && $1<200) ) && $3>-10)  ) ) print }' dump.lammpstrj > data_test")
    os.system("python2 ~/opt/script/CalcQValue.py 2lhc_part.pdb data_test test")
    os.system("tail -n 50 test | awk '{s+=$1} END {print s/50}' > gaa")
    os.system("tail -n 50 test | sort | tail -n 1 > ga_highest")

    os.system("python2 ~/opt/script/CalcQValue.py 2lhc.pdb dump.lammpstrj test2")
    os.system("tail -n 50 test2 | sort | tail -n 1 > gb_highest")

    os.system("tail -n +40 final.pdb | ghead -n -15 > final_chosen.pdb")

    os.system("python2 ~/opt/small_script/CalcQValueFrom2Pdb.py 2lhc_part.pdb final_chosen.pdb > ga")
    os.system("python2 ~/opt/small_script/CalcQValueFrom2Pdb.py 2lhd.pdb final.pdb > gb")

    os.system("cat ga_highest >> ../../results/ga_highest.dat")
    os.system("cat gb_highest >> ../../results/gb_highest.dat")
    os.system("cat gaa >> ../../results/gaa.dat")
    os.system("cat ga >> ../../results/ga.dat")
    os.system("cat gb >> ../../results/gb.dat")
# os.system("> results/cross_q")
# for i in range(n):
#     for j in range(0, i):
#         os.system("echo '"+str(i)+" "+str(j)+" \c' >> results/cross_q")
#         os.system("python2 ~/opt/CalcQValue.py analysis/"+str(i)+"/final \
#         analysis/"+str(j)+"/final.txt >> results/cross_q ")
# exec(open("config.py").read())
# n = number_of_run
# steps = simulation_steps

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
