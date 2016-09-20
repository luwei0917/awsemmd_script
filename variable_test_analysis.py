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
import glob
# from run_parameter import *
parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        do see the difference variable make \
        run simulation")

# parser.add_argument("template", help="the name of template file")
args = parser.parse_args()

# protein_name = args.template.split('_', 1)[-1].strip('/')
# protein_name = args.template.strip('/')
# simulation_steps = 4 * 10**6
# warm_up_steps = 10 * 10**5


vmd = "/Applications/VMD\ 1.9.2.app/Contents/MacOS/startup.command"

# folder_list = open('folder_list', 'w')
# n = 5
# membrane_k = [1, 2, 3]
# rg_cylindrical_spring_constants = [1, 0.1]
#
# for MemK in membrane_k:
#     add_force_strengths = range(-3-MemK, -1-MemK)
#     pre_pre_folder_name = "MemK"+str(MemK)
#     for ForceStrength in add_force_strengths:
#         pre_folder_name = pre_pre_folder_name+"ForceStrength"+str(ForceStrength)
#         for SpringConstant in rg_cylindrical_spring_constants:
#             for i in range(n):
#                 # simulation set up
#                 folder_name = pre_folder_name+"SpringConstant"+str(SpringConstant)+"_"+str(i)+"/"
#                 os.system("mkdir -p "+folder_name)
#                 folder_list.write(folder_name+"\n")
folder_list = glob.glob("MemK1ForceStrength-3SpringConstant1_*")
#folder_list = [line.rstrip('/\n') for line in open('folder_list')]

# folder_list = ["MemK2ForceStrength-4SpringConstant0.1_0"]
os.system("mkdir -p MyResults")
for folder_name in folder_list:
    os.chdir(folder_name)
    # os.system("vmd -e memmbrane_show.tcl")
    os.system("cp ~/opt/plot_scripts/2xov_movie_screenshot.tcl .")
    os.system(vmd+" -e 2xov_movie_screenshot.tcl")
    print(folder_name)
    os.system("cp frame1000.tga ../MyResults/frame"+folder_name+"_1000.tga")
    #os.system("cp frame450.tga ../Results/frame"+folder_name+"_450.tga")
    # os.system("movie.py "+protein_name)
    os.chdir("..")
