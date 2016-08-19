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
        do see the difference variable make \
        run simulation")

parser.add_argument("template", help="the name of template file")
args = parser.parse_args()

vmd = "/Applications/VMD\ 1.9.2.app/Contents/MacOS/startup.command"
# protein_name = args.template.split('_', 1)[-1].strip('/')
protein_name = args.template.strip('/')
simulation_steps = 4 * 10**6
warm_up_steps = 10 * 10**5

seed(datetime.now())

rg_cylindrical_spring_constants = [1, 0.1, 0.01, 0.001]


for SpringConstant in rg_cylindrical_spring_constants:
    # simulation set up
    folder_name = "SpringConstant"+str(SpringConstant)+"/"
    print("---"+folder_name+"---\n")
    os.chdir(folder_name)
    os.system("movie.py "+protein_name)
    os.system(vmd + " -e membrane_show.tcl")
    os.chdir("..")

# print("hello world")
