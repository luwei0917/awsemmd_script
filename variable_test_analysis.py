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

# protein_name = args.template.split('_', 1)[-1].strip('/')
protein_name = args.template.strip('/')
simulation_steps = 4 * 10**6
warm_up_steps = 10 * 10**5

seed(datetime.now())
vmd = "/Applications/VMD\ 1.9.2.app/Contents/MacOS/startup.command"

folder_list = [line.rstrip('\n') for line in open('folder_list')]
for folder_name in folder_list:
    os.chdir(folder_name)
    print(folder_name)
    # os.system("pwd")
    os.system("movie.py "+protein_name)
    os.system(vmd+" -e membrane_show.tcl")
    os.chdir("..")
