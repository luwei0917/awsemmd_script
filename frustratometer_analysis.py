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
from myPersonalFunctions import *
import glob
# Useful codes
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# tr " " "\n"
# sed 1d
# sort -u -k 3
# sed -e 's/+T//'
# awk '$5=-$5' data
mypath = os.environ["PATH"]
os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
my_env = os.environ.copy()

parser = argparse.ArgumentParser(description="This is my playground for current project")
# parser.add_argument("protein", help="the name of protein")
# parser.add_argument("template", help="the name of template file")
parser.add_argument("-t", "--test", help="test ", action="store_true", default=False)
parser.add_argument("-n", "--number", type=int, default=10, help="number of run")
parser.add_argument("-d", "--debug", action="store_true", default=False)
args = parser.parse_args()


if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

if(args.test):
    structures = ["ga", "gb"]
    for structure in structures:
        folder_list = ["ga77", "gb77"]
        for folder in folder_list:
            do("cp -r gagb "+structure + "_" + folder)
            cd(structure+"_" + folder)
            do("cp {}.seq gagb.seq".format(folder))
            do("cp data.{} data.gagb".format(structure))
            mode_list = ["single", "configuration", "mutational"]
            for mode in mode_list:
                do("~/bin/lmp_serial_old < {}.in".format(mode))
                do("mv tertiary_frustration.dat tertiary_frustration_{}.dat".format(mode))
                do("mv tertiary_frustration.tcl tertiary_frustration_{}.tcl".format(mode))
            cd("..")
