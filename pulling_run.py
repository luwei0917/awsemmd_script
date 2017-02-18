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
import numpy
# Useful codes
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# tr " " "\n"
# sed 1d
# sort -u -k 3
# sed -e 's/+T//'
mypath = os.environ["PATH"]
os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
my_env = os.environ.copy()

parser = argparse.ArgumentParser(description="This is my playground for current project")

# parser.add_argument("protein", help="the name of protein")
# parser.add_argument("template", help="the name of template file")
parser.add_argument("-t", "--test", help="test ", action="store_true", default=False)
parser.add_argument("--pulling", action="store_true", default=False)
parser.add_argument("--pulling2", action="store_true", default=False)
parser.add_argument("--qnqc", action="store_true", default=False)
parser.add_argument("--mutation", action="store_true", default=False)
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("--protein", default="2xov")
parser.add_argument("--dimension", type=int, default=1)
args = parser.parse_args()


if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

if(args.test):
    n = 20
    for i in range(n):
        cd(str(i))
        do("run.py -i 2xov -o")
        cd("..")
# def test():
#     # folder = "T_300_D_130"
#     with open("folder_list") as ins:
#         for line in ins:
#             target = line.strip('\n')
#             cd(target)
#             name_list = ["a206g", "l155a"]
#             for name in name_list:
#                 for i in range(2):
#                     my_from = "simulation/{0}".format(str(i))
#                     my_to = name
#                     cmd = "rsync -a --exclude='restart.*' --exclude='slurm-*' --exclude='movie*' --exclude='q*' --exclude='x.*' {} {}".format(my_from, my_to)
#                     do(cmd)
#                     cd(name+"/{0}".format(str(i)))
#                     do("cp ~/opt/pulling/2xov_{}_rerun.in 2xov_rerun.in".format(name))
#                     do("cp ~/opt/pulling/2xov_{}.seq .".format(name))
#                     do("cp ~/opt/pulling/rerun.slurm .")
#                     cmd = "sbatch rerun.slurm"
#                     do(cmd)
#                     cd("../..")
#             cd("..")
#     # script = "tail -n+2 cv-200-400-10.dat | sort -r -k 2 | head -n1"
#     # result = subprocess.check_output(script, shell=True).decode("utf-8").split()[0]
#     # print(result)
# if(args.test):
#     test()
