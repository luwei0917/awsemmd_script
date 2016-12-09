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
parser.add_argument("--gagb", help="Project gagb", action="store_true", default=False)
parser.add_argument("-f", "--freeEnergy", help="free energy calculation", action="store_true", default=False)
parser.add_argument("--test", help="test ", action="store_true", default=False)

args = parser.parse_args()


def test():
    print("don't show me")
    for i in range(20):
        os.chdir(str(i))
        # os.system("gg.py --qnqc")
        os.system("paste qn qc > qnqc")
        os.chdir("..")
if(args.test):
    test()

gagb_free_energy_config = "\n\
ndim = {} \n\
biasing_variable_column = 3\n\
energy_column = 6                  \n\
pmf_variable_column_1 = {}          \n\
pmf_variable_column_2 = {}          \n\
nbins1 = 20                        \n\
nbins2 = 20                        \n\
start_temperature = 200            \n\
end_temperature = 400              \n\
temperature_increment = 10         \n\
N_samples = 4000                   \n\
"
ndim = 1
pmf_variable_column_1 = 3
pmf_variable_column_2 = 2


def gagb_freeEnergy_calculation():
    print("gagb_freeEnergy_calculation")
    # wham one d on ga, gb
    # name_list = ["ga", "gb"]
    name_list = ["test"]
    for name in name_list:
        folder = "one_d_"+name
        os.system("mkdir -p "+folder)
        os.chdir(folder)
        os.system("make_metadata.py --gagb")
        config = open('config.py', 'w')
        this_config = gagb_free_energy_config.format(ndim, pmf_variable_column_2, pmf_variable_column_1)
        config.write(this_config)
        # config.write("protein_name = '%s'\nnumber_of_run = %d\nsimulation_steps = %d\n\
        # warm_up_steps = %d\n" % (protein_name, n, simulation_steps, warm_up_steps))
        config.close()
        os.system("python2 ~/opt/compute-pmf.py")
        os.system("myplot.py --gagb")
        os.chdir("..")
    # wham two d

if(args.gagb and args.freeEnergy):
    gagb_freeEnergy_calculation()
