#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import imp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import datetime
import pickle
matplotlib.style.use('fivethirtyeight')


# print(plt.style.available)
# mypath = os.environ["PATH"]
# os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
# my_env = os.environ.copy()

parser = argparse.ArgumentParser(
    description="Plot my graphs quickly")

parser.add_argument("-s", "--save", action="store_true", default=False)
parser.add_argument("-r", "--reproduce", default=None)
parser.add_argument("-t", "--test", action="store_true", default=False)
parser.add_argument("-m", "--mode", default="pulling")
parser.add_argument("-d", "--debug", action="store_true", default=False)

args = parser.parse_args()


if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir


if args.reproduce is not None:
    print("Reproducing!")
    with open(args.reproduce, "rb") as f:
        args = pickle.load(f)
        print(args)
        args.save = False

if(args.test):
    print("Hello Test World")
    force_list = np.arange(1,2.5,0.1)
    array = []
    cwd = os.getcwd()
    print(cwd)
    for force in force_list:
        folder = "2d_2_force_" + str(force)
        cd(folder)
        do("pulling_plotcontour.py")
        do("cp test.png ../results/{}.png".format(force))
        cd(cwd)

if(args.save):
    # print(os.getcwd())
    # print(args)
    print("Saving")
    # print(datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
    with open("args"+datetime.datetime.now().strftime("%Y%m%d-%H%M"), "wb") as f:
        pickle.dump(args, f)
    os.system("cp ~/opt/plot.py plot_{}.py".format(datetime.datetime.now().strftime("%Y%m%d-%H%M")))

if(args.mode == 1):
    folder = "wham"
    do("mkdir " + folder)
    cd(folder)
    do("make_metadata.py -k 200 -a 2 -t")
    do("pulling_analysis.py -f -m 8")
    do("sbatch freeEnergy.slurm")
