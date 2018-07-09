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

parser = argparse.ArgumentParser(
    description="Prepare the data for run and analysis. \
                Codes here only need run once")

# parser.add_argument("protein", help="the name of protein")
# parser.add_argument("template", help="the name of template file")
parser.add_argument("-t", "--test", help="test ", action="store_true", default=False)
parser.add_argument("-d", "--debug", action="store_true", default=False)
args = parser.parse_args()


if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

if(args.test):
    for i in range(20):
        cd(str(i))
        # do("tail -n+3 energy.log | awk '{print $NF - $13}' > etotal")
        do("tail -n+3 energy.log | awk '{print $13}' > etotal")
        do("paste etotal qw > qw_etotal")
        cd("..")
    with open("data", "w") as out:
        out.write("step, qw, run, energy\n")
        for i in range(20):
            print(i)
            with open(str(i)+"/qw_etotal") as f:
                step = 0
                for line in f:
                    step += 1
                    energy, qw = line.split()
                    out.write("{}, {}, run_{}, {}\n".format(step, qw, i, energy))
                    # out.write(str(n)+", "+qw+", run_"+str(i)+", "+energy+"\n")
