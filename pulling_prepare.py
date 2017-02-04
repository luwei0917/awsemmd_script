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
parser.add_argument("--distance", action="store_true", default=False)
parser.add_argument("--replace", action="store_true", default=False)
parser.add_argument("--make_metadata", action="store_true", default=False)
parser.add_argument("-m", "--mode", type=int, default=1)
args = parser.parse_args()


if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

# compute distance by "read dump file"
# if(args.test):
#     for i in range(40):
#         print(i)
#         cd(str(i))
#         do("read_dump_file.py")
#         cd("..")


def replace(TARGET, FROM, TO):
    do("sed -i.bak 's/{}/{}/g' {}".format(FROM,TO,TARGET))


def extract_data():
    do("tail -n+3 energy.log | awk '{print $NF}' > etotal")
    do("head -n 6000 etotal | tail -n 2000 > etotal_half")
    do("head -n 6000 qn | tail -n 2000 > qn_half")
    do("head -n 6000 qc | tail -n 2000 > qc_half")
    do("paste qn qc etotal | tail -n 4000 > data")
    do("paste qn_half qc_half etotal_half > halfdata")


if(args.test):
    for i in range(40):
        print(i)
        cd(str(i))
        do("cp ~/pulling/qo.slurm .")
        do("sbatch qo.slurm")
        # extract_data()
        cd("..")

if(args.make_metadata):
    kconstant = 1000   # double the k constant
    temp = 350
    q0 = 0.12
    metadata = open("metadatafile", "w")
    for i in range(40):
        q = q0 + i*0.02
        # target = "../simulation/350/" + str(i) + "/halfdata {} {} {:.2f}\n".format(temp, kconstant, q)
        target = "../simulation/350/" + str(i) + "/data {} {} {:.2f}\n".format(temp, kconstant, q)
        metadata.write(target)
    metadata.close()

if(args.replace):
    target = "2xov.in"
    replace(target, "TEMPERATURE", "300")
    replace(target, "RANDOM", str(randint(1, 10**6)))

if(args.distance):
    do("read_dump_file.py")
