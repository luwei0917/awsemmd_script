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


parser = argparse.ArgumentParser(description="This is my playground for current project")
# parser.add_argument("protein", help="the name of protein")
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
    do("rm average_f")
    structures = ["ga", "gb"]
    for structure in structures:
        # sequence_list = ["ga77", "gb77"]
        sequence_list = ["gb77", "gb88b", "gb91", "gb95", "gb", "ga", "ga95", "ga91", "ga88", "ga77"]
        for sequence in sequence_list:
            folder = structure+"_" + sequence
            cd(folder)
            do("awk '{ total += $19 } END { print total/NR }' tertiary_frustration_configuration.dat > average_f")
            # do("printf {}".format(folder))
            # do("echo -n sdf")
            do("printf '{} ' >> ../average_f".format(folder))
            do("cat average_f >> ../average_f")
            cd("..")
