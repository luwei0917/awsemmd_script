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
    structures = ["ga", "gb"]
    for structure in structures:
        sequence_list = ["ga77", "gb77"]
        # sequence_list = ["gb77", "gb88b", "gb91", "gb95", "gb", "ga", "ga95", "ga91", "ga88", "ga77"]
        for sequence in sequence_list:
            folder = structure+"_" + sequence
            do("cp -r gagb "+folder)
            cd(folder)
            do("cp {}.seq gagb.seq".format(sequence))
            do("cp data.{} data.gagb".format(structure))
            mode_list = ["single", "configuration", "mutational"]
            for mode in mode_list:
                do("~/bin/lmp_serial_old < {}.in".format(mode))
                do("mv tertiary_frustration.dat tertiary_frustration_{}.dat".format(mode))
                do("mv tertiary_frustration.tcl tertiary_frustration_{}.tcl".format(mode))
                do("awk '{ total += $19 } END { print total/NR }' tertiary_frustration_configuration.dat > average_f")
                do("echo '{}' >> ../average_f".format(sequence))
                do("cat average_f >> ../average_f")
            cd("..")
