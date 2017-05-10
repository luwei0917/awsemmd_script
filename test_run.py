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
from time import sleep

# from run_parameter import *
parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    automatic copy the template file, \
    run simulation and analysis")

parser.add_argument("input", help="the name of template file")
parser.add_argument("-n", "--number", type=int, default=20,
                    help="Number of simulation run")
parser.add_argument("-t", "--test", help="test mode",
                    action="store_true")
parser.add_argument("-i", "--inplace", help="change in this folder",
                    action="store_true", default=False)
parser.add_argument("-m", "--mode",
                    help="mode 2 is dependence run",
                    type=int, default=1)
parser.add_argument("-d", "--debug", action="store_true", default=False)
args = parser.parse_args()
# TODO:
# add clean command.
# test analysis, and fullfill missing anaylsis.

# protein_name = args.template.split('_', 1)[-1].strip('/')
n = args.number


if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir


seed(datetime.now())
do(  # replace RANDOM with a radnom number
    "sed -i.bak 's/RANDOM/'" +
    str(randint(1, 10**6)) +
    "'/g' "+args.input)
do("~/bin/lmp_serial < "+args.input)


# print("hello world")
