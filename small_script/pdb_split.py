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
import numpy as np
import fileinput
# from run_parameter import *
parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    do see the difference variable make \
    run simulation")

parser.add_argument("input", help="the name of template file")
args = parser.parse_args()
inputFile = args.input

a = open("A.pdb", "w")
b = open("B.pdb", "w")
with open(inputFile, "r") as f:
    for line in f:
        try:
            chain = line[21]
            if chain == "A":
                a.write(line)
            elif chain == "B":
                b.write(line)
        except:
            pass
        # individual = line.split()
        # # print(len(individual))
        # try:
        #     if(individual[4] == "A"):
        #         a.write(line)
        #     elif(individual[4] == "B"):
        #         b.write(line)
        # except:
        #     pass

# print("hello")
