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
# from run_parameter import *
parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        make metadata")
args = parser.parse_args()

temp_list = range(350, 400, 100)
# q0_list = range(10, 100, 10)

# config = open('config.py', 'w')
# config.write("number_of_run = %d\nsimulation_steps = %d\n\
# warm_up_steps = %d\n" % (n, simulation_steps, warm_up_steps))
# config.close()
kconstant = 200
n = 50
q_bias_step = 0.02

metadata = open("metadatafile", "w")
for temp in temp_list:
    q0 = 0
    for i in range(n):
        q0 += q_bias_step
        metadata.write("../simulation/%.0f/%d/halfdata.dat %.0f %d %.2f\n" % (temp, i, temp, kconstant, q0))

metadata.close()
sys.exit(0)


# print("hello world")
