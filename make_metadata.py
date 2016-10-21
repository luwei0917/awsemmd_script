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

temp_list = range(200, 500, 100)
q0_list = range(10, 100, 10)

# config = open('config.py', 'w')
# config.write("number_of_run = %d\nsimulation_steps = %d\n\
# warm_up_steps = %d\n" % (n, simulation_steps, warm_up_steps))
# config.close()
kconstant = 200
metadata = open("metadatafile", "w")
for temp in temp_list:
    for q0_percent in q0_list:
        q0 = q0_percent/100.0
        if(q0_percent % 2 == 0):
            metadata.write("simulation/%.0f/%.1f/wham.dat %.0f %d %.2f\n" % (temp, q0, temp, kconstant, q0))
        else:
            metadata.write("simulation/%.0f/%.2f/wham.dat %.0f %d %.2f\n" % (temp, q0, temp, kconstant, q0))

metadata.close()
sys.exit(0)


# print("hello world")
