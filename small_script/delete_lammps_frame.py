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
import subprocess
import glob
import re

# compute cross Q for every pdb pair in one folder
# parser = argparse.ArgumentParser(description="Compute cross q")
# parser.add_argument("-m", "--mode",
#                     type=int, default=1)

# args = parser.parse_args()

def shrinkage(n=552, shrink_size=10):
    bashCommand = "wc dump.lammpstrj"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    line_number = int(output.decode("utf-8").split()[0])
    print(line_number)
    # print(line_number/552)
    # number of atom = 543
    n = 552
    with open("small.lammpstrj", "w") as out:
        with open("dump.lammpstrj", "r") as f:
            for i, line in enumerate(f):
                if (i // n) % shrink_size == 0:
                    out.write(line)
