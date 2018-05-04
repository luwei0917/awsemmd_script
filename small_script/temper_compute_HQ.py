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
import numpy as np
# compute cross Q for every pdb pair in one folder
parser = argparse.ArgumentParser(description="Compute cross q")
parser.add_argument("-m", "--mode",
                    type=int, default=-1)
parser.add_argument("-t", "--targetMode",
                    type=int, default=0)
parser.add_argument("protein", default="2xov", help="the name of protein")
args = parser.parse_args()
# if args.mode == 1:
#

do = os.system
cd = os.chdir

name = args.protein
for i in range(12):
    cmd = f"python3 ~/opt/small_script/find_distance.py {name} -m {i} -t {args.targetMode}"
    do(cmd)