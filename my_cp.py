#!/usr/bin/env python3
# @Author: Wei Lu <weilu>
# @Date:   04-Dec-2017
# @Email:  wl45@rice.edu
# @Last modified by:   weilu
# @Last modified time: 04-Dec-2017
# @Copyright: Free

import os
import sys
import argparse
import platform
from datetime import datetime
from myPersonalFunctions import *
import glob
import numpy as np
import datetime
import pandas as pd
import re
# Useful codes
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# tr " " "\n"
# sed 1d
# sort -u -k 3
# sed -e 's/+T//'
parser = argparse.ArgumentParser(description="Analysis code, need run multiple times")

parser.add_argument("source", help="copy from inlcude /")
parser.add_argument("to", help="target place")
parser.add_argument("-m", "--mode", type=int, default=0)

args = parser.parse_args()

# find directory1 -type f -size +500k -exec cp -nv {} directory2/ \;
my_from = args.source
my_to = args.to
if args.mode == 0:
    cmd = f"rsync -a --exclude='*.pkl' --exclude='*.dat' --exclude='*.out' {my_from} {my_to}"
    print(cmd)
    os.system(cmd)
