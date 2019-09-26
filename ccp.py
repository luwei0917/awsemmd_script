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
import glob
from time import sleep
import fileinput
import numpy as np
import pandas as pd
import subprocess


# Useful codes
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# tr " " "\n"
# sed 1d
# sort -u -k 3
# sed -e 's/+T//'
# import re
# numbers = re.compile(r'(\d+)')
# def numericalSort(value):
#     parts = numbers.split(value)
#     parts[1::2] = map(int, parts[1::2])
#     return parts
# mypath = os.environ["PATH"]
# os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
# my_env = os.environ.copy()

parser = argparse.ArgumentParser(description="ccp.py from to")
parser.add_argument("source", help="copy from")
parser.add_argument("to", help="copy to")
parser.add_argument("-m", "--mode", type=int, default=0)
args = parser.parse_args()

# find directory1 -type f -size +500k -exec cp -nv {} directory2/ \;
my_from = args.source
my_to = args.to
if args.mode == 0:
    cmd = f"rsync -a --exclude='*.pkl' --exclude='*.dat' --exclude='*.out' {my_from} {my_to}"
    print(cmd)
    os.system(cmd)
if args.mode == 1:
    cmd = f"rsync -a --exclude='fraglib' {my_from} {my_to}"
    print(cmd)
    os.system(cmd)
