#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import imp

parser = argparse.ArgumentParser(
        description="This is a python3 script to scancel my jobs")

parser.add_argument("start", type=int, help="starting from")
parser.add_argument("-n", type=int, default=1, help="next n jobs")
args = parser.parse_args()

n = args.n
start = args.start
for i in range(n):
    os.system("scancel "+str(start+i))
