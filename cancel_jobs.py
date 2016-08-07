#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import imp

parser = argparse.ArgumentParser(
        description="This is a python3 script to scancel my jobs")

parser.add_argument("from", type=int, help="starting from")
parser.add_argument("N", type=int, help="next n jobs")
args = parser.parse_args()

n = args.N
start = args.from
for i in range(n):
    print(start+i)# os.system("scancel "+str(start+i))
