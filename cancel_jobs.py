#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import imp

parser = argparse.ArgumentParser(description="This is a python3 script to scancel my jobs")

parser.add_argument("start", type=int, help="starting from")
parser.add_argument("-n", type=int, default=1, help="next n jobs")
parser.add_argument("-a", "--all", action="store_true", default=False)
args = parser.parse_args()

def getFromTerminal(CMD):
    return subprocess.Popen(CMD,stdout=subprocess.PIPE,shell=True).communicate()[0].decode()

def scancel_jobs_in_folder():
    cmd = "find . -name 'slurm-*' | rev | awk -F'[-.]' '{print $2}' | rev"
    lines = getFromTerminal(cmd).splitlines()
    for line in lines:
        print(line)
        os.system("scancel " + line)

if args.all:
    print("I will cancel every one under this folder.")
    scancel_jobs_in_folder()
else:
    n = args.n
    start = args.start
    for i in range(n):
        os.system("scancel "+str(start+i))
