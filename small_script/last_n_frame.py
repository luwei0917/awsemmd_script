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

parser = argparse.ArgumentParser(description="This is my aawsem project headquarter")
parser.add_argument("-t", "--test", help="Test run", action="store_true", default=False)
parser.add_argument("-m", "--mode",
                    type=int, default=1)
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("protein", help="the name of protein file")
parser.add_argument("--total", type=int, default=-1)
parser.add_argument("-n", "--number",
                    type=int, default=20)
args = parser.parse_args()

try:
    protein_name,_ = args.protein.split('.')
except:
    try:
        protein_name,_ = args.protein.split('/')
    except:
        protein_name = args.protein
        print("ATTENSION! protein name {}\n Correct?".format(protein_name))

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

def file_width(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])
if args.mode == 1:
    # print(file_width("wham.dat"))
    if args.total == -1:
        total = file_width("wham.dat")
    else:
        total = args.total

    n = args.number
    do("mkdir -p frames")
    for i in range(n):
        frame = total - n + i
        do("BuildAllAtomsFromLammps_seq.py dump.lammpstrj frames/{0} ../{1}.seq {2}".format(frame, protein_name, frame))


if args.mode == 2:
    # print(file_width("wham.dat"))
    n = args.number
    do("mkdir -p frames")
    for i in range(n):
        # frame = total - n + i
        frame = 1000 + i*35
        do("BuildAllAtomsFromLammps_seq.py dump.lammpstrj frames/{0} {1}.seq {2}".format(i, protein_name, frame))
