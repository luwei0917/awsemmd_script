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
from myFunctions_helper import *
# compute cross Q for every pdb pair in one folder
parser = argparse.ArgumentParser(description="Compute cross q")
parser.add_argument("protein", default="2xov", help="the name of protein")
parser.add_argument("-m", "--mode",
                    type=int, default=-1)
parser.add_argument("-t", "--targetMode",
                    type=int, default=0)
args = parser.parse_args()
# if args.mode == 1:
#

if args.protein == "2xov":
    if args.targetMode == 0:
        helices_list = [(3,23), (56,77), (80, 101), (109, 126), (135, 150), (159, 178)]
elif args.protein == "tmhc2":
    if args.targetMode == 0:
        helices_list = [(1,38), (43,75), (79, 117), (122, 154), (1, 100), (100, 155)]

if args.mode == -1:
    dumpFile = "dump.lammpstrj"
    outFile = "z_complete.dat"
else:
    dumpFile = f"dump.lammpstrj.{args.mode}"
    outFile = f"z_complete_{args.mode}.dat"

with open(outFile, "w") as f:
    a = read_lammps(dumpFile)
    f.write("z_average, abs_z_average, z_h1, z_h2, z_h3, z_h4, z_h5, z_h6\n")
    for atoms in a:
        b = np.array(atoms)
        z = b.mean(axis=0)[2]
        f.write(str(z)+ ", ")
        z = np.abs(b).mean(axis=0)[2]
        f.write(str(z)+ ", ")           
        for count, (i,j) in enumerate(helices_list):
            z = np.mean(b[i:j], axis=0)[2]
            if count == 5:
                f.write(str(z))
            else:
                f.write(str(z)+ ", ")
        f.write("\n")
