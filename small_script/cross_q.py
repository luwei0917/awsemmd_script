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
parser = argparse.ArgumentParser(description="Compute cross q")
parser.add_argument("-m", "--mode",
                    type=int, default=1)
args = parser.parse_args()
if args.mode == 1:
    files = glob.glob("*.pdb")
    n = len(files)
    with open("cross_q", "w") as f:
        for name_i in files:
            print(name_i)
            for name_j in files:
                script = "python2 ~/opt/small_script/CalcQValueFrom2Pdb.py {} {}".format(name_i, name_j)
                result = subprocess.check_output(script, shell=True).decode("utf-8")
                # print(result.strip('\n'))
                # f.write(result.strip('\n'))
                i = re.findall('\d+', name_i)[0]
                j = re.findall('\d+', name_j)[0]
                f.write("{}, {}, {}".format(i, j, result))

if args.mode == 2:
    cross_q = open('cross_q', 'w')
    n = 20
    for i in range(n):
        print(i)
        for j in range(n):
            script = "python2 ~/opt/small_script/CalcQValueFrom2Pdb.py {}.pdb {}.pdb".format(i, j)
            # x = subprocess.check_output()
            result = subprocess.check_output(script, shell=True).decode("utf-8")
            # print(result.strip('\n'))
            cross_q.write(result.strip('\n'))
            cross_q.write(" ")
        cross_q.write("\n")
if args.mode == 3:
    cross_q = open('matrix', 'w')
    n = 20
    for i in range(n):
        print(i)
        for j in range(n):
            script = "python2 ~/opt/small_script/CalcQValueFrom2Pdb.py {}.pdb {}.pdb".format(i, j)
            # x = subprocess.check_output()
            result = subprocess.check_output(script, shell=True).decode("utf-8")
            # print(result.strip('\n'))
            cross_q.write(result)
            # cross_q.write(" ")
        # cross_q.write("\n")
