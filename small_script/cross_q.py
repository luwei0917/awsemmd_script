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
                    type=int, default=5)
args = parser.parse_args()
if args.mode == 5:
    # files = [32, 41, 47, 48, 57, 58, 6, 63, 69, 72, 82, 96]
    files = [0,1,10,11,12,13,14,15,16,17,18,19,2]
    n = len(files)
    cross_q = open('matrix', 'w')
    with open("cross_q", "w") as f:
        for name_i in files:
            print(name_i)
            for name_j in files:
                script = "python2 ~/opt/small_script/CalcQValueFrom2Pdb.py {} {}".format(name_i, name_j)
                result = subprocess.check_output(script, shell=True).decode("utf-8")
                # print(result.strip('\n'))
                # f.write(result.strip('\n'))
                i = name_i
                j = name_j
                f.write("{}, {}, {}".format(i, j, result))
                cross_q.write(result)

if args.mode == 1:
    files = glob.glob("*.pdb")
    n = len(files)
    cross_q = open('matrix', 'w')
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
                cross_q.write(result)

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
if args.mode == 4:
    files = glob.glob("*.pdb")
    n = len(files)
    with open("matrix", "w") as f:
        for name_i in files:
            print(name_i)
            for name_j in files:
                script = "python2 ~/opt/small_script/CalcQValueFrom2Pdb.py {} {}".format(name_i, name_j)
                result = subprocess.check_output(script, shell=True).decode("utf-8")
                # print(result.strip('\n'))
                # f.write(result.strip('\n'))
                i = re.findall('\d+', name_i)[0]
                j = re.findall('\d+', name_j)[0]
                # f.write("{}, {}, {}".format(i, j, result))
                f.write(result)
