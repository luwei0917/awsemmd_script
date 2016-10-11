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
import matplotlib.pyplot as plt
import numpy as np
# Useful codes
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# tr " " "\n"
# sed 1d
# sort -u -k 3
# sed -e 's/+T//'
mypath = os.environ["PATH"]
os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
my_env = os.environ.copy()

parser = argparse.ArgumentParser(
        description="Do what I want quickly")
parser.add_argument("file", help="the name of file")
parser.add_argument("file2", help="the name of file")
args = parser.parse_args()

vtotal = []
with open(args.file) as f:
    print(f.readline())
    for line in f:
        data = line.strip('\n').split(sep=" ")
        data = [float(i) for i in data]
        vtotal += [data[-3]]
    # print(vtotal)
v2 = []
with open(args.file2) as f:
    print(f.readline())
    f.readline()
    for line in f:
        data = line.strip('\n').split(sep="\t")
        data = [float(i) for i in data]
        v2 += [data[-1]]
    # print(v2)
# plt.plot(vtotal)
# plt.plot(v2)
plt.plot(np.array(vtotal)-np.array(v2))
plt.show()
# sum =0
# for i in data[1:-3]:
#     sum += i
# print(sum+(144.18+9.369))
