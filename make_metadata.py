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
# from run_parameter import *
parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    make metadata")

parser.add_argument("-k", "--kconstant", type=int, default=600)
parser.add_argument("-q", "--qStart", type=float, default=0.0)
parser.add_argument("-n", "--number", type=int, default=20)

parser.add_argument("-m", "--mode", type=int, default=1)
parser.add_argument("-a", "--additionalMode", type=int, default=1)

parser.add_argument("-t", "--test", help="test ", action="store_true", default=False)
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("-s", "--save", action="store_true", default=False)
parser.add_argument("-r", "--reproduce", default=None)
args = parser.parse_args()

if(args.reproduce):
    print("Reproducing!")
    with open(args.reproduce, "rb") as f:
        args = pickle.load(f)
        print(args)


if(args.save):
    print(os.getcwd())
    print("Saving")
    args.save = False
    with open("args"+datetime.datetime.now().strftime("%m%d-%H%M"), "wb") as f:
        pickle.dump(args, f)

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir



if(args.test):
    kconstant = args.kconstant * 2   # double the k constant
    q0 = args.qStart
    metadata = open("metadatafile", "w")
    if(args.mode == 1):
        qStep = 0.05
        temp_list = [135, 160, 185, 210]
        # temp_list = [160]
        # temp_list = ['300', '200', '250']
        if(args.additionalMode == 1):
            pre_fix = "first_2000_"
        else:
            pre_fix = ""
    for i in range(args.number):
        q = q0 + i * qStep
        name = pre_fix + str(i)
        for temp in temp_list:
            target = "../data/{}/".format(temp) + name + " {} {} {:.2f}\n".format(temp, kconstant, q)
            # target = "../data/t_{}/small_".format(temp) + str(i) + " {} {} {:.2f}\n".format(temp, kconstant, q)
            # target = "../data/{}/".format(temp) + name + " {} {} {:.2f}\n".format(temp, kconstant, q)
            metadata.write(target)
    metadata.close()
