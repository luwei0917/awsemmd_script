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
parser.add_argument("-p", "--pick", default=None)
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


if(args.pick):
    # print("Hello")
    metadata = open("metadatafile", "w")
    # folder_list = [4, 23, 27, 29, 30, 35, 36, 37, 38]
    # folder_list = [41, 42, 45, 91, 95, 7, 33]
    # force 0.3, dis 55-65
    # folder_list = [0, 2, 4, 6, 7, 9, 12, 13, 16, 17, 18, 25, 27, 30, 32, 33, 39]
    # force 0.3, dis 45-55
    # folder_list = [13, 17, 23, 24, 6, 8]
    # dis 170-185, force 0.3
    folder_list = args.pick.split(",")
    folder_list = [int(i) for i in folder_list]
    print(folder_list)
    cwd = os.getcwd()
    print(cwd)
    if(args.mode == 1):
        for i in folder_list:
            for j in range(0, 7):
                target = cwd + "/../wt_2/simulation/{0}/{1}/data\n".format(i, j)
                metadata.write(target)

    if(args.mode == 2):
        for i in folder_list:
            for j in [0, 1, 2]:
                if(i > 39):
                    ii = i - 40
                    target = cwd + "/../wt_2/simulation/{0}/{1}/data\n".format(ii, j)
                else:
                    target = cwd + "/../wt/simulation/{0}/{1}/data\n".format(i, j)
                metadata.write(target)
