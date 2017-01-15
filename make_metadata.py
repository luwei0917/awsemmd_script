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

parser.add_argument("--qnqc", help="for all calculate q of n terminal and q of c terminal ", action="store_true", default=False)
parser.add_argument("--qnqc2", help="for all calculate q of n terminal and q of c terminal ", action="store_true", default=False)
parser.add_argument("--pulling", help="for all calculate q of n terminal and q of c terminal ", action="store_true", default=False)
parser.add_argument("--gagb", help="for all calculate q of n terminal and q of c terminal ", action="store_true", default=False)
parser.add_argument("-m", "--mode", type=int, default=1)
parser.add_argument("--server", action="store_true", default=False)
args = parser.parse_args()


def gagb():
    print("GaGb")
    # os.system("cp folder_list .")
    kconstant = 112   # double the k constant
    temp = 350
    q0 = 0.12
    metadata = open("metadatafile", "w")
    for i in range(40):
        q = q0 + i*0.02
        target = "../simulation/350/" + str(i) + "/halfdata.dat {} {} {:.2f}\n".format(temp, kconstant, q)
        metadata.write(target)
    metadata.close()

if(args.gagb):
    gagb()


def qnqc():
    print("QnQc")
    # os.system("cp folder_list .")
    kconstant = 0.01

    metadata = open("metadatafile", "w")
    with open('folder_list', 'r') as ins:
        for line in ins:
            target = line.strip('\n')
            temp = target.split("_")[1]
            x = target.split("_")[3]
            # print(temp)
            if(args.mode == 1):
                t1 = "../" + target + "/simulation/0/halfdata {} {} -{}\n".format(temp, kconstant, x)
                metadata.write(t1)
            elif(args.mode == 2):
                t2 = "../" + target + "/simulation/1/halfdata {} {} -{}\n".format(temp, kconstant, x)
                metadata.write(t2)
    metadata.close()
if(args.qnqc):
    qnqc()


if(args.pulling):
    print("2xov pulling")
    # os.system("cp folder_list .")
    kconstant = 0.01

    metadata = open("metadatafile", "w")
    with open('folder_list', 'r') as ins:
        for line in ins:
            target = line.strip(' \n')
            temp = target.split("_")[1]
            x = target.split("_")[3]
            # print(temp)
            if(args.server):
                if(args.mode == 1):
                    t1 = "/scratch/wl45/freeEnergy_2xov/pullingDistance/" + target + "/simulation/0/halfdata {} {} {}\n".format(temp, kconstant, x)
                    metadata.write(t1)
                elif(args.mode == 2):
                    t2 = "/scratch/wl45/freeEnergy_2xov/pullingDistance/" + target + "/simulation/1/halfdata {} {} {}\n".format(temp, kconstant, x)
                    metadata.write(t2)
            else:
                if(args.mode == 1):
                    t1 = "/Users/weilu/Research/server/freeEnergy_2xov/pullingDistance/" + target + "/simulation/0/halfdata {} {} {}\n".format(temp, kconstant, x)
                    metadata.write(t1)
                elif(args.mode == 2):
                    t2 = "/Users/weilu/Research/server/freeEnergy_2xov/pullingDistance/" + target + "/simulation/1/halfdata {} {} {}\n".format(temp, kconstant, x)
                    metadata.write(t2)
    metadata.close()



def qnqc2():
    print("QnQc2")
    # os.system("cp folder_list .")
    kconstant = 0.02

    metadata = open("metadatafile", "w")
    with open('folder_list', 'r') as ins:
        for line in ins:
            target = line.strip('\n')
            temp = target.split("_")[1]
            x = target.split("_")[3]
            # print(temp)
            if(args.mode == 1):
                t1 = "../" + target + "/simulation/0/halfdata {} {} {}\n".format(temp, kconstant, x)
                metadata.write(t1)
            elif(args.mode == 2):
                t2 = "../" + target + "/simulation/1/halfdata {} {} {}\n".format(temp, kconstant, x)
                metadata.write(t2)
    metadata.close()
if(args.qnqc2):
    qnqc2()

# def qnqc2():
#     print("QnQc2")
#     # os.system("cp folder_list .")
#     kconstant = 400   # double the k constant
#     temp = 350
#     q0 = 0.12
#     metadata = open("metadatafile", "w")
#     for i in range(40):
#         q = q0 + i*0.02
#         target = "../simulation/350/" + str(i) + "/halfdata.dat {} {} {:.2f}\n".format(temp, kconstant, q)
#         metadata.write(target)
#     metadata.close()
#
# if(args.qnqc2):
#     qnqc2()
# temp_list = range(350, 400, 100)
# # q0_list = range(10, 100, 10)
#
# # config = open('config.py', 'w')
# # config.write("number_of_run = %d\nsimulation_steps = %d\n\
# # warm_up_steps = %d\n" % (n, simulation_steps, warm_up_steps))
# # config.close()
# kconstant = 200
# n = 50
# q_bias_step = 0.02
#
# metadata = open("metadatafile", "w")
# for temp in temp_list:
#     q0 = 0
#     for i in range(n):
#         q0 += q_bias_step
#         metadata.write("../simulation/%.0f/%d/halfdata.dat %.0f %d %.2f\n" % (temp, i, temp, kconstant, q0))
#
# metadata.close()
# sys.exit(0)


# print("hello world")
