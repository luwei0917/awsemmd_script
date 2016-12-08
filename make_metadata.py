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

args = parser.parse_args()


def qnqc():
    print("QnQc")
    # os.system("cp folder_list .")
    kconstant = 0.05

    metadata = open("metadatafile", "w")
    with open('folder_list', 'r') as ins:
        for line in ins:
            target = line.strip('\n')
            temp = target.split("_")[1]
            x = target.split("_")[3]
            # print(temp)
            # t1 = "../" + target + "/simulation/0/halfdata {} {} -{}\n".format(temp, kconstant, x)
            t2 = "../" + target + "/simulation/1/halfdata {} {} {}\n".format(temp, kconstant, x)
            # metadata.write(t1)
            metadata.write(t2)
    metadata.close()
if(args.qnqc):
    qnqc()
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
