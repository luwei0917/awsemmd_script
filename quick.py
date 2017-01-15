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
import numpy as np
from myPersonalFunctions import *
# import matplotlib.pyplot as plt

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

parser.add_argument("--qnqc", help="for all calculate q of n terminal and q of c terminal ", action="store_true", default=False)
parser.add_argument("--qnqc2", help="step2", action="store_true", default=False)
parser.add_argument("--wham", help="compute wham ", action="store_true", default=False)
args = parser.parse_args()


def wham():
    array = []
    cwd = os.getcwd()
    os.system("mkdir -p wham_analysis")
    os.system("rm wham_analysis/data")
    with open('folder_list', 'r') as ins:
        for line in ins:
            target = line.strip('\n')
            temp = target.split("_")[1]
            x = target.split("_")[3]
            if(temp == "300" and x != "50"):
                t1 = "../"+target + "/simulation/0"
                t2 = "../"+target + "/simulation/1"
                array.append(t1)
            # array.append(t2)
    for i in array:
        os.chdir(i)
        os.system("pwd")
        folder_name = os.getcwd()
        os.system("cat halfdata >> "+cwd+"/wham_analysis/data")
        os.chdir(cwd)

    # os.system("cp ~/opt/wham_analysis/*.m wham_analysis/")
    os.chdir("wham_analysis")
    os.system(" awk '{print $1}' data > qn")
    os.system(" awk '{print $2}' data > qc")
    os.system(" awk '{print $3}' data > qc2")
    os.system(" awk '{print $4}' data > e_total")
    os.system(" awk '{print $5}' data > p_total")
    # os.system("~/opt/script/wham/fused_calc_cv.sc wham_analysis/ 2xov 20 300 250 350 10 50 200 0.05 1")


if(args.wham):
    wham()


# def qnqc():
#     array = []
#     cwd = os.getcwd()
#     print(cwd)
#     with open('folder_list', 'r') as ins:
#         for line in ins:
#             target = line.strip('\n')
#             t1 = target + "/simulation/0"
#             t2 = target + "/simulation/1"
#             array.append(t1)
#             array.append(t2)
#     for i in array:
#         os.chdir(i)
#         os.system("pwd")
#         folder_name = os.getcwd()
#         os.system("cp ~/opt/pulling/qnqc.slurm .")
#         os.system("sbatch qnqc.slurm")
#         os.system("cp ~/opt/pulling/2xov_rerun.in .")
#         os.system("cp ~/opt/pulling/rerun.slurm .")
#         os.system("sbatch rerun.slurm")
#
#         os.chdir(cwd)
def qnqc():
    n = 20
    cwd = os.getcwd()
    print(cwd)
    for i in range(n):
        os.chdir(str(i))
        os.system("pwd")
        folder_name = os.getcwd()
        os.system("cp ~/opt/pulling/qnqc.slurm .")
        os.system("sbatch qnqc.slurm")

        os.chdir(cwd)
if(args.qnqc):
    qnqc()


def qnqc2():
    array = []
    cwd = os.getcwd()
    print(cwd)
    with open('folder_list', 'r') as ins:
        for line in ins:
            target = line.strip('\n')
            t1 = target + "/simulation/0"
            t2 = target + "/simulation/1"
            array.append(t1)
            array.append(t2)
    for i in array:
        os.chdir(i)
        os.system("pwd")
        folder_name = os.getcwd()
        os.system("tail -n+3 energy.log > energy")
        os.system("head -n 4000 energy > energy_all")
        os.system("tail -n 2000 energy_all > energy_half")
        os.system("awk '{print $17}' energy_half > etotal")
        os.system("sed '/^#/ d' x.colvars.traj > test")
        os.system("awk 'NR % 10 == 1 ' test > x")
        os.system("head -n 4000 x > x_all")
        os.system("tail -n 2000 x_all > x_half")
        os.system("awk '{print $2}' x_half > myx_half")
        os.system("paste etotal myx_half > halfdata")
        os.system("paste qn_half qc_half qc2_half etotal myx_half > halfdata")
        # os.system("pwd")
        os.chdir(cwd)

if(args.qnqc2):
    qnqc2()
# parser.add_argument("file", help="the name of file")
# parser.add_argument("file2", help="the name of file")
# vtotal = []
# with open(args.file) as f:
#     print(f.readline())
#     for line in f:
#         data = line.strip('\n').split(sep=" ")
#         data = [float(i) for i in data]
#         vtotal += [data[-3]]
#     # print(vtotal)
# v2 = []
# with open(args.file2) as f:
#     print(f.readline())
#     f.readline()
#     for line in f:
#         data = line.strip('\n').split(sep="\t")
#         data = [float(i) for i in data]
#         v2 += [data[-1]]
#     # print(v2)
# # plt.plot(vtotal)
# # plt.plot(v2)
# plt.plot(np.array(vtotal)-np.array(v2))
# plt.show()
# # sum =0
# # for i in data[1:-3]:
# #     sum += i
# # print(sum+(144.18+9.369))
