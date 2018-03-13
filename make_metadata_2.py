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
import glob
import numpy as np
# from run_parameter import *
parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    make metadata")

parser.add_argument("-k", "--kconstant", type=float, default=0.02)
parser.add_argument("-t", "--temperature", type=float, default=350)
parser.add_argument("-q", "--qStart", type=float, default=0.0)
parser.add_argument("-n", "--number", type=int, default=20)
parser.add_argument("--submode", type=int, default=-1)
parser.add_argument("-m", "--mode", type=int, default=0)
parser.add_argument("-a", "--additionalMode", type=int, default=1)
parser.add_argument("-p", "--pick", default=None)
# parser.add_argument("-t", "--test", help="test ", action="store_true", default=False)
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

if args.mode == 22:
    print("Distance biased 1D Free Energy, mode {}".format(args.mode))
    print("multiple temp")
    cwd = os.getcwd()
    files = glob.glob("../data/*")
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        for oneFile in files:
            tmp = oneFile.split("/")[-1].replace('.dat', '')
            t = tmp.split("_")[1]
            dis = tmp.split("_")[3]
            # print(tmp)
            if int(float(dis)) > 150:
                continue
            if int(t) == 450:
                target = "../{} {} {} {}\n".format(oneFile, t, kconstant, dis)
                out.write(target)
if args.mode == 21:
    print("Distance biased 1D Free Energy, mode {}".format(args.mode))
    print("multiple temp")
    cwd = os.getcwd()
    files = glob.glob("../data/*")
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        for oneFile in files:
            tmp = oneFile.split("/")[-1].replace('.dat', '')
            t = tmp.split("_")[1]
            dis = tmp.split("_")[3]
            # print(tmp)
            if int(float(dis)) > 150:
                continue
            if int(t)< 750:
                target = "../{} {} {} {}\n".format(oneFile, t, kconstant, dis)
                out.write(target)

if args.mode == 20:
    print("Distance biased 1D Free Energy, mode {}".format(args.mode))
    print("multiple temp")
    cwd = os.getcwd()
    files = glob.glob("../data/*")
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        for oneFile in files:
            tmp = oneFile.split("/")[-1].replace('.dat', '')
            t = tmp.split("_")[1]
            dis = tmp.split("_")[3]
            # print(tmp)
            if int(float(dis)) > 150:
                continue
            if int(t)<= 450:
                target = "../{} {} {} {}\n".format(oneFile, t, kconstant, dis)
                out.write(target)

if args.mode == 19:
    print("Distance biased 1D Free Energy, mode {}".format(args.mode))
    print("multiple temp")
    cwd = os.getcwd()
    files = glob.glob("../data/*")
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        for oneFile in files:
            tmp = oneFile.split("/")[-1].replace('.dat', '')
            t = tmp.split("_")[1]
            dis = tmp.split("_")[3]
            # print(tmp)
            if int(float(dis)) > 150:
                continue
            if int(t)== 500:
                target = "../{} {} {} {}\n".format(oneFile, t, kconstant, dis)
                out.write(target)

if args.mode == 18:
    print("Distance biased 1D Free Energy, mode {}".format(args.mode))
    print("multiple temp")
    cwd = os.getcwd()
    files = glob.glob("../data/*")
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        for oneFile in files:
            tmp = oneFile.split("/")[-1].replace('.dat', '')
            t = tmp.split("_")[1]
            dis = tmp.split("_")[3]
            # print(tmp)
            if int(float(dis)) > 150:
                continue
            if int(t)< 750 and int(t) > 390:
                target = "../{} {} {} {}\n".format(oneFile, t, kconstant, dis)
                out.write(target)


if args.mode == 17:
    print("Distance biased 1D Free Energy, mode {}".format(args.mode))
    print("multiple temp")
    cwd = os.getcwd()
    files = glob.glob("../data/*")
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        for oneFile in files:
            tmp = oneFile.split("/")[-1].replace('.dat', '')
            t = tmp.split("_")[1]
            dis = tmp.split("_")[3]
            # print(tmp)
            if int(t) < 500 or int(t) > 600:
                continue
            if int(float(dis)) % 2 == 1:
                continue
            target = "../{} {} {} {}\n".format(oneFile, t, kconstant, dis)
            out.write(target)

if args.mode == 16:
    print("Distance biased 1D Free Energy, mode {}".format(args.mode))
    print("multiple temp")
    cwd = os.getcwd()
    files = glob.glob("../data/*")
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        for oneFile in files:
            tmp = oneFile.split("/")[-1].replace('.dat', '')
            t = tmp.split("_")[1]
            dis = tmp.split("_")[3]
            # print(tmp)
            if int(t) < 400 or int(t) > 500:
                continue
            if int(float(dis)) % 2 == 1:
                continue
            target = "../{} {} {} {}\n".format(oneFile, t, kconstant, dis)
            out.write(target)

if args.mode == 15:
    print("Distance biased 1D Free Energy, mode {}".format(args.mode))
    print("multiple temp")
    cwd = os.getcwd()
    files = glob.glob("../data/*")
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        for oneFile in files:
            tmp = oneFile.split("/")[-1].replace('.dat', '')
            t = tmp.split("_")[1]
            dis = tmp.split("_")[3]
            # print(tmp)
            if int(t) < 500:
                continue
            target = "../{} {} {} {}\n".format(oneFile, t, kconstant, dis)
            out.write(target)


if args.mode == 14:
    print("Distance biased 1D Free Energy, mode {}".format(args.mode))
    print("multiple temp")
    cwd = os.getcwd()
    dis_list = np.linspace(30, 150, 121)
    temp_list = [500, 550, 600, 650]
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        for dis in dis_list:
            for temp in temp_list:
                q = dis
                if int(float(q)) % 2 == 0:
                    target = "../data/dis_{}_temp_{}.dat {} {} {}\n".format(dis, int(temp), temp, kconstant, q)
                    out.write(target)

if args.mode == 13:
    print("Distance biased 1D Free Energy, mode {}".format(args.mode))
    print("multiple temp")
    cwd = os.getcwd()
    dis_list = np.linspace(30, 150, 121)
    temp_list = [450, 500, 550, 600]
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        for dis in dis_list:
            for temp in temp_list:
                q = dis
                if int(float(q)) % 3 == 0:
                    target = "../data/dis_{}_temp_{}.dat {} {} {}\n".format(dis, int(temp), temp, kconstant, q)
                    out.write(target)

if args.mode == 12:
    print("Distance biased 1D Free Energy, mode {}".format(args.mode))
    print("multiple temp")
    cwd = os.getcwd()
    dis_list = np.linspace(30, 180, 151)
    temp_list = [450, 500, 550, 600]
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        for dis in dis_list:
            for temp in temp_list:
                q = dis
                if int(float(q)) % 3 == 0:
                    target = "../../../data/dis_{}_temp_{}.dat {} {} {}\n".format(dis, int(temp), temp, kconstant, q)
                    out.write(target)

if args.mode == 9 or args.mode == 10 or args.mode == 11:
    print("Distance biased 1D Free Energy, mode {}".format(args.mode))
    print("multiple temp")
    cwd = os.getcwd()
    # print(files)
    # temp = args.temperature
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        # temp_list = [350, 400, 450, 500]
        if args.mode == 10:
            temp_list = [args.temperature]
        if args.mode == 9:
            temp_list = [350, 400, 450, 500, 550, 600]
        if args.mode == 11:
            temp_list = [450, 500, 550, 600]
        files = glob.glob("../simulation/dis_*")
        for oneFile in files:
            for temp in temp_list:
                q = oneFile.split("_")[-1]
                if args.submode == 4:
                    if int(float(q)) % 5 == 0:
                        target = cwd + "/" + oneFile + "/0/t{}_new.dat {} {} {}\n".format(int(temp), temp, kconstant, q)
                        out.write(target)
                if args.submode == 3:
                    if float(q) < 130:
                        target = cwd + "/" + oneFile + "/0/t{}_new.dat {} {} {}\n".format(int(temp), temp, kconstant, q)
                        out.write(target)
                if args.submode == 1:
                    if float(q) > 131:
                        target = cwd + "/" + oneFile + "/0/t{}_new.dat {} {} {}\n".format(int(temp), temp, kconstant, q)
                        out.write(target)
                if args.submode == 2:
                    target = cwd + "/" + oneFile + "/0/t{}_new.dat {} {} {}\n".format(int(temp), temp, kconstant, q)
                    out.write(target)
# if args.mode == 10:
#     print("Distance biased 1D Free Energy, mode {}".format(args.mode))
#     print("One temp")
#     cwd = os.getcwd()
#     # print(files)
#     temp = args.temperature
#     kconstant = args.kconstant
#     with open("metadatafile", "w") as out:
#         files = glob.glob("../simulation/dis_*".format(temp))
#         for oneFile in files:
#             q = oneFile.split("_")[-1]
#             target = cwd + "/" + oneFile + "/0/t{}.dat {} {} {}\n".format(int(temp), temp, kconstant, q)
#             out.write(target)
#
# if args.mode == 9:
#     print("Distance biased 1D Free Energy, mode {}".format(args.mode))
#     print("multiple temp")
#     cwd = os.getcwd()
#     # print(files)
#     # temp = args.temperature
#     kconstant = args.kconstant
#     with open("metadatafile", "w") as out:
#         # temp_list = [350, 400, 450, 500]
#         temp_list = [350, 400, 450, 500, 550, 600]
#         files = glob.glob("../simulation/dis_*")
#         for oneFile in files:
#             for temp in temp_list:
#                 q = oneFile.split("_")[-1]
#                 target = cwd + "/" + oneFile + "/0/t{}.dat {} {} {}\n".format(temp, temp, kconstant, q)
#                 out.write(target)

if args.mode == 8:
    print("Distance biased 1D Free Energy, mode {}".format(args.mode))
    cwd = os.getcwd()
    # print(files)
    # temp = args.temperature
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        temp_list = [220, 230]
        for temp in temp_list:
            files = glob.glob("../temp_{}/simulation/dis_*".format(temp))
            for oneFile in files:
                q = oneFile.split("_")[-1]
                target = cwd + "/" + oneFile + "/0/data {} {} {}\n".format(temp, kconstant, q)
                out.write(target)


if args.mode == 7:
    print("Distance biased 1D Free Energy, mode {}".format(args.mode))
    cwd = os.getcwd()
    # print(files)
    temp = args.temperature
    kconstant = args.kconstant

    with open("metadatafile", "w") as out:
        temp_list = [400, 450, 500]
        for temp in temp_list:
            files = glob.glob("../temp_{}/simulation/dis_*".format(temp))
            for oneFile in files:
                q = oneFile.split("_")[-1]
                target = cwd + "/" + oneFile + "/0/data {} {} {}\n".format(temp, kconstant, q)
                out.write(target)


if args.mode == 6:
    print("Distance biased 1D Free Energy, mode {}".format(args.mode))
    files = glob.glob("../data/*")
    cwd = os.getcwd()
    # print(files)
    temp = args.temperature
    kconstant = 0.0
    with open("metadatafile", "w") as out:
        for oneFile in files:
            q = 0.0
            target = cwd + "/" + oneFile + "/all_data {} {} {}\n".format(temp, kconstant, q)
            out.write(target)

if args.mode == 5:
    print("Distance biased 1D Free Energy, mode {}".format(args.mode))
    files = glob.glob("simulation/*")
    cwd = os.getcwd()
    print(files)
    temp = args.temperature
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        for oneFile in files:
            q = int(oneFile.split("/")[-1])*0.02
            target = cwd + "/" + oneFile + "/0/data {} {} {}\n".format(temp, kconstant, q)
            out.write(target)


if args.mode == 1:
    print("Distance biased 1D Free Energy, mode 1")
    files = glob.glob("../simulation/dis_*")
    cwd = os.getcwd()
    # print(files)
    temp = args.temperature
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        for oneFile in files:
            q = oneFile.split("_")[-1]
            target = cwd + "/" + oneFile + "/0/data {} {} {}\n".format(temp, kconstant, q)
            out.write(target)

if args.mode == 2:
    print("Distance biased 1D Free Energy, mode 2")
    files = glob.glob("../simulation/dis_*")
    # print(files)
    temp = args.temperature
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        for oneFile in files:
            q = oneFile.split("_")[-1]
            if float(q) < 100:
                target = "../" + oneFile + "/0/data {} {} {}\n".format(temp, kconstant, q)
                out.write(target)

if args.mode == 3:
    print("Distance biased 1D Free Energy, mode 3")
    files = glob.glob("../simulation/dis_*")
    cwd = os.getcwd()
    # print(files)
    temp = args.temperature
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        for oneFile in files:
            q = oneFile.split("_")[-1]
            target = cwd + "/" + oneFile + "/{}/data {} {} {}\n".format(args.submode, temp, kconstant, q)
            out.write(target)

if args.mode == 4:
    print("Distance biased 1D Free Energy, mode 4")
    files = glob.glob("simulation/dis_*")
    cwd = os.getcwd()
    # print(files)
    temp = args.temperature
    kconstant = args.kconstant
    with open("metadatafile", "w") as out:
        for oneFile in files:
            q = float(oneFile.split("_")[-1])/400
            target = cwd + "/" + oneFile + "/0/data {} {} {}\n".format(temp, kconstant, q)
            out.write(target)


# if(args.test):
#     kconstant = args.kconstant * 2   # double the k constant
#     q0 = args.qStart
#     metadata = open("metadatafile", "w")
#     if(args.mode == 1):
#         qStep = 0.05
#         temp_list = [135, 160, 185, 210]
#         # temp_list = [160]
#         # temp_list = ['300', '200', '250']
#         if(args.additionalMode == 1):
#             pre_fix = "first_2000_"
#         else:
#             pre_fix = ""
#         for i in range(args.number):
#             q = q0 + i * qStep
#             name = pre_fix + str(i)
#             for temp in temp_list:
#                 target = "../data/{}/".format(temp) + name + " {} {} {:.2f}\n".format(temp, kconstant, q)
#                 # target = "../data/t_{}/small_".format(temp) + str(i) + " {} {} {:.2f}\n".format(temp, kconstant, q)
#                 # target = "../data/{}/".format(temp) + name + " {} {} {:.2f}\n".format(temp, kconstant, q)
#                 metadata.write(target)
#         metadata.close()
    # See pulling prepare
    # if(args.mode == 2):
    #     print("2xov pulling")
    #     # os.system("cp folder_list .")
    #     kconstant = 0.04
    #     metadata = open("metadatafile", "w")
    #     with open('../folder_list', 'r') as ins:
    #         for line in ins:
    #             target = line.strip(' \n')
    #             temp = target.split("_")[1]
    #             x = target.split("_")[3]
    #             # print(temp)
    #             cwd = os.getcwd()
    #             t1 = cwd + "../" + target + "/halfdata {} {} {}\n".format(temp, kconstant, x)
    #             metadata.write(t1)
    #     metadata.close()

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
            for j in range(0, 13):
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