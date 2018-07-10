#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import imp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import datetime
import pickle
matplotlib.style.use('fivethirtyeight')


# print(plt.style.available)
# mypath = os.environ["PATH"]
# os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
# my_env = os.environ.copy()

parser = argparse.ArgumentParser(
    description="Plot my graphs quickly")
# parser.add_argument("data", help="the name of data file")
parser.add_argument("--qnqc", help="for all calculate q of n terminal and q of c terminal ", action="store_true", default=False)
parser.add_argument("--qnqc_pull", help="for all calculate q of n terminal and q of c terminal ", action="store_true", default=False)

# parser.add_argument("--qnqc2", help="for all calculate q of n terminal and q of c terminal ", action="store_true", default=False)
parser.add_argument("--gagb", help="for all calculate q of n terminal and q of c terminal ", action="store_true", default=False)
parser.add_argument("--compare", help="for all calculate q of n terminal and q of c terminal ", action="store_true", default=False)
parser.add_argument("--all_temp", type=int, default=0)
parser.add_argument("outname", nargs='?', help="output filename", default="test.png")
parser.add_argument("--temperature", type=int, default=350,
                    help="temperature")
parser.add_argument("-n", "--number", type=int, default=10,
                    help="number")
parser.add_argument("--minor", type=int, default=1,
                    help="minor control")
parser.add_argument("-s", "--save", action="store_true", default=False)
parser.add_argument("-r", "--reproduce", default=None)
parser.add_argument("-t", "--test", action="store_true", default=False)
parser.add_argument("-m", "--mode", default="gagb")
parser.add_argument("-d", "--debug", action="store_true", default=False)
args = parser.parse_args()


if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir


if args.reproduce is not None:
    print("Reproducing!")
    with open(args.reproduce, "rb") as f:
        args = pickle.load(f)
        print(args)
        args.save = False

if(args.test):
    print("Hello Test World")
    output = args.outname
    temp = args.temperature
    ax = plt.subplot(1, 1, 1)
    # name_list = ["gb77", "gb88b", "gb91", "gb95", "gb", "ga", "ga95", "ga91", "ga88", "ga77"]
    name_list = ["gb77", "gb88b", "gb91", "gb95", "gb"]
    # name_list = ["ga", "ga95", "ga91", "ga88", "ga77"]
    # temp_list = range(300,330,10)
    for gagb in name_list:
        name = gagb+'-'+str(temp)+'.dat'
        data = pd.read_table(name, sep='\s+', comment='#', names=["bin","bin_center_1","f","df","e","s"])
        # print(data)
        data["f"] = data["f"] - data["f"].iloc[-1]
        # data.plot(ax=ax, x='bin_center_1', y='f', linewidth=5.0)
        data.plot(ax=ax, x='bin_center_1', y='f',xlim=(-1, 1), label="\n"+str(gagb))
    ax.set_xlabel("Q of gb")
    ax.set_ylabel("free energy(kT)")
    # ax.set_title("Force at 0.7 Kcal/mole-Angstrom")
    ax.legend_.remove()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.4))
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.gcf().subplots_adjust(left=0.15)
    plt.gcf().subplots_adjust(right=0.80)
    fig = plt.gcf()
    fig.savefig(output)
    os.system("open " + output)

if(args.save):
    # print(os.getcwd())
    # print(args)
    print("Saving")
    # print(datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
    with open("args"+datetime.datetime.now().strftime("%Y%m%d-%H%M"), "wb") as f:
        pickle.dump(args, f)
    os.system("cp ~/opt/gagb_plot.py plot_{}.py".format(datetime.datetime.now().strftime("%Y%m%d-%H%M")))


if(args.all_temp > 0):
    print("Hello World GaGb1")
    output = args.outname
    temp = args.temperature
    ax = plt.subplot(1, 1, 1)
    if(args.all_temp == 1):
        temp_list = range(310,400,20)
    elif(args.all_temp == 2):
        temp_list = [temp]
    # temp_list = range(300,330,10)

    if(args.mode == "gagb"):
        for temp in temp_list:
            name = 'pmf-'+str(temp)+'.dat'
            data = pd.read_table(name, sep='\s+', comment='#', names=["bin","bin_center_1","f","df","e","s"])
            # print(data)
            # data.plot(ax=ax, x='bin_center_1', y='f', linewidth=5.0)
            data.plot(ax=ax, x='bin_center_1', y='f',xlim=(-1, 1), label="T= \n"+str(temp))
        ax.set_xlabel("Q of gb")
    if(args.mode == "pulling"):
        for temp in temp_list:
            name = 'pmf-'+str(temp)+'.dat'
            data = pd.read_table(name, sep='\s+', comment='#', names=["bin","bin_center_1","f","df","e","s"])
            # print(data)
            # data.plot(ax=ax, x='bin_center_1', y='f', linewidth=5.0)
            data.plot(ax=ax, x='bin_center_1', y='f',xlim=(0, 400), label="T= \n"+str(temp))
        ax.set_xlabel("Distance(Å)")
    ax.set_ylabel("free energy(kT)")
    # ax.set_title("Force at 0.7 Kcal/mole-Angstrom")
    ax.legend_.remove()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.6))
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.gcf().subplots_adjust(left=0.15)
    plt.gcf().subplots_adjust(right=0.80)
    fig = plt.gcf()
    fig.savefig(output)
    os.system("open " + output)


if(args.compare):
    print("Hello World gagb_compare")
    output = args.outname
    temp = args.temperature
    ax = plt.subplot(1, 1, 1)
    # ax = plt.figure()
    ax.set_title("gagb")

    # name = 'pmf-'+str(temp)+'.dat'
    # target_list = ["gb77", "gb88b", "gb91", "gb95", "gb", "ga", "ga95", "ga88", "ga77"]
    target_list = ["gb77", "gb88b", "gb91", "gb95", "gb"]
    for target in target_list:
        name = target+"-"+str(temp)+'.dat'
        with open(name) as f:
            x = []
            y = []
            next(f)
            first = next(f)
            norm = 0
            for line in f:
                # print(line.split())
                bin_num, bin_center, freeEnergy, *rest = line.split()
                print(bin_num, bin_center, freeEnergy)
                x += [float(bin_center)]
                y += [float(freeEnergy)]
                if(bin_num == '15'):
                    norm = float(freeEnergy)
            for i in range(len(y)):
                y[i] -= norm
            plt.plot(x,y, label=target)
        # data = pd.read_table(name, sep='\s+', comment='#', names=["bin","bin_center_1","f","df","e","s"])
        # print(data["f"].iloc[-1])
        # print("!!")
        # print(data[data["bin_center_1"]==0.582])
        # # data["f"] = data["f"] - data[data['bin_center_1'] == "0.582"]["f"]
        # # if(target == "ga95"):
        # #     print(data["f"])
        # data.plot(ax=ax, x='bin_center_1', y='f', label=target, linewidth=3.0)
    # name = 'gb-'+str(temp)+'.dat'
    # data2 = pd.read_table(name, sep='\s+', comment='#', names=["bin","bin_center_1","f","df","e","s"])
    # data2.plot(ax=ax, x='bin_center_1', y='f', label="gb")
    plt.ylabel("free energy is units of kT")
    plt.xlabel("qo based on ga's sturcture")

    # ax.set_legend_bgcolor('white')

    ax.legend(loc='center left', bbox_to_anchor=(0.97, 0.8))
    # # ax.set_axis_bgcolor('white')
    # # ax.legend.set_facecolor('white')
    # # legend = plt.legend(frameon = 1)
    # # frame = legend.get_frame()
    # # frame.set_facecolor('green')
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.gcf().subplots_adjust(left=0.15)
    plt.gcf().subplots_adjust(right=0.80)
    fig = plt.gcf()
    fig.savefig(output)
    # fig.savefig(output, transparent=True)
    os.system("open " + output)
    # data.show()


def gagb():
    print("Hello World GaGb")
    output = args.outname
    temp = args.temperature

    name = 'pmf-'+str(temp)+'.dat'
    data = pd.read_table(name, sep='\s+', comment='#', names=["bin","bin_center_1","f","df","e","s"])
    print(data)
    data.plot(x='bin_center_1', y='f')
    fig = plt.gcf()
    fig.savefig(output)
    os.system("open " + output)
    # data.show()
if(args.gagb):
    gagb()
# def qnqc():
#     print("Hello World")
#     n = 2
#     f, axarr = plt.subplots(4, 4, sharex='col', sharey='row')

#     for i in range(n):
#         name = str(i) + "/qnqc"
#         data = pd.read_table(name, header=None, names=["qn", "qc"])

#         plt.tick_params(
#             axis='x',
#             which='both',
#             bottom='off',
#             top='off',
#             labelbottom='off'
#         )
#         data.plot(ax=axarr[i, i], legend=False)
#     # print(data)
#     # data.plot()
#     # plt.plot(data[1])
#     # plt.show()
#     # ts = pd.Series(np.random.randn(1000), index=pd.date_range('1/1/2000', periods=1000))
#     # df = pd.DataFrame(np.random.randn(1000, 4), index=ts.index, columns=list('ABCD'))
#     # df = df.cumsum()
#     # fig = plt.figure()

#     # plt.subplot(4, 4, 13)
#     # plt.plot(range(12))
#     fig = plt.gcf()
#     fig.savefig('figure.pdf')
#     os.system("open figure.pdf")
#     # data.show()


def qnqc():
    print("Hello World")
    n = args.number
    for i in range(n):
        name = str(i) + "/qnqc"
        data = pd.read_table(name, header=None, names=["qn", "qc"])
        ax = plt.subplot(4, 5, i + 1)
        plt.tick_params(
            axis='x',
            which='both',
            bottom='off',
            top='off',
            labelbottom='off'
        )
        plt.tick_params(
            axis='y',
            which='both',
            right='off',
        )
        data.plot(ax=ax, legend=False, ylim=(0, 1))
    # print(data)
    # data.plot()
    # plt.plot(data[1])
    # plt.show()
    # ts = pd.Series(np.random.randn(1000), index=pd.date_range('1/1/2000', periods=1000))
    # df = pd.DataFrame(np.random.randn(1000, 4), index=ts.index, columns=list('ABCD'))
    # df = df.cumsum()
    # fig = plt.figure()
    # plt.subplot(4, 4, 13)
    # plt.plot(range(12))
    plt.suptitle('Red is Qn, Blue is Qc')
    fig = plt.gcf()
    output = args.outname
    fig.savefig(output)
    os.system("open " + output)
    # fig.savefig('figure.pdf')
    # os.system("open figure.pdf")
    # data.show()
if(args.qnqc):
    qnqc()


def qnqc_pull():
    print("Hello World")
    n = args.number
    for i in range(n):
        if(args.minor == 1):
            name = "" + str(i) + "/addforce.dat"
        elif(args.minor == 2):
            name = "rerun_" + str(i) + "/addforce.dat"
        elif(args.minor == 3):
            name = "restart_" + str(i) + "/addforce.dat"
            name2 = "rerun_" + str(i) + "/addforce.dat"
            data2 = pd.read_table(name2, sep='\s+')
        data = pd.read_table(name, sep='\s+')
        # print(data)
        # data = pd.read_table(name, header=None, names=["qn", "qc"])
        ax = plt.subplot(4, 5, i + 1)
        plt.tick_params(
            axis='x',
            which='both',
            bottom='off',
            top='off',
            labelbottom='off'
        )
        plt.tick_params(
            axis='y',
            which='both',
            right='off',
        )
        if(args.minor == 1):
            data["addedFroce"] = -data["addedFroce"]
            data.plot(ax=ax, legend=False, x="position", y="addedFroce", xlim=(0, 400), ylim=(0, 2))
        elif(args.minor == 2):
            data["addedFroce"] = data["addedFroce"]/2.0
            data.plot(ax=ax, legend=False, x="Distance", y="addedFroce", xlim=(0, 400), ylim=(0, 2))
        elif(args.minor == 3):
            data["addedFroce"] = -data["addedFroce"]
            data.plot(ax=ax, legend=False, x="Distance", y="addedFroce", xlim=(0, 600), ylim=(0, 2))
            # data2["addedFroce"] = -data2["addedFroce"]
            data2.plot(ax=ax, legend=False, x="Distance", y="addedFroce", xlim=(0, 600), ylim=(0, 2))

    # print(data)
    # data.plot()
    # plt.plot(data[1])
    # plt.show()
    # ts = pd.Series(np.random.randn(1000), index=pd.date_range('1/1/2000', periods=1000))
    # df = pd.DataFrame(np.random.randn(1000, 4), index=ts.index, columns=list('ABCD'))
    # df = df.cumsum()
    # fig = plt.figure()
    # plt.subplot(4, 4, 13)
    # plt.plot(range(12))
    # plt.suptitle('Red is Qn, Blue is Qc')
    fig = plt.gcf()
    output = args.outname
    fig.savefig(output)
    os.system("open " + output)
    # fig.savefig('figure.pdf')
    # os.system("open figure.pdf")
    # data.show()
if(args.qnqc_pull):
    qnqc_pull()
# exec(open("config.py").read())
# # print(n, x, y, type(y))
# n = number_of_run
# steps = simulation_steps
# # protein_name
#
# # print(n, steps)
# # sys.exit(0)
#
# os.system("mkdir -p results")
# os.system("cp ~/opt/plot_scripts/qw_all.plt .")
# os.system("gnuplot -e 'number_of_run={}' qw_all.plt".format(n-1))
# for i in range(n):
#     print(i)
#     # analysis
#     os.chdir("analysis/"+str(i))
#     os.system("cp ~/opt/plot_scripts/*.plt .")
#     os.system("gnuplot qw.plt")
#     os.system("mv qw.pdf ../../results/qw_{0}.pdf".format(str(i)))
#     os.chdir("../..")
# protein_name = args.template.strip('/')
# os.system("cp ~/opt/plot_scripts/free_energy.plt .")
# os.system("gnuplot free_energy.plt ")
# os.system("open free_energy.pdf")
#
# print("Hello World GaGb")
# output = args.outname
# temp = args.temperature
# ax = plt.subplot(1, 1, 1)
#
# name = 'pmf-'+str(temp)+'.dat'
# data = pd.read_table(name, sep='\s+', comment='#', names=["bin","bin_center_1","f","df","e","s"])
# print(data)
# data.plot(ax=ax, x='bin_center_1', y='f', label=str(temp))
# ax.set_xlabel("Q of gb")
# fig = plt.gcf()
# fig.savefig(output)
# os.system("open " + output)
