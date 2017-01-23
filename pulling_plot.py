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
parser.add_argument("--gagb_compare", help="for all calculate q of n terminal and q of c terminal ", action="store_true", default=False)
parser.add_argument("--all_temp", type=int, default=0)
parser.add_argument("outname", nargs='?', help="output filename", default="test.png")
parser.add_argument("--temperature", type=int, default=340,
                    help="temperature")
parser.add_argument("-n", "--number", type=int, default=10,
                    help="number")
parser.add_argument("--minor", type=int, default=1,
                    help="minor control")
parser.add_argument("--one", action="store_true", default=False)
parser.add_argument("--multi", action="store_true", default=False)
parser.add_argument("-s", "--save", action="store_true", default=False)
parser.add_argument("-r", "--reproduce", default=None)
parser.add_argument("-t", "--test", action="store_true", default=False)
parser.add_argument("-m", "--mode", default="pulling")
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
    force_list = np.arange(1,2.5,0.1)
    array = []
    cwd = os.getcwd()
    print(cwd)
    for force in force_list:
        folder = "2d_2_force_" + str(force)
        cd(folder)
        do("pulling_plotcontour.py")
        do("cp test.png ../results/{}.png".format(force))
        cd(cwd)

if(args.save):
    # print(os.getcwd())
    # print(args)
    print("Saving")
    # print(datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
    with open("args"+datetime.datetime.now().strftime("%Y%m%d-%H%M"), "wb") as f:
        pickle.dump(args, f)
    os.system("cp ~/opt/plot.py plot_{}.py".format(datetime.datetime.now().strftime("%Y%m%d-%H%M")))

if(args.one):
    output = args.outname
    temp = args.temperature
    ax = plt.subplot(1, 1, 1)

    name = 'pmf-'+str(temp)+'.dat'
    data = pd.read_table(name, sep='\s+', comment='#', names=["bin","bin_center_1","f","df","e","s"])
    # print(data)
    # data.plot(ax=ax, x='bin_center_1', y='f', linewidth=5.0)
    data.plot(ax=ax, x='bin_center_1', y='f',xlim=(0, 200), label="T= \n"+str(temp))
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

if(args.multi):
    output = args.outname
    temp = args.temperature
    ax = plt.subplot(1, 1, 1)

    sequence = "a206g"
    force = 1.0
    i = 1
    force_list = [1.6, 1.8, 2.0]
    for force in force_list:
        target = str(i)+"_force_"+str(force)+"_"+sequence+"/"
        # target = "1_force_1.0_a206g/"
        name = target + 'pmf-'+str(temp)+'.dat'
        data = pd.read_table(name, sep='\s+', comment='#', names=["bin","bin_center_1","f","df","e","s"])
        # print(data)
        # data.plot(ax=ax, x='bin_center_1', y='f', linewidth=5.0)
        data.plot(ax=ax, x='bin_center_1', y='f',xlim=(0, 200), label=str(force))
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
            data.plot(ax=ax, x='bin_center_1', y='f',xlim=(0, 1), label="T= \n"+str(temp))
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
