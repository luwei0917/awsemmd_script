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
matplotlib.style.use('ggplot')
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
parser.add_argument("outname", help="output filename")
parser.add_argument("-t", "--temperature", type=int, default=330,
                    help="temperature")
parser.add_argument("-n", "--number", type=int, default=10,
                    help="number")
parser.add_argument("-m", "--minor", type=int, default=1,
                    help="minor control")
args = parser.parse_args()
# protein_name = args.template.strip('/')
# os.system("cp ~/opt/plot_scripts/free_energy.plt .")
# os.system("gnuplot free_energy.plt ")
# os.system("open free_energy.pdf")
#
def gagb_compare():
    print("Hello World gagb_compare")
    output = args.outname
    temp = args.temperature
    ax = plt.subplot(1, 1, 1)
    # ax = plt.figure()
    ax.set_title("ga77")

    # name = 'pmf-'+str(temp)+'.dat'
    name = "ga77gb_pmf-330.dat"
    data = pd.read_table(name, sep='\s+', comment='#', names=["bin","bin_center_1","f","df","e","s"])
    print(data)
    data.plot(ax=ax, x='bin_center_1', y='f', label="ga")
    name = "gagb_pmf-330.dat"
    data2 = pd.read_table(name, sep='\s+', comment='#', names=["bin","bin_center_1","f","df","e","s"])
    print(data2)
    data2.plot(ax=ax, x='bin_center_1', y='f', label="gb")
    plt.ylabel("free energy is units of kT")
    plt.xlabel("qo based on gb's sturcture")
    fig = plt.gcf()
    fig.savefig(output)
    os.system("open " + output)
    # data.show()
if(args.gagb_compare):
    gagb_compare()

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
            data.plot(ax=ax, legend=False, x="Distance", y="addedFroce", xlim=(0, 400), ylim=(0, 2))
        elif(args.minor == 2):
            data["addedFroce"] = data["addedFroce"]/2.0
            data.plot(ax=ax, legend=False, x="Distance", y="addedFroce", xlim=(0, 400), ylim=(0, 2))

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
