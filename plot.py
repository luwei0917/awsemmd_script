#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import imp

# mypath = os.environ["PATH"]
# os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
# my_env = os.environ.copy()

parser = argparse.ArgumentParser(
    description="Plot my graphs quickly")

parser.add_argument("--qnqc", help="for all calculate q of n terminal and q of c terminal ", action="store_true", default=False)
parser.add_argument("--qnqc2", help="for all calculate q of n terminal and q of c terminal ", action="store_true", default=False)

parser.add_argument("-t", "--temperature", type=int, default=400,
                    help="temperature")
args = parser.parse_args()
# protein_name = args.template.strip('/')
# os.system("cp ~/opt/plot_scripts/free_energy.plt .")
# os.system("gnuplot free_energy.plt ")
# os.system("open free_energy.pdf")
#


def qnqc2():

if(args.qnqc2):
    qnqc2()


def qnqc():
    os.system("cp ~/opt/pulling/freeEnergy.gp .")
    os.system("cp ~/opt/pulling/oneFreeEnergy.gp .")
    os.system("gnuplot freeEnergy.gp")
    os.system("gnuplot oneFreeEnergy.gp")
    os.system("open freeEnergy.pdf ")
    os.system("open oneFreeEnergy.pdf")
if(args.qnqc):
    qnqc()
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
