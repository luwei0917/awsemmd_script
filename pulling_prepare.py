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
from myPersonalFunctions import *
import glob
import numpy

parser = argparse.ArgumentParser(
    description="Prepare the data for run and analysis. \
                Codes here only need run once")

# parser.add_argument("protein", help="the name of protein")
# parser.add_argument("template", help="the name of template file")
parser.add_argument("-t", "--test", help="test ", action="store_true", default=False)
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("--distance", action="store_true", default=False)
parser.add_argument("--replace", action="store_true", default=False)
parser.add_argument("--summary", action="store_true", default=False)
parser.add_argument("--make_metadata", action="store_true", default=False)
parser.add_argument("--qnqc", action="store_true", default=False)
parser.add_argument("--data", action="store_true", default=False)
parser.add_argument("--continue_run", action="store_true", default=False)
parser.add_argument("-m", "--mode", type=int, default=1)
parser.add_argument("--submode", type=int, default=0)
args = parser.parse_args()


if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

# compute distance by "read dump file"
# if(args.test):
#     for i in range(40):
#         print(i)
#         cd(str(i))
#         do("read_dump_file.py")
#         cd("..")


def replace(TARGET, FROM, TO):
    do("sed -i.bak 's/{}/{}/g' {}".format(FROM,TO,TARGET))


def extract_data():
    # do("tail -n+3 energy.log | awk '{print $NF}' > etotal")
    do("awk '{print $17}' energy.dat  > etotal")
    do("head -n 6000 etotal | tail -n 2000 > etotal_half")
    do("head -n 6000 qn | tail -n 2000 > qn_half")
    do("head -n 6000 qc | tail -n 2000 > qc_half")
    do("head -n 6000 qo | tail -n 2000 > qo_half")
    do("paste qn qc etotal qo | tail -n 4000 > data")
    do("paste qn_half qc_half etotal_half qo_half > halfdata")
if args.mode == 8:
    print("Reorganize data, mode {}".format(args.mode))
    files = glob.glob("*")
    for one in files:
        cd(one)
        cd(str(args.submode))
        do("paste wham.dat lipid.dat | awk '{print $5+$7}' | sed 1d > e.dat")
        do("awk '{print $2}' addforce.dat |  sed 's/,$//' | sed 1d > d.dat")
        do("awk '{print $2}' wham.dat |  sed 's/,$//' | sed 1d > qw.dat")
        do("awk '{print $2/400}' addforce.dat |  sed 's/,$//' | sed 1d > d_2.dat")
        # do("paste membrane.dat rg.dat sum.dat d.dat | tail -n 600 > data")
        do("paste e.dat d.dat qw.dat d_2.dat | tail -n 2000 > data")
        cd("../..")

if args.mode == 7:
    print("Reorganize data, mode {}".format(args.mode))
    files = glob.glob("*")
    for one in files:
        cd(one)
        cd("0")
        do("paste wham.dat pressure.dat | awk '{print $5+$7}' | sed 1d > e.dat")
        do("awk '{print $2}' addforce.dat |  sed 's/,$//' | sed 1d > d.dat")
        do("awk '{print $2}' wham.dat |  sed 's/,$//' | sed 1d > qw.dat")
        do("awk '{print $2/400}' addforce.dat |  sed 's/,$//' | sed 1d > d_2.dat")
        # do("paste membrane.dat rg.dat sum.dat d.dat | tail -n 600 > data")
        do("paste e.dat d.dat qw.dat d_2.dat | tail -n 2000 > data")
        cd("../..")

if args.mode == 6:
    print("Reorganize data, mode {}".format(args.mode))
    files = glob.glob("*")
    for one in files:
        cd(one)
        do("cat cluster.dat >> ../cluster.dat")
        do("cat all_data >> ../all_data")
        # do("paste tmp tmp_2 > all_data")
        cd("..")

if args.mode == 5:
    print("Reorganize data, mode {}".format(args.mode))
    files = glob.glob("*")
    for one in files:
        cd(one)
        do("cat states.csv | sed 1d > cluster.dat")
        do("head -n 10000 data  > all_data")
        # do("paste tmp tmp_2 > all_data")
        cd("..")
if args.mode == 4:
    print("Reorganize data, mode {}".format(args.mode))
    files = glob.glob("*")
    for one in files:
        cd(one)
        cd("0")
        do("awk '{print $5}' wham.dat |  sed 's/,$//' | sed 1d > e.dat")
        do("awk '{print $2}' distance.dat |  sed 's/,$//' | sed 1,2d > d.dat")
        do("awk '{print $2}' wham.dat |  sed 's/,$//' | sed 1d > qw.dat")
        do("awk '{print $2}' addforce.dat |  sed 's/,$//' | sed 1d > d_2.dat")
        # do("paste membrane.dat rg.dat sum.dat d.dat | tail -n 600 > data")
        # do("paste e.dat d.dat qw.dat d_2.dat | tail -n 3000 > data")
        location = "/scratch/wl45/aug_2017/freeEnergy/data/{}".format(one)
        do("mkdir -p "+location)
        do("paste e.dat d.dat qw.dat d_2.dat > data")
        do("cp data " + location + "/")
        cd("../..")

if args.mode == 3:
    print("Reorganize data, mode {}".format(args.mode))
    files = glob.glob("*")
    for one in files:
        cd(one)
        cd("0")
        do("awk '{print $5}' wham.dat |  sed 's/,$//' | sed 1d > e.dat")
        do("awk '{print $2}' addforce.dat |  sed 's/,$//' | sed 1d > d.dat")
        do("awk '{print $2}' wham.dat |  sed 's/,$//' | sed 1d > qw.dat")
        do("awk '{print $2/400}' addforce.dat |  sed 's/,$//' | sed 1d > d_2.dat")
        # do("paste membrane.dat rg.dat sum.dat d.dat | tail -n 600 > data")
        do("paste e.dat d.dat qw.dat d_2.dat | tail -n 2000 > data")
        cd("../..")

if args.mode == 1:
    print("Reorganize data, mode 1")
    files = glob.glob("*")
    for one in files:
        cd(one)
        cd("0")
        do("awk '{print $5}' wham.dat |  sed 's/,$//' | sed 1d > e.dat")
        do("awk '{print $13}' energy.dat |  sed 's/,$//' | sed 1d > rg.dat")
        do("awk '{print $9}' energy.dat |  sed 's/,$//' | sed 1d > membrane.dat")
        do("paste membrane.dat rg.dat | awk '{print $1+$2}' > sum.dat")
        do("awk '{print $2}' addforce.dat |  sed 's/,$//' | sed 1d > d.dat")
        do("awk '{print $2}' wham.dat |  sed 's/,$//' | sed 1d > qw.dat")
        do("awk '{print $2/400}' addforce.dat |  sed 's/,$//' | sed 1d > d_2.dat")
        # do("paste membrane.dat rg.dat sum.dat d.dat | tail -n 600 > data")
        do("paste e.dat d.dat qw.dat d_2.dat | tail -n 600 > data")
        cd("../..")
if args.mode == 2:
    print("Reorganize data, mode 2")
    files = glob.glob("*")
    for one in files:
        cd(one)
        cd("0")
        do("awk '{print $5}' wham.dat |  sed 's/,$//' | sed 1d > e.dat")
        do("awk '{print $13}' energy.dat |  sed 's/,$//' | sed 1d > rg.dat")
        do("awk '{print $9}' energy.dat |  sed 's/,$//' | sed 1d > membrane.dat")
        do("paste membrane.dat rg.dat | awk '{print $1+$2}' > sum.dat")
        do("awk '{print $2}' addforce.dat |  sed 's/,$//' | sed 1d > d.dat")
        do("paste membrane.dat rg.dat sum.dat d.dat | tail -n 600 > data")
        # do("paste e.dat d.dat | tail -n 600 > data")
        cd("../..")
