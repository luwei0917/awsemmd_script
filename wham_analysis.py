#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import imp

mypath = os.environ["PATH"]
os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
my_env = os.environ.copy()

parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        automatically prepare for wham analysis")

# parser.add_argument("protein", help="the name of protein")
parser.add_argument("--nick", action="store_true", default=False)
parser.add_argument("--zheng", action="store_true", default=False)
parser.add_argument("-m", "--mode", type=int, default=1)
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("--run", action="store_true", default=False)
args = parser.parse_args()

# protein_name = args.protein.strip('/')


if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

if(args.run):
    with open("run.sh", "w") as f:
        f.write("~/opt/script/wham/fused_calc_cv.sc . 2xov 40 300 150 250 10 30 200 0.12 0.9\n")
if(args.zheng):
    # do("mkdir wham")
    # cd("wham")
    kconstant = 400   # double the k constant
    q0 = 0.12
    metadata = open("metadatafile", "w")
    for i in range(40):
        q = q0 + i*0.02
        temp_list = [300]
        for temp in temp_list:
            target = "../data/" + str(i) + "/tinydata {} {} {:.2f}\n".format(temp, kconstant, q)
            metadata.write(target)
    metadata.close()

# cmd = "~/opt/script/wham/fused_calc_cv.sc {} {} 40 300 100 200 10 50 200 0.12 0.90"
# location = "test"
# cmd = cmd.format(location, protein_name)
# # os.system("~/opt/script/wham/fused_calc_cv.sc top_he_t400_q100/ top7 50 400 350 450 5 50 100 0 0.98")
# os.system(cmd)
