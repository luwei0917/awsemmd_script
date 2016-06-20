#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess

mypath = os.environ["PATH"]
os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
my_env = os.environ.copy()

parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        automatically analysis the simulation")

parser.add_argument("template", help="the name of template file")
parser.add_argument("-n", "--number", type=int, default=20,
                    help="Number of simulation run")
args = parser.parse_args()

list_of_max_q = []

n = args.number
protein_name = args.template.strip('/')

for i in range(n):
    # move
    os.chdir("analysis/"+str(i))
    os.system("cp final.png ../../results/final_"+str(i)+".png")
    os.chdir("../..")
