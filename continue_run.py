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
        continues previous simulation")

parser.add_argument("template", help="the name of template file")
parser.add_argument("-n", "--number", type=int, default=20,
                    help="Number of simulation run")
parser.add_argument("-m", "--movie", help="generate the movie",
                    action="store_true")
parser.add_argument("-p", "--plotOnly", help="only generate the plot",
                    action="store_true")
args = parser.parse_args()

list_of_max_q = []

n = args.number
protein_name = args.template.strip('/')
