#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import myPersonalFunctions
import fileinput
from small_script.myFunctions import *

parser = argparse.ArgumentParser(
    description="The goal of this python3 code is to automatically create \
    zim file. Written by Wei Lu."
)
parser.add_argument("-f", "--fastaFile", default="crystal_structure.fasta", help="The name of the fasta file")



args = parser.parse_args()

seq = ""
with open(args.fastaFile, "r") as f:
    for line in f:
        if line[0] == ">":
            pass
        else:
            # print(line)
            seq += line.strip()
# create_zim(f"crystal.seq")
print(seq, len(seq))

data = read_hydrophobicity_scale(seq, isNew=False)
z = data["DGwoct"].values
np.savetxt("zim", z, fmt="%.2f")