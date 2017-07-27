#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
from Bio.PDB import *
import numpy as np



parser = argparse.ArgumentParser(description="This is used to compute the contact")
parser.add_argument("-m", "--mode", help="Test run", type=int, default=0)
parser.add_argument("-d", "--debug", action="store_true", default=False)
args = parser.parse_args()

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir


def find_seq_and_save(name, getFrom, getTo, title):
    with open("casp11.seq", "r") as f:
        seq = ""
        for line in f:
            if line[0] == ">":
                if line[1:6] == name:
                    firstLine = line
                    seq = ""
            else:
                seq += line.strip()
    fileName = "fasta/" + name + ".fasta"
    with open(fileName, "w") as out:
        out.write(title + "\n")
        out.write(seq[(getFrom-1):getTo] + "\n")


if args.mode == 1:
    with open("chosen.csv", "r") as f:
        next(f)
        for line in f:
            tmp = line.split(",")
            Name = tmp[0]
            length = tmp[2]
            getFrom = int(tmp[5].split("-")[0])
            getTo = int(tmp[5].split("-")[1])
            pdb = tmp[6]
            print(Name, length, pdb)
            title = ">" + Name + "_" + tmp[1] + ", " + pdb + ", " + length + " residues"
            do("mv all_pdbs/{0} pdbs/".format(Name+"-"+tmp[1]+".pdb"))
            # find_seq_and_save(Name, getFrom, getTo, title)
            # do("hhblits -i fasta/{0}.fasta  -n 1 -d pdb70_from_mmcif_25May17/pdb70 -o results/{0}.hhr".format(Name))
