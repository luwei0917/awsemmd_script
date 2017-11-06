#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
from Bio.PDB import *
import numpy as np



parser = argparse.ArgumentParser(description="This is used to compute the contact")
parser.add_argument("-m", "--mode", help="Test run", type=int, default=-1)
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("protein", help="the name of protein file")
args = parser.parse_args()

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir


if args.mode == 2:
    protein_name = args.protein.split('.')[0]
    raptorFile = "raptor." + protein_name + ".dat"
    with open("raptor.csv", "w") as out:
        out.write("i, j, Probability\n")
        with open(raptorFile, "r") as f:
            for line in f:
                oneLine = line.split()
                if len(oneLine) == 5:
                    i = int(oneLine[0])
                    j = int(oneLine[1])
                    p = float(oneLine[4])
                    out_line = str(i) + ", " + str(j) + ", " + str(p) + "\n"
                    out.write(out_line)
    # with open("raptor.csv", "w") as out:
    #     out.write("i, j, Distance\n")

if args.mode == 1:
    parser = PDBParser()
    target = args.protein
    structure = parser.get_structure('target', target)
    # Assume only one chain
    model = structure[0]
    # chain = model['A']
    chain = model.child_list[0]
    with open("contact.csv", "w") as out:
        out.write("i, j, Distance\n")
        for i_residue in chain:
            is_regular_res = i_residue.has_id('CA')
            # print(i_residue, is_regular_res)
            if i_residue.get_resname() == "GLY":
                res_i = i_residue['CA']
            else:
                res_i = i_residue['CB']
            for j_residue in chain:
                if j_residue.get_resname() == "GLY":
                    res_j = j_residue['CA']
                else:
                    res_j = j_residue['CB']
                distance = res_i-res_j
                i = i_residue.id[1]
                j = j_residue.id[1]
                out_line = str(i) + ", " + str(j) + ", " + str(distance) + "\n"
                out.write(out_line)
