#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import myPersonalFunctions

parser = argparse.ArgumentParser(description="take two pdb, align them, and show them")
parser.add_argument("p1", help="first pdb")
parser.add_argument("p2", help="second pdb")
parser.add_argument("-p", "--plot", action="store_true", default=False, help="Plot the result")
args = parser.parse_args()

try:
    protein1_name,_ = args.p1.split('.')
except:
    try:
        protein1_name,_ = args.p1.split('/')
    except:
        protein1_name = args.p1
        print("ATTENSION! protein1 name {}\n Correct?".format(protein1_name))
try:
    protein2_name,_ = args.p2.split('.')
except:
    try:
        protein2_name,_ = args.p2.split('/')
    except:
        protein2_name = args.p2
        print("ATTENSION! protein2 name {}\n Correct?".format(protein2_name))


do = os.system
cd = os.chdir


do(f"~/opt/TMalign/TMalign {protein1_name}.pdb {protein2_name}.pdb -o result")
do("cp result_all_atm result_all_atm.pdb")
do("grep 'TM-score=' result > tmscore.dat")
# Seq_ID=n_identical/n_aligned
do("cat tmscore.dat")
with open("tmscore.dat", "r") as f:
    for line in f:
        aligned_length,rmsd,tmscore,seqid = line.split(",")
        aligned_length = int(aligned_length.split("=")[1])
        rmsd = float(rmsd.split("=")[1])
        tmscore = float(tmscore.split("=")[1])
        seqid = float(seqid.split("=")[1])
        print("aligned_length, rmsd, tmscore, seqid")
        print(aligned_length, rmsd, tmscore, seqid)
if(args.plot):
    do("pymol ~/opt/plot_scripts/tmalign_all.pml")

