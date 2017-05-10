#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import myPersonalFunctions

parser = argparse.ArgumentParser(description="This is my aawsem project headquarter")
parser.add_argument("-t", "--test", help="Test run", action="store_true", default=False)
parser.add_argument("-m", "--mode",
                    type=int, default=1)
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("protein", help="the name of protein file")
parser.add_argument("-p", "--plot", action="store_true", default=False, help="Plot the result")
parser.add_argument("--casp", action="store_true", default=False)
parser.add_argument("--step", type=int, default=8000,
                    help="Which step to show")
parser.add_argument("--tmalign", action="store_true", default=False)
args = parser.parse_args()

try:
    protein_name,_ = args.protein.split('.')
except:
    try:
        protein_name,_ = args.protein.split('/')
    except:
        protein_name = args.protein
        print("ATTENSION! protein name {}\n Correct?".format(protein_name))
if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir


if(args.tmalign):
    do("")


if(args.test):
    do("BuildAllAtomsFromLammps_seq.py dump.lammpstrj awsem {}.seq 8000".format(protein_name))
    do("~/opt/TMalign/TMalign awsem.pdb {}.pdb -o result".format(protein_name))
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

if(args.casp):
    if(args.mode == 1):         # Only calculate the last frame
        do("BuildAllAtomsFromLammps_seq.py dump.lammpstrj awsem {}.seq 8000".format(protein_name))
        do("~/opt/TMalign/TMalign awsem.pdb ~/opt/AAWSEM/casp11/targets/{}.pdb -o result".format(protein_name))
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
    if(args.mode == 2):                 # Compute for dump structure every stride number
        stride = 1000
        start = 0
        tmp = start
        end = 8000
        do("rm tmscore_all_tmp.dat")
        do("rm -r aligned_structures")
        do("mkdir aligned_structures")
        while tmp <= end:
            do("python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj awsem {}.seq {}".format(protein_name, tmp))
            do("~/opt/TMalign/TMalign awsem.pdb ~/opt/AAWSEM/casp11/targets/{}.pdb -o result".format(protein_name))
            do("grep 'TM-score=' result >> tmscore_all_tmp.dat")
            do("cp result_all_atm aligned_structures/{}_{}".format(protein_name, tmp))
            # Seq_ID=n_identical/n_aligned
            tmp += stride
        with open("tmscore_all.dat", "w") as out:
            out.write("Step, Aligned_length, RMSD, TMscore, Seqid\n")
            with open("tmscore_all_tmp.dat", "r") as f:
                tmp = start
                for line in f:
                    aligned_length,rmsd,tmscore,seqid = line.split(",")
                    aligned_length = int(aligned_length.split("=")[1])
                    rmsd = float(rmsd.split("=")[1])
                    tmscore = float(tmscore.split("=")[1])
                    seqid = float(seqid.split("=")[1])
                    out.write("{}, {}, {}, {}, {}\n".format(tmp, aligned_length, rmsd, tmscore, seqid))
                    tmp += stride

if(args.plot):
    do("cp ~/opt/plot_scripts/tmalign_all_tmp.pml .")
    do(  # replace NAME with target file name
        "sed -i.bak 's/NAME/" +
        str(args.protein) +
        "/g' tmalign_all_tmp.pml")
    do("pymol tmalign_all_tmp.pml")
