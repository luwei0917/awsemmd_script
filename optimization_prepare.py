#!/usr/bin/env python3
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
# import imp
import glob
from time import sleep
import fileinput
import numpy as np
import pandas as pd
import subprocess
from small_script.myFunctions import *
from collections import defaultdict

parser = argparse.ArgumentParser(description="Prepare the optimization folder")
parser.add_argument("DatabaseFolder", help="your database folder")
parser.add_argument("OptimizationFolder", help="your optimization folder")
args = parser.parse_args()

# if args.test:
#     do = print
# else:
#     do = os.system
with open('log_optimization_preparation.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')


def do(cmd, get=False, show=True):
    if get:
        out = subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True).communicate()[0].decode()
        if show:
            print(out, end="")
        return out
    else:
        return subprocess.Popen(cmd, shell=True).wait()
cd = os.chdir

def get_pdbList(pdbFolderList):
    # filter out those has incomplete peptide.txt
    _all = []
    # generate decoys
    protein_list = []
    for pdbFolder in pdbFolderList:
        # pdbName = pdbFolder.split("_")[-1]
        pdbName = pdbFolder.split("/")[-1]
        source = pdbFolder + f"/*.pdb"
        p_list = glob.glob(source)
        # print(p_list, source)
        # assert len(p_list) == 1
        p = p_list[0]
        chain_seq, seq = getChainsAndSeq(p)

        # print(pdbName, i, len(seq), len(chain_seq))
        decoy_list = []
        if not os.path.exists(f"{pdbFolder}/peptide.txt"):
            print(pdbFolder, "not exist")
            continue
        with open(f"{pdbFolder}/peptide.txt") as f:
            for line in f:
                pep = line.strip()
                pep_len = len(pep)
                # assert len(pep) == 9
                for c in list(set(chain_seq)):
                    if chain_seq.count(c) == pep_len:
                        first_c = chain_seq.find(c)
                        a = list(seq)
                        a[first_c:first_c+pep_len] = pep
                        assert len(set(chain_seq[first_c:first_c+pep_len])) == 1
                        decoy = "".join(a)
                        decoy_list.append(decoy)
        _all.append([pdbName, len(decoy_list)])
    a = pd.DataFrame(_all, columns=["Name", "Length"])
    if len(a.query("Length != 1000")) > 0:
        print(a.query("Length < 1000"))
    pdbList = list(a.query("Length == 1000")["Name"])
    return pdbList

def mycp(source, target):
    os.system(f"cp {source} {target}")

do = os.system

def getSeq(fileLocation):
    p = PDBParser()
    s = p.get_structure("test", fileLocation)
    seq = ""
    residues = list(s.get_residues())
    for residue in residues:
        res_id = residue.get_id()[0]
        if res_id==' ':
            residue_name = residue.get_resname()
            seq += three_to_one(residue_name)
    return seq

# get chains nad seq
def getChainsAndSeq(fileLocation):
    # fileLocation = "/Users/weilu/Research/examples/optimization/optimization/Structure_Ensemble/1.pdb"
    p = PDBParser()
    pdb = p.get_structure("test", fileLocation)
    residues = list(pdb.get_residues())
    seq = ""
    chains = ""
    for residue in residues:
        res_id = residue.get_id()[0]
        chain = residue.get_full_id()[2]
        if res_id==' ':
            residue_name = residue.get_resname()
            seq += three_to_one(residue_name)
            chains += chain
    return chains, seq

def create_folders_and_copy_pdbs_and_setup_decoys(pre, pdbFolderList, pdbList, n_decoys=-1):
    # n_decoys = -1 means infer from decoy_list size.
    do(f"mkdir -p {pre}")
    do(f"mkdir -p {pre}/database/dompdb")
    do(f"mkdir -p {pre}/database/S20_seq")
    do(f"mkdir -p {pre}/optimization/decoys/shuffle")
    do(f"mkdir -p {pre}/phis")
    do(f"mkdir -p {pre}/optimization/proteins_name_list")
    do(f"mkdir -p {pre}/optimization/slurms")
    do(f"mkdir -p {pre}/optimization/gammas")
    do(f"mkdir -p {pre}/optimization/outs")
    # generate decoys
    protein_list = []
    for pdbFolder in pdbFolderList:
        pdbName = pdbFolder.split("/")[-1]
        if pdbName not in pdbList:
            continue
        source = pdbFolder + f"/*.pdb"
        p_list = glob.glob(source)
        # print(p_list, source)
    #     assert len(p_list) == 1
        p = p_list[0]
        chain_seq, seq = getChainsAndSeq(p)

        # print(pdbName, i, len(seq), len(chain_seq))
        decoy_list = []
        f = open(f"{pdbFolder}/native.txt")
        native = f.readlines()[0].strip()
        pep_len = len(native)
        # assert len(pep) == 9
        for c in list(set(chain_seq)):
            if chain_seq.count(c) == pep_len:
                first_c = chain_seq.find(c)
                a = list(seq)
                a[first_c:first_c+pep_len] = native
                assert len(set(chain_seq[first_c:first_c+pep_len])) == 1
                native_seq = "".join(a)
        with open(f"{pdbFolder}/peptide.txt") as f:
            for line in f:
                pep = line.strip()
                pep_len = len(pep)
                # assert len(pep) == 9
                for c in list(set(chain_seq)):
                    if chain_seq.count(c) == pep_len:
                        first_c = chain_seq.find(c)
                        a = list(seq)
                        a[first_c:first_c+pep_len] = pep
                        assert len(set(chain_seq[first_c:first_c+pep_len])) == 1
                        decoy = "".join(a)
                        decoy_list.append(decoy)
        # two branch. one: all become one.
        if len(decoy_list) == 0:
            continue
        if n_decoys == -1:
            n_decoys == len(decoy_list)
        # assert len(decoy_list) == 1000
        for i in range(10):
            # idx = i*22 + 1
            idx = i*9 + 1
            p = os.path.dirname(p_list[0]) + f"/{idx}.pdb"
            opt_pdbName = f"{pdbName}_{idx}"
            fileLocation = f"{pre}/optimization/decoys/shuffle/{opt_pdbName}.decoys"
            # print("1", pdbName, len(decoy_list), pep_len)
            with open(fileLocation, "w") as out:
                # if len(decoy_list) != 1000:
                #     print("wrong")
                #     decoy_list = random.choices(decoy_list, k=1000)
                for i_decoy in range(n_decoys):
                    decoy = decoy_list[i_decoy]
                    out.write(decoy+"\n")
            protein_list.append(opt_pdbName)


            target = f"{pre}/database/dompdb/{opt_pdbName}.pdb"
            ## move native pdbs to dompdb
            mycp(p, target)
            ## move native seq to S20_seq
    #         seq = getSeq(target)
            seq = native_seq
            fileLocation = f"{pre}/database/S20_seq/{opt_pdbName}.seq"
            with open(fileLocation, "w") as out:
                out.write(seq+"\n")
    return protein_list

folder = args.DatabaseFolder
pdbFolderList = glob.glob(f"{folder}/*")
print(pdbFolderList)
pdbList = get_pdbList(pdbFolderList)
print(pdbList)
pre = args.OptimizationFolder
# create_folders_and_copy_pdbs_and_setup_decoys(pre, pdbFolderList, pdbList, n_decoys=-1)
protein_list = create_folders_and_copy_pdbs_and_setup_decoys(pre, pdbFolderList, pdbList, n_decoys=10)



## write protein_list
fileLocation = f"{pre}/optimization/protein_list"
with open(fileLocation, "w") as out:
    for pdbName in protein_list:
        out.write(f"{pdbName}\n")

phi_list_contact = '''\
phi_pairwise_contact_well 4.5 6.5 5.0 10
phi_density_mediated_contact_well 6.5 9.5 5.0 10 2.6 7.0
phi_burial_well 4.0
'''
phi_list = phi_list_contact
with open(f"{pre}/optimization/phi_list.txt", "w") as out:
    out.write(phi_list)


