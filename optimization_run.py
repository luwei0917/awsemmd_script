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

parser = argparse.ArgumentParser(description="Compute phis under the optimization folder")
# parser.add_argument("DatabaseFolder", help="your database folder")
# parser.add_argument("OptimizationFolder", help="your optimization folder")
parser.add_argument("-m", "--mode", type=int, default=0)
parser.add_argument("-l", "--label", type=str, default="label")
args = parser.parse_args()

# if args.test:
#     do = print
# else:
#     do = os.system
with open('log_optimization_run.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')


# shuffle iter0.
# ligand optimization.
from pyCodeLib import *
import warnings
warnings.filterwarnings('ignore')

scavenge_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=commons
#SBATCH --partition=scavenge
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=04:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH -o outs/slurm-%j.out
echo "My job ran on:"
echo $SLURM_NODELIST
srun {}\n'''

base_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH -o outs/slurm-%j.out
echo "My job ran on:"
echo $SLURM_NODELIST
srun {}\n'''


n_decoys = 10
separateDecoysNum = -1
# template = base_slurm
template = scavenge_slurm

def do(cmd, get=False, show=True):
    if get:
        out = subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True).communicate()[0].decode()
        if show:
            print(out, end="")
        return out
    else:
        return subprocess.Popen(cmd, shell=True).wait()
cd = os.chdir

def slurmRun(slurmFileName, cmd, template=scavenge_slurm, memory=1):
    with open(slurmFileName, "w") as out:
        out.write(template.format(cmd))
        # out.write(scavenge_slurm.format(f"python3 ~/opt/compute_phis.py -m 1 proteins_name_list/proteins_name_list_{name}.txt"))
    replace(slurmFileName, "#SBATCH --mem-per-cpu=1G", f"#SBATCH --mem-per-cpu={memory}G")
    a = getFromTerminal(f"sbatch {slurmFileName}")
    jobId = a.split(" ")[-1].strip()
    return jobId


if args.mode == 1:
    # time.sleep(36000)
    with open("protein_list") as f:
        content = f.readlines()
    pos = 0
    i = 0
    n = len(content)
    # n = 100  # for testing
    while pos < n:
        with open(f"proteins_name_list/proteins_name_list_{i}.txt", "w") as out:
            for ii in range(1):
                if pos < n:
                    out.write(content[pos])
                pos += 1
            i += 1
    print(i)
    n = i
    i = 0
    jobIdList = []
    for i in range(n):
        proteins = f"proteins_name_list/proteins_name_list_{i}.txt"
        # generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")
        jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/optimization_run.py -m 222 -l {proteins}", template=template)
        # jobId = slurmRun(f"slurms/run_{i}.slurm", f"python3 ~/opt/compute_phis.py -m 0 proteins_name_list/proteins_name_list_{i}.txt")
        jobIdList.append(jobId)
        do(f"cat {proteins} >> iter0_complete.txt")
    # exit()
    waitForJobs(jobIdList, sleepInterval=300)
    with open(f"slurms/run_on_scavenge.slurm", "w") as out:
        out.write(base_slurm.format(f"python3 ~/opt/optimization_run.py -m 4"))
    replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
    do(f"sbatch slurms/run_on_scavenge.slurm")
if args.mode == 22:
    protein = args.label
    do(f"python3 ~/opt/compute_phis.py -m 0 {protein}")
    # do("python3 ~/opt/compute_phis.py -m 0 test_protein")
if args.mode == 222:
    proteins = args.label
    evaluate_phis_over_training_set_for_native_structures_Wei(proteins, "phi_list.txt",
                decoy_method='shuffle', max_decoys=1e+10, tm_only=False, num_processors=1, separateDecoysNum=separateDecoysNum)
if args.mode == 3:
    with open(f"slurms/run_on_scavenge.slurm", "w") as out:
        out.write(scavenge_slurm.format(f"python3 ~/opt/optimization_run.py -m 4"))
    replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
    do(f"sbatch slurms/run_on_scavenge.slurm")
if args.mode == 4:
    # complete_proteins = "iter0.txt"
    complete_proteins = "protein_list"
    # A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_parallel(complete_proteins, "phi_list.txt", decoy_method='shuffle',
    #                                 num_decoys=1000, noise_filtering=True, jackhmmer=False, subset=None, read=2)
    A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='shuffle',
                                    num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False)
# if args.mode == 44:
#     with open(f"slurms/run_on_scavenge.slurm", "w") as out:
#         out.write(base_slurm.format(f"python3 ~/opt/gg_server.py -d jun10 -m 4"))
#     replace(f"slurms/run_on_scavenge.slurm", "#SBATCH --mem-per-cpu=1G", "#SBATCH --mem-per-cpu=60G")
#     do(f"sbatch slurms/run_on_scavenge.slurm")