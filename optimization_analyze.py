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

parser = argparse.ArgumentParser(description="Compute gammas under the optimization folder")
# parser.add_argument("DatabaseFolder", help="your database folder")
# parser.add_argument("OptimizationFolder", help="your optimization folder")
parser.add_argument("name", help="name of gamma")
parser.add_argument("-c", "--constant", type=float, default=0)
parser.add_argument("--proteinList", type=str, default="protein_list")
parser.add_argument("--gammaFile", type=str, default="/home/wl45/opt/parameters/original_gamma")

# parser.add_argument("-l", "--label", type=str, default="label")
args = parser.parse_args()

# if args.test:
#     do = print
# else:
#     do = os.system
with open('log_optimization_analyze.txt', 'a') as f:
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


# n_decoys = 10
# separateDecoysNum = -1
# template = base_slurm
template = scavenge_slurm
save_gamma_pre = "saved_gammas"

trial_name = args.name
pre = "gammas/"
# pp = f"protein_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0"
# pp = f"protein_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0phi_debye_huckel_well0"
# complete_proteins = "protein_list"
complete_proteins = args.proteinList
phi_list = read_phi_list("phi_list.txt")
training_set = read_column_from_file(complete_proteins, 1)
total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string(
    phi_list, training_set)
# print(full_parameters_string)
pp = complete_proteins + "_" + full_parameters_string
# original gamma is used for compute c, make sure decoy energy is the same for new gamma and original gamma.
# gamma_file_name = "/home/wl45/opt/parameters/original_gamma"
gamma_file_name = args.gammaFile
original_gamma = np.loadtxt(gamma_file_name)
a = list(original_gamma)
# print(total_phis, num_phis)
if total_phis == 691:
    a.append(1)
    original_gamma = np.array(a)
if total_phis == 630:
    original_gamma = original_gamma[:630]

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

def get_filtered_gamma(pre, cutoff, pp):
    # pp = "cath-dataset-nonredundant-S20Clean_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0"
    # pp = "proteins_name_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0"

    A_name = pp + "_A"
    B_name = pp + "_B"
    B_filtered_name = pp + "_B_filtered"
    P_name = pp + "_P"
    Gamma_name = pp + "_gamma"
    Gamma_filtered_name = pp + "_gamma_filtered"
    Lamb_name = pp + "_lamb"
    Lamb_filtered_name = pp + "_lamb_filtered"

    A = np.loadtxt(pre+A_name)
    B = np.loadtxt(pre+B_name)
    B_filtered = np.loadtxt(pre+B_filtered_name, dtype=complex, converters={
                               0: lambda s: complex(s.decode().replace('+-', '-'))})
    Gamma = np.loadtxt(pre+Gamma_name)
    Gamma_filtered = np.loadtxt(pre+Gamma_filtered_name, dtype=complex, converters={
                               0: lambda s: complex(s.decode().replace('+-', '-'))})
    Lamb = np.loadtxt(pre+Lamb_name, dtype=complex, converters={
                               0: lambda s: complex(s.decode().replace('+-', '-'))})
    Lamb_filtered = np.loadtxt(pre+Lamb_filtered_name, dtype=complex, converters={
                               0: lambda s: complex(s.decode().replace('+-', '-'))})

    half_B_name = pp + "_half_B"
    half_B = np.loadtxt(pre+half_B_name)
    other_half_B_name = pp + "_other_half_B"
    other_half_B = np.loadtxt(pre+other_half_B_name)
    std_half_B_name = pp + "_std_half_B"
    std_half_B = np.loadtxt(pre+std_half_B_name)


    # pre = "/Users/weilu/Research/server/april_2019/"
    location = pre + f"../../phis/{pp}_phi_decoy_summary.txt"
    A_prime = np.loadtxt(location)

    lamb, P = np.linalg.eig(B)
    lamb, P = sort_eigenvalues_and_eigenvectors(lamb, P)
    filtered_lamb = np.copy(lamb)
    cutoff_mode = cutoff
    filtered_B_inv, filtered_lamb, P = get_filtered_B_inv_lambda_and_P(filtered_lamb,
                                                                       cutoff_mode, P)
    filtered_gamma = np.dot(filtered_B_inv, A)
    filtered_B = np.linalg.inv(filtered_B_inv)

    return A, A_prime, filtered_gamma, filtered_B_inv






# we want to impose additional contraint so that A' * gamma = constnat.(-562.23)
# cutoff_list = [100, 200, 300, 400, 500, 600]
# cutoff_list += [10, 20, 30, 40, 50, 80]
cutoff_list = list(np.arange(100, total_phis, 100))
cutoff_list += [630, 650, 670, 690]
print("cutoff_list: ", cutoff_list)
do("mkdir -p saved_gammas")
for cutoff_i in cutoff_list:
    A, A_prime, filtered_gamma, filtered_B_inv = get_filtered_gamma(pre, cutoff_i, pp)
    # c = np.dot(A_prime, original_gamma_deybe)
    if args.constant == 0.0:
        c = np.dot(A_prime, original_gamma)
    else:
        c = args.constant
    print("A' gamma = constant:", c)
    B_inv = filtered_B_inv
    lambda_2 = (A_prime.dot(B_inv).dot(A) - c) / (A_prime.dot(B_inv).dot(A_prime))
    gamma_new = B_inv.dot(A-A_prime*lambda_2)
    # impose A'gamma
    # save_gamma_pre = "/Users/weilu/Research/server/sep_2019/saved_gammas/"

    np.savetxt(f"{save_gamma_pre}/{trial_name}_cutoff{cutoff_i}_impose_Aprime_constraint", gamma_new)

    name = f"{save_gamma_pre}/{trial_name}_cutoff{cutoff_i}_impose_Aprime_constraint"
    cmd = f"convert_to_simulation_format.py {name} Oct26_{name}"
    do(cmd)
