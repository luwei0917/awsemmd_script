#!/usr/bin/env python3
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
import imp
import glob
from time import sleep
import fileinput
import numpy as np
import pandas as pd
from small_script.variable_test import variable_test
from small_script.variable_test2 import variable_test2
import subprocess
from small_script.myFunctions import compute_theta_for_each_helix
from small_script.myFunctions import *
from collections import defaultdict


parser = argparse.ArgumentParser(description="This is my playground for current project")
parser.add_argument("projectName", help="name of the folder your simulation will run in")
parser.add_argument("gamma", help="pre name of your gamma")
parser.add_argument("-k", "--kfrag", type=float, default=0.01)
parser.add_argument("-m", "--mode", type=int, default=0)
parser.add_argument("--singleMemory", action="store_true", default=False)
parser.add_argument("-q", "--qbias", action="store_true", default=False)
args = parser.parse_args()


with open('generate_decoys.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')

def do(cmd):
    subprocess.Popen(cmd, shell=True).wait()
# do = os.system
cd = os.chdir



dataset = {"old":"1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "),
            "new":"1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "),
            "test":["t089", "t120", "t251", "top7", "1ubq", "t0766", "t0778", "t0782", "t0792", "t0803", "t0815", "t0833", "t0842", "t0844"]}
dataset["combined"] = dataset["old"] + dataset["new"]
def returnSteps(p):
    if p in "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "):
        steps = 80
    elif p in "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR".split(", "):
        steps = 40
    elif p in ["1MBA", "2FHA"]:
        steps = 30
    return steps

simulation_location = args.projectName
# gamma = args.gamma + "gamma.dat"
# burial = args.gamma + "burial_gamma.dat"
gamma = args.gamma.replace("burial_gamma", "gamma")
burial = args.gamma


if args.mode == 0:
    data = pd.read_csv("/scratch/wl45/may_2019/protein_info.csv", index_col=0)
    k_rg = 1
    pdb_list = dataset["combined"]

    for p in pdb_list:
        name = p.lower()[:4]

        print(name)

        do(f"mkdir -p {simulation_location}/{name}")
        pre = f"{simulation_location}/{name}/{name}"
        # do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
        do(f"cp -r all_simulations/{name}/ {simulation_location}/")
        # change the k frag.
        replace(f"{pre}/fix_backbone_coeff.data", "0.01", args.kfrag)
        if args.singleMemory:
            replace(f"{pre}/fix_backbone_coeff.data", "frags.mem", "single_frags.mem")

        rg = data.query(f"Protein == '{name}'")["Rg"].values[0]
        replace(f"{pre}/{name}_multi.in", "#FIXRG", f"fix rg alpha_carbons spring/rg {k_rg} {rg}")
        # replace(f"{pre}/{name}_multi.in", "minimize", "#minimize")
        if args.qbias:
            replace(f"{pre}/{name}_multi.in", "#FIXBIAS", "fix               qbias alpha_carbons qbias fix_qbias_coeff.data\\nfix_modify        qbias energy no\\nvariable          biasinge equal f_qbias\\n")
            do(f"cp ~/opt/fix_qbias_coeff.data {pre}/")
            qbias = "--bias 1"
        else:
            qbias = ""
        # replace(f"{pre}/fix_backbone_coeff.data", "\[Dssp_Hdrgn\]-", "\[Dssp_Hdrgn\]")
        # replace(f"{pre}/fix_backbone_coeff.data", "\[P_AP\]-", "\[P_AP\]")

        # replace(f"{pre}/fix_backbone_coeff.data", "\[Water\]", "\[Water\]-")
        # replace(f"{pre}/fix_backbone_coeff.data", "\[Burial\]", "\[Burial\]-")

        do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
        do(f"cp {gamma} {pre}/gamma.dat")
        # print(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
        do(f"cp {pre}/burial_gamma.dat {pre}/original_burial_gamma.dat")
        do(f"cp {burial} {pre}/burial_gamma.dat")
    cd(simulation_location)
    for p in pdb_list:
        name = p.lower()[:4]

        # steps = 40
        # steps = 2
        # steps = 80
        steps = returnSteps(p)

        cd(name)
        do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1 {qbias}")
        # do(f"run.py -n 20 {name} --commons 2 -s {steps} --runs 1")
        # do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 2 --start crystal")
        # do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
        cd("..")
    cd("..")


if args.mode == 1:
    d = pd.read_csv("seq_info.csv", index_col=0)
    pdb_list = d.query("length < 150 and index % 2 == 0")["protein"].tolist()

    for p in pdb_list:
        name = p
        print(name)

        do(f"mkdir -p {simulation_location}/{name}")
        pre = f"{simulation_location}/{name}/{name}"
        # do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
        do(f"cp -r all_simulations/{name}/ {simulation_location}/")
        # change the k frag.
        replace(f"{pre}/fix_backbone_coeff.data", "0.01", args.kfrag)
        if args.singleMemory:
            replace(f"{pre}/fix_backbone_coeff.data", "frags.mem", "single_frags.mem")
        # if args.mode == 0:
        #     rg = data.query(f"Protein == '{name}'")["Rg"].values[0]
        #     replace(f"{pre}/{name}_multi.in", "#FIXRG", f"fix rg alpha_carbons spring/rg {k_rg} {rg}")
        # replace(f"{pre}/{name}_multi.in", "minimize", "#minimize")

        # replace(f"{pre}/fix_backbone_coeff.data", "\[Dssp_Hdrgn\]-", "\[Dssp_Hdrgn\]")
        # replace(f"{pre}/fix_backbone_coeff.data", "\[P_AP\]-", "\[P_AP\]")

        # replace(f"{pre}/fix_backbone_coeff.data", "\[Water\]", "\[Water\]-")
        # replace(f"{pre}/fix_backbone_coeff.data", "\[Burial\]", "\[Burial\]-")

        do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
        do(f"cp {gamma} {pre}/gamma.dat")
        # print(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
        do(f"cp {pre}/burial_gamma.dat {pre}/original_burial_gamma.dat")
        do(f"cp {burial} {pre}/burial_gamma.dat")
    cd(simulation_location)
    for p in pdb_list:
        name = p
        # steps = 40
        # steps = 2
        # steps = 80
        steps = 40
        cd(name)
        do(f"run.py -n 10 {name} --commons 2 -s {steps} --runs 1")
        # do(f"run.py -n 20 {name} --commons 2 -s {steps} --runs 1")
        # do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 2 --start crystal")
        # do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
        cd("..")
    cd("..")


# if args.mode == 2:
#     d = pd.read_csv("seq_info.csv", index_col=0)
#     pdb_list = d.query("length < 150 and index % 2 == 0")["protein"].tolist()

#     for p in pdb_list:
#         name = p
#         print(name)
#         do(f"mkdir -p {simulation_location}/{name}")
#         pre = f"{simulation_location}/{name}/{name}"
#         # do(f"cp -r all_simulations/{name}/{name} {simulation_location}/{name}/")
#         do(f"cp -r all_simulations/{name}/ {simulation_location}/")
#         # change the k frag.
#         replace(f"{pre}/fix_backbone_coeff.data", "0.01", args.kfrag)
#         replace(f"{pre}/fix_backbone_coeff.data", "frags.mem", "single_frags.mem")
#         # if args.mode == 0:
#         #     rg = data.query(f"Protein == '{name}'")["Rg"].values[0]
#         #     replace(f"{pre}/{name}_multi.in", "#FIXRG", f"fix rg alpha_carbons spring/rg {k_rg} {rg}")
#         # replace(f"{pre}/{name}_multi.in", "minimize", "#minimize")

#         # replace(f"{pre}/fix_backbone_coeff.data", "\[Dssp_Hdrgn\]-", "\[Dssp_Hdrgn\]")
#         # replace(f"{pre}/fix_backbone_coeff.data", "\[P_AP\]-", "\[P_AP\]")

#         # replace(f"{pre}/fix_backbone_coeff.data", "\[Water\]", "\[Water\]-")
#         # replace(f"{pre}/fix_backbone_coeff.data", "\[Burial\]", "\[Burial\]-")

#         do(f"cp {pre}/gamma.dat {pre}/original_gamma.dat")
#         do(f"cp {gamma} {pre}/gamma.dat")
#         # print(f"cp {gammaSource}/iteration_gamma.dat {pre}/gamma.dat")
#         do(f"cp {pre}/burial_gamma.dat {pre}/original_burial_gamma.dat")
#         do(f"cp {burial} {pre}/burial_gamma.dat")
#     cd(simulation_location)
#     for p in pdb_list:
#         name = p
#         # steps = 40
#         # steps = 2
#         # steps = 80
#         steps = 40
#         cd(name)
#         do(f"run.py -n 10 {name} --commons 2 -s {steps} --runs 1")
#         # do(f"run.py -n 20 {name} --commons 2 -s {steps} --runs 1")
#         # do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 2 --start crystal")
#         # do(f"run.py -n 30 {name} --commons 2 -s {steps} --runs 1")
#         cd("..")
#     cd("..")