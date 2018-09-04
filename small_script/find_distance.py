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
import subprocess
import glob
import re
import numpy as np
# compute cross Q for every pdb pair in one folder
parser = argparse.ArgumentParser(description="Compute cross q")
parser.add_argument("-m", "--mode",
                    type=int, default=-1)
parser.add_argument("-t", "--targetMode",
                    type=int, default=0)
parser.add_argument("protein", default="2xov", help="the name of protein")
args = parser.parse_args()
# if args.mode == 1:
#

if args.protein == "2xov":
    if args.targetMode == 0:
        target_residue_i = 1
        target_residue_j = 181
    elif args.targetMode == 1:
        target_residue_i = 136
        target_residue_j = 181
    elif args.targetMode == 2:
        # Helix 3,4
        target_residue_i = 80
        target_residue_j = 131
    elif args.targetMode == 3:
        # Helix 1,2
        target_residue_i = 1
        target_residue_j = 80
elif args.protein == "tmhc2":
    if args.targetMode == 0:
        target_residue_i = 1
        target_residue_j = 156
nFrame = 0
found = False
if args.mode == -1:
    lammps_file = "dump.lammpstrj"
else:
    lammps_file = f"dump.lammpstrj.{args.mode}"
target_xi_list = []
target_yi_list = []
target_zi_list = []

target_xj_list = []
target_yj_list = []
target_zj_list = []
with open(lammps_file, "r") as lfile:
    for line in lfile:
        l = line.strip()
        if l[:5]=="ITEM:":
            item = l[6:]
        else:
            if item == "TIMESTEP":
                step = int(l)
                atoms = []
                atoms2 = []
                atoms3 = []
                box = []
                A = []
                nFrame = nFrame + 1
            elif item == "NUMBER OF ATOMS":
                    n_atoms = int(l)
                    xyz_count = 0
            elif item[:10] == "BOX BOUNDS":
                # I center x, y not z
                if xyz_count <= 1:
                    xyz_count += 1
                    box.append(l)
                    l = l.split()
                    # A.append([float(l[0]), float(l[1])])
                    l_left = (float(l[0]) - float(l[1]))/2.0
                    l_right = (float(l[1]) - float(l[0]))/2.0
                    A.append([l_left, l_right])
                    # print l_right - l_left
                else:
                    xyz_count = 0
                    box.append(l)
                    l = l.split()
                    A.append([float(l[0]), float(l[1])])
                    # l_left = (float(l[0]) - float(l[1]))/2.0
                    # l_right = (float(l[1]) - float(l[0]))/2.0
                    # A.append([l_left, l_right])
                    # print l_right - l_left
            elif item[:5] == "ATOMS":
                l = l.split()
                i_atom = int(l[0])
                x = float(l[2])
                y = float(l[3])
                z = float(l[4])
                x = (A[0][1] - A[0][0])*x + A[0][0]
                y = (A[1][1] - A[1][0])*y + A[1][0]
                z = (A[2][1] - A[2][0])*z + A[2][0]
                # C alpha distance
                if i_atom == target_residue_i*3-2:
                    target_xi_list.append(x)
                    target_yi_list.append(y)
                    target_zi_list.append(z)
                if i_atom == target_residue_j*3-2:
                    target_xj_list.append(x)
                    target_yj_list.append(y)
                    target_zj_list.append(z)
                # desc = atom_desc[l[1]]
                # atom = Atom(i_atom, atom_type[l[1]], l[1], x, y, z, desc)
                # atoms.append(atom)
xi = np.array(target_xi_list)
yi = np.array(target_yi_list)
zi = np.array(target_zi_list)
xj = np.array(target_xj_list)
yj = np.array(target_yj_list)
zj = np.array(target_zj_list)
dis = np.sqrt((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2)

if args.targetMode == 0:
    extra = ""
elif args.targetMode == 1:
    extra = "_h56"
elif args.targetMode == 2:
    extra = "_h34"
elif args.targetMode == 3:
    extra = "_h12"
if args.mode == -1:
    disName = f"distance{extra}.dat"
else:
    disName = f"distance{extra}_{args.mode}.dat"

with open(disName, "w") as f:
    if args.targetMode == 0:
        f.write("Steps, DisReal\n")
    elif args.targetMode == 1:
        f.write("Steps, Dis_h56\n")
    elif args.targetMode == 2:
        f.write("Steps, Dis_h34\n")
    elif args.targetMode == 3:
        f.write("Steps, Dis_h12\n")
    for i, d in enumerate(dis):
        f.write("{}, {}\n".format(i*4000, d))
