#!/usr/bin/env python3
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
from time import sleep
from Bio.PDB.PDBParser import PDBParser
import numpy as np
from scipy.interpolate import interp1d

def get_r(all_res, i, j):
    res1 = all_res[i]
    res2 = all_res[j]
    try:
        cb1 = res1["CB"]
    except:
        cb1 = res1["CA"]
    try:
        cb2 = res2["CB"]
    except:
        cb2 = res2["CA"]
    r = cb1 - cb2
    return r

def compute_machine_learning_energy(s, index_array, interaction_array, num_of_points=100):
    x = [0.0, 2.0, 3.5, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.25, 12.75, 13.25, 13.75, 14.25, 14.75, 15.25, 15.75, 16.25, 16.75, 17.25, 17.75, 18.25, 18.75, 19.25, 19.75]
    
    all_res = list(s.get_residues())
    r_max = max(x)
    r_min = min(x)
    dr = (r_max-r_min)/(num_of_points-1)

    max_r_index_1 = num_of_points - 2
    energy = 0.0

    for index_pair, interaction in zip(index_array, interaction_array):
        i,j = index_pair
        r = get_r(all_res, i, j)
        # print(i,j,interaction, r)
        r_index_1 = int(min(max_r_index_1, np.floor(r/dr)))
        r_index_2 = r_index_1 + 1
        r_1 = r_min + dr * r_index_1
        r_2 = r_min + dr * r_index_2
        v1 = interaction[r_index_1]
        v2 = interaction[r_index_2]
        e = ((v2-v1)*r+v1*r_2-v2*r_1)/(r_2-r_1)
        energy += e
    return energy


# from run_parameter import *
parser = argparse.ArgumentParser(
    description="compute the energy")

parser.add_argument("protein", help="The name of the protein")
parser.add_argument("--dataFile", default="dist.npz", help="pair info input")
parser.add_argument("--saved_file", default="ml_data.npz", help="saved pair info input")
parser.add_argument("--UseSavedFile", action="store_true", default=False)

args = parser.parse_args()

parser = PDBParser()
pdbFile = args.protein
s = parser.get_structure("X", pdbFile)
all_res = list(s.get_residues())

x = [0.0, 2.0, 3.5, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.25, 12.75, 13.25, 13.75, 14.25, 14.75, 15.25, 15.75, 16.25, 16.75, 17.25, 17.75, 18.25, 18.75, 19.25, 19.75]
num_of_points = 100

UseSavedFile = args.UseSavedFile
saved_file = args.saved_file
dataFile = args.dataFile

if UseSavedFile and os.path.isfile(saved_file):
    data = np.load(saved_file)
    index_array = data["index_array"]
    interaction_array = data["interaction_array"]
else:
    # spline fit
    a = np.load(dataFile)
    distspline = a['distspline']

    n = distspline.shape[0]
    interaction_list = []
    index_list = []
    xnew = np.linspace(min(x), max(x), num=num_of_points, endpoint=True)
    for i in range(n):
        for j in range(i+1, n):
            if np.alltrue(distspline[i][j] == 0):
                continue
            y = distspline[i][j]
            f = interp1d(x, y)
            ynew = f(xnew)
            interaction_list.append(ynew)
            index_list.append([i, j])
    index_array = np.array(index_list)
    interaction_array = np.array(interaction_list)
    np.savez(saved_file, index_array=index_array, interaction_array=interaction_array)

energy = compute_machine_learning_energy(s, index_array, interaction_array)
print(energy)
