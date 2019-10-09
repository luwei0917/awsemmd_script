import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
import imp
import numpy as np
import fileinput
from itertools import product
import pandas as pd
from scipy.interpolate import griddata
from scipy.interpolate import interp2d
import seaborn as sns
from os import listdir

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import griddata
import matplotlib as mpl
# sys.path.insert(0,'..')
# from notebookFunctions import *
# from .. import notebookFunctions

from Bio.PDB.PDBParser import PDBParser
from pyCodeLib import *



parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    automatic copy the template file, \
    run simulation and analysis")

parser.add_argument("protein", help="The name of the protein")
parser.add_argument("-l", "--label", type=str, default="/Users/weilu/opt/ff_contact/1r69/")
parser.add_argument("-f", "--familyFold", action="store_true", default=False)
parser.add_argument("-m", "--membrane", action="store_true", default=False)
args = parser.parse_args()

code = {"GLY" : "G", "ALA" : "A", "LEU" : "L", "ILE" : "I",
        "ARG" : "R", "LYS" : "K", "MET" : "M", "CYS" : "C",
        "TYR" : "Y", "THR" : "T", "PRO" : "P", "SER" : "S",
        "TRP" : "W", "ASP" : "D", "GLU" : "E", "ASN" : "N",
        "GLN" : "Q", "PHE" : "F", "HIS" : "H", "VAL" : "V",
        "M3L" : "K", "MSE" : "M", "CAS" : "C"}
gamma_se_map_1_letter = {   'A': 0,  'R': 1,  'N': 2,  'D': 3,  'C': 4,
                            'Q': 5,  'E': 6,  'G': 7,  'H': 8,  'I': 9,
                            'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14,
                            'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19}


def read_gamma(gammaFile):
    data = np.loadtxt(gammaFile)
    gamma_direct = data[:210]
    gamma_mediated = data[210:]
    return gamma_direct, gamma_mediated
from pathlib import Path
import platform

if(platform.system() == 'Darwin'):
    pre = str(Path.home()) + "/openmmawsem/parameters/"
elif(platform.system() == 'Linux'):
    pre = "/projects/pw8/wl45/openawsem/parameters/"
else:
    print("system unkown")


# pre = "/Users/weilu/Research/server/sep_2019/saved_gammas/no_mix"
pre = "/Users/weilu/Research/server/sep_2019/saved_gammas/Sep28_trial_5_cutoff500_impose_Aprime_constraint"
if args.membrane:
    gamma_direct, gamma_mediated = read_gamma(f"{pre}/membrane_gamma.dat")
else:
    gamma_direct, gamma_mediated = read_gamma(f"{pre}/gamma.dat")
burial_gamma = np.loadtxt(f"{pre}/burial_gamma.dat").T

nwell = 1
gamma_ijm = np.zeros((nwell, 20, 20))
water_gamma_ijm = np.zeros((nwell, 20, 20))
protein_gamma_ijm = np.zeros((nwell, 20, 20))
m = 0
count = 0
for i in range(20):
    for j in range(i, 20):
        gamma_ijm[m][i][j] = gamma_direct[count][0]
        gamma_ijm[m][j][i] = gamma_direct[count][0]
        count += 1
count = 0
for i in range(20):
    for j in range(i, 20):
        water_gamma_ijm[m][i][j] = gamma_mediated[count][1]
        water_gamma_ijm[m][j][i] = gamma_mediated[count][1]
        count += 1
count = 0
for i in range(20):
    for j in range(i, 20):
        protein_gamma_ijm[m][i][j] = gamma_mediated[count][0]
        protein_gamma_ijm[m][j][i] = gamma_mediated[count][0]
        count += 1

# def compute_chi(data):
#     ca_all = data.query("type == 'CA'")[["x","y","z"]].values
#     cb_all = data.query("type == 'CB'")[["x","y","z"]].values
#     c_all = data.query("type == 'C'")[["x","y","z"]].values
#     n_all = data.query("type == 'N'")[["x","y","z"]].values
#     print(len(ca_all), len(cb_all), len(c_all), len(n_all))
#     energy = 0
#     for i in range(len(n_all)):
#         ca = ca_all[i]
#         cb = cb_all[i]
#         c = c_all[i]
#         n = n_all[i]
#         chi0 = -0.83
#         k_chi = 20*4.184
#         r_ca_cb = cb-ca
#         r_c_ca = ca-c
#         r_ca_n = n-ca
#         norm_r_ca_cb = np.sum(r_ca_cb**2)**0.5
#         norm_r_c_ca = np.sum(r_c_ca**2)**0.5
#         norm_r_ca_n = np.sum(r_ca_n**2)**0.5
#         a = np.cross(-r_c_ca,r_ca_n)/norm_r_c_ca/norm_r_ca_n
#         chi = np.dot(a,r_ca_cb)/norm_r_ca_cb
#         dchi = chi - chi0
#         energy += k_chi*dchi*dchi
#     return energy

def compute_chi(data):
    res_list = get_res_list(structure)
    energy = 0
    for res1globalindex, res1 in enumerate(res_list):
        ca = ca_all[i]
        cb = cb_all[i]
        c = c_all[i]
        n = n_all[i]
        chi0 = -0.83
        k_chi = 20*4.184
        r_ca_cb = cb-ca
        r_c_ca = ca-c
        r_ca_n = n-ca
        norm_r_ca_cb = np.sum(r_ca_cb**2)**0.5
        norm_r_c_ca = np.sum(r_c_ca**2)**0.5
        norm_r_ca_n = np.sum(r_ca_n**2)**0.5
        a = np.cross(-r_c_ca,r_ca_n)/norm_r_c_ca/norm_r_ca_n
        chi = np.dot(a,r_ca_cb)/norm_r_ca_cb
        dchi = chi - chi0
        energy += k_chi*dchi*dchi
    return energy


def compute_debye_huckel(data):
    res_list = get_res_list(structure)
    k_dh = 4.15
    debye_huckel = 0
    k_screening = 1.0
    screening_length = 10  # (in the unit of A)
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        for res2globalindex, res2 in enumerate(res_list):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
            if res2globalindex > res1globalindex:
                res1Name = three_to_one(res1.get_resname())
                res2Name = three_to_one(res2.get_resname())
                charge_1 = 0
                charge_2 = 0
                if res1Name == "R" or res1Name == "K":
                    charge_1 = 1
                if res1Name == "D" or res1Name == "E":
                    charge_1 = -1
                if res2Name == "R" or res2Name == "K":
                    charge_2 = 1
                if res2Name == "D" or res2Name == "E":
                    charge_2 = -1
                if charge_1 * charge_2 != 0:
                    r = get_interaction_distance(res1, res2)
                    debye_huckel += charge_1*charge_2/r*math.exp(-k_screening*r/screening_length)
    debye_huckel *= k_dh
    return debye_huckel

input_pdb_filename = "/Users/weilu/Research/server_backup/jan_2019/compute_energy/12asA00"
def compute_mediated(structure, kappa=5.0):
    res_list = get_res_list(structure)
    neighbor_list = get_neighbor_list(structure)
    sequence = get_sequence_from_structure(structure)
    cb_density = calculate_cb_density(res_list, neighbor_list)
    r_min = 6.5
    r_max = 9.5
    # kappa = 5.0
    min_seq_sep = 10
    density_threshold = 2.6
    density_kappa = 7.0
    # phi_mediated_contact_well = np.zeros((2, 20,20))
    v_mediated = 0
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        rho_i = cb_density[res1globalindex]
        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max+2.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)
            rho_j = cb_density[res2globalindex]
            # if 1 is wrong. because B20 will interact with A1 twice.
            if_1 = res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex)
            # if 2 is the correct one. should be used.
            if_2 = res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex)
            if_3 = res2globalindex - res1globalindex >= min_seq_sep
            # if if_1 and not if_2:
            #     print("true 1, false 2",res2globalindex, res1globalindex, res2chain, res1chain, res2index, res1index)
            # if not if_1 and if_2:
            #     print("false 1, true 2", res2globalindex, res1globalindex, res2chain, res1chain, res2index, res1index)
            # if if_3 and not if_1:
            #     print("true 3, false 1",res2globalindex, res1globalindex, res2chain, res1chain, res2index, res1index)
            # if not if_3 and if_1:
            #     print("false 3, true 1, if3 stricker than if1,",res2globalindex, res1globalindex, res2chain, res1chain, res2index, res1index)
            if if_3 and not if_2:
                print("true 3, false 2",res2globalindex, res1globalindex, res2chain, res1chain, res2index, res1index)
            if not if_3 and if_2:
                print("false 3, true 2, if3 stricker than if2,",res2globalindex, res1globalindex, res2chain, res1chain, res2index, res1index)
            # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
            # if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain):
            # if res2globalindex - res1globalindex >= min_seq_sep:
            if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                res1type = get_res_type(res_list, res1)
                res2type = get_res_type(res_list, res2)
                rij = get_interaction_distance(res1, res2)
                _pij_protein = prot_water_switchFunc_sigmaProt(
                    rho_i, rho_j, density_threshold, density_kappa) * protein_gamma_ijm[0][res1type][res2type]
                _pij_water = prot_water_switchFunc_sigmaWater(
                    rho_i, rho_j, density_threshold, density_kappa) * water_gamma_ijm[0][res1type][res2type]
                v_mediated += (_pij_protein + _pij_water) * interaction_well(rij, r_min, r_max, kappa)
    return v_mediated

input_pdb_filename = "/Users/weilu/Research/server_backup/jan_2019/compute_energy/12asA00"
def compute_direct(structure, kappa=5.0):
    res_list = get_res_list(structure)
    neighbor_list = get_neighbor_list(structure)
    sequence = get_sequence_from_structure(structure)
    r_min = 4.5
    r_max = 6.5
    # kappa = 5
    min_seq_sep = 10
    # phi_pairwise_contact_well = np.zeros((20,20))
    v_direct = 0
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        # print(get_interaction_atom(res1).get_vector()[2], type(get_interaction_atom(res1).get_vector()[2]))
        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max+2.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)

            # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
            # if res2globalindex - res1globalindex >= min_seq_sep:
            if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                # print(i)
                res1type = get_res_type(res_list, res1)
                res2type = get_res_type(res_list, res2)
                rij = get_interaction_distance(res1, res2)
                gamma = gamma_ijm[0][res1type][res2type]
    #             phi_pairwise_contact_well[res1type][res2type] += interaction_well(rij, r_min, r_max, kappa)
                v_direct += gamma * interaction_well(rij, r_min, r_max, kappa)
    return v_direct

def compute_burial(structure, kappa=4.0):
    res_list = get_res_list(structure)
    neighbor_list = get_neighbor_list(structure)
    sequence = get_sequence_from_structure(structure)
    cb_density = calculate_cb_density(res_list, neighbor_list)
    rho_table = [[0.0, 3.0], [3.0, 6.0], [6.0, 9.0]]
    v_burial = 0
    for i in range(3):
        for res1globalindex, res1 in enumerate(res_list):
            res1index = get_local_index(res1)
            res1chain = get_chain(res1)
            res1type = get_res_type(res_list, res1)
            res1density = cb_density[res1globalindex]
            # print res1globalindex, res1index, res1chain, res1type, res1density
            v_burial += burial_gamma[i][res1type] * interaction_well(res1density, rho_table[i][0], rho_table[i][1], kappa)
    return v_burial


input_pdb_filename = "/Users/weilu/Research/server_backup/jan_2019/compute_energy/12asA00.pdb"
def compute_direct_2(input_pdb_filename):
    _all = []
    seq = ""
    p=PDBParser()
    structure=p.get_structure("x", input_pdb_filename)
    for model in structure:
        for chain in model:
            for residue in chain:
                seq += code[residue.resname]
                if residue.resname == "GLY":
                    x,y,z = residue["CA"].get_coord()
                else:
                    x,y,z = residue["CB"].get_coord()
                _all.append([x,y,z])
    v_direct = 0
    data = np.array(_all)
    n = len(data)
    for i in range(n):
        x, y, z = data[i]
        ai = gamma_se_map_1_letter[seq[i]]
        for j in range(i+10, n):
            xj, yj, zj = data[j]
            aj = gamma_se_map_1_letter[seq[j]]
            r = ((x-xj)**2 + (y-yj)**2 + (z-zj)**2)**0.5
            gamma = gamma_ijm[0][ai][aj]
    #         gamma = 1
            v_direct += gamma * interaction_well(r, 4.5, 6.5, 5)
    #         v_direct += 1
    return v_direct


def compute_mediated_multiDensity(structure, kappa=5.0):
    res_list = get_res_list(structure)
    neighbor_list = get_neighbor_list(structure)
    sequence = get_sequence_from_structure(structure)
    cb_density = calculate_cb_density(res_list, neighbor_list)
    weight_density = calculate_cb_weight_density(res_list, neighbor_list)
    r_min = 6.5
    r_max = 9.5
    # kappa = 5.0
    min_seq_sep = 10
    density_threshold = 2.6
    weight_density_threshold = 3.0
    density_kappa = 7.0
    # phi_mediated_contact_well = np.zeros((2, 20,20))
    v_mediated = 0
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        rho_i = cb_density[res1globalindex]
        rho_i_weight = weight_density[res1globalindex]
        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max+2.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)
            rho_j = cb_density[res2globalindex]
            rho_j_weight = weight_density[res2globalindex]
            # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
            if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                res1type = get_res_type(res_list, res1)
                res2type = get_res_type(res_list, res2)
                protein_gamma = protein_gamma_ijm[0][res1type][res2type]
                water_gamma = water_gamma_ijm[0][res1type][res2type]

                rij = get_interaction_distance(res1, res2)
                _pij_protein = prot_water_switchFunc_sigmaProt(
                    rho_i, rho_j, density_threshold, density_kappa) * protein_gamma
                _pij_water = prot_water_switchFunc_sigmaWater(
                    rho_i, rho_j, density_threshold, density_kappa) * water_gamma
                v_mediated += (_pij_protein + _pij_water) * interaction_well(rij, r_min, r_max, kappa)

                heavy_gamma = protein_gamma_ijm[0][res1type][res2type]
                light_gamma = water_gamma_ijm[0][res1type][res2type]
                _pij_heavy = prot_water_switchFunc_sigmaProt(
                    rho_i_weight, rho_j_weight, weight_density_threshold, density_kappa) * heavy_gamma
                _pij_light = prot_water_switchFunc_sigmaWater(
                    rho_i_weight, rho_j_weight, weight_density_threshold, density_kappa) * light_gamma
                v_mediated += (_pij_heavy + _pij_light) * interaction_well(rij, r_min, r_max, kappa)
    return v_mediated

def compute_burial_multiDensity(structure, kappa=4.0):
    res_list = get_res_list(structure)
    neighbor_list = get_neighbor_list(structure)
    sequence = get_sequence_from_structure(structure)
    cb_density = calculate_cb_density(res_list, neighbor_list)
    weight_density = calculate_cb_weight_density(res_list, neighbor_list)
    rho_table = [[0.0, 3.0], [3.0, 6.0], [6.0, 9.0]]
    weight_rho_table = [[0.0, 4.0], [4.0, 8.0], [8.0, 28.0]]
    v_burial = 0
    for i in range(3):
        for res1globalindex, res1 in enumerate(res_list):
            res1index = get_local_index(res1)
            res1chain = get_chain(res1)
            res1type = get_res_type(res_list, res1)
            res1density = cb_density[res1globalindex]
            res1weight = weight_density[res1globalindex]
            # print res1globalindex, res1index, res1chain, res1type, res1density
            v_burial += burial_gamma[i][res1type] * interaction_well(res1density, rho_table[i][0], rho_table[i][1], kappa)
            v_burial += burial_gamma[i][res1type] * interaction_well(res1weight, weight_rho_table[i][0], weight_rho_table[i][1], kappa)
    return v_burial

def compute_direct_family_fold(structure, f_direct, kappa=5.0):
    res_list = get_res_list(structure)
    neighbor_list = get_neighbor_list(structure)
    sequence = get_sequence_from_structure(structure)
    r_min = 4.5
    r_max = 6.5
    # kappa = 5
    min_seq_sep = 10
    # phi_pairwise_contact_well = np.zeros((20,20))
    v_direct = 0
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        # print(get_interaction_atom(res1).get_vector()[2], type(get_interaction_atom(res1).get_vector()[2]))
        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max+2.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)

            # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
            if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                res1type = get_res_type(res_list, res1)
                res2type = get_res_type(res_list, res2)
                rij = get_interaction_distance(res1, res2)
                # gamma = gamma_ijm[0][res1type][res2type]
                gamma = f_direct[res1globalindex][res2globalindex]
    #             phi_pairwise_contact_well[res1type][res2type] += interaction_well(rij, r_min, r_max, kappa)
                v_direct += gamma * interaction_well(rij, r_min, r_max, kappa)
    return v_direct


def compute_mediated_family_fold(structure, f_water, f_protein, kappa=5.0):
    res_list = get_res_list(structure)
    neighbor_list = get_neighbor_list(structure)
    sequence = get_sequence_from_structure(structure)
    cb_density = calculate_cb_density(res_list, neighbor_list)
    r_min = 6.5
    r_max = 9.5
    # kappa = 5.0
    min_seq_sep = 10
    density_threshold = 2.6
    density_kappa = 7.0
    # phi_mediated_contact_well = np.zeros((2, 20,20))
    v_mediated = 0
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        rho_i = cb_density[res1globalindex]
        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max+2.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)
            rho_j = cb_density[res2globalindex]
            # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
            if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                res1type = get_res_type(res_list, res1)
                res2type = get_res_type(res_list, res2)
                rij = get_interaction_distance(res1, res2)
                # protein_gamma = protein_gamma_ijm[0][res1type][res2type]
                # water_gamma = water_gamma_ijm[0][res1type][res2type]
                protein_gamma = f_protein[res1globalindex][res2globalindex]
                water_gamma = f_water[res1globalindex][res2globalindex]
                _pij_protein = prot_water_switchFunc_sigmaProt(
                    rho_i, rho_j, density_threshold, density_kappa) * protein_gamma
                _pij_water = prot_water_switchFunc_sigmaWater(
                    rho_i, rho_j, density_threshold, density_kappa) * water_gamma
                v_mediated += (_pij_protein + _pij_water) * interaction_well(rij, r_min, r_max, kappa)
    return v_mediated

def compute_burial_family_fold(structure, f_burial, kappa=4.0):
    res_list = get_res_list(structure)
    neighbor_list = get_neighbor_list(structure)
    sequence = get_sequence_from_structure(structure)
    cb_density = calculate_cb_density(res_list, neighbor_list)
    rho_table = [[0.0, 3.0], [3.0, 6.0], [6.0, 9.0]]
    v_burial = 0
    for i in range(3):
        for res1globalindex, res1 in enumerate(res_list):
            res1index = get_local_index(res1)
            res1chain = get_chain(res1)
            res1type = get_res_type(res_list, res1)
            res1density = cb_density[res1globalindex]
            # print res1globalindex, res1index, res1chain, res1type, res1density
            # b_gamma = burial_gamma[i][res1type]
            b_gamma = f_burial[res1globalindex][i]
            v_burial += b_gamma * interaction_well(res1density, rho_table[i][0], rho_table[i][1], kappa)
    return v_burial

def get_pre_and_post(res_list, index):
    n = len(res_list)
    if index == 0:
        return res_list[0], res_list[1]
    elif index == n - 1:
        return res_list[index-1], res_list[index]
    else:
        return res_list[index-1], res_list[index+1]

def compute_direct_multiLetter(structure, kappa=5.0):
    # gamma_ij_multiLetter = np.zeros((4, 4, 20, 20))
    gamma_ij_multiLetter = np.zeros((80, 80))
    for i in range(4):
        for j in range(4):
            for ii in range(20):
                for jj in range(20):
                    gamma_ij_multiLetter[i*20+ii][j*20+jj] = gamma_ijm[0][ii][jj]
    res_list = get_res_list(structure)
    neighbor_list = get_neighbor_list(structure)
    sequence = get_sequence_from_structure(structure)
    r_min = 4.5
    r_max = 6.5
    # kappa = 5
    min_seq_sep = 10
    # phi_pairwise_contact_well = np.zeros((20,20))
    v_direct = 0
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        # print(get_interaction_atom(res1).get_vector()[2], type(get_interaction_atom(res1).get_vector()[2]))
        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max+2.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)

            # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
            if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                res1type = get_res_type(res_list, res1)
                res2type = get_res_type(res_list, res2)

                res1_pre, res1_post = get_pre_and_post(res_list, res1globalindex)
                res2_pre, res2_post = get_pre_and_post(res_list, res2globalindex)
                res1_neighbor_type = get_neighbor_res_type(res1_pre, res1_post)
                res2_neighbor_type = get_neighbor_res_type(res2_pre, res2_post)

                rij = get_interaction_distance(res1, res2)

                gamma = gamma_ij_multiLetter[res1_neighbor_type*20+res1type][res2_neighbor_type*20+res2type]
    #             phi_pairwise_contact_well[res1type][res2type] += interaction_well(rij, r_min, r_max, kappa)
                v_direct += gamma * interaction_well(rij, r_min, r_max, kappa)
    return v_direct

def compute_mediated_multiLetter(structure, kappa=5.0):
    # protein_gamma_ij_multiLetter = np.zeros((4, 4, 20, 20))
    # water_gamma_ij_multiLetter = np.zeros((4, 4, 20, 20))
    # for i in range(4):
    #     for j in range(4):
    #         protein_gamma_ij_multiLetter[i][j] = protein_gamma_ijm[0]
    #         water_gamma_ij_multiLetter[i][j] = water_gamma_ijm[0]
    protein_gamma_ij_multiLetter = np.zeros((80, 80))
    water_gamma_ij_multiLetter = np.zeros((80, 80))
    for i in range(4):
        for j in range(4):
            for ii in range(20):
                for jj in range(20):
                    protein_gamma_ij_multiLetter[i*20+ii][j*20+jj] = protein_gamma_ijm[0][ii][jj]
                    water_gamma_ij_multiLetter[i*20+ii][j*20+jj] = water_gamma_ijm[0][ii][jj]
    res_list = get_res_list(structure)
    neighbor_list = get_neighbor_list(structure)
    sequence = get_sequence_from_structure(structure)
    cb_density = calculate_cb_density(res_list, neighbor_list)
    r_min = 6.5
    r_max = 9.5
    # kappa = 5.0
    min_seq_sep = 10
    density_threshold = 2.6
    density_kappa = 7.0
    # phi_mediated_contact_well = np.zeros((2, 20,20))
    v_mediated = 0
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        rho_i = cb_density[res1globalindex]
        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max+2.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)
            rho_j = cb_density[res2globalindex]
            # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
            if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                res1type = get_res_type(res_list, res1)
                res2type = get_res_type(res_list, res2)
                rij = get_interaction_distance(res1, res2)

                res1_pre, res1_post = get_pre_and_post(res_list, res1globalindex)
                res2_pre, res2_post = get_pre_and_post(res_list, res2globalindex)
                res1_neighbor_type = get_neighbor_res_type(res1_pre, res1_post)
                res2_neighbor_type = get_neighbor_res_type(res2_pre, res2_post)
                gamma_p = protein_gamma_ij_multiLetter[res1_neighbor_type*20+res1type][res2_neighbor_type*20+res2type]
                gamma_w = water_gamma_ij_multiLetter[res1_neighbor_type*20+res1type][res2_neighbor_type*20+res2type]

                _pij_protein = prot_water_switchFunc_sigmaProt(
                    rho_i, rho_j, density_threshold, density_kappa) * gamma_p

                _pij_water = prot_water_switchFunc_sigmaWater(
                    rho_i, rho_j, density_threshold, density_kappa) * gamma_w
                v_mediated += (_pij_protein + _pij_water) * interaction_well(rij, r_min, r_max, kappa)
    return v_mediated

def compute_burial_multiLetter(structure, kappa=4.0):
    burial_gamma_multiLetter = np.zeros((4, 3, 20))
    for i in range(4):
        burial_gamma_multiLetter[i] = burial_gamma
    res_list = get_res_list(structure)
    neighbor_list = get_neighbor_list(structure)
    sequence = get_sequence_from_structure(structure)
    cb_density = calculate_cb_density(res_list, neighbor_list)
    rho_table = [[0.0, 3.0], [3.0, 6.0], [6.0, 9.0]]
    v_burial = 0
    for i in range(3):
        for res1globalindex, res1 in enumerate(res_list):
            res1index = get_local_index(res1)
            res1chain = get_chain(res1)
            res1type = get_res_type(res_list, res1)

            res1_pre, res1_post = get_pre_and_post(res_list, res1globalindex)
            res1_neighbor_type = get_neighbor_res_type(res1_pre, res1_post)

            res1density = cb_density[res1globalindex]
            # print res1globalindex, res1index, res1chain, res1type, res1density
            v_burial += burial_gamma_multiLetter[res1_neighbor_type][i][res1type] * interaction_well(res1density, rho_table[i][0], rho_table[i][1], kappa)
    return v_burial

# def compute_single_helix_orientation(structure):
#     res_list = get_res_list(structure)
#     for res1globalindex, res1 in enumerate(res_list):
#         for res2globalindex, res2 in enumerate(res_list):

def read_beta_parameters():
    ### directly copied from Nick Schafer's
    # os.chdir(parameter_directory)
    in_anti_HB = open("anti_HB", 'r').readlines()
    in_anti_NHB = open("anti_NHB", 'r').readlines()
    in_para_HB = open("para_HB", 'r').readlines()
    in_para_one = open("para_one", 'r').readlines()
    in_anti_one = open("anti_one", 'r').readlines()

    p_par = np.zeros((20))
    p_anti = np.zeros((20))
    p_antihb = np.zeros((20,20,2))
    p_antinhb = np.zeros((20,20,2))
    p_parhb = np.zeros((20,20,2))

    for i in range(20):
        p_par[i] = float(in_para_one[i].strip())
        p_anti[i] = float(in_anti_one[i].strip())
        for j in range(20):
            p_antihb[i][j][0] = float(in_anti_HB[i].strip().split()[j])
            p_antinhb[i][j][0] = float(in_anti_NHB[i].strip().split()[j])
            p_parhb[i][j][0] = float(in_para_HB[i].strip().split()[j])

    for i in range(20):
        for j in range(20):
            p_antihb[i][j][1] = float(in_anti_HB[i+21].strip().split()[j])
            p_antinhb[i][j][1] = float(in_anti_NHB[i+21].strip().split()[j])
            p_parhb[i][j][1] = float(in_para_HB[i+21].strip().split()[j])
    return p_par, p_anti, p_antihb, p_antinhb, p_parhb


pdb = (args.protein).split(".")[0]
structure = parse_pdb(pdb)
e_debye_huckel = compute_debye_huckel(structure)
print("Debye Huckel", e_debye_huckel)
# e_chi = compute_chi(structure)
# print("Chi", e_chi)
e_mediated = compute_mediated(structure)
e_direct = compute_direct(structure)
e_burial = compute_burial(structure)

# print("Mediated, Direct, Mediated+Direct, Burial, Mediated+Direct+Burial")
name_list = "Protein, Mediated, Direct, Mediated+Direct, Burial, Mediated+Direct+Burial".split(",")
out_line = " ".join([f"{i:<20}" for i in name_list])

print(out_line)
# print(e_mediated, e_direct, e_mediated + e_direct, e_burial, e_mediated + e_direct + e_burial)

energy_out_list = [-e_mediated, -e_direct, -(e_mediated+e_direct), -e_burial, -(e_mediated + e_direct + e_burial)]
out_line = " ".join([f"{pdb:<20}"] + ["{0:<20.3f}".format(i) for i in energy_out_list])
print(out_line)

# print(compute_direct_2("8ab.pdb"))
exit()
if args.familyFold:
    pre = args.label
    # pre = "/Users/weilu/Research/server/may_2019/family_fold/ff_contact/1r69/"
    f_direct = np.loadtxt(f"{pre}/direct.dat")
    f_water = np.loadtxt(f"{pre}/water.dat")
    f_protein = np.loadtxt(f"{pre}/protein.dat")
    f_burial = np.loadtxt(f"{pre}/burial.dat")

    e_direct_ff = compute_direct_family_fold(structure, f_direct, kappa=5.0)
    e_mediated_ff = compute_mediated_family_fold(structure, f_water, f_protein)
    e_burial_ff = compute_burial_family_fold(structure, f_burial)
    energy_out_list = [-e_direct_ff, -e_mediated_ff, -e_burial_ff, -(e_direct_ff+e_mediated_ff), -(e_direct_ff+e_mediated_ff+e_burial_ff)]
    out_line = " ".join([f"{pdb:<20}"] + ["{0:<20.3f}".format(i) for i in energy_out_list])
    print(out_line)


e_direct_multiLetter = compute_direct_multiLetter(structure)
print(e_direct_multiLetter)
e_mediated_multiLetter = compute_mediated_multiLetter(structure)
print(e_mediated_multiLetter)
e_burial_multiLetter = compute_burial_multiLetter(structure)
print(e_burial_multiLetter)
# kappa = 10.0
# kappa_list = [5.0, 10.0]
# for kappa in kappa_list:
#     e_mediated = compute_mediated(structure, kappa=kappa)
#     e_direct = compute_direct(structure, kappa=kappa)
#     e_burial = compute_burial(structure)
#     print(f"kappa: {kappa}", e_mediated, e_direct, e_mediated + e_direct)

print("multiDensity")
e_mediated_multiLetter = compute_mediated_multiDensity(structure)
print(e_mediated_multiLetter)
e_burial_multiLetter = compute_burial_multiDensity(structure)
print(e_burial_multiLetter)