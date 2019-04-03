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
gamma_direct, gamma_mediated = read_gamma("/Users/weilu/openmmawsem/parameters/gamma.dat")

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

def compute_chi(data):
    ca_all = data.query("type == 'CA'")[["x","y","z"]].values
    cb_all = data.query("type == 'CB'")[["x","y","z"]].values
    c_all = data.query("type == 'C'")[["x","y","z"]].values
    n_all = data.query("type == 'N'")[["x","y","z"]].values
    print(len(ca_all), len(cb_all), len(c_all), len(n_all))
    energy = 0
    for i in range(len(n_all)):
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
            if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
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
    kappa = 5
    min_seq_sep = 10
    # phi_pairwise_contact_well = np.zeros((20,20))
    v_direct = 0
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max+2.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)
            if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                res1type = get_res_type(res_list, res1)
                res2type = get_res_type(res_list, res2)
                rij = get_interaction_distance(res1, res2)
                gamma = gamma_ijm[0][res1type][res2type]
    #             phi_pairwise_contact_well[res1type][res2type] += interaction_well(rij, r_min, r_max, kappa)
                v_direct += gamma * interaction_well(rij, r_min, r_max, kappa)
    return v_direct

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


def read_beta_parameters():
    ### directly copied from Nick Schafer's
    #os.chdir(parameter_directory)
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
# e_mediated = compute_mediated(structure)
# e_direct = compute_direct(structure)
# print(e_mediated, e_direct, e_mediated + e_direct)
# kappa = 10.0
kappa_list = [5.0, 10.0]
for kappa in kappa_list:
    e_mediated = compute_mediated(structure, kappa=kappa)
    e_direct = compute_direct(structure, kappa=kappa)
    print(f"kappa: {kappa}", e_mediated, e_direct, e_mediated + e_direct)
