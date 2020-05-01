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

def change_gamma_format(gamma_direct, gamma_mediated):
    nwell = 2
    gamma_ijm = np.zeros((nwell, 20, 20))
    water_gamma_ijm = np.zeros((nwell, 20, 20))
    protein_gamma_ijm = np.zeros((nwell, 20, 20))
    m = 0
    count = 0
    for i in range(20):
        for j in range(i, 20):
            gamma_ijm[0][i][j] = gamma_direct[count][0]
            gamma_ijm[0][j][i] = gamma_direct[count][0]
            gamma_ijm[1][i][j] = gamma_direct[count][1]
            gamma_ijm[1][j][i] = gamma_direct[count][1]
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
    return gamma_ijm, water_gamma_ijm, protein_gamma_ijm

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
    min_seq_sep = 10
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        for res2globalindex, res2 in enumerate(res_list):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
            # if res2globalindex > res1globalindex:
            if res2globalindex >= res1globalindex + min_seq_sep:
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
def phosphorylation(res1globalindex, res2globalindex, res1type, res2type, m, phosphorylated_residue_index, phosphorylated_residue_seq):
    # // Four letter classes
    # // 1) SHL: Small Hydrophilic (ALA, GLY, PRO, SER THR) or (A, G, P, S, T) or {0, 7, 14, 15, 16}
    # // 2) AHL: Acidic Hydrophilic (ASN, ASP, GLN, GLU) or (N, D, Q, E) or {2, 3, 5, 6}
    # // 3) BAS: Basic (ARG HIS LYS) or (R, H, K) or {1, 8, 11}
    # // 4) HPB: Hydrophobic (CYS, ILE, LEU, MET, PHE, TRP, TYR, VAL) or (C, I, L, M, F, W, Y, V)  or {4, 9, 10, 12, 13, 17, 18, 19}
    bb_four_letter_map = [1, 3, 2, 2, 4, 2, 2, 1, 3, 4, 4, 3, 4, 4, 1, 1, 1, 4, 4, 4]
    k_hypercharge = 1
    if (res1globalindex+1) in phosphorylated_residue_index:
        # print(res1globalindex, res2globalindex, k_hypercharge)
        idx = phosphorylated_residue_index.index(res1globalindex+1)

        if bb_four_letter_map[res2type] == 1:
            k_hypercharge = m
        elif bb_four_letter_map[res2type] == 2 or bb_four_letter_map[res2type] == 3:
            k_hypercharge = m*m
        else:
            k_hypercharge = 1
        res1type = res_type_map[phosphorylated_residue_seq[idx]]
    if (res2globalindex+1) in phosphorylated_residue_index:
        # print(res1globalindex, res2globalindex, k_hypercharge)
        idx = phosphorylated_residue_index.index(res2globalindex+1)
        if bb_four_letter_map[res1type] == 1:
            k_hypercharge = m
        elif bb_four_letter_map[res1type] == 2 or bb_four_letter_map[res1type] == 3:
            k_hypercharge = m*m
        else:
            k_hypercharge = 1
        res2type = res_type_map[phosphorylated_residue_seq[idx]]

    return k_hypercharge, res1type, res2type

def compute_mediated(structure, protein_gamma_ijm, water_gamma_ijm, kappa=5.0, hasPhosphorylation=False, fixWellCenter=True):
    if hasPhosphorylation:
        import configparser
        config = configparser.ConfigParser()
        config.read("phosphorylation.dat")
        m = eval(config['phosphorylation']['m'])
        phosphorylated_residue_index = eval(config['phosphorylation']['phosphorylated_residue_index'])
        phosphorylated_residue_seq = eval(config['phosphorylation']['phosphorylated_residue_seq'])
        # print(m, phosphorylated_residue_index, phosphorylated_residue_seq)
        # print(res_type_map['E'])
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
    if not fixWellCenter:
        a = pd.read_csv("/Users/weilu/opt/parameters/side_chain/cbd_cbd_real_contact_symmetric.csv")
        cb_density = calculate_cb_density_wellCenter(res_list, neighbor_list, a)
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
                res1type_old = res1type
                res2type_old = res2type
                if hasPhosphorylation:
                    k_hypercharge, res1type, res2type = phosphorylation(res1globalindex, res2globalindex, res1type, res2type, m, phosphorylated_residue_index, phosphorylated_residue_seq)
                else:
                    k_hypercharge = 1

                protein_gamma = protein_gamma_ijm[0][res1type][res2type]*k_hypercharge
                water_gamma = water_gamma_ijm[0][res1type][res2type]*k_hypercharge
                if k_hypercharge != 1:
                    print(res1globalindex, res2globalindex, res1type_old, res2type_old, res1type, res2type, protein_gamma_ijm[0][res1type_old][res2type_old], water_gamma_ijm[0][res1type_old][res2type_old], protein_gamma, water_gamma, k_hypercharge)

                _pij_protein = prot_water_switchFunc_sigmaProt(
                    rho_i, rho_j, density_threshold, density_kappa) * protein_gamma
                _pij_water = prot_water_switchFunc_sigmaWater(
                    rho_i, rho_j, density_threshold, density_kappa) * water_gamma

                if not fixWellCenter:
                    res1_name = res1.get_resname()
                    res2_name = res2.get_resname()
                    if res1_name == "GLY" or res2_name == "GLY":
                        r_min_res1_res2 = 6.5
                        r_max_res1_res2 = 9.5
                    else:
                        b = a.query(f"ResName1=='{res1_name}' and ResName2=='{res2_name}'")
                        if len(b) == 0:
                            b = a.query(f"ResName1=='{res2_name}' and ResName2=='{res1_name}'")
                        try:
                            r_min_res1_res2 = float(b["r_max"]) + 1.5
                            r_max_res1_res2 = float(b["r_max"]) + 4.5
                        except:
                            print(b)
                    # r_min_res1_res2 = 6.5
                    # r_max_res1_res2 = 9.5
                    v_mediated += (_pij_protein + _pij_water) * interaction_well(rij, r_min_res1_res2, r_max_res1_res2, kappa)
                else:
                    v_mediated += (_pij_protein + _pij_water) * interaction_well(rij, r_min, r_max, kappa)
    return v_mediated



input_pdb_filename = "/Users/weilu/Research/server_backup/jan_2019/compute_energy/12asA00"
def compute_direct(structure, gamma_ijm, kappa=5.0, hasPhosphorylation=False, r_min=2.5, fixWellCenter=True, environment=False):
    if hasPhosphorylation:
        import configparser
        config = configparser.ConfigParser()
        config.read("phosphorylation.dat")
        m = eval(config['phosphorylation']['m'])
        phosphorylated_residue_index = eval(config['phosphorylation']['phosphorylated_residue_index'])
        phosphorylated_residue_seq = eval(config['phosphorylation']['phosphorylated_residue_seq'])
        # print(m, phosphorylated_residue_index, phosphorylated_residue_seq)
        # print(res_type_map['E'])

    res_list = get_res_list(structure)
    neighbor_list = get_neighbor_list(structure)
    sequence = get_sequence_from_structure(structure)

    if environment:
        isH = {}
        isP = {}
        for i in range(20):
            isH[dindex_to_1[i]] = res_type_map_HP[dindex_to_1[i]]
            isP[dindex_to_1[i]] = 1 - res_type_map_HP[dindex_to_1[i]]
        cbd_info = pd.read_csv("/Users/weilu/opt/parameters/side_chain/cbd_cbd_real_contact_symmetric.csv")
        density_H = calculate_property_density_with_cbd_info(res_list, neighbor_list, isH, cbd_info).round(3)
        density_P = calculate_property_density_with_cbd_info(res_list, neighbor_list, isP, cbd_info).round(3)
        # print(density_H)
        # print(density_P)
        # print(isH, isP)
        density_kappa = 1
        d_HP0 = 0
    # r_min = 4.5
    r_max = 6.5
    # kappa = 5
    min_seq_sep = 10
    # phi_pairwise_contact_well = np.zeros((20,20))
    v_direct = 0
    if not fixWellCenter:
        a = pd.read_csv("/Users/weilu/opt/parameters/side_chain/cbd_cbd_real_contact_symmetric.csv")
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
                if hasPhosphorylation:
                    k_hypercharge, res1type, res2type = phosphorylation(res1globalindex, res2globalindex, res1type, res2type, m, phosphorylated_residue_index, phosphorylated_residue_seq)
                else:
                    k_hypercharge = 1
                gamma = gamma_ijm[0][res1type][res2type] * k_hypercharge
    #             phi_pairwise_contact_well[res1type][res2type] += interaction_well(rij, r_min, r_max, kappa)
                if not fixWellCenter:
                    res1_name = res1.get_resname()
                    res2_name = res2.get_resname()
                    if res1_name == "GLY" or res2_name == "GLY":
                        r_min_res1_res2 = 2.5
                        r_max_res1_res2 = 6.5
                    else:
                        b = a.query(f"ResName1=='{res1_name}' and ResName2=='{res2_name}'")
                        if len(b) == 0:
                            b = a.query(f"ResName1=='{res2_name}' and ResName2=='{res1_name}'")
                        try:
                            r_min_res1_res2 = float(b["r_min"]) - 0.5
                            r_max_res1_res2 = float(b["r_max"]) + 1.5
                        except:
                            print(b)
                    # r_min_res1_res2 = 2.5
                    # r_max_res1_res2 = 6.5
                else:
                    r_max_res1_res2 = r_min
                    r_max_res1_res2 = r_max
                if environment:
                    d_H_i = density_H[res1globalindex]
                    d_P_i = density_P[res1globalindex]
                    d_H_j = density_H[res2globalindex]
                    d_P_j = density_P[res2globalindex]
                    d_H = d_H_i + d_H_j
                    d_P = d_P_i + d_P_j
                    sigma_H = 0.5 * np.tanh(density_kappa * (d_H - d_P - d_HP0)) + 0.5
                    sigma_P = 1 - sigma_H
                    gamma_H = gamma_ijm[0][res1type][res2type]
                    gamma_P = gamma_ijm[1][res1type][res2type]
                    theta = interaction_well(rij, r_min_res1_res2, r_max_res1_res2, kappa)
                    v_direct += (gamma_H * sigma_H + gamma_P * sigma_P) * theta
                else:
                    v_direct += gamma * interaction_well(rij, r_min_res1_res2, r_max_res1_res2, kappa)
    return v_direct


def compute_burial(structure, burial_gamma, kappa=4.0, hasPhosphorylation=False):
    if hasPhosphorylation:
        import configparser
        config = configparser.ConfigParser()
        config.read("phosphorylation.dat")
        m = eval(config['phosphorylation']['m'])
        phosphorylated_residue_index = eval(config['phosphorylation']['phosphorylated_residue_index'])
        phosphorylated_residue_seq = eval(config['phosphorylation']['phosphorylated_residue_seq'])
        print(m, phosphorylated_residue_index, phosphorylated_residue_seq)
        print(res_type_map['E'])
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
            if hasPhosphorylation and (res1globalindex+1) in phosphorylated_residue_index:
                idx = phosphorylated_residue_index.index(res1globalindex+1)
                res1type = res_type_map[phosphorylated_residue_seq[idx]]
            gamma = burial_gamma[i][res1type]
            # print res1globalindex, res1index, res1chain, res1type, res1density
            v_burial += gamma * interaction_well(res1density, rho_table[i][0], rho_table[i][1], kappa)
    return v_burial

def read_hydrophobicity_scale(seq, tableLocation, isNew=False):
    seq_dataFrame = pd.DataFrame({"oneLetterCode":list(seq)})
    # HFscales = pd.read_table("~/opt/small_script/Whole_residue_HFscales.txt")
    # print(f"reading hydrophobicity scale table from {tableLocation}/Whole_residue_HFscales.txt")
    HFscales = pd.read_csv(f"{tableLocation}/Whole_residue_HFscales.txt", sep="\t")
    if not isNew:
        # Octanol Scale
        # new and old difference is at HIS.
        code = {"GLY" : "G", "ALA" : "A", "LEU" : "L", "ILE" : "I",
                "ARG+" : "R", "LYS+" : "K", "MET" : "M", "CYS" : "C",
                "TYR" : "Y", "THR" : "T", "PRO" : "P", "SER" : "S",
                "TRP" : "W", "ASP-" : "D", "GLU-" : "E", "ASN" : "N",
                "GLN" : "Q", "PHE" : "F", "HIS+" : "H", "VAL" : "V",
                "M3L" : "K", "MSE" : "M", "CAS" : "C"}
    else:
        code = {"GLY" : "G", "ALA" : "A", "LEU" : "L", "ILE" : "I",
                "ARG+" : "R", "LYS+" : "K", "MET" : "M", "CYS" : "C",
                "TYR" : "Y", "THR" : "T", "PRO" : "P", "SER" : "S",
                "TRP" : "W", "ASP-" : "D", "GLU-" : "E", "ASN" : "N",
                "GLN" : "Q", "PHE" : "F", "HIS0" : "H", "VAL" : "V",
                "M3L" : "K", "MSE" : "M", "CAS" : "C"}
    HFscales_with_oneLetterCode = HFscales.assign(oneLetterCode=HFscales.AA.str.upper().map(code)).dropna()
    data = seq_dataFrame.merge(HFscales_with_oneLetterCode, on="oneLetterCode", how="left")
    return data


def compute_membrane(structure, kappa=4.0):
    k_membrane = 1
    membrane_center = 0
    k_m = 2
    z_m = 15
    tanh = np.tanh
    res_list = get_res_list(structure)
    neighbor_list = get_neighbor_list(structure)
    # sequence = get_sequence_from_structure(structure)
    seq = [three_to_one(res.get_resname()) for res in res_list]
    sequence = "".join(seq)
    v_membrane = 0
    hydrophobicityScale_list = read_hydrophobicity_scale(sequence, "/Users/weilu/openmmawsem/helperFunctions")["DGwoct"].values
    # print(hydrophobicityScale_list)
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        res1type = get_res_type(res_list, res1)
        z = res1['CA'].get_coord()[-1]
        # print res1globalindex, res1index, res1chain, res1type, res1density
        v_membrane += k_membrane*(0.5*tanh(k_m*((z-membrane_center)+z_m))+0.5*tanh(k_m*(z_m-(z-membrane_center))))*hydrophobicityScale_list[res1globalindex]
    return v_membrane

def compute_positive_inside_rule(structure, kappa=4.0):
    k_membrane = 1
    membrane_center = 0
    k_m = 2
    z_m = 15
    tanh = np.tanh
    res_list = get_res_list(structure)
    neighbor_list = get_neighbor_list(structure)
    # sequence = get_sequence_from_structure(structure)
    seq = [three_to_one(res.get_resname()) for res in res_list]
    sequence = "".join(seq)
    v_membrane = 0
    positive_inside_residue_table = {"G":0, "A":0, "V":0, "C":0, "P":0, "L":0, "I":0, "M":0, "W":0, "F":0,
                                        "S":0, "T":0, "Y":0, "N":0, "Q":0,
                                        "K":-1, "R":-1, "H":0,
                                        "D":0, "E":0}
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        res1type = get_res_type(res_list, res1)
        try:
            z = res1['CB'].get_coord()[-1]
        except:
            z = res1['CA'].get_coord()[-1]
        # print res1globalindex, res1index, res1chain, res1type, res1density
        thickness = 15
        z_m = thickness * positive_inside_residue_table[three_to_one(res1.get_resname())]
        v_membrane += k_membrane*(z-membrane_center-z_m)**2
    v_membrane /= 100
    return v_membrane

input_pdb_filename = "/Users/weilu/Research/server_backup/jan_2019/compute_energy/12asA00.pdb"
def compute_direct_2(input_pdb_filename, gamma_ijm):
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


def compute_mediated_multiDensity(structure, protein_gamma_ijm, kappa=5.0):
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

def compute_burial_multiDensity(structure, burial_gamma, kappa=4.0):
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

def compute_direct_multiLetter(structure, gamma_ijm, kappa=5.0):
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

def compute_mediated_multiLetter(structure, protein_gamma_ijm, kappa=5.0):
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

def compute_burial_multiLetter(structure, burial_gamma_multiLetter, kappa=4.0):
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


def get_pap_gamma_APH(donor_idx, acceptor_idx, chain_i, chain_j, gamma_APH):
    # if chain_i == chain_j and abs(j-i) < 13 or abs(j-i) > 16:
    # if abs(j-i) < 13 or abs(j-i) > 16:
    # if i-j < 13 or i-j > 16:
    # if (donor_idx - acceptor_idx >= 13 and donor_idx - acceptor_idx <= 16) or chain_i != chain_j:
    if (donor_idx - acceptor_idx >= 13 and donor_idx - acceptor_idx <= 16):
        return gamma_APH
    else:
        return 0
def get_pap_gamma_AP(donor_idx, acceptor_idx, chain_i, chain_j, gamma_AP):
    if (donor_idx - acceptor_idx >= 17) or chain_i != chain_j:
        #     if (donor_idx - acceptor_idx >= 17):
        return gamma_AP
    else:
        return 0
def compute_pap1(structure):
    all_res = list(structure.get_residues())
    chains = [res.get_full_id()[2] for res in all_res]
    n = len(all_res)
    eta = 7
    r0 = 8
    e_aph = 0
    e_ap = 0
    for i in range(n-4):
        for j in range(4, n):
            chain_i = chains[i]
            chain_j = chains[j]
            gamma_APH = get_pap_gamma_APH(j, i, chain_i, chain_j, 1)
            gamma_AP = get_pap_gamma_AP(j, i, chain_i, chain_j, 0.4)
            dis = all_res[i]["CA"] - all_res[j]["CA"]
            dis_p4 = all_res[i+4]["CA"] - all_res[j-4]["CA"]
            rho1 = 0.5*(1+np.tanh(eta*(r0-dis)))
            rho2 = 0.5*(1+np.tanh(eta*(r0-dis_p4)))
            e_aph += -gamma_APH * rho1 * rho2
            e_ap += -gamma_AP * rho1 * rho2
        #     if show > 1e-6:
        #         print(i, j, show, rho1, rho2)
        #     if i == 0:
        #         print(i, j, show, rho1, rho2)
        # break
    # print(a_)
    return e_aph, e_ap
def get_pap_gamma_P(donor_idx, acceptor_idx, chain_i, chain_j, gamma_P):
    if (donor_idx - acceptor_idx >= 9) or chain_i != chain_j:
        return gamma_P
    else:
        return 0

def compute_pap2(structure):
    all_res = list(structure.get_residues())
    chains = [res.get_full_id()[2] for res in all_res]
    n = len(all_res)
    eta = 7
    r0 = 8
    e_p = 0
    for i in range(n-4):
        for j in range(n-4):
            chain_i = chains[i]
            chain_j = chains[j]
            gamma_P = get_pap_gamma_AP(j, i, chain_i, chain_j, 0.4)
            dis = all_res[i]["CA"] - all_res[j]["CA"]
            dis_p4 = all_res[i+4]["CA"] - all_res[j+4]["CA"]
            rho1 = 0.5*(1+np.tanh(eta*(r0-dis)))
            rho2 = 0.5*(1+np.tanh(eta*(r0-dis_p4)))
            e_p += -gamma_P * rho1 * rho2
        #     if show > 1e-6:
        #         print(i, j, show, rho1, rho2)
        #     if i == 0:
        #         print(i, j, show, rho1, rho2)
        # break
    # print(a_)
    return e_p

def dis(a, b):
    return ((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)**0.5

def compute_side_chain_energy_for_x(x, means, precisions_chol, log_det, weights):
    n_features = 3
    n_components, _ = means.shape

    mean_dot_precisions_chol = np.zeros((3,3))
    log_prob = np.zeros(3)
    for i in range(n_components):
        mean_dot_precisions_chol[i] = np.dot(means[i], precisions_chol[i])
        y = np.dot(x, precisions_chol[i]) - mean_dot_precisions_chol[i]
        log_prob[i] = np.sum(np.square(y))

    log_gaussian_prob = -.5 * (n_features * np.log(2 * np.pi) + log_prob) + log_det
    c = np.max(log_gaussian_prob + np.log(weights))
    score = np.log(np.sum(np.exp(log_gaussian_prob + np.log(weights) - c))) + c
    kt = 1
    E_side_chain = -score*kt
    # print(E_side_chain)
    return E_side_chain

def read_fasta(fastaFile):
    seq = ""
    with open(fastaFile, "r") as f:
        for line in f:
            if line[0] == ">":
                pass
            else:
                # print(line)
                seq += line.strip()
    return seq

def compute_side_chain_energy(structure, seq):
    E_side_chain_energy = 0
    # parser = PDBParser()
    # pdbFile = "/Users/weilu/Research/server/feb_2020/compare_side_chain_with_and_without/native/256_cbd_submode_7_debug/crystal_structure.pdb"
    # fastaFile = "/Users/weilu/Research/server/feb_2020/compare_side_chain_with_and_without/native/256_cbd_submode_7_debug/crystal_structure.fasta"
    # structure = parser.get_structure("x", pdbFile)
    print(seq)

    means_dic = {}
    precisions_chol_dic = {}
    log_det_dic = {}
    weights_dic = {}
    res_type_list = ['GLY', 'ALA', 'VAL', 'CYS', 'PRO', 'LEU', 'ILE', 'MET', 'TRP', 'PHE', 'SER', 'THR', 'TYR', 'GLN', 'ASN', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU']
    for res_type in res_type_list:
        if res_type == "GLY":
            continue

        means = np.loadtxt(f"/Users/weilu/opt/parameters/side_chain/{res_type}_means.txt")
        precisions_chol = np.loadtxt(f"/Users/weilu/opt/parameters/side_chain/{res_type}_precisions_chol.txt").reshape(3,3,3)
        log_det = np.loadtxt(f"/Users/weilu/opt/parameters/side_chain/{res_type}_log_det.txt")
        weights = np.loadtxt(f"/Users/weilu/opt/parameters/side_chain/{res_type}_weights.txt")
        means_dic[res_type] = means

        precisions_chol_dic[res_type] = precisions_chol
        log_det_dic[res_type] = log_det
        weights_dic[res_type] = weights

    for res in structure.get_residues():
        if res.get_full_id()[1] != 0:
            continue
        # x_com = get_side_chain_center_of_mass(res)
        # resname = res.resname
        resname = one_to_three(seq[res.id[1]-1])
        if resname == "GLY":
            continue
        try:
            n = res["N"].get_coord()
            ca = res["CA"].get_coord()
            c = res["C"].get_coord()
        except:
            continue
        x_com = res["CB"].get_coord()
        x = np.array([dis(x_com, n), dis(x_com, ca), dis(x_com, c)])
        r_ca_com = dis(x_com, ca)
    #     resname = "TYR"
        if resname == "GLY":
            side_chain_energy = 0
        else:
            side_chain_energy = compute_side_chain_energy_for_x(x, means_dic[resname],
                                                                precisions_chol_dic[resname],
                                                                log_det_dic[resname],
                                                                weights_dic[resname])
        if abs(side_chain_energy) > 10:
            print(res.id[1], resname, x_com, x, round(side_chain_energy,3), round(r_ca_com,3))
        # print(res.id[1], resname, x_com, round(side_chain_energy,3), round(r_ca_com,3))
        E_side_chain_energy += side_chain_energy
    return E_side_chain_energy

def get_side_chain_center_of_mass(atoms):
    # ensure complete first
    total = np.array([0., 0., 0.])
    total_mass = 0
    for atom in atoms:
        if atom.get_name() in ["N", "CA", "C", "O", "OXT"]:
            continue
        if atom.element == "H":
            continue
        total += atom.mass * atom.get_coord()
        total_mass += atom.mass
        # print(atom.get_name(), atom.get_coord())
    x_com = total / total_mass
    return x_com

def compute_side_chain_exclude_volume_energy(structure, fileLocation='./cbd_cbd_real_contact_symmetric.csv'):
    gamma_se_map_1_letter = {   'A': 0,  'R': 1,  'N': 2,  'D': 3,  'C': 4,
                                'Q': 5,  'E': 6,  'G': 7,  'H': 8,  'I': 9,
                                'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14,
                                'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19}

    r_min_table = np.zeros((20,20))
    r_max_table = np.zeros((20,20))
    # fileLocation = '/Users/weilu/Research/server/mar_2020/cmd_cmd_exclude_volume/cbd_cbd_real_contact_symmetric.csv'
    df = pd.read_csv(fileLocation)
    for i, line in df.iterrows():
        res1 = line["ResName1"]
        res2 = line["ResName2"]
        r_min_table[gamma_se_map_1_letter[three_to_one(res1)]][gamma_se_map_1_letter[three_to_one(res2)]] = line["r_min"]
        r_min_table[gamma_se_map_1_letter[three_to_one(res2)]][gamma_se_map_1_letter[three_to_one(res1)]] = line["r_min"]
        r_max_table[gamma_se_map_1_letter[three_to_one(res1)]][gamma_se_map_1_letter[three_to_one(res2)]] = line["r_max"]
        r_max_table[gamma_se_map_1_letter[three_to_one(res2)]][gamma_se_map_1_letter[three_to_one(res1)]] = line["r_max"]

    all_res = get_res_list(structure)
    n = len(all_res)
    e = 0
    for i in range(n):
        for j in range(i+1, n):
            res1 = all_res[i]
            res2 = all_res[j]
            resname1 = res1.resname
            resname2 = res2.resname
            if resname1 == "GLY" or resname2 == "GLY":
                continue
            cbd_1 = get_side_chain_center_of_mass(res1.get_atoms())
            cbd_2 = get_side_chain_center_of_mass(res2.get_atoms())
            r = dis(cbd_1, cbd_2)
            r_max = r_max_table[gamma_se_map_1_letter[three_to_one(resname1)]][gamma_se_map_1_letter[three_to_one(resname2)]]
            r_min = r_min_table[gamma_se_map_1_letter[three_to_one(resname1)]][gamma_se_map_1_letter[three_to_one(resname2)]]
            if r_max - r_min < 0.1:
                print(res1, res2, r_max, r_min)
            e += np.heaviside(r_max-r, 0)*((r-r_max)/(r_max-r_min))**2
        print(res1, cbd_1)
    return e
