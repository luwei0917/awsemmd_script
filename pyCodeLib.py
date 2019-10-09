import time
import subprocess
import os
import math
import sys
import functools
import itertools
import numpy as np
from scipy import stats
import random
from functools import partial
import datetime
import pandas as pd


# MYHOME = "/Users/weilu/Research/optimization/"
# MYHOME = "/scratch/wl45/oct_2018/03_week/optimization/"
# MYHOME = "/scratch/wl45/jan_2019/optimization/"
# MYHOME = "optimization/"
# For matplot
import matplotlib.pyplot as plt
from cycler import cycler
# plt.rc('axes', prop_cycle=(cycler('color', [
#        'k', 'r', 'g', 'b', 'y', '#8c564b', '#7f7f7f', '#9467bd', 'c', 'm', '#ff7f0e', '#1f77b4'])))

# # Lable font size;
# plt.rcParams.update({'font.size': 12})

# # Tight lay out
# plt.tight_layout()

# For Biopython
from Bio.PDB import *
from Bio.PDB.Polypeptide import one_to_three, three_to_one


##########

# sys.path.append('/opt/home/xl23/script/python3')


###### Code for for loop #########################
def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step


def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step
#############################################

def say_hello():
    print("Hi you geys")


def new_function(words_to_print):
    print("%s" % words_to_print)


###########################################################################
# This script will help split the given string into letters;
#
# Written by Xingcheng Lin, 03/08/2017;
###########################################################################


def splitString(inputString):

    outputList = list(inputString)

    return outputList


###########################################################################
# Optimization of AWSEM gamma potential;
#
# Borrowed from the codes from Nick Schafer, 08/30/2018
###########################################################################


# structures_directory = "database/dompdb/"
structures_directory = "../database/dompdb/"
phis_directory = "../phis/"
decoys_root_directory = "decoys/"
aligments_root_directory = "alignments/"
tms_directory = "/opt/home/xl23/Working/Levine/jason/optimization/awsem/tms/"
gammas_directory = "gammas/"
sub_gammas_directory = "sub_gamma/"
# These three functions provide a way of calling a function multiple times (that run independently)
# on a certain number of processors so that a new function call starts when a processor becomes available


def call_independent_functions_on_n_processors(function, arguments_lists, num_processors):
    from multiprocessing import Pool
    pool = Pool(int(num_processors))
    results = pool.map(universal_worker, pool_args(function, *arguments_lists))


def universal_worker(input_pair):
    function, args = input_pair
    return function(*args)


def pool_args(function, *args):
    return list(zip(itertools.repeat(function), list(zip(*args))))

##############################################################################################################


res_type_map_letters = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
                        'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

inverse_res_type_map = dict(list(zip(list(range(20)), res_type_map_letters)))

hydrophobicity_letters = ['R', 'K', 'N', 'Q', 'D', 'E', 'H', 'Y',
                          'W', 'S', 'T', 'G', 'P', 'A', 'M', 'C', 'F', 'L', 'V', 'I']

hydrophobicity_map = dict(list(zip(hydrophobicity_letters, list(range(20)))))

res_type_map = {
    'A': 0,
    'C': 4,
    'D': 3,
    'E': 6,
    'F': 13,
    'G': 7,
    'H': 8,
    'I': 9,
    'K': 11,
    'L': 10,
    'M': 12,
    'N': 2,
    'P': 14,
    'Q': 5,
    'R': 1,
    'S': 15,
    'T': 16,
    'V': 19,
    'W': 17,
    'Y': 18
}

res_type_map_HP = {
    'C': 0,
    'M': 0,
    'F': 0,
    'I': 0,
    'L': 0,
    'V': 0,
    'W': 0,
    'Y': 0,
    'A': 1,
    'H': 1,
    'T': 1,
    'G': 1,
    'P': 1,
    'D': 1,
    'E': 1,
    'N': 1,
    'Q': 1,
    'R': 1,
    'K': 1,
    'S': 1
}

def parse_pdb(pdb_id):
    parser = PDBParser()
    return parser.get_structure(pdb_id, "%s.pdb" % pdb_id)


def interaction_well(r, r_min, r_max, kappa):
    return 0.5 * (np.tanh(kappa * (r - r_min)) * np.tanh(kappa * (r_max - r))) + 0.5

def interaction_well_2(r, r_min, r_max, kappa):
    return 0.5*np.tanh(kappa * (r - r_min)) + 0.5*np.tanh(kappa * (r_max - r))
# Switching function dictating the difference between the protein and water-mediated environment


def prot_water_switchFunc_sigmaWater(rho_i, rho_j, rho_0, kappa):
    return 0.25 * (1 - np.tanh(kappa * (rho_i - rho_0))) * (1 - np.tanh(kappa * (rho_j - rho_0)))


def prot_water_switchFunc_sigmaProt(rho_i, rho_j, rho_0, kappa):
    return 1 - prot_water_switchFunc_sigmaWater(rho_i, rho_j, rho_0, kappa)


def get_local_index(residue):
    return residue.get_id()[1]


def get_global_index(residue_list, residue):
    return residue_list.index(residue)


def get_chain(residue):
    return residue.get_parent().get_id()


def get_neighbors_within_radius(neighbor_list, residue, radius):
    return neighbor_list.search(get_interaction_atom(residue).get_coord(), radius, level='R')


def get_interaction_atom(residue):
    try:
        if residue.resname == "GLY":
            res = residue['CA']
            return res
        else:
            res = residue['CB']
            return res
    except:
        # print(residue)
        # print("----------Use CA instead---------------")
        # probably because mutation.
        try:
            res = residue['CA']
            return res
            # raise
        except:
            # print("no CA found, work around is to just use any atom")
            # return list(residue.get_atoms())[0]
            print("need debug", residue)
            raise


def get_interaction_distance(res1, res2):
    return get_interaction_atom(res1) - get_interaction_atom(res2)


def get_res_type(res_list, residue):
    return res_type_map[three_to_one(residue.get_resname())]

def get_neighbor_res_type(res_pre, res_post):
    table = np.zeros((2,2))
    table[0][0] = 0
    table[0][1] = 1
    table[1][0] = 2
    table[1][1] = 3
    r1 = res_type_map_HP[three_to_one(res_pre.get_resname())]
    r2 = res_type_map_HP[three_to_one(res_post.get_resname())]
    return int(table[r1][r2])

def get_res_type_HP(res_list, residue):
    return res_type_map_HP[three_to_one(residue.get_resname())]


def get_res_list(structure, tm_only=False):
    pdb_id = structure.get_id().split('/')[-1]
    res_list = Selection.unfold_entities(structure, 'R')

    # Get all residues from a structure
    res_list = [residue for residue in res_list if not is_hetero(residue)]

    if tm_only:
        tm = read_column_from_file(os.path.join(
            tms_directory, pdb_id + '.tm'), 1)
        res_list = [residue for i, residue in enumerate(
            res_list) if tm[i] == '2']

    return res_list


def get_protein_name(structure):
    return structure.get_id().split('/')[-1].split('.')[0]


def get_atom_list(structure):
    atom_list = Selection.unfold_entities(structure, 'A')  # A for atoms
    return atom_list


def get_sequence_from_structure(structure):
    sequence = ""
    ppb = PPBuilder(radius=10.0)
    for pp in ppb.build_peptides(structure, aa_only=False):
        sequence += '%s\n' % pp.get_sequence()
    return sequence.replace('\n', '')


def get_neighbor_list(structure, tm_only=False):
    protein = get_protein_name(structure)
    res_list = get_res_list(structure)
    atom_list = [a for a in get_atom_list(
        structure) if not is_hetero(a.get_parent())]
    if tm_only:
        tm = read_column_from_file(os.path.join(
            tms_directory, protein + '.tm'), 1)
        atom_list = [a for a in atom_list if tm[get_global_index(
            res_list, a.get_parent())] == '2']

    neighbor_list = NeighborSearch(atom_list)
    return neighbor_list


def get_number_of_lines_in_file(filename):
    try:
        return open(filename, 'r').read().count("\n")
    except IOError:
        return 0


def read_decoy_sequences(sequence_file_name):
    sequences = []
    with open(sequence_file_name, "r") as sequence_file:
        for line in sequence_file:
            line = line.strip()
            sequences.append(line)
    return sequences

def read_decoy_structures(structure_file_name):
    structures = []
    with open(structure_file_name, "r") as structure_file:
        for line in structure_file:
            line = line.strip()
            s = parse_pdb(os.path.join(line))
            structures.append(s)
    return structures

def read_decoy_structures_andQ(structure_file_name):
    if structure_file_name[-3:] == "pkl":
        a = pd.read_pickle(structure_file_name)
        structures = a["structure"].tolist()
        Qs = a["Qw"].tolist()
        return structures, Qs

    structures = []
    Qs = []
    with open(structure_file_name, "r") as structure_file:
        for line in structure_file:
            line, Q = line.strip().split()
            s = parse_pdb(os.path.join(line))
            structures.append(s)
            Qs.append(Q)
    return structures, Qs


def is_hetero(residue):
    if residue.id[0] != ' ':
        return True
    else:
        return False

def plot_contact_well(gammas, invert_sign=True, fix_colorbar=True, inferBound=False,
                        vmin=-0.3, vmax=0.3, fix_confidence_colorbar=True, confidence_vmin=0,
                        confidence_vmax=1.0, plot_confidence=False, confidence_lower=None, confidence_upper=None):
    size = 20
    interaction_matrix = np.zeros((size, size))
    i_content = 0
    for i in range(size):
        for j in range(i, size):
            index1 = hydrophobicity_map[inverse_res_type_map[i]]
            index2 = hydrophobicity_map[inverse_res_type_map[j]]
            interaction_matrix[index1][index2] = gammas[i_content]
            interaction_matrix[index2][index1] = gammas[i_content]
            i_content += 1

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # The minus sign is here to be consistent with the way AWSEM thinks about gammas
    if invert_sign:
        interaction_matrix *= -1

    if inferBound:
        vmin = np.min(interaction_matrix)
        vmax = np.max(interaction_matrix)

    if fix_colorbar:
        cax = ax.pcolor(interaction_matrix, vmin=vmin,
                        vmax=vmax, cmap="bwr")
    else:
        cax = ax.pcolor(interaction_matrix, cmap="RdBu_r")
    fig.colorbar(cax)

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(interaction_matrix.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(interaction_matrix.shape[1]) + 0.5, minor=False)

    ax.set_xticklabels(hydrophobicity_letters)
    ax.set_yticklabels(hydrophobicity_letters)

    if plot_confidence:
        confidence_interval_size = confidence_upper - confidence_lower
        confidence_matrix = np.zeros((size, size))
        i_content = 0
        for i in range(size):
            for j in range(i, size):
                index1 = hydrophobicity_map[inverse_res_type_map[i]]
                index2 = hydrophobicity_map[inverse_res_type_map[j]]
                confidence_matrix[index1][index2] = confidence_interval_size[i_content]
                confidence_matrix[index2][index1] = confidence_interval_size[i_content]
                i_content += 1

        fig = plt.figure()
        ax = fig.add_subplot(111)
        if fix_confidence_colorbar:
            cax = ax.pcolor(confidence_matrix, vmin=confidence_vmin,
                            vmax=confidence_vmax, cmap="RdBu_r")
        else:
            cax = ax.pcolor(confidence_matrix, cmap="RdBu_r")
        fig.colorbar(cax)

        # put the major ticks at the middle of each cell
        ax.set_yticks(np.arange(confidence_matrix.shape[0]) + 0.5, minor=False)
        ax.set_xticks(np.arange(confidence_matrix.shape[1]) + 0.5, minor=False)

        ax.set_xticklabels(hydrophobicity_letters)
        ax.set_yticklabels(hydrophobicity_letters)

    # plt.savefig('direct_contact.pdf')
    plt.show()


def plot_phi_pairwise_contact_well(gammas, invert_sign=True, fix_colorbar=True, vmin=-0.3, vmax=0.3, fix_confidence_colorbar=True, confidence_vmin=0, confidence_vmax=1.0, plot_confidence=False, confidence_lower=None, confidence_upper=None):
    size = 20
    interaction_matrix = np.zeros((size, size))
    i_content = 0
    for i in range(size):
        for j in range(i, size):
            index1 = hydrophobicity_map[inverse_res_type_map[i]]
            index2 = hydrophobicity_map[inverse_res_type_map[j]]
            interaction_matrix[index1][index2] = gammas[i_content]
            interaction_matrix[index2][index1] = gammas[i_content]
            i_content += 1

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # The minus sign is here to be consistent with the way AWSEM thinks about gammas
    if invert_sign:
        interaction_matrix *= -1
    if fix_colorbar:
        cax = ax.pcolor(interaction_matrix, vmin=vmin,
                        vmax=vmax, cmap="bwr")
    else:
        cax = ax.pcolor(interaction_matrix, cmap="RdBu_r")
    fig.colorbar(cax)

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(interaction_matrix.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(interaction_matrix.shape[1]) + 0.5, minor=False)

    ax.set_xticklabels(hydrophobicity_letters)
    ax.set_yticklabels(hydrophobicity_letters)

    if plot_confidence:
        confidence_interval_size = confidence_upper - confidence_lower
        confidence_matrix = np.zeros((size, size))
        i_content = 0
        for i in range(size):
            for j in range(i, size):
                index1 = hydrophobicity_map[inverse_res_type_map[i]]
                index2 = hydrophobicity_map[inverse_res_type_map[j]]
                confidence_matrix[index1][index2] = confidence_interval_size[i_content]
                confidence_matrix[index2][index1] = confidence_interval_size[i_content]
                i_content += 1

        fig = plt.figure()
        ax = fig.add_subplot(111)
        if fix_confidence_colorbar:
            cax = ax.pcolor(confidence_matrix, vmin=confidence_vmin,
                            vmax=confidence_vmax, cmap="RdBu_r")
        else:
            cax = ax.pcolor(confidence_matrix, cmap="RdBu_r")
        fig.colorbar(cax)

        # put the major ticks at the middle of each cell
        ax.set_yticks(np.arange(confidence_matrix.shape[0]) + 0.5, minor=False)
        ax.set_xticks(np.arange(confidence_matrix.shape[1]) + 0.5, minor=False)

        ax.set_xticklabels(hydrophobicity_letters)
        ax.set_yticklabels(hydrophobicity_letters)

    plt.savefig('direct_contact.pdf')
    plt.show()


def plot_phi_protein_mediated_contact_well(gammas, invert_sign=True, fix_colorbar=True, vmin=-0.3, vmax=0.3, fix_confidence_colorbar=True, confidence_vmin=0, confidence_vmax=1.0, plot_confidence=False, confidence_lower=None, confidence_upper=None):
    size = 20
    interaction_matrix = np.zeros((size, size))
    i_content = 0
    for i in range(size):
        for j in range(i, size):
            index1 = hydrophobicity_map[inverse_res_type_map[i]]
            index2 = hydrophobicity_map[inverse_res_type_map[j]]
            interaction_matrix[index1][index2] = gammas[i_content]
            interaction_matrix[index2][index1] = gammas[i_content]
            i_content += 1

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # The minus sign is here to be consistent with the way AWSEM thinks about gammas
    if invert_sign:
        interaction_matrix *= -1
    if fix_colorbar:
        cax = ax.pcolor(interaction_matrix, vmin=vmin,
                        vmax=vmax, cmap="bwr")
    else:
        cax = ax.pcolor(interaction_matrix, cmap="RdBu_r")
    fig.colorbar(cax)

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(interaction_matrix.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(interaction_matrix.shape[1]) + 0.5, minor=False)

    ax.set_xticklabels(hydrophobicity_letters)
    ax.set_yticklabels(hydrophobicity_letters)

    if plot_confidence:
        confidence_interval_size = confidence_upper - confidence_lower
        confidence_matrix = np.zeros((size, size))
        i_content = 0
        for i in range(size):
            for j in range(i, size):
                index1 = hydrophobicity_map[inverse_res_type_map[i]]
                index2 = hydrophobicity_map[inverse_res_type_map[j]]
                confidence_matrix[index1][index2] = confidence_interval_size[i_content]
                confidence_matrix[index2][index1] = confidence_interval_size[i_content]
                i_content += 1

        fig = plt.figure()
        ax = fig.add_subplot(111)
        if fix_confidence_colorbar:
            cax = ax.pcolor(confidence_matrix, vmin=confidence_vmin,
                            vmax=confidence_vmax, cmap="RdBu_r")
        else:
            cax = ax.pcolor(confidence_matrix, cmap="RdBu_r")
        fig.colorbar(cax)

        # put the major ticks at the middle of each cell
        ax.set_yticks(np.arange(confidence_matrix.shape[0]) + 0.5, minor=False)
        ax.set_xticks(np.arange(confidence_matrix.shape[1]) + 0.5, minor=False)

        ax.set_xticklabels(hydrophobicity_letters)
        ax.set_yticklabels(hydrophobicity_letters)

    plt.savefig('protein_mediated.pdf')
    plt.show()


def plot_phi_water_mediated_contact_well(gammas, invert_sign=True, fix_colorbar=True, vmin=-0.3, vmax=0.3, fix_confidence_colorbar=True, confidence_vmin=0, confidence_vmax=1.0, plot_confidence=False, confidence_lower=None, confidence_upper=None):
    size = 20
    interaction_matrix = np.zeros((size, size))
    i_content = 0
    for i in range(size):
        for j in range(i, size):
            index1 = hydrophobicity_map[inverse_res_type_map[i]]
            index2 = hydrophobicity_map[inverse_res_type_map[j]]
            interaction_matrix[index1][index2] = gammas[i_content]
            interaction_matrix[index2][index1] = gammas[i_content]
            i_content += 1

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # The minus sign is here to be consistent with the way AWSEM thinks about gammas
    if invert_sign:
        interaction_matrix *= -1
    if fix_colorbar:
        cax = ax.pcolor(interaction_matrix, vmin=vmin,
                        vmax=vmax, cmap="bwr")
    else:
        cax = ax.pcolor(interaction_matrix, cmap="RdBu_r")
    fig.colorbar(cax)

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(interaction_matrix.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(interaction_matrix.shape[1]) + 0.5, minor=False)

    ax.set_xticklabels(hydrophobicity_letters)
    ax.set_yticklabels(hydrophobicity_letters)

    if plot_confidence:
        confidence_interval_size = confidence_upper - confidence_lower
        confidence_matrix = np.zeros((size, size))
        i_content = 0
        for i in range(size):
            for j in range(i, size):
                index1 = hydrophobicity_map[inverse_res_type_map[i]]
                index2 = hydrophobicity_map[inverse_res_type_map[j]]
                confidence_matrix[index1][index2] = confidence_interval_size[i_content]
                confidence_matrix[index2][index1] = confidence_interval_size[i_content]
                i_content += 1

        fig = plt.figure()
        ax = fig.add_subplot(111)
        if fix_confidence_colorbar:
            cax = ax.pcolor(confidence_matrix, vmin=confidence_vmin,
                            vmax=confidence_vmax, cmap="RdBu_r")
        else:
            cax = ax.pcolor(confidence_matrix, cmap="RdBu_r")
        fig.colorbar(cax)

        # put the major ticks at the middle of each cell
        ax.set_yticks(np.arange(confidence_matrix.shape[0]) + 0.5, minor=False)
        ax.set_xticks(np.arange(confidence_matrix.shape[1]) + 0.5, minor=False)

        ax.set_xticklabels(hydrophobicity_letters)
        ax.set_yticklabels(hydrophobicity_letters)

    plt.savefig('water_mediated.pdf')
    plt.show()


def read_all_gammas(phi_list_file_name, training_set_file, training_decoy_method, gamma_file_name=None, noise_filtering=True, read_confidence=False, bootstrapping_confidence=95, bootstrapping_iterations=1000, read_averaged_gammas=False, read_original_phis=False):
    phi_list = read_phi_list(phi_list_file_name)
    training_set = read_column_from_file(training_set_file, 1)

    # If we need to read in the cases where the decoy structures are explicitly provided, we need to change the name correspondingly;
    if read_original_phis == "decoy" and training_decoy_method == "TCR_modeling":
        total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string_decoy_structures_provided(
            phi_list, training_set)
    else:
        total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string(
            phi_list, training_set)
    if gamma_file_name == None:
        if noise_filtering:
            gamma_file_name = os.path.join(gammas_directory, "%s_%s_gamma_filtered" % (
                training_set_file.split('/')[-1].split('.')[0], full_parameters_string))
        elif read_averaged_gammas:
            gamma_file_name = os.path.join(gammas_directory, "%s_%s_gamma_averaged" % (
                training_set_file.split('/')[-1].split('.')[0], full_parameters_string))
        elif read_original_phis == "native":
            gamma_file_name = os.path.join(phis_directory, "%s_%s_phi_native_summary.txt" % (
                training_set_file.split('/')[-1].split('.')[0], full_parameters_string))
        elif read_original_phis == "decoy":
            gamma_file_name = os.path.join(phis_directory, "%s_%s_phi_decoy_summary.txt" % (
                training_set_file.split('/')[-1].split('.')[0], full_parameters_string))
        else:
            gamma_file_name = os.path.join(gammas_directory, "%s_%s_gamma" % (
                training_set_file.split('/')[-1].split('.')[0], full_parameters_string))
    if read_confidence:
        confidence_lower = np.loadtxt(os.path.join(gammas_directory, "%s_%s_confidence_lower_%d_%d" % (training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string, bootstrapping_confidence, bootstrapping_iterations)))
        confidence_upper = np.loadtxt(os.path.join(gammas_directory, "%s_%s_confidence_upper_%d_%d" % (training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string, bootstrapping_confidence, bootstrapping_iterations)))

    # Need to filter out the complex number if in the "filtered" mode;
    if noise_filtering:
        print(gamma_file_name)
        gamma = np.loadtxt(gamma_file_name, dtype=complex, converters={
                           0: lambda s: complex(s.decode().replace('+-', '-'))})
    else:
        gamma = np.loadtxt(gamma_file_name)

    individual_gammas = []
    individual_confidence_lower = []
    individual_confidence_upper = []
    for i, (phi, parameters) in enumerate(phi_list):
        individual_gammas.append(gamma[0:num_phis[i]])
        np.savetxt(gammas_directory + 'individual_gamma' +
                   str(i) + '.txt', individual_gammas[i], fmt='%1.5f')
        gamma = gamma[num_phis[i]:]
        if read_confidence:
            individual_confidence_lower.append(confidence_lower[0:num_phis[i]])
            confidence_lower = confidence_lower[num_phis[i]:]
            individual_confidence_upper.append(confidence_upper[0:num_phis[i]])
            confidence_upper = confidence_upper[num_phis[i]:]

    if read_confidence:
        return individual_gammas, individual_confidence_lower, individual_confidence_upper
    else:
        return individual_gammas


def plot_all_gammas(phi_list_file_name, individual_gammas, vmin=-0.3, vmax=0.3, invert_sign=False, gammas_to_plot=None, plot_confidence=False, individual_confidence_lower=None, individual_confidence_upper=None):
    phi_list = read_phi_list(phi_list_file_name)
    if gammas_to_plot == None:
        gammas_to_plot = list(range(len(phi_list)))
    for i, (phi, parameters) in enumerate(phi_list):
        if not i in gammas_to_plot:
            continue
        plot_phi = globals()['plot_' + phi]
        print(phi, parameters)
        if plot_confidence:
            plot_phi(phi, parameters, individual_gammas[i], plot_confidence=plot_confidence,
                     confidence_lower=individual_confidence_lower[i], confidence_upper=individual_confidence_upper[i])
        else:
            plot_phi(individual_gammas[i], vmin=vmin,
                     vmax=vmax, invert_sign=invert_sign)

'''
def phi_pairwise_contact_well(res_list_tmonly, res_list_entire, neighbor_list, parameter_list, TCRmodeling=False):

    r_min, r_max, kappa, min_seq_sep = parameter_list
    r_min = float(r_min)
    r_max = float(r_max)
    kappa = float(kappa)
    min_seq_sep = int(min_seq_sep)
    phi_pairwise_contact_well = np.zeros((20, 20))
    for res1globalindex, res1 in enumerate(res_list_entire):

        res1index = get_local_index(res1)
        res1chain = get_chain(res1)

        # For TCR modeling, we only need the sequence in the peptide;
        if TCRmodeling:

            if (res1 in res_list_tmonly):
                for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max + 2.0):
                    res2index = get_local_index(res2)
                    res2chain = get_chain(res2)
                    res2globalindex = get_global_index(res_list_entire, res2)
                    # Here, we strictly consider only between the peptide and the TCR (chain C & D):
                    # Res1 through tm_only, is already in chain E, we only need to control the res2 to
                    # be in chian C or D;
                    if (res2chain == 'C' or res2chain == 'D'):
                        res1type = get_res_type(res_list_entire, res1)
                        res2type = get_res_type(res_list_entire, res2)
                        rij = get_interaction_distance(res1, res2)
                        phi_pairwise_contact_well[res1type][res2type] += interaction_well(
                            rij, r_min, r_max, kappa)
                        if not res1type == res2type:
                            phi_pairwise_contact_well[res2type][res1type] += interaction_well(
                                rij, r_min, r_max, kappa)
            else:
                continue

        else:

            for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max + 2.0):
                res2index = get_local_index(res2)
                res2chain = get_chain(res2)
                res2globalindex = get_global_index(res_list_entire, res2)
                if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                    res1type = get_res_type(res_list_entire, res1)
                    res2type = get_res_type(res_list_entire, res2)
                    rij = get_interaction_distance(res1, res2)
                    phi_pairwise_contact_well[res1type][res2type] += interaction_well(
                        rij, r_min, r_max, kappa)
                    if not res1type == res2type:
                        phi_pairwise_contact_well[res2type][res1type] += interaction_well(
                            rij, r_min, r_max, kappa)

    phis_to_return = []
    for i in range(20):
        for j in range(i, 20):
            phis_to_return.append(phi_pairwise_contact_well[i][j])
    return phis_to_return
'''
def phi_pairwise_contact_well(res_list, neighbor_list, parameter_list):
    r_min, r_max, kappa, min_seq_sep = parameter_list
    r_min = float(r_min)
    r_max = float(r_max)
    kappa = float(kappa)
    min_seq_sep = int(min_seq_sep)
    phi_pairwise_contact_well = np.zeros((20,20))
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max+2.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)
            # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
            if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                res1type = get_res_type(res_list, res1)
                res2type = get_res_type(res_list, res2)
                rij = get_interaction_distance(res1, res2)
                phi_pairwise_contact_well[res1type][res2type] += interaction_well(rij, r_min, r_max, kappa)
                if not res1type == res2type:
                    phi_pairwise_contact_well[res2type][res1type] += interaction_well(rij, r_min, r_max, kappa)

    phis_to_return = []
    for i in range(20):
        for j in range(i, 20):
            phis_to_return.append(phi_pairwise_contact_well[i][j])
    return phis_to_return


def calculate_cb_density(res_list, neighbor_list, min_seq_sep=2):
    num_residues = len(res_list)
    density = np.zeros(num_residues)
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        for res2 in get_neighbors_within_radius(neighbor_list, res1, 9.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)
            if abs(res2index - res1index) >= min_seq_sep or (res1chain != res2chain):
                rij = get_interaction_distance(res1, res2)
                density[res1globalindex] += interaction_well(rij, 4.5, 6.5, 5)
    return density

def calculate_cb_weight_density(res_list, neighbor_list, min_seq_sep=2):
    weight_info = pd.read_csv("~/opt/parameters/amino_acid_side_chain_weight", comment="#", sep="\s+")
    weight_info["normalized_weight"] = weight_info["weight"] /(weight_info["weight"].min())
    num_residues = len(res_list)
    density = np.zeros(num_residues)
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        for res2 in get_neighbors_within_radius(neighbor_list, res1, 9.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)
            if abs(res2index - res1index) >= min_seq_sep or (res1chain != res2chain):
                rij = get_interaction_distance(res1, res2)
                res2type = three_to_one(res2.get_resname())
                try:
                    weight = float(weight_info.query(f"oneLetterCode == '{res2type}'")["normalized_weight"])
                except:
                    print(res2)
                    print(res2type)
                    print(weight_info.query(f"oneLetterCode == '{res2type}'")["normalized_weight"])
                density[res1globalindex] += weight*interaction_well(rij, 4.5, 6.5, 5)
    return density

def phi_burial_well(res_list, neighbor_list, parameter_list):
    kappa = parameter_list[0]
    kappa = float(kappa)

    cb_density = calculate_cb_density(res_list, neighbor_list)

    phi_burial = np.zeros((3, 20))
    rho_table = [[0.0, 3.0], [3.0, 6.0], [6.0, 9.0]]
    for i in range(3):
        for res1globalindex, res1 in enumerate(res_list):
            res1index = get_local_index(res1)
            res1chain = get_chain(res1)
            res1type = get_res_type(res_list, res1)
            res1density = cb_density[res1globalindex]
            # print res1globalindex, res1index, res1chain, res1type, res1density
            burial_i = interaction_well(res1density, rho_table[i][0], rho_table[i][1], kappa)
            phi_burial[i][res1type] += burial_i

    phis_to_return = []
    for i in range(3):
        for j in range(20):
            phis_to_return.append(phi_burial[i][j])
    return phis_to_return

def phi_debye_huckel_well(res_list, neighbor_list, parameter_list):
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
    return [debye_huckel]

def phi_burial_multiDensity_well(res_list, neighbor_list, parameter_list):
    kappa = parameter_list[0]
    kappa = float(kappa)

    cb_density = calculate_cb_density(res_list, neighbor_list)
    weight_density = calculate_cb_weight_density(res_list, neighbor_list)

    phi_burial = np.zeros((3, 20))
    rho_table = [[0.0, 3.0], [3.0, 6.0], [6.0, 9.0]]

    weight_phi_burial = np.zeros((3, 20))
    weight_rho_table = [[0.0, 4.0], [4.0, 8.0], [8.0, 28.0]]
    for i in range(3):
        for res1globalindex, res1 in enumerate(res_list):
            res1index = get_local_index(res1)
            res1chain = get_chain(res1)
            res1type = get_res_type(res_list, res1)
            res1density = cb_density[res1globalindex]
            # print res1globalindex, res1index, res1chain, res1type, res1density
            burial_i = interaction_well(res1density, rho_table[i][0], rho_table[i][1], kappa)
            phi_burial[i][res1type] += burial_i

            weight_res1density = weight_density[res1globalindex]
            burial_i = interaction_well(weight_res1density, weight_rho_table[i][0], weight_rho_table[i][1], kappa)
            weight_phi_burial[i][res1type] += burial_i

    phis_to_return = []
    for i in range(3):
        for j in range(20):
            phis_to_return.append(phi_burial[i][j])
    for i in range(3):
        for j in range(20):
            phis_to_return.append(weight_phi_burial[i][j])
    return phis_to_return

def phi_density_mediated_contact_well(res_list, neighbor_list, parameter_list):
    r_min, r_max, kappa, min_seq_sep, density_threshold, density_kappa = parameter_list
    cb_density = calculate_cb_density(res_list, neighbor_list)
    r_min = float(r_min)
    r_max = float(r_max)
    kappa = float(kappa)
    min_seq_sep = int(min_seq_sep)
    density_threshold = float(density_threshold)
    density_kappa = float(density_kappa)
    phi_mediated_contact_well = np.zeros((2, 20,20))
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
                _pij_protein = prot_water_switchFunc_sigmaProt(
                    rho_i, rho_j, density_threshold, density_kappa) * interaction_well(rij, r_min, r_max, kappa)
                _pij_water = prot_water_switchFunc_sigmaWater(
                    rho_i, rho_j, density_threshold, density_kappa) * interaction_well(rij, r_min, r_max, kappa)
                phi_mediated_contact_well[0][res1type][res2type] += _pij_protein
                phi_mediated_contact_well[1][res1type][res2type] += _pij_water
                if not res1type == res2type:
                    phi_mediated_contact_well[0][res2type][res1type] += _pij_protein
                    phi_mediated_contact_well[1][res2type][res1type] += _pij_water

    phis_to_return = []
    for i in range(2):
        for j in range(20):
            for k in range(j, 20):
                phis_to_return.append(phi_mediated_contact_well[i][j][k])
    return phis_to_return

def phi_contact_hybrid_well(res_list, neighbor_list, parameter_list):
    r_min, r_max, kappa, min_seq_sep, density_threshold, density_kappa = parameter_list
    cb_density = calculate_cb_density(res_list, neighbor_list)
    r_min = float(r_min)
    r_max = float(r_max)
    kappa = float(kappa)
    min_seq_sep = int(min_seq_sep)
    density_threshold = float(density_threshold)
    density_kappa = float(density_kappa)
    phi_mediated_contact_well = np.zeros((2,3,20,20))

    eta_switching = 1
    z_m = 15
    r_min_I = 4.5
    r_max_I = 6.5

    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        rho_i = cb_density[res1globalindex]
        z1 = get_z_position(res1)
        alphaMembrane1 = interaction_well_2(z1, -z_m, z_m, eta_switching)

        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max+2.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)
            rho_j = cb_density[res2globalindex]
            z2 = get_z_position(res2)
            alphaMembrane2 = interaction_well_2(z1, -z_m, z_m, eta_switching)

            # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
            if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                res1type = get_res_type(res_list, res1)
                res2type = get_res_type(res_list, res2)
                rij = get_interaction_distance(res1, res2)

                alpha_water = (1 - alphaMembrane1 * alphaMembrane2)
                alpha_membrane = alphaMembrane1 * alphaMembrane2

                direct = interaction_well(rij, r_min_I, r_max_I, kappa)
                _pij_high_density = prot_water_switchFunc_sigmaProt(
                    rho_i, rho_j, density_threshold, density_kappa) * interaction_well(rij, r_min, r_max, kappa)
                _pij_low_density = prot_water_switchFunc_sigmaWater(
                    rho_i, rho_j, density_threshold, density_kappa) * interaction_well(rij, r_min, r_max, kappa)

                phi_mediated_contact_well[0][0][res1type][res2type] += alpha_water * direct
                phi_mediated_contact_well[0][1][res1type][res2type] += alpha_water * _pij_high_density
                phi_mediated_contact_well[0][2][res1type][res2type] += alpha_water * _pij_low_density
                if not res1type == res2type:
                    phi_mediated_contact_well[0][0][res2type][res1type] += alpha_water * direct
                    phi_mediated_contact_well[0][1][res2type][res1type] += alpha_water * _pij_high_density
                    phi_mediated_contact_well[0][2][res2type][res1type] += alpha_water * _pij_low_density

                phi_mediated_contact_well[1][0][res1type][res2type] += alpha_membrane * direct
                phi_mediated_contact_well[1][1][res1type][res2type] += alpha_membrane * _pij_high_density
                phi_mediated_contact_well[1][2][res1type][res2type] += alpha_membrane * _pij_low_density
                if not res1type == res2type:
                    phi_mediated_contact_well[1][0][res2type][res1type] += alpha_membrane * direct
                    phi_mediated_contact_well[1][1][res2type][res1type] += alpha_membrane * _pij_high_density
                    phi_mediated_contact_well[1][2][res2type][res1type] += alpha_membrane * _pij_low_density

    phis_to_return = []
    for ii in range(2):
        for jj in range(3):
            for i in range(20):
                for j in range(i, 20):
                    phis_to_return.append(phi_mediated_contact_well[ii][jj][i][j])

    kappa = 4.0
    phi_burial = np.zeros((2, 3, 20))
    rho_table = [[0.0, 3.0], [3.0, 6.0], [6.0, 9.0]]
    for i in range(3):
        for res1globalindex, res1 in enumerate(res_list):
            res1index = get_local_index(res1)
            res1chain = get_chain(res1)
            res1type = get_res_type(res_list, res1)
            res1density = cb_density[res1globalindex]
            z1 = get_z_position(res1)
            alphaMembrane1 = interaction_well_2(z1, -z_m, z_m, eta_switching)
            # print res1globalindex, res1index, res1chain, res1type, res1density
            burial_i = interaction_well(res1density, rho_table[i][0], rho_table[i][1], kappa)
            phi_burial[0][i][res1type] += burial_i * (1 - alphaMembrane1)
            phi_burial[1][i][res1type] += burial_i * alphaMembrane1

    for ii in range(2):
        for i in range(3):
            for j in range(20):
                phis_to_return.append(phi_burial[ii][i][j])

    return phis_to_return


def phi_contact_membrane_well(res_list, neighbor_list, parameter_list):
    r_min, r_max, kappa, min_seq_sep, density_threshold, density_kappa = parameter_list
    cb_density = calculate_cb_density(res_list, neighbor_list)
    r_min = float(r_min)
    r_max = float(r_max)
    kappa = float(kappa)
    min_seq_sep = int(min_seq_sep)
    density_threshold = float(density_threshold)
    density_kappa = float(density_kappa)
    phi_mediated_contact_well = np.zeros((2,3,20,20))

    eta_switching = 1
    z_m = 15
    r_min_I = 4.5
    r_max_I = 6.5

    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        rho_i = cb_density[res1globalindex]
        # z1 = get_z_position(res1)
        # alphaMembrane1 = interaction_well_2(z1, -z_m, z_m, eta_switching)

        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max+2.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)
            rho_j = cb_density[res2globalindex]
            # z2 = get_z_position(res2)
            # alphaMembrane2 = interaction_well_2(z1, -z_m, z_m, eta_switching)

            # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
            if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                res1type = get_res_type(res_list, res1)
                res2type = get_res_type(res_list, res2)
                rij = get_interaction_distance(res1, res2)

                # alpha_water = (1 - alphaMembrane1 * alphaMembrane2)
                # alpha_membrane = alphaMembrane1 * alphaMembrane2
                alpha_membrane = 1

                direct = interaction_well(rij, r_min_I, r_max_I, kappa)
                _pij_high_density = prot_water_switchFunc_sigmaProt(
                    rho_i, rho_j, density_threshold, density_kappa) * interaction_well(rij, r_min, r_max, kappa)
                _pij_low_density = prot_water_switchFunc_sigmaWater(
                    rho_i, rho_j, density_threshold, density_kappa) * interaction_well(rij, r_min, r_max, kappa)

                phi_mediated_contact_well[1][0][res1type][res2type] += alpha_membrane * direct
                phi_mediated_contact_well[1][1][res1type][res2type] += alpha_membrane * _pij_high_density
                phi_mediated_contact_well[1][2][res1type][res2type] += alpha_membrane * _pij_low_density
                if not res1type == res2type:
                    phi_mediated_contact_well[1][0][res2type][res1type] += alpha_membrane * direct
                    phi_mediated_contact_well[1][1][res2type][res1type] += alpha_membrane * _pij_high_density
                    phi_mediated_contact_well[1][2][res2type][res1type] += alpha_membrane * _pij_low_density

    phis_to_return = []
    ii = 1
    for jj in range(3):
        for i in range(20):
            for j in range(i, 20):
                phis_to_return.append(phi_mediated_contact_well[ii][jj][i][j])

    kappa = 4.0
    phi_burial = np.zeros((2, 3, 20))
    rho_table = [[0.0, 3.0], [3.0, 6.0], [6.0, 9.0]]
    for i in range(3):
        for res1globalindex, res1 in enumerate(res_list):
            res1index = get_local_index(res1)
            res1chain = get_chain(res1)
            res1type = get_res_type(res_list, res1)
            res1density = cb_density[res1globalindex]
            # z1 = get_z_position(res1)
            # alphaMembrane1 = interaction_well_2(z1, -z_m, z_m, eta_switching)
            # print res1globalindex, res1index, res1chain, res1type, res1density
            burial_i = interaction_well(res1density, rho_table[i][0], rho_table[i][1], kappa)
            # phi_burial[0][i][res1type] += burial_i * (1 - alphaMembrane1)
            phi_burial[1][i][res1type] += burial_i * 1

    ii = 1
    for i in range(3):
        for j in range(20):
            phis_to_return.append(phi_burial[ii][i][j])

    return phis_to_return


def phi_contact_only_membrane_well(res_list, neighbor_list, parameter_list):
    r_min, r_max, kappa, min_seq_sep, density_threshold, density_kappa = parameter_list
    cb_density = calculate_cb_density(res_list, neighbor_list)
    r_min = float(r_min)
    r_max = float(r_max)
    kappa = float(kappa)
    min_seq_sep = int(min_seq_sep)
    density_threshold = float(density_threshold)
    density_kappa = float(density_kappa)
    phi_mediated_contact_well = np.zeros((2,3,20,20))

    eta_switching = 1
    z_m = 15
    r_min_I = 4.5
    r_max_I = 6.5

    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        rho_i = cb_density[res1globalindex]
        # z1 = get_z_position(res1)
        # alphaMembrane1 = interaction_well_2(z1, -z_m, z_m, eta_switching)

        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max+2.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)
            rho_j = cb_density[res2globalindex]
            # z2 = get_z_position(res2)
            # alphaMembrane2 = interaction_well_2(z1, -z_m, z_m, eta_switching)

            # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
            if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                res1type = get_res_type(res_list, res1)
                res2type = get_res_type(res_list, res2)
                rij = get_interaction_distance(res1, res2)

                # alpha_water = (1 - alphaMembrane1 * alphaMembrane2)
                # alpha_membrane = alphaMembrane1 * alphaMembrane2
                alpha_membrane = 1

                direct = interaction_well(rij, r_min_I, r_max_I, kappa)
                _pij_high_density = prot_water_switchFunc_sigmaProt(
                    rho_i, rho_j, density_threshold, density_kappa) * interaction_well(rij, r_min, r_max, kappa)
                _pij_low_density = prot_water_switchFunc_sigmaWater(
                    rho_i, rho_j, density_threshold, density_kappa) * interaction_well(rij, r_min, r_max, kappa)

                phi_mediated_contact_well[1][0][res1type][res2type] += alpha_membrane * direct
                phi_mediated_contact_well[1][1][res1type][res2type] += alpha_membrane * _pij_high_density
                phi_mediated_contact_well[1][2][res1type][res2type] += alpha_membrane * _pij_low_density
                if not res1type == res2type:
                    phi_mediated_contact_well[1][0][res2type][res1type] += alpha_membrane * direct
                    phi_mediated_contact_well[1][1][res2type][res1type] += alpha_membrane * _pij_high_density
                    phi_mediated_contact_well[1][2][res2type][res1type] += alpha_membrane * _pij_low_density

    phis_to_return = []
    ii = 1
    for jj in range(3):
        for i in range(20):
            for j in range(i, 20):
                phis_to_return.append(phi_mediated_contact_well[ii][jj][i][j])

    # kappa = 4.0
    # phi_burial = np.zeros((2, 3, 20))
    # rho_table = [[0.0, 3.0], [3.0, 6.0], [6.0, 9.0]]
    # for i in range(3):
    #     for res1globalindex, res1 in enumerate(res_list):
    #         res1index = get_local_index(res1)
    #         res1chain = get_chain(res1)
    #         res1type = get_res_type(res_list, res1)
    #         res1density = cb_density[res1globalindex]
    #         # z1 = get_z_position(res1)
    #         # alphaMembrane1 = interaction_well_2(z1, -z_m, z_m, eta_switching)
    #         # print res1globalindex, res1index, res1chain, res1type, res1density
    #         burial_i = interaction_well(res1density, rho_table[i][0], rho_table[i][1], kappa)
    #         # phi_burial[0][i][res1type] += burial_i * (1 - alphaMembrane1)
    #         phi_burial[1][i][res1type] += burial_i * 1

    # ii = 1
    # for i in range(3):
    #     for j in range(20):
    #         phis_to_return.append(phi_burial[ii][i][j])

    return phis_to_return

def get_pre_and_post(res_list, index):
    n = len(res_list)
    if index == 0:
        return res_list[0], res_list[1]
    elif index == n - 1:
        return res_list[index-1], res_list[index]
    else:
        return res_list[index-1], res_list[index+1]


def phi_pairwise_contact_multiLetter_well(res_list, neighbor_list, parameter_list):
    r_min, r_max, kappa, min_seq_sep = parameter_list
    r_min = float(r_min)
    r_max = float(r_max)
    kappa = float(kappa)
    min_seq_sep = int(min_seq_sep)
    # phi_pairwise_contact_multiLetter_well = np.zeros((4,4,20,20))
    phi_pairwise_contact_multiLetter_well = np.zeros((80,80))
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)

        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max+2.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)

            # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
            if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                res1type = get_res_type(res_list, res1)
                res2type = get_res_type(res_list, res2)
                rij = get_interaction_distance(res1, res2)

                res1_pre, res1_post = get_pre_and_post(res_list, res1globalindex)
                res2_pre, res2_post = get_pre_and_post(res_list, res2globalindex)
                res1_neighbor_type = get_neighbor_res_type(res1_pre, res1_post)
                res2_neighbor_type = get_neighbor_res_type(res2_pre, res2_post)
                res1_with_neighbor_type = res1_neighbor_type*20+res1type
                res2_with_neighbor_type = res2_neighbor_type*20+res2type

                phi_pairwise_contact_multiLetter_well[res1_with_neighbor_type][res2_with_neighbor_type] += interaction_well(rij, r_min, r_max, kappa)
                if not res1_with_neighbor_type == res2_with_neighbor_type:
                    phi_pairwise_contact_multiLetter_well[res2_with_neighbor_type][res1_with_neighbor_type] += interaction_well(rij, r_min, r_max, kappa)

    phis_to_return = []
    for i in range(80):
        for j in range(i, 80):
            phis_to_return.append(phi_pairwise_contact_multiLetter_well[i][j])
    return phis_to_return


def phi_density_mediated_contact_multiLetter_well(res_list, neighbor_list, parameter_list):
    r_min, r_max, kappa, min_seq_sep, density_threshold, density_kappa = parameter_list
    cb_density = calculate_cb_density(res_list, neighbor_list)
    r_min = float(r_min)
    r_max = float(r_max)
    kappa = float(kappa)
    min_seq_sep = int(min_seq_sep)
    density_threshold = float(density_threshold)
    density_kappa = float(density_kappa)
    # phi_mediated_contact_multiLetter_well = np.zeros((4,4, 2, 20,20))
    phi_mediated_contact_multiLetter_well = np.zeros((2, 80, 80))
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
                res1_with_neighbor_type = res1_neighbor_type*20+res1type
                res2_with_neighbor_type = res2_neighbor_type*20+res2type

                _pij_protein = prot_water_switchFunc_sigmaProt(
                    rho_i, rho_j, density_threshold, density_kappa) * interaction_well(rij, r_min, r_max, kappa)
                _pij_water = prot_water_switchFunc_sigmaWater(
                    rho_i, rho_j, density_threshold, density_kappa) * interaction_well(rij, r_min, r_max, kappa)
                phi_mediated_contact_multiLetter_well[0][res1_with_neighbor_type][res2_with_neighbor_type] += _pij_protein
                phi_mediated_contact_multiLetter_well[1][res1_with_neighbor_type][res2_with_neighbor_type] += _pij_water
                if not res1_with_neighbor_type == res2_with_neighbor_type:
                    phi_mediated_contact_multiLetter_well[0][res2_with_neighbor_type][res1_with_neighbor_type] += _pij_protein
                    phi_mediated_contact_multiLetter_well[1][res2_with_neighbor_type][res1_with_neighbor_type] += _pij_water

    phis_to_return = []
    for i in range(2):
        for j in range(80):
            for k in range(j, 80):
                phis_to_return.append(phi_mediated_contact_multiLetter_well[i][j][k])
    return phis_to_return

def phi_burial_multiLetter_well(res_list, neighbor_list, parameter_list):
    kappa = parameter_list[0]
    kappa = float(kappa)

    cb_density = calculate_cb_density(res_list, neighbor_list)

    # phi_burial_multiLetter = np.zeros((4, 3, 20))
    phi_burial_multiLetter = np.zeros((3, 80))
    rho_table = [[0.0, 3.0], [3.0, 6.0], [6.0, 9.0]]
    for i in range(3):
        for res1globalindex, res1 in enumerate(res_list):
            res1index = get_local_index(res1)
            res1chain = get_chain(res1)
            res1type = get_res_type(res_list, res1)
            res1density = cb_density[res1globalindex]

            res1_pre, res1_post = get_pre_and_post(res_list, res1globalindex)
            res1_neighbor_type = get_neighbor_res_type(res1_pre, res1_post)
            res1_with_neighbor_type = res1_neighbor_type*20+res1type
            # print res1globalindex, res1index, res1chain, res1type, res1density
            burial_i = interaction_well(res1density, rho_table[i][0], rho_table[i][1], kappa)
            # phi_burial_multiLetter[res1_neighbor_type][i][res1type] += burial_i
            phi_burial_multiLetter[i][res1_with_neighbor_type] += burial_i

    phis_to_return = []
    # for ii in range(4):
    #     for i in range(3):
    #         for j in range(20):
    #             phis_to_return.append(phi_burial_multiLetter[ii][i][j])

    for i in range(3):
        for j in range(80):
            phis_to_return.append(phi_burial_multiLetter[i][j])

    return phis_to_return

def read_gamma(gammaFile):
    data = np.loadtxt(gammaFile)
    gamma_direct = data[:210]
    gamma_mediated = data[210:]
    return gamma_direct, gamma_mediated
def get_gammas(gammaFile="./gamma.dat", memGammaFile="./membrane_gamma_rescaled.dat"):
    gamma_direct, gamma_mediated = read_gamma(gammaFile)
    nwell = 2
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
    if memGammaFile:
        mem_gamma_direct, mem_gamma_mediated = read_gamma(memGammaFile)
        m = 1  # membrane environment
        count = 0
        for i in range(20):
            for j in range(i, 20):
                gamma_ijm[m][i][j] = mem_gamma_direct[count][0]
                gamma_ijm[m][j][i] = mem_gamma_direct[count][0]
                count += 1
        count = 0
        for i in range(20):
            for j in range(i, 20):
                water_gamma_ijm[m][i][j] = mem_gamma_mediated[count][1]
                water_gamma_ijm[m][j][i] = mem_gamma_mediated[count][1]
                count += 1
        count = 0
        for i in range(20):
            for j in range(i, 20):
                protein_gamma_ijm[m][i][j] = mem_gamma_mediated[count][0]
                protein_gamma_ijm[m][j][i] = mem_gamma_mediated[count][0]
                count += 1
    return gamma_ijm, water_gamma_ijm, protein_gamma_ijm


def get_burial_gammas(burialGammaFile="./burial_gamma.dat"):
    return np.loadtxt(burialGammaFile)


def get_z_position(res1):
    return get_interaction_atom(res1).get_vector()[2]


def phi_normalize_relative_k(res_list, neighbor_list, parameter_list):
    phi_relative_k_well = np.zeros(11)
    cb_density = calculate_cb_density(res_list, neighbor_list)
    gamma_ijm, water_gamma_ijm, protein_gamma_ijm = get_gammas(memGammaFile=None)
    burial_gamma = get_burial_gammas()

    r_min_I = 4.5
    r_max_I = 6.5
    kappa = 5.0
    min_seq_sep = 10

    # Direct Term
    phi_pairwise_contact_well = np.zeros((2,2))
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max_I+2.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)
            # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
            if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                res1type = get_res_type(res_list, res1)
                res2type = get_res_type(res_list, res2)
                res1type_HP = get_res_type_HP(res_list, res1)
                res2type_HP = get_res_type_HP(res_list, res2)
                rij = get_interaction_distance(res1, res2)
                _direct = gamma_ijm[0][res1type][res2type] * interaction_well(rij, r_min_I, r_max_I, kappa)
                phi_pairwise_contact_well[res1type_HP][res2type_HP] += _direct
                if not res2type_HP == res1type_HP:
                    phi_pairwise_contact_well[res2type_HP][res1type_HP] += _direct
    # Mediated Term
    r_min = 6.5
    r_max = 9.5
    density_threshold = 2.6
    density_kappa = 7.0
    phi_mediated_contact_well = np.zeros((2,2,2))
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
                res1type_HP = get_res_type_HP(res_list, res1)
                res2type_HP = get_res_type_HP(res_list, res2)
                rij = get_interaction_distance(res1, res2)
                _pij_protein = prot_water_switchFunc_sigmaProt(
                    rho_i, rho_j, density_threshold, density_kappa) * protein_gamma_ijm[0][res1type][res2type]
                _pij_water = prot_water_switchFunc_sigmaWater(
                    rho_i, rho_j, density_threshold, density_kappa) * water_gamma_ijm[0][res1type][res2type]
                _pij_protein *= interaction_well(rij, r_min, r_max, kappa)
                _pij_water *= interaction_well(rij, r_min, r_max, kappa)
                phi_mediated_contact_well[0][res1type_HP][res2type_HP] += _pij_protein
                phi_mediated_contact_well[1][res1type_HP][res2type_HP] += _pij_water
                if not res1type_HP == res2type_HP:
                    phi_mediated_contact_well[0][res2type_HP][res1type_HP] += _pij_protein
                    phi_mediated_contact_well[1][res2type_HP][res1type_HP] += _pij_water
    # Burial Term
    phi_burial = np.zeros((2))
    burial_kappa = 4.0
    rho_table = [[0.0, 3.0], [3.0, 6.0], [6.0, 9.0]]
    for i in range(3):
        for res1globalindex, res1 in enumerate(res_list):
            res1index = get_local_index(res1)
            res1chain = get_chain(res1)
            res1type = get_res_type(res_list, res1)
            res1type_HP = get_res_type_HP(res_list, res1)
            res1density = cb_density[res1globalindex]
            # print res1globalindex, res1index, res1chain, res1type, res1density
            burial_i = interaction_well(res1density, rho_table[i][0], rho_table[i][1], burial_kappa)
            phi_burial[res1type_HP] += burial_i * burial_gamma[res1type][i]

    phis_to_return = []
    for i in range(2):
        for j in range(i, 2):
            phis_to_return.append(-phi_pairwise_contact_well[i][j])
    for i in range(2):
        for j in range(2):
            for k in range(j, 2):
                phis_to_return.append(-phi_mediated_contact_well[i][j][k])
    for i in range(2):
        phis_to_return.append(-phi_burial[i])
    return phis_to_return

def phi_relative_k_well(res_list, neighbor_list, parameter_list, z_m_high=None, z_m_low=None):
    phi_relative_k_well = np.zeros(2)
    cb_density = calculate_cb_density(res_list, neighbor_list)

    gamma_ijm, water_gamma_ijm, protein_gamma_ijm = get_gammas()
    r_min = 6.5
    r_max = 9.5

    r_min_I = 4.5
    r_max_I = 6.5
    kappa = 5.0
    min_seq_sep = 10
    # min_sequence_separation_mem = 13
    eta_switching = 10
    z_m = 15
    if z_m_high is None:
        z_m_high = z_m
    if z_m_low is None:
        z_m_low = -z_m
    density_threshold = 2.6
    density_kappa = 7.0
    # print("z_m_high", z_m_high)
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        rho_i = cb_density[res1globalindex]
        z1 = get_z_position(res1)
        # print("res1globalindex", res1globalindex, "z1", z1)
        alphaMembrane1 = interaction_well_2(z1, z_m_low, z_m_high, eta_switching)
        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max+2.0):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)
            z2 = get_z_position(res2)
            alphaMembrane2 = interaction_well_2(z1, z_m_low, z_m_high, eta_switching)
            rho_j = cb_density[res2globalindex]
            # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
            if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
                res1type = get_res_type(res_list, res1)
                res2type = get_res_type(res_list, res2)
                rij = get_interaction_distance(res1, res2)
                _pij_protein = prot_water_switchFunc_sigmaProt(
                    rho_i, rho_j, density_threshold, density_kappa) * protein_gamma_ijm[0][res1type][res2type]
                _pij_water = prot_water_switchFunc_sigmaWater(
                    rho_i, rho_j, density_threshold, density_kappa) * water_gamma_ijm[0][res1type][res2type]
                water_part = (_pij_protein + _pij_water) * interaction_well(rij, r_min, r_max, kappa)
                # direct term
                gamma = gamma_ijm[0][res1type][res2type]
                water_part += gamma * interaction_well(rij, r_min_I, r_max_I, kappa)
                _pij_protein = prot_water_switchFunc_sigmaProt(
                    rho_i, rho_j, density_threshold, density_kappa) * protein_gamma_ijm[1][res1type][res2type]
                _pij_water = prot_water_switchFunc_sigmaWater(
                    rho_i, rho_j, density_threshold, density_kappa) * water_gamma_ijm[1][res1type][res2type]
                membrane_part = (_pij_protein + _pij_water) * interaction_well(rij, r_min, r_max, kappa)
                gamma = gamma_ijm[1][res1type][res2type]
                membrane_part += gamma * interaction_well(rij, r_min_I, r_max_I, kappa)
                phi_relative_k_well[0] += (1 - alphaMembrane1 * alphaMembrane2) * water_part
                # phi_relative_k_well[1] += (1 - alphaMembrane1 * alphaMembrane2) * membrane_part
                phi_relative_k_well[1] += alphaMembrane1 * alphaMembrane2 * membrane_part
    phis_to_return = []
    for i in range(2):
        phis_to_return.append(phi_relative_k_well[i])
    return phis_to_return



six_letter_code_combinations_index = list(range(210))
six_letter_code_combinations = ['000004', '000013', '000022', '000031', '000040', '000103', '000112', '000121', '000130', '000202', '000211', '000220', '000301', '000310', '000400', '001003', '001012', '001021', '001030', '001102', '001111', '001120', '001201', '001210', '001300', '002002', '002011', '002020', '002101', '002110', '002200', '003001', '003010', '003100', '004000', '010003', '010012', '010021', '010030', '010102', '010111', '010120', '010201', '010210', '010300', '011002', '011011', '011020', '011101', '011110', '011200', '012001', '012010', '012100', '013000', '020002', '020011', '020020', '020101', '020110', '020200', '021001', '021010', '021100', '022000', '030001', '030010', '030100', '031000', '040000', '100003', '100012', '100021', '100030', '100102', '100111', '100120', '100201', '100210', '100300', '101002', '101011', '101020', '101101', '101110', '101200', '102001', '102010', '102100', '103000', '110002', '110011', '110020', '110101', '110110', '110200', '111001', '111010', '111100', '112000', '120001', '120010', '120100', '121000', '130000', '200002', '200011', '200020', '200101', '200110', '200200', '201001', '201010', '201100', '202000', '210001', '210010', '210100', '211000', '220000', '300001', '300010', '300100', '301000', '310000', '400000']
six_letter_code_combinations_map = dict(list(zip(six_letter_code_combinations, six_letter_code_combinations_index)))
inverse_six_letter_code_combinations_map = dict(list(zip(six_letter_code_combinations_index, six_letter_code_combinations)))
# (G)(CAST)(P)(IMLV)(NHQEDRK)(FYW)
six_letter_code_letters = {
    'I':3,
    'M':3,
    'L':3,
    'V':3,
    'F':5,
    'Y':5,
    'W':5,
    'G':0,
    'P':2,
    'C':1,
    'A':1,
    'S':1,
    'T':1,
    'N':4,
    'H':4,
    'Q':4,
    'E':4,
    'D':4,
    'R':4,
    'K':4
}

def count_general_amino_acid_combinations(number_of_amino_acid_types, number_of_sites):
    number_of_valid_codes = 0
    possible_codes = list(set(itertools.combinations_with_replacement(list(range(number_of_sites+1))*number_of_amino_acid_types, number_of_amino_acid_types)))
    for code in possible_codes:
        if sum(list(code)) == number_of_sites:
            print(code)
            number_of_valid_codes += 1
    print(number_of_valid_codes)

def calculate_six_letter_index(restypelist):
    six_letter_list = [0]*6
    for restype in restypelist:
        six_letter_list[six_letter_code_letters[restype]] += 1
    six_letter_list = [str(i) for i in six_letter_list]
    six_letter_string = ''.join(six_letter_list)
    return six_letter_code_combinations_map[six_letter_string], six_letter_string

def str2bool(string):
    if string == 'True':
        return True
    elif string == 'False':
        return False
    else:
        return None

def get_ca_ca_distance(res1, res2):
    ca1 = res1['CA']
    ca2 = res2['CA']
    return ca1 - ca2

def get_cb_cb_distance(res1, res2):
    return get_interaction_atom(res1) - get_interaction_atom(res2)

def get_protein_name_from_res_list(res_list):
    return res_list[0].get_parent().get_parent().get_parent().get_id().split('/')[-1].split('.')[0]


def phi_six_letter_four_body_helix_docking(res_list, neighbor_list, parameter_list):
    r_min, r_max, kappa, min_seq_sep, output_helix_docking_locations = parameter_list
    output_helix_docking_locations = str2bool(output_helix_docking_locations)
    helix_docking_locations_directory = "./data/"
    protein = get_protein_name_from_res_list(res_list)
    fileName = "%s%s-helix_docking_locations.dat" % (helix_docking_locations_directory, protein)
    import os
    write = False
    if output_helix_docking_locations and not os.path.exists(fileName):
        write = True
        output_helix_docking_locations_file = open(fileName, 'w')
    r_min = float(r_min)
    r_max = float(r_max)
    kappa = float(kappa)
    min_seq_sep = int(min_seq_sep)
    phis_to_return = [0]*126
    i_ip4_pairs = []
    # Find all (i, i+4) pairs
    for res_i_globalindex, res_i in enumerate(res_list):
        res_i_index = get_local_index(res_i)
        res_i_chain = get_chain(res_i)
        try:
            res_ip4 = res_list[res_i_globalindex+4]
        except IndexError:
            continue
        res_ip4_index = get_local_index(res_ip4)
        res_ip4_chain = get_chain(res_ip4)
        if res_ip4_index == res_i_index+4 and res_ip4_chain == res_i_chain:
            i_ip4_pairs.append([res_i, res_ip4])
    # Loop over pairs
    # print("Protein Size", len(res_list), "Pair Size", len(i_ip4_pairs))
    for i, res_pair1 in enumerate(i_ip4_pairs):
        for j, res_pair2 in enumerate(i_ip4_pairs):
            res_pair1_index_1 = get_local_index(res_pair1[0])
            res_pair1_index_2 = get_local_index(res_pair1[1])
            res_pair2_index_1 = get_local_index(res_pair2[0])
            res_pair2_index_2 = get_local_index(res_pair2[1])
            res_pair1_chain = get_chain(res_pair1[0])
            res_pair2_chain = get_chain(res_pair2[0])
            res_pair1_globalindex = get_global_index(res_list, res_pair1[0])
            res_pair2_globalindex = get_global_index(res_list, res_pair2[0])
            if (res_pair1_chain != res_pair2_chain and res_pair2_globalindex > res_pair1_globalindex) or \
                (res_pair2_index_1-res_pair1_index_1 >= min_seq_sep and \
                 res_pair2_index_2-res_pair1_index_1 >= min_seq_sep):
                # calculate six letter indices
                res_pair1_type_1 = three_to_one(res_pair1[0].get_resname())
                res_pair1_type_2 = three_to_one(res_pair1[1].get_resname())
                res_pair2_type_1 = three_to_one(res_pair2[0].get_resname())
                res_pair2_type_2 = three_to_one(res_pair2[1].get_resname())
                phi_index, six_letter_string = calculate_six_letter_index([res_pair1_type_1, res_pair1_type_2, res_pair2_type_1, res_pair2_type_2])
                # Calculate parallel phi

                d00 = get_cb_cb_distance(res_pair1[0], res_pair2[0])
                d11 = get_cb_cb_distance(res_pair1[1], res_pair2[1])
                d01 = get_cb_cb_distance(res_pair1[0], res_pair2[1])
                d10 = get_cb_cb_distance(res_pair1[1], res_pair2[0])
                parallel_phi = interaction_well(d00, r_min, r_max, kappa)*\
                               interaction_well(d11, r_min, r_max, kappa)
                # Calculate anti-parallel phi
                anti_parallel_phi = interaction_well(d01, r_min, r_max, kappa)*\
                                    interaction_well(d10, r_min, r_max, kappa)
                # Add the total into the correct index
                total_phi = parallel_phi + anti_parallel_phi
                if total_phi > 0.5 and write:
                    output_helix_docking_locations_file.write("%s %d %d %s %d %d %s%s%s%s %s %f" % (res_pair1_chain, res_pair1_index_1, res_pair1_index_2, res_pair2_chain, res_pair2_index_1, res_pair2_index_2, res_pair1_type_1, res_pair1_type_2, res_pair2_type_1, res_pair2_type_2, six_letter_string, total_phi))
                    output_helix_docking_locations_file.write(" %f %f %f %f\n" % (d00, d11, d01, d10))
                phis_to_return[phi_index] += total_phi
    # Return phis
    return phis_to_return


def make_helix_docking_contacts_tcl(input_file_name, output_file_name):
    output_file = open(output_file_name, 'w')
    output_file.write("mol modstyle 0 0 NewCartoon 0.300000 10.000000 4.100000 0\n")
    output_file.write("mol modmaterial 0 0 Transparent\n")
    output_file.write("mol modcolor 0 0 ColorID 0\n")
    output_file.write("mol modcolor 0 0 ColorID 2\n")
    output_file.write("mol color ColorID 2\n")
    output_file.write("mol representation NewCartoon 0.300000 10.000000 4.100000 0\n")
    output_file.write("mol selection all\n")
    output_file.write("mol material Transparent\n")
    for i, line in enumerate(open(input_file_name, 'r')):
        line = line.split()
        if len(line) == 0:
            continue
        chain1, res1, res2, chain2, res3, res4, aa, code_string = line
        output_file.write("mol addrep 0\n")
        output_file.write("mol modstyle %d 0 VDW 1.000000 12.000000\n" % (i+1))
        output_file.write("mol modcolor %d 0 ColorID 4\n" % (i+1))
        output_file.write("mol modselect %d 0 resid %s %s %s %s\n" % (i+1, res1, res2, res3, res4))
        output_file.write("mol color ColorID 4\n")
        output_file.write("mol material Opaque\n")

# def phi_six_letter_four_body_helix_docking(res_list, neighbor_list, parameter_list):
#     r_min, r_max, kappa, min_seq_sep, output_helix_docking_locations = parameter_list
#     output_helix_docking_locations = str2bool(output_helix_docking_locations)
#     helix_docking_locations_directory = "./"
#     if output_helix_docking_locations:
#         protein = get_protein_name_from_res_list(res_list)
#         output_helix_docking_locations_file = open("%s%s-helix_docking_locations.dat" % (helix_docking_locations_directory, protein), 'w')
#     r_min = float(r_min)
#     r_max = float(r_max)
#     kappa = float(kappa)
#     min_seq_sep = int(min_seq_sep)
#     phis_to_return = [0]*126
#     i_ip4_pairs = []
#     # Find all (i, i+4) pairs
#     for res_i_globalindex, res_i in enumerate(res_list):
#         res_i_index = get_local_index(res_i)
#         res_i_chain = get_chain(res_i)
#         try:
#             res_ip4 = res_list[res_i_globalindex+4]
#         except IndexError:
#             continue
#         res_ip4_index = get_local_index(res_ip4)
#         res_ip4_chain = get_chain(res_ip4)
#         if res_ip4_index == res_i_index+4 and res_ip4_chain == res_i_chain:
#             i_ip4_pairs.append([res_i, res_ip4])
#     # Loop over pairs
#     for i, res_pair1 in enumerate(i_ip4_pairs):
#         for j, res_pair2 in enumerate(i_ip4_pairs):
#             res_pair1_index_1 = get_local_index(res_pair1[0])
#             res_pair1_index_2 = get_local_index(res_pair1[1])
#             res_pair2_index_1 = get_local_index(res_pair2[0])
#             res_pair2_index_2 = get_local_index(res_pair2[1])
#             res_pair1_chain = get_chain(res_pair1[0])
#             res_pair2_chain = get_chain(res_pair2[0])
#             res_pair1_globalindex = get_global_index(res_list, res_pair1[0])
#             res_pair2_globalindex = get_global_index(res_list, res_pair2[0])
#             if (res_pair1_chain != res_pair2_chain and res_pair2_globalindex > res_pair1_globalindex) or \
#                 (res_pair2_index_1-res_pair1_index_1 >= min_seq_sep and \
#                  res_pair2_index_2-res_pair1_index_1 >= min_seq_sep):
#                 # calculate six letter indices
#                 res_pair1_type_1 = three_to_one(res_pair1[0].get_resname())
#                 res_pair1_type_2 = three_to_one(res_pair1[1].get_resname())
#                 res_pair2_type_1 = three_to_one(res_pair2[0].get_resname())
#                 res_pair2_type_2 = three_to_one(res_pair2[1].get_resname())
#                 phi_index, six_letter_string = calculate_six_letter_index([res_pair1_type_1, res_pair1_type_2, res_pair2_type_1, res_pair2_type_2])
#                 # Calculate parallel phi
#                 parallel_phi = interaction_well(get_ca_ca_distance(res_pair1[0], res_pair2[0]), r_min, r_max, kappa)*\
#                                interaction_well(get_ca_ca_distance(res_pair1[1], res_pair2[1]), r_min, r_max, kappa)
#                 # Calculate anti-parallel phi
#                 anti_parallel_phi = interaction_well(get_ca_ca_distance(res_pair1[0], res_pair2[1]), r_min, r_max, kappa)*\
#                                     interaction_well(get_ca_ca_distance(res_pair1[1], res_pair2[0]), r_min, r_max, kappa)
#                 # Add the total into the correct index
#                 total_phi = parallel_phi + anti_parallel_phi
#                 if total_phi > 0.8 and output_helix_docking_locations:
#                     output_helix_docking_locations_file.write("%s %d %d %s %d %d %s%s%s%s %s\n" % (res_pair1_chain, res_pair1_index_1, res_pair1_index_2, res_pair2_chain, res_pair2_index_1, res_pair2_index_2, res_pair1_type_1, res_pair1_type_2, res_pair2_type_1, res_pair2_type_2, six_letter_string))
#                 phis_to_return[phi_index] += total_phi
#     # Return phis
#     return phis_to_return
# def phi_protein_mediated_contact_well(res_list_tmonly, res_list_entire, neighbor_list, parameter_list, TCRmodeling=False):
#     r_min, r_max, kappa, min_seq_sep, density_threshold, density_kappa = parameter_list
#     cb_density = calculate_cb_density(res_list_entire, neighbor_list)
#     r_min = float(r_min)
#     r_max = float(r_max)
#     kappa = float(kappa)
#     min_seq_sep = int(min_seq_sep)
#     density_threshold = float(density_threshold)
#     density_kappa = float(density_kappa)
#     phi_mediated_contact_well = np.zeros((20, 20))
#     for res1globalindex, res1 in enumerate(res_list_entire):
#         res1index = get_local_index(res1)
#         res1chain = get_chain(res1)

#         # For TCR modeling, we only need the sequence in the peptide;
#         if TCRmodeling:

#             if (res1 in res_list_tmonly):

#                 for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max + 2.0):
#                     res2index = get_local_index(res2)
#                     res2chain = get_chain(res2)
#                     res2globalindex = get_global_index(res_list_entire, res2)

#                     # Get the density parameters of i an j
#                     rho_i = cb_density[res1globalindex]
#                     rho_j = cb_density[res2globalindex]
#                     rho_0 = density_threshold
#                     kappa = density_kappa

#                     # Here, we strictly consider only between the peptide and the TCR (chain C & D):
#                     # Res1 through tm_only, is already in chain E, we only need to control the res2 to
#                     # be in chian C or D;
#                     if (res2chain == 'C' or res2chain == 'D'):
#                         res1type = get_res_type(res_list_entire, res1)
#                         res2type = get_res_type(res_list_entire, res2)
#                         rij = get_interaction_distance(res1, res2)
#                         phi_mediated_contact_well[res1type][res2type] += prot_water_switchFunc_sigmaProt(
#                             rho_i, rho_j, rho_0, kappa) * interaction_well(rij, r_min, r_max, kappa)
#                         if not res1type == res2type:
#                             phi_mediated_contact_well[res2type][res1type] += prot_water_switchFunc_sigmaProt(
#                                 rho_j, rho_i, rho_0, kappa) * interaction_well(rij, r_min, r_max, kappa)
#             else:
#                 continue

#         else:

#             for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max + 2.0):
#                 res2index = get_local_index(res2)
#                 res2chain = get_chain(res2)
#                 res2globalindex = get_global_index(res_list_entire, res2)

#                 # Get the density parameters of i an j
#                 rho_i = cb_density[res1globalindex]
#                 rho_j = cb_density[res2globalindex]
#                 rho_0 = density_threshold
#                 kappa = density_kappa

#                 # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
#                 if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
#                     res1type = get_res_type(res_list_entire, res1)
#                     res2type = get_res_type(res_list_entire, res2)
#                     rij = get_interaction_distance(res1, res2)
#                     phi_mediated_contact_well[res1type][res2type] += prot_water_switchFunc_sigmaProt(
#                         rho_i, rho_j, rho_0, kappa) * interaction_well(rij, r_min, r_max, kappa)
#                     if not res1type == res2type:
#                         phi_mediated_contact_well[res2type][res1type] += prot_water_switchFunc_sigmaProt(
#                             rho_j, rho_i, rho_0, kappa) * interaction_well(rij, r_min, r_max, kappa)

#     phis_to_return = []
#     for i in range(20):
#         for j in range(i, 20):
#             phis_to_return.append(phi_mediated_contact_well[i][j])
#     return phis_to_return


# def phi_water_mediated_contact_well(res_list_tmonly, res_list_entire, neighbor_list, parameter_list, TCRmodeling=False):
#     r_min, r_max, kappa, min_seq_sep, density_threshold, density_kappa = parameter_list
#     cb_density = calculate_cb_density(res_list_entire, neighbor_list)
#     r_min = float(r_min)
#     r_max = float(r_max)
#     kappa = float(kappa)
#     min_seq_sep = int(min_seq_sep)
#     density_threshold = float(density_threshold)
#     density_kappa = float(density_kappa)
#     phi_mediated_contact_well = np.zeros((20, 20))
#     for res1globalindex, res1 in enumerate(res_list_entire):
#         res1index = get_local_index(res1)
#         res1chain = get_chain(res1)

#         # For TCR modeling, we only need the sequence in the peptide;
#         if TCRmodeling:

#             if (res1 in res_list_tmonly):

#                 for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max + 2.0):
#                     res2index = get_local_index(res2)
#                     res2chain = get_chain(res2)
#                     res2globalindex = get_global_index(res_list_entire, res2)

#                     # Get the density parameters of i an j
#                     rho_i = cb_density[res1globalindex]
#                     rho_j = cb_density[res2globalindex]
#                     rho_0 = density_threshold
#                     kappa = density_kappa

#                     # Here, we strictly consider only between the peptide and the TCR (chain C & D):
#                     # Res1 through tm_only, is already in chain E, we only need to control the res2 to
#                     # be in chian C or D;
#                     if (res2chain == 'C' or res2chain == 'D'):
#                         res1type = get_res_type(res_list_entire, res1)
#                         res2type = get_res_type(res_list_entire, res2)
#                         rij = get_interaction_distance(res1, res2)
#                         phi_mediated_contact_well[res1type][res2type] += prot_water_switchFunc_sigmaWater(
#                             rho_i, rho_j, rho_0, kappa) * interaction_well(rij, r_min, r_max, kappa)
#                         if not res1type == res2type:
#                             phi_mediated_contact_well[res2type][res1type] += prot_water_switchFunc_sigmaWater(
#                                 rho_j, rho_i, rho_0, kappa) * interaction_well(rij, r_min, r_max, kappa)

#             else:
#                 continue

#         else:

#             for res2 in get_neighbors_within_radius(neighbor_list, res1, r_max + 2.0):
#                 res2index = get_local_index(res2)
#                 res2chain = get_chain(res2)
#                 res2globalindex = get_global_index(res_list_entire, res2)

#                 # Get the density parameters of i an j
#                 rho_i = cb_density[res1globalindex]
#                 rho_j = cb_density[res2globalindex]
#                 rho_0 = density_threshold
#                 kappa = density_kappa

#                 # if res2index - res1index >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
#                 if res2globalindex - res1globalindex >= min_seq_sep or (res1chain != res2chain and res2globalindex > res1globalindex):
#                     res1type = get_res_type(res_list_entire, res1)
#                     res2type = get_res_type(res_list_entire, res2)
#                     rij = get_interaction_distance(res1, res2)
#                     phi_mediated_contact_well[res1type][res2type] += prot_water_switchFunc_sigmaWater(
#                         rho_i, rho_j, rho_0, kappa) * interaction_well(rij, r_min, r_max, kappa)
#                     if not res1type == res2type:
#                         phi_mediated_contact_well[res2type][res1type] += prot_water_switchFunc_sigmaWater(
#                             rho_j, rho_i, rho_0, kappa) * interaction_well(rij, r_min, r_max, kappa)

#     phis_to_return = []
#     for i in range(20):
#         for j in range(i, 20):
#             phis_to_return.append(phi_mediated_contact_well[i][j])
#     return phis_to_return


# Evaluate for those that decoy strucutres are provided;
native_structures_directory = "/opt/home/xl23/Working/Levine/jason/optimization/awsem/native_structures_pdbs_with_virtual_cbs"
decoy_structures_directory = "/opt/home/xl23/Working/Levine/jason/optimization/awsem/decoy_structures_pdbs_with_virtual_cbs"


def evaluate_phis_over_training_set_for_native_structures(native_training_set_file, phi_list_file_name, decoy_method, max_decoys, tm_only=False, num_processors=1, TCRmodeling=False):
    phi_list = read_phi_list(phi_list_file_name)
    print(phi_list)
    training_set = read_column_from_file(native_training_set_file, 1)
    print(training_set)

    function_to_evaluate = functools.partial(
        evaluate_phis_for_native_protein, tm_only=tm_only, TCRmodeling=TCRmodeling)
    arguments_lists = [training_set, [phi_list] * len(training_set), [
        decoy_method] * len(training_set), [max_decoys] * len(training_set)]
    call_independent_functions_on_n_processors(
        function_to_evaluate, arguments_lists, num_processors)
    # for protein in training_set:
    #   evaluate_phis_for_protein(protein, phi_list, decoy_method, max_decoys, tm_only=tm_only)


def evaluate_phis_over_training_set_for_decoy_structures(decoy_training_set_file, phi_list_file_name, decoy_method, max_decoys, tm_only=False, num_processors=1, TCRmodeling=False):
    phi_list = read_phi_list(phi_list_file_name)
    print(phi_list)
    training_set = read_column_from_file(decoy_training_set_file, 1)
    print(training_set)

    function_to_evaluate = functools.partial(
        evaluate_phis_for_decoy_protein, tm_only=tm_only, TCRmodeling=TCRmodeling)
    arguments_lists = [training_set, [phi_list] * len(training_set), [
        decoy_method] * len(training_set), [max_decoys] * len(training_set)]
    call_independent_functions_on_n_processors(
        function_to_evaluate, arguments_lists, num_processors)
    # for protein in training_set:
    #   evaluate_phis_for_protein(protein, phi_list, decoy_method, max_decoys, tm_only=tm_only)


def evaluate_phis_for_native_protein(protein, phi_list, decoy_method, max_decoys, tm_only=False, TCRmodeling=False):
    print(protein)
    structure = parse_pdb(os.path.join(native_structures_directory, protein))

    # Two lists of res_list, one for the peptide (selected by the .tm file), one for the entire list
    res_list_tmonly = get_res_list(structure, tm_only=True)
    res_list_entire = get_res_list(structure, tm_only=False)
    # Here, we are going to take every residues close to the pMHC peptide, so there is no restriction (tm_only) on what is going to be taken;
    neighbor_list = get_neighbor_list(structure, tm_only=False)

    sequence = get_sequence_from_structure(structure)

    for phi, parameters in phi_list:
        phi = globals()[phi]
        parameters_string = get_parameters_string(parameters)
        # check to see if the decoys are already generated
        number_of_lines_in_file = get_number_of_lines_in_file(os.path.join(
            phis_directory, "%s_%s_native_%s" % (phi.__name__, protein, parameters_string)))
        if not number_of_lines_in_file >= 1:
            output_file = open(os.path.join(phis_directory, "%s_%s_native_%s" % (
                phi.__name__, protein, parameters_string)), 'w')
            phis_to_write = phi(res_list_tmonly, res_list_entire,
                                neighbor_list, parameters, TCRmodeling=TCRmodeling)
            output_file.write(str(phis_to_write).strip(
                '[]').replace(',', '') + '\n')
            output_file.close()


def evaluate_phis_for_decoy_protein(protein, phi_list, decoy_method, max_decoys, tm_only=False, TCRmodeling=False):
    print(protein)
    structure = parse_pdb(os.path.join(decoy_structures_directory, protein))

    # Two lists of res_list, one for the peptide (selected by the .tm file), one for the entire list
    res_list_tmonly = get_res_list(structure, tm_only=True)
    res_list_entire = get_res_list(structure, tm_only=False)
    # Here, we are going to take every residues close to the pMHC peptide, so there is no restriction (tm_only) on what is going to be taken;
    neighbor_list = get_neighbor_list(structure, tm_only=False)

    sequence = get_sequence_from_structure(structure)

    for phi, parameters in phi_list:
        phi = globals()[phi]
        parameters_string = get_parameters_string(parameters)
        # check to see if the decoys are already generated
        number_of_lines_in_file = get_number_of_lines_in_file(os.path.join(
            phis_directory, "%s_%s_decoy_%s" % (phi.__name__, protein, parameters_string)))
        if not number_of_lines_in_file >= 1:
            output_file = open(os.path.join(phis_directory, "%s_%s_decoy_%s" % (
                phi.__name__, protein, parameters_string)), 'w')
            phis_to_write = phi(res_list_tmonly, res_list_entire,
                                neighbor_list, parameters, TCRmodeling=TCRmodeling)
            output_file.write(str(phis_to_write).strip(
                '[]').replace(',', '') + '\n')
            output_file.close()


def evaluate_phis_over_training_set_for_native_structures_Wei(training_set_file, phi_list_file_name, decoy_method, max_decoys, tm_only=False, multi_seq=False, sampleK=1000, num_processors=1, **kwargs):
    phi_list = read_phi_list(phi_list_file_name)
    print(phi_list)
    training_set = read_column_from_file(training_set_file, 1)
    print(training_set)
    function_to_evaluate = partial(evaluate_phis_for_protein_Wei, tm_only=tm_only, multi_seq=multi_seq, sampleK=sampleK, **kwargs)
    arguments_lists = [training_set, [phi_list]*len(training_set), [decoy_method]*len(training_set), [max_decoys]*len(training_set)]
    call_independent_functions_on_n_processors(function_to_evaluate, arguments_lists, num_processors)
    # for protein in training_set:
    #   evaluate_phis_for_protein(protein, phi_list, decoy_method, max_decoys, tm_only=tm_only)

def evaluate_phis_over_training_set_for_decoy_structures_Wei(training_set_file, phi_list_file_name, decoy_method, max_decoys, tm_only=False, num_processors=1, **kwargs):
    phi_list = read_phi_list(phi_list_file_name)
    print(phi_list)
    training_set = read_column_from_file(training_set_file, 1)
    print(training_set)
    function_to_evaluate = partial(evaluate_phis_for_decoy_protein_Wei, tm_only=tm_only, **kwargs)
    arguments_lists = [training_set, [phi_list]*len(training_set), [decoy_method]*len(training_set), [max_decoys]*len(training_set)]
    call_independent_functions_on_n_processors(function_to_evaluate, arguments_lists, num_processors)
# def evaluate_phis_over_training_set(training_set_file, phi_list_file_name, decoy_method, max_decoys, tm_only=False, num_processors=1, TCRmodeling=False):
#     phi_list = read_phi_list(phi_list_file_name)
#     print(phi_list)
#     training_set = read_column_from_file(training_set_file, 1)
#     print(training_set)

#     function_to_evaluate = functools.partial(
#         evaluate_phis_for_protein, tm_only=tm_only, TCRmodeling=TCRmodeling)
#     arguments_lists = [training_set, [phi_list] * len(training_set), [
#         decoy_method] * len(training_set), [max_decoys] * len(training_set)]
#     call_independent_functions_on_n_processors(
#         function_to_evaluate, arguments_lists, num_processors)
#     # for protein in training_set:
#     #   evaluate_phis_for_protein(protein, phi_list, decoy_method, max_decoys, tm_only=tm_only)

def evaluate_phis_for_protein_Wei(protein, phi_list, decoy_method, max_decoys, multi_seq=False, tm_only=False, sampleK=1000, separateDecoysNum=-1):
        print(protein)
        structure = parse_pdb(os.path.join(structures_directory,protein))
        res_list = get_res_list(structure, tm_only=tm_only)
        neighbor_list = get_neighbor_list(structure, tm_only=tm_only)
        sequence = get_sequence_from_structure(structure)
        for phi, parameters in phi_list:
            phiF = globals()[phi]
            parameters_string = get_parameters_string(parameters)
            mutate_whole_sequence(res_list, sequence)
            print(phi, parameters, parameters_string)
            # print("Naitve seq", get_sequence_from_structure(structure))
            # check to see if the decoys are already generated
            # number_of_lines_in_file = get_number_of_lines_in_file(os.path.join(phis_directory, "%s_%s_native_%s" % (phiF.__name__, protein, parameters_string)))
            # if not number_of_lines_in_file >= 1:
            output_file = open(os.path.join(phis_directory, "%s_%s_native_%s" % (phiF.__name__, protein, parameters_string)), 'w')
            if multi_seq is False:
                phis_to_write = phiF(res_list, neighbor_list, parameters)
                output_file.write(str(phis_to_write).strip('[]').replace(',', '')+'\n')
            if multi_seq is True:
                mutli_sequences = read_decoy_sequences(os.path.join(aligments_root_directory, "%s_filtered_0.05.seqs" % (protein)))
                size = len(mutli_sequences)
                if size < sampleK:
                    mutli_sequences = mutli_sequences * (sampleK//size+1)
                    print(protein, len(mutli_sequences), size)
                mutli_sequences = random.sample(mutli_sequences, sampleK)
                for i_decoy, decoy_sequence in enumerate(mutli_sequences):
                    mutate_whole_sequence(res_list, decoy_sequence)
                    phis_to_write = phiF(res_list, neighbor_list, parameters)
                    output_file.write(str(phis_to_write).strip('[]').replace(',', ' ')+'\n')
            output_file.close()
            # number_of_lines_in_file = get_number_of_lines_in_file(os.path.join(phis_directory, "%s_%s_decoys_%s_%s" % (phiF.__name__, protein, decoy_method, parameters_string)))
            # if not number_of_lines_in_file >= max_decoys:
            output_file = open(os.path.join(phis_directory, "%s_%s_decoys_%s_%s" % (phiF.__name__, protein, decoy_method, parameters_string)), 'w')
            if separateDecoysNum == -1:
                decoy_sequences = read_decoy_sequences(os.path.join(decoys_root_directory, "%s/%s.decoys" % (decoy_method, protein)))
                for i_decoy, decoy_sequence in enumerate(decoy_sequences):
                    if i_decoy >= max_decoys:
                        break
                    mutate_whole_sequence(res_list, decoy_sequence)
                    phis_to_write = phiF(res_list, neighbor_list, parameters)
                    output_file.write(str(phis_to_write).strip('[]').replace(',', ' ')+'\n')
                output_file.close()
            elif separateDecoysNum >= 1:
                for ii in range(separateDecoysNum):
                    proteinName = f"{protein}_{ii}"
                    decoy_sequences = read_decoy_sequences(os.path.join(decoys_root_directory, "%s/%s.decoys" % (decoy_method, proteinName)))
                    for i_decoy, decoy_sequence in enumerate(decoy_sequences):
                        if i_decoy >= max_decoys:
                            break
                        mutate_whole_sequence(res_list, decoy_sequence)
                        phis_to_write = phiF(res_list, neighbor_list, parameters)
                        output_file.write(str(phis_to_write).strip('[]').replace(',', ' ')+'\n')
                output_file.close()
            else:
                print("must be at least 1 or -1.")
def evaluate_phis_for_decoy_protein_Wei(protein, phi_list, decoy_method, max_decoys, tm_only=False, withBiased=False, mode=0, pickle=True):
        print(protein, withBiased)
        if mode == 1:
            baseProtein = protein.split("_")[0]
        elif mode == 0:
            baseProtein = protein
        structure = parse_pdb(os.path.join(structures_directory, baseProtein))
        res_list = get_res_list(structure)
        neighbor_list = get_neighbor_list(structure)
        sequence = get_sequence_from_structure(structure)
        for phi, parameters in phi_list:
            phiF = globals()[phi]
            parameters_string = get_parameters_string(parameters)
            print(phi, parameters, parameters_string)
            mutate_whole_sequence(res_list, sequence)
            # check to see if the decoys are already generated
            # number_of_lines_in_file = get_number_of_lines_in_file(os.path.join(phis_directory, "%s_%s_native_%s" % (phiF.__name__, protein, parameters_string)))
            # if not number_of_lines_in_file >= 1:
            output_file = open(os.path.join(phis_directory, "%s_%s_native_%s" % (phiF.__name__, protein, parameters_string)), 'w')
            phis_to_write = phiF(res_list, neighbor_list, parameters)
            output_file.write(str(phis_to_write).strip('[]').replace(',', '')+'\n')
            output_file.close()
            # number_of_lines_in_file = get_number_of_lines_in_file(os.path.join(phis_directory, "%s_%s_decoys_%s_%s" % (phiF.__name__, protein, decoy_method, parameters_string)))
            # if not number_of_lines_in_file >= max_decoys:
            output_file = open(os.path.join(phis_directory, "%s_%s_decoys_%s_%s" % (phiF.__name__, protein, decoy_method, parameters_string)), 'w')
            # decoy_sequences = read_decoy_sequences(os.path.join(decoys_root_directory, "%s/%s.decoys" % (decoy_method, protein)))
            if withBiased:
                if pickle:
                    decoy_structures, Qs = read_decoy_structures_andQ(os.path.join(decoys_root_directory, "%s/%s.pkl" % (decoy_method, protein)))
                else:
                    decoy_structures, Qs = read_decoy_structures_andQ(os.path.join(decoys_root_directory, "%s/%s.decoys" % (decoy_method, protein)))
                output_file_Q = open(os.path.join(phis_directory, "%s_%s_decoysQ_%s_%s" % (phiF.__name__, protein, decoy_method, parameters_string)), 'w')
                for i_decoy, Q in enumerate(Qs):
                    if i_decoy >= max_decoys:
                        print("i_decoy larger than max_decoys")
                    output_file_Q.write(str(Q)+"\n")
                output_file_Q.close()
            else:
                decoy_structures = read_decoy_structures(os.path.join(decoys_root_directory, "%s/%s.decoys" % (decoy_method, protein)))

            for i_decoy, decoy_structure in enumerate(decoy_structures):
                if i_decoy % 1000 == 0:
                    print(i_decoy)
                if i_decoy >= max_decoys:
                    break
                # decoy_structure = parse_pdb(os.path.join(structures_directory,protein))
                decoy_res_list = get_res_list(decoy_structure)
                decoy_neighbor_list = get_neighbor_list(decoy_structure)
                phis_to_write = phiF(decoy_res_list, decoy_neighbor_list, parameters)
                # if withBiased:
                #     # phis_to_write *= 1 - Qs[i_decoy]
                #     # phis_to_write = [x * (1 - float(Qs[i_decoy])) for x in phis_to_write]
                #     output_file_Q.write(str(Qs[i_decoy])+"\n")
                output_file.write(str(phis_to_write).strip('[]').replace(',', ' ')+'\n')
            output_file.close()
            # if withBiased:
            #     output_file_Q.close()

'''
def evaluate_phis_for_protein(protein, phi_list, decoy_method, max_decoys, tm_only=False, TCRmodeling=False):
    print(protein)
    structure = parse_pdb(os.path.join(native_structures_directory, protein))

    # Two lists of res_list, one for the peptide (selected by the .tm file), one for the entire list
    res_list_tmonly = get_res_list(structure, tm_only=True)
    res_list_entire = get_res_list(structure, tm_only=False)
    # Here, we are going to take every residues close to the pMHC peptide, so there is no restriction (tm_only) on what is going to be taken;
    neighbor_list = get_neighbor_list(structure, tm_only=False)

    sequence = get_sequence_from_structure(structure)

    for phi, parameters in phi_list:
        phi = globals()[phi]
        parameters_string = get_parameters_string(parameters)
        # check to see if the decoys are already generated
        number_of_lines_in_file = get_number_of_lines_in_file(os.path.join(
            phis_directory, "%s_%s_native_%s" % (phi.__name__, protein, parameters_string)))
        if not number_of_lines_in_file >= 1:
            output_file = open(os.path.join(phis_directory, "%s_%s_native_%s" % (
                phi.__name__, protein, parameters_string)), 'w')
            phis_to_write = phi(res_list_tmonly, res_list_entire,
                                neighbor_list, parameters, TCRmodeling=TCRmodeling)
            output_file.write(str(phis_to_write).strip(
                '[]').replace(',', '') + '\n')
            output_file.close()
        number_of_lines_in_file = get_number_of_lines_in_file(os.path.join(
            phis_directory, "%s_%s_decoys_%s_%s" % (phi.__name__, protein, decoy_method, parameters_string)))
        if not number_of_lines_in_file >= max_decoys:
            output_file = open(os.path.join(phis_directory, "%s_%s_decoys_%s_%s" % (
                phi.__name__, protein, decoy_method, parameters_string)), 'w')
            decoy_sequences = read_decoy_sequences(os.path.join(
                decoys_root_directory, "%s/%s.decoys" % (decoy_method, protein)))
            for i_decoy, decoy_sequence in enumerate(decoy_sequences):
                if i_decoy >= max_decoys:
                    break
                mutate_whole_sequence(res_list_entire, decoy_sequence)
                phis_to_write = phi(res_list_tmonly, res_list_entire,
                                    neighbor_list, parameters, TCRmodeling=TCRmodeling)
                output_file.write(str(phis_to_write).strip(
                    '[]').replace(',', ' ') + '\n')
            output_file.close()
'''

def mutate_whole_sequence(res_list, new_sequence):
    for i in range(len(res_list)):
        try:
            res_list[i].resname = one_to_three(new_sequence[i])
        except:
            pass
            # print(f"Skip {i}th, {new_sequence[i]}, using the old residue name")
    return res_list


def read_column_from_file(file_name, column, header_comment_syntax="#", num_header_lines=0, column_delimiter=''):
    list_to_return = []
    for i, line in enumerate(open(file_name, 'r')):
        line = line.strip('\n')
        if len(line) == 0:
            continue
        if i < num_header_lines:
            continue
        if line[0] == header_comment_syntax:
            continue
        if column_delimiter == '':
            line = line.split()
        else:
            line = line.split(column_delimiter)
        list_to_return.append(line[int(column) - 1])
    return list_to_return


def read_phi_list(phi_list_file_name, header_comment_syntax="#", num_header_lines=0, column_delimiter=' '):
    input_file = open(phi_list_file_name, 'r')
    phi_list = []
    for i, line in enumerate(input_file):
        line = line.strip()
        if len(line) == 0:
            continue
        if i < num_header_lines:
            continue
        if line[0] == header_comment_syntax:
            continue
        line = line.split(column_delimiter, 1)
        try:
            parameters = line[1].split()
        except IndexError:
            parameters = []
        phi_list.append([line[0], parameters])
    return phi_list


def get_parameters_string(parameters):
    parameter_string = ""
    for parameter in parameters:
        parameter_string += parameter
        if not parameters.index(parameter) + 1 == len(parameters):
            parameter_string += '_'
    return parameter_string


def generate_decoy_sequences(proteins_list_file_name, methods=['shuffle', 'cyclic'], num_decoys=[1000, 1000], databaseLocation=".", separateDecoysNum=-1):
    protein_list = read_column_from_file(proteins_list_file_name, 1)
    os.chdir(decoys_root_directory)
    for i, method in enumerate(methods):
        if not os.path.exists(method):
            os.makedirs(method)
        os.chdir(method)
        print(os.getcwd())
        for protein in protein_list:
            print(method, protein)
            if separateDecoysNum == -1:
                output_file = open("%s.decoys" % protein, 'w')
                for j in range(num_decoys[i]):
                    output_file.write(generate_decoy_sequence(
                        protein, method=method, degree=j, databaseLocation=databaseLocation) + '\n')
                output_file.close()
            elif separateDecoysNum >= 1:
                for ii in range(separateDecoysNum):
                    outName = f"{protein}_{ii}.decoys"
                    output_file = open(outName, 'w')
                    oneFileDecoysN = int(num_decoys[i]/separateDecoysNum)
                    for j in range(oneFileDecoysN):
                        output_file.write(generate_decoy_sequence(
                            protein, method=method, degree=j, databaseLocation=databaseLocation) + '\n')
                    output_file.close()
            else:
                print("must be at least 1 or -1.")
        os.chdir('..')
    os.chdir("..")


# def generate_decoy_structures(proteins_list_file_name, methods=['lammps'], num_decoys=[1000], databaseLocation="."):
#     protein_list = read_column_from_file(proteins_list_file_name, 1)
#     os.chdir(decoys_root_directory)
#     for i, method in enumerate(methods):
#         if not os.path.exists(method):
#             os.makedirs(method)
#         os.chdir(method)
#         print(os.getcwd())
#         for protein in protein_list:
#             print(method, protein)
#             output_file = open("%s.decoys" % protein, 'w')
#             for j in range(num_decoys[i]):
#                 output_file.write(generate_decoy_structure(
#                     protein, method=method, degree=j, databaseLocation=databaseLocation) + '\n')
#             output_file.close()
#         os.chdir('..')
#     os.chdir("..")


def shuffle_string(string):
    list_string = list(string)
    random.shuffle(list_string)
    return ''.join(list_string)

def cyclically_permute_string(string, degree):
    string = list(string)
    n = len(string)
    return ''.join([string[i - (degree % n)] for i in range(n)])

def get_sublist(list_name, indices_list):
    return [list_name[i] for i in indices_list]

def get_sublist_complement(list_name, indices_list):
    return [list_name[i] for i in range(len(list_name)) if i not in indices_list]

membrane_database_root = os.getenv('MEMBRANE_DATABASE_ROOT', 'C:\\Users\\Dell\\research\\databases\\membrane_protein_database\\')
tm_root_directory = os.path.join(membrane_database_root, 'cleaned_database', 'tms')

def generate_decoy_sequence(protein, method='TCR_randomization', degree=None, databaseLocation="."):

    sequences_root_directory = os.path.join(databaseLocation,"database/S20_seq/")

    with open("%s%s.seq" % (sequences_root_directory, protein), "r") as sequence_file:
        native_sequence = sequence_file.read().replace('\n', '')

    if method == 'shuffle':
        return shuffle_string(native_sequence)
    elif method == 'cyclic':
        if degree == None:
            print("Must specify degree with method cyclic")
            sys.exit()
        return cyclically_permute_string(native_sequence, degree)
    elif method == 'constrained_shuffle' or method == 'constrained_cyclic':
        native_sequence = list(native_sequence)
        tm = read_column_from_file(os.path.join(
            tm_root_directory, protein + '.tm'), 1)
        outside_indices = [i for i, x in enumerate(tm) if x == '1' or x == '3']
        outside_list = get_sublist(native_sequence, outside_indices)
        outside_string = ''.join(outside_list)
        membrane_indices = [i for i, x in enumerate(tm) if x == '2']
        membrane_list = get_sublist(native_sequence, membrane_indices)
        membrane_string = ''.join(membrane_list)
        new_sequence = []
        if method == 'constrained_shuffle':
            outside_string = shuffle_string(outside_string)
            membrane_string = shuffle_string(membrane_string)
        elif method == 'constrained_cyclic':
            if degree == None:
                print("Must specify degree with method cyclic")
                sys.exit()
            outside_string = cyclically_permute_string(outside_string, degree)
            membrane_string = cyclically_permute_string(
                membrane_string, degree)
        outside_list = list(outside_string)
        membrane_list = list(membrane_string)
        for region in tm:
            if region == '1' or region == '3':
                new_sequence.append(outside_list.pop(0))
            if region == '2':
                new_sequence.append(membrane_list.pop(0))
        return ''.join(new_sequence)
    elif method == 'TCR_randomization':
        sequence_toberandomized = list(native_sequence)
        gBinder_sequences = open(
            "gBinder_sequences.txt", 'r').read().splitlines()
        resids_toberandomized = open(
            "randomize_position_file.txt", 'r').readline().split(' ')
        for resid_toberandomized in resids_toberandomized:
            # Residue to be randomly mutated to (randomly selected from 20 amino acids);
            resAbbr = random.choice(list(
                ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]))

            # Replace the corresponding amino acid with 20 possibilities;
            sequence_toberandomized[int(resid_toberandomized) - 1] = resAbbr

        newsequence = ''.join(sequence_toberandomized)
        if newsequence in gBinder_sequences:
            return
        else:
            return newsequence

# def generate_decoy_structure(protein, method='lammps', degree=None, databaseLocation="."):
#     structures_root_directory = os.path.join(databaseLocation,"database/S20_seq/")
#     with open("%s%s.seq" % (structures_root_directory, protein), "r") as structure_file:
#         native_structure = structure_file.read().replace('\n', '')
#     if method == 'shuffle':
#         return shuffle_string(native_structure)

def get_total_phis_and_parameter_string(phi_list, training_set, mode=0, simulation_location_list=None, simulation_location_list_dic=None, multiSeq=False, **kwargs):
    full_parameters_string = ""
    # Find out how many total phi_i there are
    total_phis = 0
    num_phis = []
    for phi_and_parameters in phi_list:
        phi = phi_and_parameters[0]
        full_parameters_string += phi
        parameters = phi_and_parameters[1]
        i_phi = phi_list.index(phi_and_parameters)
        parameters_string = get_parameters_string(parameters)
        full_parameters_string += parameters_string
        for i_protein, protein in enumerate(training_set):
            if i_protein > 0:
                break
            if mode == 1:
                protein = protein + f"_{simulation_location_list[0]}"
            elif mode == 2:
                simulation_location_list = simulation_location_list_dic[protein]
                protein = protein + f"_{simulation_location_list[0]}"

            if multiSeq:
                input_file = np.loadtxt(os.path.join(phis_directory, "%s_%s_native_%s" % (
                                        phi, protein, parameters_string)))
                _, n = input_file.shape
                num_phis.append(n)
                total_phis += n
            else:
                input_file = open(os.path.join(phis_directory, "%s_%s_native_%s" % (
                    phi, protein, parameters_string)), 'r')
                for line in input_file:
                    line = line.strip().split()
                    num_phis.append(len(line))
                    total_phis += len(line)
    return total_phis, full_parameters_string, num_phis


def get_total_phis_and_parameter_string_decoy_structures_provided(phi_list, training_set):
    full_parameters_string = ""
    # Find out how many total phi_i there are
    total_phis = 0
    num_phis = []
    for phi_and_parameters in phi_list:
        phi = phi_and_parameters[0]
        full_parameters_string += phi
        parameters = phi_and_parameters[1]
        i_phi = phi_list.index(phi_and_parameters)
        parameters_string = get_parameters_string(parameters)
        full_parameters_string += parameters_string
        for i_protein, protein in enumerate(training_set):
            if i_protein > 0:
                break
            input_file = open(os.path.join(phis_directory, "%s_%s_decoys_lammps_%s" % (
                phi, protein, parameters_string)), 'r')
            for line in input_file:
                line = line.strip().split()
                num_phis.append(len(line))
                total_phis += len(line)
    return total_phis, full_parameters_string, num_phis


def read_decoy_phi_structures_provided(protein, phi_list, total_phis, jackhmmer=False):
    phi_decoy = np.zeros(total_phis)
    i_phi = 0
    for phi_and_parameters in phi_list:
        phi = phi_and_parameters[0]
        parameters = phi_and_parameters[1]
        parameters_string = get_parameters_string(parameters)
        if jackhmmer:
            input_file = open(os.path.join(jackhmmer_phis_directory, "%s_%s_decoys_lammps_%s" % (
                phi, protein, parameters_string)), 'r')
        else:
            input_file = open(os.path.join(phis_directory, "%s_%s_decoys_lammps_%s" % (
                phi, protein, parameters_string)), 'r')

        for line in input_file:
            line = line.strip().split()
            for i_value, value_i in enumerate(line):
                phi_decoy[i_phi] = float(line[i_value])
                i_phi += 1
    return phi_decoy


def read_native_phi(protein, phi_list, total_phis, jackhmmer=False, mode=0, simulation_location_list=None, simulation_location_list_dic=None, multiSeq=False):
    i_phi = 0
    if mode == 1:
        protein = protein + f"_{simulation_location_list[0]}"
    elif mode == 2:
        simulation_location_list = simulation_location_list_dic[protein]
        protein = protein + f"_{simulation_location_list[0]}"
    if multiSeq:
        # print("MultiSeq")
        phi_native = np.zeros(total_phis)
        for phi_and_parameters in phi_list:
            phi = phi_and_parameters[0]
            parameters = phi_and_parameters[1]
            parameters_string = get_parameters_string(parameters)

            input_file = np.loadtxt(os.path.join(phis_directory, "%s_%s_native_%s" % (
                                    phi, protein, parameters_string)))
            average_input = np.average(input_file, axis=0)
            n = len(average_input)
            phi_native[i_phi:i_phi+n] = average_input
            i_phi += n
    else:
        phi_native = np.zeros(total_phis)
        for phi_and_parameters in phi_list:
            phi = phi_and_parameters[0]
            parameters = phi_and_parameters[1]
            parameters_string = get_parameters_string(parameters)
            if jackhmmer:
                input_file = open(os.path.join(jackhmmer_phis_directory, "%s_%s_native_%s" % (
                    phi, protein, parameters_string)), 'r')
            else:
                input_file = open(os.path.join(phis_directory, "%s_%s_native_%s" % (
                    phi, protein, parameters_string)), 'r')

            for line in input_file:
                # print(line, len(line))

                line = line.strip().split()
                # print(len(line))
                for i_value, value_i in enumerate(line):
                    try:
                        phi_native[i_phi] = float(line[i_value])
                        i_phi += 1
                    except:
                        print("ERROR", phi, protein, parameters_string)
                        print(line, len(line))
                        exit()
    return phi_native


def read_decoy_phis(protein, phi_list, total_phis, num_phis, num_decoys, decoy_method, jackhmmer=False, mode=0, simulation_location_list=None, simulation_location_list_dic=None, **kwargs):
    if mode == 0:
        protein_list = [protein]
    elif mode == 1:
        protein_list = [protein+f"_{x}" for x in simulation_location_list]
        num_decoys = int(num_decoys/len(protein_list))
    elif mode == 2:
        simulation_location_list = simulation_location_list_dic[protein]
        protein_list = [protein+f"_{x}" for x in simulation_location_list]
        num_decoys = int(num_decoys/len(protein_list))
    all_phis = []
    for protein in protein_list:
        # print(protein)
        phi_i_decoy = np.zeros((num_decoys, total_phis))
        for i_phi_function, phi_and_parameters in enumerate(phi_list):
            # print(protein, i_phi_function, phi_and_parameters)
            phi = phi_and_parameters[0]
            parameters = phi_and_parameters[1]
            i_phi = phi_list.index(phi_and_parameters)
            parameters_string = get_parameters_string(parameters)
            if jackhmmer:
                input_file = open(os.path.join(jackhmmer_phis_directory, "%s_%s_decoys_%s" % (
                    phi, protein, parameters_string)), 'r')
            else:
                input_file = open(os.path.join(phis_directory, "%s_%s_decoys_%s_%s" % (
                    phi, protein, decoy_method, parameters_string)), 'r')

            for i_decoy, line in enumerate(input_file):
                if i_decoy >= num_decoys:
                    break
                first_phi = np.cumsum(num_phis)[
                    i_phi_function] - num_phis[i_phi_function]
                i_phi = first_phi
                line = line.strip().split()
                for i_value, value_i in enumerate(line):
                    phi_i_decoy[i_decoy][i_phi] = float(line[i_value])
                    i_phi += 1
        all_phis.append(phi_i_decoy)
    return np.concatenate(all_phis)

def read_decoyQ_phis(protein, phi_list, total_phis, num_phis, num_decoys, decoy_method, jackhmmer=False, mode=0, simulation_location_list=None, simulation_location_list_dic=None, **kwargs):
    if mode == 0:
        protein_list = [protein]
    elif mode == 1:
        protein_list = [protein+f"_{x}" for x in simulation_location_list]
        num_decoys = int(num_decoys/len(protein_list))
    elif mode == 2:
        simulation_location_list = simulation_location_list_dic[protein]
        protein_list = [protein+f"_{x}" for x in simulation_location_list]
        num_decoys = int(num_decoys/len(protein_list))
    all_phis = []
    for protein in protein_list:
        # print(protein)
        phi_i_decoy = np.zeros((num_decoys, 1))
        for i_phi_function, phi_and_parameters in enumerate(phi_list):
            # print(protein, i_phi_function, phi_and_parameters)
            phi = phi_and_parameters[0]
            parameters = phi_and_parameters[1]
            i_phi = phi_list.index(phi_and_parameters)
            parameters_string = get_parameters_string(parameters)
            # try:
            input_file = open(os.path.join(phis_directory, "%s_%s_decoysQ_%s_%s" % (
                phi, protein, decoy_method, parameters_string)), 'r')
            # except:
            #     input_file = open(os.path.join(phis_directory, "%s_%s_decoysQ_%s_%s" % (
            #         "phi_normalize_relative_k", protein, "lammps", "1")), 'r')
            for i_decoy, line in enumerate(input_file):
                # if i_decoy >= num_decoys:
                #     break
                first_phi = np.cumsum(num_phis)[
                    i_phi_function] - num_phis[i_phi_function]
                i_phi = first_phi
                line = line.strip()
                Q_i = float(line)
                phi_i_decoy[i_decoy] = Q_i
                # i_phi += 1
            break
        all_phis.append(phi_i_decoy)
    return np.concatenate(all_phis)

def calculate_A_and_B(average_phi_decoy, phi_native, total_phis, num_decoys, phi_i_decoy):
    A = average_phi_decoy - phi_native

    half_B = np.zeros((total_phis, total_phis))
    std_half_B = np.zeros((total_phis, total_phis))
    for i in range(total_phis):
        for j in range(total_phis):
            phi_i_phi_j_values = np.zeros(num_decoys)
            for i_decoy in range(num_decoys):
                phi_i_phi_j_i_decoy = phi_i_decoy[i_decoy][i] * \
                    phi_i_decoy[i_decoy][j]
                phi_i_phi_j_values[i_decoy] = phi_i_phi_j_i_decoy
                half_B[i][j] += phi_i_phi_j_i_decoy
            half_B[i][j] /= float(num_decoys)
            std_half_B[i][j] = np.std(phi_i_phi_j_values)

    other_half_B = np.zeros((total_phis, total_phis))
    for i in range(total_phis):
        for j in range(total_phis):
            other_half_B[i][j] += average_phi_decoy[i] * average_phi_decoy[j]

    B = half_B - other_half_B

    return A, B, half_B, other_half_B, std_half_B

def calculate_A_and_B_wei(average_phi_decoy, phi_native, all_phis, reduceMemory=True):
    print("calculate_A_and_B")
    os.system("echo 'calculate_A_and_B' >> log")
    print(datetime.datetime.now())
    A = average_phi_decoy - phi_native
    size_of_training_set, num_decoys, total_phis = all_phis.shape
    half_B = np.zeros((total_phis, total_phis))
    std_half_B = np.zeros((total_phis, total_phis))
    other_half_B = np.zeros((total_phis, total_phis))
    if reduceMemory:
        for p in range(size_of_training_set):
            os.system(f"echo '{p}' >> log")
            phis_i = all_phis[p].reshape(num_decoys, total_phis, 1)
            for j in range(total_phis):
                phis_j = phis_i[:, j].reshape(num_decoys, 1, 1)
                half_B[j] += np.average(phis_i * phis_j, axis=0).reshape(total_phis)
                std_half_B[j] += np.std(phis_i * phis_j, axis=0).reshape(total_phis)
    else:
        for p in range(size_of_training_set):
            os.system(f"echo '{p}' >> log")
            phis_i = all_phis[p].reshape(num_decoys, total_phis, 1)
            phis_j = all_phis[p].reshape(num_decoys, 1, total_phis)
            half_B += np.average(phis_i * phis_j, axis=0)
            std_half_B += np.std(phis_i * phis_j, axis=0)
    half_B /= size_of_training_set
    std_half_B /= size_of_training_set

    for p in range(size_of_training_set):
        average_phi = np.average(all_phis[p], axis=0)
        other_half_B += average_phi.reshape(total_phis, 1) * average_phi.reshape(1, total_phis)
    other_half_B /= size_of_training_set
    print("End")
    print(datetime.datetime.now())
    os.system("echo 'End' >> log")
    # for i in range(total_phis):
    #     for j in range(total_phis):
    #         for protein_i in range(size_of_training_set):
    #             protein_i_phi_i_phi_j_values = np.zeros(num_decoys)
    #             protein_i_half_B = np.zeros((total_phis, total_phis))
    #             protein_i_std_half_B = np.zeros((total_phis, total_phis))
    #             for i_decoy in range(num_decoys):
    #                 phi_i_phi_j = all_phis[protein_i][i_decoy][i] * \
    #                     all_phis[protein_i][i_decoy][j]
    #                 protein_i_phi_i_phi_j_values[i_decoy] = phi_i_phi_j
    #                 protein_i_half_B[i][j] += phi_i_phi_j
    #             protein_i_half_B[i][j] /= float(num_decoys)
    #             protein_i_std_half_B[i][j] = np.std(protein_i_phi_i_phi_j_values)
    #             half_B[i][j] += protein_i_half_B[i][j]
    #             std_half_B[i][j] += protein_i_std_half_B[i][j]
    # half_B /= size_of_training_set
    # std_half_B /= size_of_training_set

    # for i in range(total_phis):
    #     for j in range(total_phis):
    #         other_half_B[i][j] += average_phi_decoy[i] * average_phi_decoy[j]

    B = half_B - other_half_B

    return A, B, half_B, other_half_B, std_half_B

# def calculate_A_B_and_gamma_decoy_structures_provided(native_training_set_file, decoy_training_set_file, phi_list_file_name, decoy_method, num_decoys, noise_filtering=True, jackhmmer=False):
#     phi_list = read_phi_list(phi_list_file_name)
#     native_training_set = read_column_from_file(native_training_set_file, 1)

#     # Find out how many total phi_i there are and get full parameter string
#     total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string(
#         phi_list, native_training_set)

#     phi_native_i_protein = np.zeros((len(native_training_set), total_phis))
#     for i_protein, protein in enumerate(native_training_set):
#         phi_native_i_protein[i_protein] = read_native_phi(
#             protein, phi_list, total_phis, jackhmmer=jackhmmer)

#     phi_native = np.average(phi_native_i_protein, axis=0)

#     # Output to a file;
#     file_prefix = "%s%s_%s" % (phis_directory, native_training_set_file.split(
#         '/')[-1].split('.')[0], full_parameters_string)
#     phi_summary_file_name = file_prefix + '_phi_native_summary.txt'
#     np.savetxt(phi_summary_file_name, phi_native, fmt='%1.5f')

#     # For decoys, we will use the same function for native sets to do the calculation;
#     decoy_training_set = read_column_from_file(decoy_training_set_file, 1)

#     # Find out how many total phi_i there are and get full parameter string
#     total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string_decoy_structures_provided(
#         phi_list, decoy_training_set)

#     phi_decoy_i_protein = np.zeros((len(decoy_training_set), total_phis))
#     for i_protein, protein in enumerate(decoy_training_set):
#         phi_decoy_i_protein[i_protein] = read_decoy_phi_structures_provided(
#             protein, phi_list, total_phis, jackhmmer=jackhmmer)

#     phi_decoy = np.average(phi_decoy_i_protein, axis=0)

#     # Output to a file;
#     file_prefix = "%s%s_%s" % (phis_directory, decoy_training_set_file.split(
#         '/')[-1].split('.')[0], full_parameters_string)
#     phi_summary_file_name = file_prefix + '_phi_decoy_summary.txt'
#     np.savetxt(phi_summary_file_name, phi_decoy, fmt='%1.5f')

#     A, B, half_B, other_half_B, std_half_B = calculate_A_and_B(
#         phi_decoy, phi_native, total_phis, num_decoys, phi_decoy_i_protein)

#     gamma = np.dot(np.linalg.pinv(B), A)

#     # write gamma file
#     file_prefix = "%s%s_%s" % (gammas_directory, native_training_set_file.split(
#         '/')[-1].split('.')[0], full_parameters_string)

#     gamma_file_name = file_prefix + '_gamma'
# #    gamma_file = open(gamma_file_name, 'w')
#     np.savetxt(gamma_file_name, gamma, '%1.5f')

#     A_file_name = file_prefix + '_A'
# #    A_file = open(A_file_name, 'w')
#     np.savetxt(A_file_name, A, fmt='%1.5f')

#     B_file_name = file_prefix + '_B'
# #    B_file = open(B_file_name, 'w')
#     np.savetxt(B_file_name, B, fmt='%1.5f')
#     #open("%s%s_%s_gamma.dat" % (gammas_directory, training_set_file.split('/')[-1].split('.')[0], full_parameters_string), 'w').write(str(gamma).strip('[]').replace('\n', ' '))

#     if noise_filtering:
#         filtered_gamma, filtered_B, filtered_lamb, P, lamb = get_filtered_gamma_B_lamb_P_and_lamb(
#             A, B, half_B, other_half_B, std_half_B, total_phis, num_decoys)
#         # gamma_file_name = "%sfiltered_%s_%s_gamma.dat" % (gammas_directory, training_set_file.split('/')[-1].split('.')[0], full_parameters_string)
#         # gamma_file = open(gamma_file_name, 'w')
#         filtered_gamma_file_name = file_prefix + '_gamma_filtered'
# #        filtered_gamma_file = open(filtered_gamma_file_name, 'w')
#         np.savetxt(filtered_gamma_file_name, filtered_gamma, fmt='%1.5f')

#         filtered_B_file_name = file_prefix + '_B_filtered'
# #        filtered_B_file = open(filtered_B_file_name, 'w')
#         np.savetxt(filtered_B_file_name, filtered_B, fmt='%1.5f')

#         filtered_lamb_file_name = file_prefix + '_lamb_filtered'
# #        filtered_lamb_file = open(filtered_lamb_file_name, 'w')
#         np.savetxt(filtered_lamb_file_name, filtered_lamb, fmt='%1.5f')

#         P_file_name = file_prefix + '_P'
#         P_file = open(P_file_name, 'w')
#         np.savetxt(P_file, P, fmt='%1.5f')

#         lamb_file_name = file_prefix + '_lamb'
# #        lamb_file = open(lamb_file_name, 'w')
#         np.savetxt(lamb_file_name, lamb, fmt='%1.5f')

#         # open("%sfiltered_%s_%s_gamma.dat" % (gammas_directory, training_set_file.split('/')[-1].split('.')[0], full_parameters_string), 'w').write(str(filtered_gamma).strip('[]').replace('\n', ' '))

#     if noise_filtering:
#         return
#         # return A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb
#     else:
#         return A, B, gamma
def calculate_A_B_and_gamma_parallel(training_set_file, phi_list_file_name, decoy_method, num_decoys,
                                        noise_filtering=True, jackhmmer=False, read=0, subset=None,
                                        subset_index=-1):
    phi_list = read_phi_list(phi_list_file_name)
    training_set = read_column_from_file(training_set_file, 1)
    print(len(training_set))
    # Find out how many total phi_i there are and get full parameter string
    total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string(
        phi_list, training_set)
    print(num_phis)
    if read == 0:  # only compute when I need A
        pass
    elif read == 1:
        print("reading native")
        os.system("echo 'Reading native' >> log")
        file_prefix = "%s%s_%s" % (phis_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)
        phi_summary_file_name = file_prefix + '_phi_native_summary.txt'
        phi_native = np.loadtxt(phi_summary_file_name)
    elif read == 2:
        phi_native_i_protein = np.zeros((len(training_set), total_phis))
        for i_protein, protein in enumerate(training_set):
            phi_native_i_protein[i_protein] = read_native_phi(
                protein, phi_list, total_phis, jackhmmer=jackhmmer)

        phi_native = np.average(phi_native_i_protein, axis=0)

        # Output to a file;
        file_prefix = "%s%s_%s" % (phis_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)
        phi_summary_file_name = file_prefix + '_phi_native_summary.txt'
        np.savetxt(phi_summary_file_name, phi_native, fmt='%1.5f')

    if read == 0:
        pass
    elif read == 1:
        print("Reading phi decoy")
        os.system("echo 'Reading phi decoy' >> log")
        # Output to a file;
        file_prefix = "%s%s_%s" % (phis_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)
        phi_summary_file_name = file_prefix + '_phi_decoy_summary.txt'
        average_phi_decoy = np.loadtxt(phi_summary_file_name)
        # phi_all_summary_file_name = file_prefix + '_phi_decoy_all_summary.txt'
        # phi_i_protein_i_decoy = np.zeros(
        #     (len(training_set), num_decoys, total_phis))
        # # phi_i_protein_i_decoy = np.loadtxt(phi_all_summary_file_name)
        # phi_i_protein_i_decoy = np.reshape(phi_i_protein_i_decoy, (len(training_set), num_decoys, total_phis))
    elif read == 2:
        phi_i_protein_i_decoy = np.zeros(
            (len(training_set), num_decoys, total_phis))

        for i_protein, protein in enumerate(training_set):
            phi_i_protein_i_decoy[i_protein] = read_decoy_phis(
                protein, phi_list, total_phis, num_phis, num_decoys, decoy_method, jackhmmer=jackhmmer)

        # The phi_i decoy is constructed as the union of all decoys of all proteins in the training set;
        phi_i_decoy_reshaped = np.reshape(phi_i_protein_i_decoy,
                                            (len(training_set) * num_decoys, total_phis))
        average_phi_decoy = np.average(phi_i_decoy_reshaped, axis=0)

        # Output to a file;
        file_prefix = "%s%s_%s" % (phis_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)
        phi_summary_file_name = file_prefix + '_phi_decoy_summary.txt'
        np.savetxt(phi_summary_file_name, average_phi_decoy, fmt='%1.5f')
        exit()
        # phi_all_summary_file_name = file_prefix + '_phi_decoy_all_summary.txt'
        # np.savetxt(phi_all_summary_file_name, phi_i_decoy_reshaped, fmt='%1.5f')
    # A, B, half_B, other_half_B, std_half_B = calculate_A_and_B(
    #     average_phi_decoy, phi_native, total_phis, num_decoys, phi_i_decoy_reshaped)
    print("done reading")
    os.system("echo 'Done reading' >> log")
    if read == 0:
        pass
    else:
        A = average_phi_decoy - phi_native
    half_B = np.zeros((total_phis, total_phis))
    std_half_B = np.zeros((total_phis, total_phis))
    other_half_B = np.zeros((total_phis, total_phis))
    if subset:
        for i_protein, protein in enumerate(subset):
            os.system(f"echo '{i_protein}, {protein}' >> log")
            print(f"{i_protein}, {protein}")
            one_half_B = np.zeros((total_phis, total_phis))
            one_std_half_B = np.zeros((total_phis, total_phis))
            one_other_half_B = np.zeros((total_phis, total_phis))

            one_phi_i_protein_i_decoy = read_decoy_phis(
                protein, phi_list, total_phis, num_phis, num_decoys, decoy_method, jackhmmer=jackhmmer)

            # phis_i = one_phi_i_protein_i_decoy.reshape(num_decoys, total_phis, 1)
            # phis_j = one_phi_i_protein_i_decoy.reshape(num_decoys, 1, total_phis)
            # one_half_B = np.average(phis_i * phis_j, axis=0)
            # one_std_half_B = np.std(phis_i * phis_j, axis=0)
            phis_i = one_phi_i_protein_i_decoy.reshape(num_decoys, total_phis, 1)
            for j in range(total_phis):
                phis_j = phis_i[:, j].reshape(num_decoys, 1, 1)
                one_half_B[j] += np.average(phis_i * phis_j, axis=0).reshape(total_phis)
                one_std_half_B[j] += np.std(phis_i * phis_j, axis=0).reshape(total_phis)

            average_phi = np.average(one_phi_i_protein_i_decoy, axis=0)
            one_other_half_B += average_phi.reshape(total_phis, 1) * average_phi.reshape(1, total_phis)
            # write gamma file
            file_prefix = "subset%s%s_%s" % (gammas_directory, training_set_file.split(
                '/')[-1].split('.')[0], full_parameters_string)
            half_B_file_name = file_prefix + f'_half_B_{protein}'
            np.savetxt(half_B_file_name, one_half_B, '%1.5f')

            other_half_B_file_name = file_prefix + f'_other_half_B_{protein}'
            np.savetxt(other_half_B_file_name, one_other_half_B, '%1.5f')

            std_half_B_file_name = file_prefix + f'_std_half_B_{protein}'
            np.savetxt(std_half_B_file_name, one_std_half_B, '%1.5f')

            half_B += one_half_B
            std_half_B += one_std_half_B
            other_half_B += one_other_half_B
        # write gamma file
        file_prefix = "subset%s%s_%s" % (gammas_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)
        # print("------")
        # print(file_prefix)
        # print("------")
        half_B_file_name = file_prefix + f'_half_B_{subset_index}'
        np.savetxt(half_B_file_name, half_B, '%1.5f')

        other_half_B_file_name = file_prefix + f'_other_half_B_{subset_index}'
        np.savetxt(other_half_B_file_name, other_half_B, '%1.5f')

        std_half_B_file_name = file_prefix + f'_std_half_B_{subset_index}'
        np.savetxt(std_half_B_file_name, std_half_B, '%1.5f')
        exit()
    else:
        # for i_protein, protein in enumerate(training_set):
        # for i_protein in range(475):
        for ii in range(712):
            # os.system(f"echo 'Reading {i_protein}, {protein}' >> log")
            i_protein = f"proteins_name_list_{ii}"
            os.system(f"echo 'Reading {i_protein}' >> log")

            print(i_protein)
            file_prefix = "subset%s%s_%s" % (gammas_directory, training_set_file.split(
                '/')[-1].split('.')[0], full_parameters_string)

            half_B_file_name = file_prefix + f'_half_B_{i_protein}'
            one_half_B = np.loadtxt(half_B_file_name)

            std_half_B_file_name = file_prefix + f'_std_half_B_{i_protein}'
            one_std_half_B = np.loadtxt(std_half_B_file_name)

            other_half_B_file_name = file_prefix + f'_other_half_B_{i_protein}'
            one_other_half_B = np.loadtxt(other_half_B_file_name)

            half_B += one_half_B
            std_half_B += one_std_half_B
            other_half_B += one_other_half_B
    size_of_training_set = len(training_set)
    half_B /= size_of_training_set
    std_half_B /= size_of_training_set
    other_half_B /= size_of_training_set
    B = half_B - other_half_B
    # A, B, half_B, other_half_B, std_half_B = calculate_A_and_B_wei_parallel(
    #     average_phi_decoy, phi_native, phi_i_protein_i_decoy)

    gamma = np.dot(np.linalg.pinv(B), A)

    # write gamma file
    file_prefix = "%s%s_%s" % (gammas_directory, training_set_file.split(
        '/')[-1].split('.')[0], full_parameters_string)

    gamma_file_name = file_prefix + '_gamma'
#    gamma_file = open(gamma_file_name, 'w')
    np.savetxt(gamma_file_name, gamma, '%1.5f')

    A_file_name = file_prefix + '_A'
#    A_file = open(A_file_name, 'w')
    np.savetxt(A_file_name, A, fmt='%1.5f')

    A_prime_file_name = file_prefix + '_A_prime'
#    A_file = open(A_file_name, 'w')
    np.savetxt(A_prime_file_name, average_phi_decoy, fmt='%1.5f')

    B_file_name = file_prefix + '_B'
#    B_file = open(B_file_name, 'w')
    np.savetxt(B_file_name, B, fmt='%1.5f')

    half_B_file_name = file_prefix + '_half_B'
#    gamma_file = open(gamma_file_name, 'w')
    np.savetxt(half_B_file_name, half_B, '%1.5f')

    other_half_B_file_name = file_prefix + '_other_half_B'
#    gamma_file = open(gamma_file_name, 'w')
    np.savetxt(other_half_B_file_name, other_half_B, '%1.5f')

    std_half_B_file_name = file_prefix + '_std_half_B'
#    gamma_file = open(gamma_file_name, 'w')
    np.savetxt(std_half_B_file_name, std_half_B, '%1.5f')

    #open("%s%s_%s_gamma.dat" % (gammas_directory, training_set_file.split('/')[-1].split('.')[0], full_parameters_string), 'w').write(str(gamma).strip('[]').replace('\n', ' '))

    if noise_filtering:
        filtered_gamma, filtered_B, filtered_lamb, P, lamb = get_filtered_gamma_B_lamb_P_and_lamb(
            A, B, half_B, other_half_B, std_half_B, total_phis, num_decoys)
        # gamma_file_name = "%sfiltered_%s_%s_gamma.dat" % (gammas_directory, training_set_file.split('/')[-1].split('.')[0], full_parameters_string)
        # gamma_file = open(gamma_file_name, 'w')
        filtered_gamma_file_name = file_prefix + '_gamma_filtered'
#        filtered_gamma_file = open(filtered_gamma_file_name, 'w')
        np.savetxt(filtered_gamma_file_name, filtered_gamma, fmt='%1.5f')

        filtered_B_file_name = file_prefix + '_B_filtered'
#        filtered_B_file = open(filtered_B_file_name, 'w')
        np.savetxt(filtered_B_file_name, filtered_B, fmt='%1.5f')

        filtered_lamb_file_name = file_prefix + '_lamb_filtered'
#        filtered_lamb_file = open(filtered_lamb_file_name, 'w')
        np.savetxt(filtered_lamb_file_name, filtered_lamb, fmt='%1.5f')

        P_file_name = file_prefix + '_P'
        # print(P)
        # P_file = open(P_file_name, 'wb')
        # np.savetxt(P_file, P, fmt='%1.5f')
        np.savetxt(P_file_name, P, fmt='%1.5f')

        lamb_file_name = file_prefix + '_lamb'
#        lamb_file = open(lamb_file_name, 'w')
        np.savetxt(lamb_file_name, lamb, fmt='%1.5f')

        # open("%sfiltered_%s_%s_gamma.dat" % (gammas_directory, training_set_file.split('/')[-1].split('.')[0], full_parameters_string), 'w').write(str(filtered_gamma).strip('[]').replace('\n', ' '))

    if noise_filtering:
        # return
        return A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb
    else:
        return A, B, gamma

def calculate_A_B_and_gamma_wl45_parallel(training_set_file, phi_list_file_name, decoy_method, num_decoys,
                                            noise_filtering=True, jackhmmer=False, read=True, withBiased=False, **kwargs):
    phi_list = read_phi_list(phi_list_file_name)
    training_set = read_column_from_file(training_set_file, 1)
    print("Size of training set", len(training_set))
    os.system(f"echo 'Size {len(training_set)}' >> log")
    # Find out how many total phi_i there are and get full parameter string
    total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string(
        phi_list, training_set, **kwargs)
    print("Size of phis", total_phis)
    if read:
        print("reading native")
        os.system("echo 'Reading native' >> log")
        file_prefix = "%s%s_%s" % (phis_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)
        phi_summary_file_name = file_prefix + '_phi_native_summary.txt'
        phi_native = np.loadtxt(phi_summary_file_name)
    else:
        phi_native_i_protein = np.zeros((len(training_set), total_phis))
        for i_protein, protein in enumerate(training_set):
            phi_native_i_protein[i_protein] = read_native_phi(
                protein, phi_list, total_phis, jackhmmer=jackhmmer, **kwargs)
        phi_native = np.average(phi_native_i_protein, axis=0)

    if read:
        print("Reading phi decoy")
        os.system("echo 'Reading phi decoy' >> log")
        # Output to a file;
        file_prefix = "%s%s_%s" % (phis_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)
        phi_summary_file_name = file_prefix + '_phi_decoy_summary.txt'
        average_phi_decoy = np.loadtxt(phi_summary_file_name)
        phi_all_summary_file_name = file_prefix + '_phi_decoy_all_summary.txt'
        phi_i_protein_i_decoy = np.zeros(
            (len(training_set), num_decoys, total_phis))
        # phi_i_protein_i_decoy = np.loadtxt(phi_all_summary_file_name)
        phi_i_protein_i_decoy = np.reshape(phi_i_protein_i_decoy, (len(training_set), num_decoys, total_phis))
    else:
        phi_i_protein_i_decoy = np.zeros(
            (len(training_set), num_decoys, total_phis))

        for i_protein, protein in enumerate(training_set):
            phi_i_protein_i_decoy[i_protein] = read_decoy_phis(
                protein, phi_list, total_phis, num_phis, num_decoys, decoy_method, jackhmmer=jackhmmer, **kwargs)
        if withBiased:
            phi_i_protein_i_decoyQ = np.zeros(
                (len(training_set), num_decoys, 1))
            for i_protein, protein in enumerate(training_set):
                phi_i_protein_i_decoyQ[i_protein] = read_decoyQ_phis(
                    protein, phi_list, total_phis, num_phis, num_decoys, decoy_method, jackhmmer=jackhmmer, **kwargs)
            phi_i_protein_i_decoy *= 1 - phi_i_protein_i_decoyQ
            normalization = np.sum(1 - phi_i_protein_i_decoyQ)
        else:
            normalization = len(training_set) * num_decoys
        # The phi_i decoy is constructed as the union of all decoys of all proteins in the training set;
        phi_i_decoy_reshaped = np.reshape(phi_i_protein_i_decoy,
                                            (len(training_set) * num_decoys, total_phis))
        # average_phi_decoy = np.average(phi_i_decoy_reshaped, axis=0)
        average_phi_decoy = np.sum(phi_i_decoy_reshaped, axis=0) / normalization
        phi_i_decoy_reshaped = None
    # A, B, half_B, other_half_B, std_half_B = calculate_A_and_B(
    #     average_phi_decoy, phi_native, total_phis, num_decoys, phi_i_decoy_reshaped)
    print("done reading")
    # os.system("echo 'Done reading' >> log")

    A = average_phi_decoy - phi_native
    size_of_training_set, num_decoys, total_phis = phi_i_protein_i_decoy.shape
    half_B = np.zeros((total_phis, total_phis))
    std_half_B = np.zeros((total_phis, total_phis))
    other_half_B = np.zeros((total_phis, total_phis))

    for p in range(size_of_training_set):
        phis_i = phi_i_protein_i_decoy[p].reshape(num_decoys, total_phis, 1)
        for j in range(total_phis):
            phis_j = phis_i[:, j].reshape(num_decoys, 1, 1)
            half_B[j] += np.average(phis_i * phis_j, axis=0).reshape(total_phis)
            std_half_B[j] += np.std(phis_i * phis_j, axis=0).reshape(total_phis)

    # half_B /= size_of_training_set
    # std_half_B /= size_of_training_set

    for p in range(size_of_training_set):
        average_phi = np.average(phi_i_protein_i_decoy[p], axis=0)
        other_half_B += average_phi.reshape(total_phis, 1) * average_phi.reshape(1, total_phis)

    # other_half_B /= size_of_training_set
    # B = half_B - other_half_B

    phi_i_protein_i_decoy = None  # free its memory

    # write gamma file
    file_prefix = "%s%s_%s" % (sub_gammas_directory, training_set_file.split(
        '/')[-1].split('.')[0], full_parameters_string)


    A_file_name = file_prefix + '_A'
#    A_file = open(A_file_name, 'w')
    np.savetxt(A_file_name, A, fmt='%1.5f')

    A_prime_file_name = file_prefix + '_A_prime'
#    A_file = open(A_file_name, 'w')
    np.savetxt(A_prime_file_name, average_phi_decoy, fmt='%1.5f')


    half_B_file_name = file_prefix + '_half_B'
#    gamma_file = open(gamma_file_name, 'w')
    np.savetxt(half_B_file_name, half_B, '%1.5f')

    other_half_B_file_name = file_prefix + '_other_half_B'
#    gamma_file = open(gamma_file_name, 'w')
    np.savetxt(other_half_B_file_name, other_half_B, '%1.5f')

    std_half_B_file_name = file_prefix + '_std_half_B'
#    gamma_file = open(gamma_file_name, 'w')
    np.savetxt(std_half_B_file_name, std_half_B, '%1.5f')

    # open("%s%s_%s_gamma.dat" % (gammas_directory, training_set_file.split('/')[-1].split('.')[0], full_parameters_string), 'w').write(str(gamma).strip('[]').replace('\n', ' '))
    return True


def calculate_A_B_and_gamma_wl45(training_set_file, phi_list_file_name, decoy_method, num_decoys,
                                    noise_filtering=True, jackhmmer=False, read=True, withBiased=False, **kwargs):
    phi_list = read_phi_list(phi_list_file_name)
    training_set = read_column_from_file(training_set_file, 1)
    print("Size of training set", len(training_set))
    os.system(f"echo 'Size {len(training_set)}' >> log")
    # Find out how many total phi_i there are and get full parameter string
    total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string(
        phi_list, training_set, **kwargs)
    print("Size of phis", total_phis)
    if read:
        print("reading native")
        os.system("echo 'Reading native' >> log")
        file_prefix = "%s%s_%s" % (phis_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)
        phi_summary_file_name = file_prefix + '_phi_native_summary.txt'
        phi_native = np.loadtxt(phi_summary_file_name)
    else:
        phi_native_i_protein = np.zeros((len(training_set), total_phis))
        for i_protein, protein in enumerate(training_set):
            phi_native_i_protein[i_protein] = read_native_phi(
                protein, phi_list, total_phis, jackhmmer=jackhmmer, **kwargs)
        phi_native = np.average(phi_native_i_protein, axis=0)

        # Output to a file;
        file_prefix = "%s%s_%s" % (phis_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)
        phi_summary_file_name = file_prefix + '_phi_native_summary.txt'
        np.savetxt(phi_summary_file_name, phi_native, fmt='%1.5f')

    if read:
        print("Reading phi decoy")
        os.system("echo 'Reading phi decoy' >> log")
        # Output to a file;
        file_prefix = "%s%s_%s" % (phis_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)
        phi_summary_file_name = file_prefix + '_phi_decoy_summary.txt'
        average_phi_decoy = np.loadtxt(phi_summary_file_name)
        phi_all_summary_file_name = file_prefix + '_phi_decoy_all_summary.txt'
        phi_i_protein_i_decoy = np.zeros(
            (len(training_set), num_decoys, total_phis))
        # phi_i_protein_i_decoy = np.loadtxt(phi_all_summary_file_name)
        phi_i_protein_i_decoy = np.reshape(phi_i_protein_i_decoy, (len(training_set), num_decoys, total_phis))
    else:
        phi_i_protein_i_decoy = np.zeros(
            (len(training_set), num_decoys, total_phis))

        for i_protein, protein in enumerate(training_set):
            phi_i_protein_i_decoy[i_protein] = read_decoy_phis(
                protein, phi_list, total_phis, num_phis, num_decoys, decoy_method, jackhmmer=jackhmmer, **kwargs)
        if withBiased:
            phi_i_protein_i_decoyQ = np.zeros(
                (len(training_set), num_decoys, 1))
            for i_protein, protein in enumerate(training_set):
                phi_i_protein_i_decoyQ[i_protein] = read_decoyQ_phis(
                    protein, phi_list, total_phis, num_phis, num_decoys, decoy_method, jackhmmer=jackhmmer, **kwargs)
            phi_i_protein_i_decoy *= 1 - phi_i_protein_i_decoyQ
            normalization = np.sum(1 - phi_i_protein_i_decoyQ)
        else:
            normalization = len(training_set) * num_decoys
        # The phi_i decoy is constructed as the union of all decoys of all proteins in the training set;
        phi_i_decoy_reshaped = np.reshape(phi_i_protein_i_decoy,
                                            (len(training_set) * num_decoys, total_phis))
        # average_phi_decoy = np.average(phi_i_decoy_reshaped, axis=0)
        # print("normalization", normalization)
        # print("phi_i_decoy_reshaped", phi_i_decoy_reshaped)

        average_phi_decoy = np.sum(phi_i_decoy_reshaped, axis=0) / normalization
        # print("average_phi_decoy", average_phi_decoy)
        # Output to a file;
        file_prefix = "%s%s_%s" % (phis_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)
        phi_summary_file_name = file_prefix + '_phi_decoy_summary.txt'
        np.savetxt(phi_summary_file_name, average_phi_decoy, fmt='%1.5f')

        phi_all_summary_file_name = file_prefix + '_phi_decoy_all_summary.txt'
        np.savetxt(phi_all_summary_file_name, phi_i_decoy_reshaped, fmt='%1.5f')
        phi_i_decoy_reshaped = None
    # A, B, half_B, other_half_B, std_half_B = calculate_A_and_B(
    #     average_phi_decoy, phi_native, total_phis, num_decoys, phi_i_decoy_reshaped)
    print("done reading")
    os.system("echo 'Done reading' >> log")
    A, B, half_B, other_half_B, std_half_B = calculate_A_and_B_wei(
        average_phi_decoy, phi_native, phi_i_protein_i_decoy)
    phi_i_protein_i_decoy = None  # free its memory
    gamma = np.dot(np.linalg.pinv(B), A)

    # write gamma file
    file_prefix = "%s%s_%s" % (gammas_directory, training_set_file.split(
        '/')[-1].split('.')[0], full_parameters_string)

    gamma_file_name = file_prefix + '_gamma'
#    gamma_file = open(gamma_file_name, 'w')
    np.savetxt(gamma_file_name, gamma, '%1.5f')

    A_file_name = file_prefix + '_A'
#    A_file = open(A_file_name, 'w')
    np.savetxt(A_file_name, A, fmt='%1.5f')

    A_prime_file_name = file_prefix + '_A_prime'
#    A_file = open(A_file_name, 'w')
    np.savetxt(A_prime_file_name, average_phi_decoy, fmt='%1.5f')

    B_file_name = file_prefix + '_B'
#    B_file = open(B_file_name, 'w')
    np.savetxt(B_file_name, B, fmt='%1.5f')

    half_B_file_name = file_prefix + '_half_B'
#    gamma_file = open(gamma_file_name, 'w')
    np.savetxt(half_B_file_name, half_B, '%1.5f')

    other_half_B_file_name = file_prefix + '_other_half_B'
#    gamma_file = open(gamma_file_name, 'w')
    np.savetxt(other_half_B_file_name, other_half_B, '%1.5f')

    std_half_B_file_name = file_prefix + '_std_half_B'
#    gamma_file = open(gamma_file_name, 'w')
    np.savetxt(std_half_B_file_name, std_half_B, '%1.5f')

    # open("%s%s_%s_gamma.dat" % (gammas_directory, training_set_file.split('/')[-1].split('.')[0], full_parameters_string), 'w').write(str(gamma).strip('[]').replace('\n', ' '))

    if noise_filtering:
        filtered_gamma, filtered_B, filtered_lamb, P, lamb = get_filtered_gamma_B_lamb_P_and_lamb(
            A, B, half_B, other_half_B, std_half_B, total_phis, num_decoys)
        # gamma_file_name = "%sfiltered_%s_%s_gamma.dat" % (gammas_directory, training_set_file.split('/')[-1].split('.')[0], full_parameters_string)
        # gamma_file = open(gamma_file_name, 'w')

        filtered_gamma_file_name = file_prefix + '_gamma_filtered'
#        filtered_gamma_file = open(filtered_gamma_file_name, 'w')
        np.savetxt(filtered_gamma_file_name, filtered_gamma, fmt='%1.5f')

        filtered_B_file_name = file_prefix + '_B_filtered'
#        filtered_B_file = open(filtered_B_file_name, 'w')
        np.savetxt(filtered_B_file_name, filtered_B, fmt='%1.5f')

        filtered_lamb_file_name = file_prefix + '_lamb_filtered'
#        filtered_lamb_file = open(filtered_lamb_file_name, 'w')
        np.savetxt(filtered_lamb_file_name, filtered_lamb, fmt='%1.5f')

        P_file_name = file_prefix + '_P'
        # print(P)
        # P_file = open(P_file_name, 'wb')
        # np.savetxt(P_file, P, fmt='%1.5f')
        np.savetxt(P_file_name, P, fmt='%1.5f')

        lamb_file_name = file_prefix + '_lamb'
#        lamb_file = open(lamb_file_name, 'w')
        np.savetxt(lamb_file_name, lamb, fmt='%1.5f')

        # open("%sfiltered_%s_%s_gamma.dat" % (gammas_directory, training_set_file.split('/')[-1].split('.')[0], full_parameters_string), 'w').write(str(filtered_gamma).strip('[]').replace('\n', ' '))

    if noise_filtering:
        # return
        return A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb
    else:
        return A, B, gamma

def calculate_A_B_and_gamma_xl23(training_set_file, phi_list_file_name, decoy_method, num_decoys,
                                    noise_filtering=True, jackhmmer=False, read=True):
    phi_list = read_phi_list(phi_list_file_name)
    training_set = read_column_from_file(training_set_file, 1)
    print(len(training_set))
    os.system(f"echo 'Size {len(training_set)}' >> log")
    # Find out how many total phi_i there are and get full parameter string
    total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string(
        phi_list, training_set)
    # print(num_phis)
    if read:
        print("reading native")
        os.system("echo 'Reading native' >> log")
        file_prefix = "%s%s_%s" % (phis_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)
        phi_summary_file_name = file_prefix + '_phi_native_summary.txt'
        phi_native = np.loadtxt(phi_summary_file_name)
    else:
        phi_native_i_protein = np.zeros((len(training_set), total_phis))
        for i_protein, protein in enumerate(training_set):
            phi_native_i_protein[i_protein] = read_native_phi(
                protein, phi_list, total_phis, jackhmmer=jackhmmer)

        phi_native = np.average(phi_native_i_protein, axis=0)

        # Output to a file;
        file_prefix = "%s%s_%s" % (phis_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)
        phi_summary_file_name = file_prefix + '_phi_native_summary.txt'
        np.savetxt(phi_summary_file_name, phi_native, fmt='%1.5f')

    if read:
        print("Reading phi decoy")
        os.system("echo 'Reading phi decoy' >> log")
        # Output to a file;
        file_prefix = "%s%s_%s" % (phis_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)
        phi_summary_file_name = file_prefix + '_phi_decoy_summary.txt'
        average_phi_decoy = np.loadtxt(phi_summary_file_name)
        phi_all_summary_file_name = file_prefix + '_phi_decoy_all_summary.txt'
        phi_i_protein_i_decoy = np.zeros(
            (len(training_set), num_decoys, total_phis))
        # phi_i_protein_i_decoy = np.loadtxt(phi_all_summary_file_name)
        phi_i_protein_i_decoy = np.reshape(phi_i_protein_i_decoy, (len(training_set), num_decoys, total_phis))
    else:
        phi_i_protein_i_decoy = np.zeros(
            (len(training_set), num_decoys, total_phis))

        for i_protein, protein in enumerate(training_set):
            phi_i_protein_i_decoy[i_protein] = read_decoy_phis(
                protein, phi_list, total_phis, num_phis, num_decoys, decoy_method, jackhmmer=jackhmmer)

        # The phi_i decoy is constructed as the union of all decoys of all proteins in the training set;
        phi_i_decoy_reshaped = np.reshape(phi_i_protein_i_decoy,
                                            (len(training_set) * num_decoys, total_phis))
        average_phi_decoy = np.average(phi_i_decoy_reshaped, axis=0)

        # Output to a file;
        file_prefix = "%s%s_%s" % (phis_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)
        phi_summary_file_name = file_prefix + '_phi_decoy_summary.txt'
        np.savetxt(phi_summary_file_name, average_phi_decoy, fmt='%1.5f')

        phi_all_summary_file_name = file_prefix + '_phi_decoy_all_summary.txt'
        np.savetxt(phi_all_summary_file_name, phi_i_decoy_reshaped, fmt='%1.5f')
        phi_i_decoy_reshaped = None
    # A, B, half_B, other_half_B, std_half_B = calculate_A_and_B(
    #     average_phi_decoy, phi_native, total_phis, num_decoys, phi_i_decoy_reshaped)
    print("done reading")
    os.system("echo 'Done reading' >> log")
    A, B, half_B, other_half_B, std_half_B = calculate_A_and_B_wei(
        average_phi_decoy, phi_native, phi_i_protein_i_decoy)
    phi_i_protein_i_decoy = None  # free its memory
    gamma = np.dot(np.linalg.pinv(B), A)

    # write gamma file
    file_prefix = "%s%s_%s" % (gammas_directory, training_set_file.split(
        '/')[-1].split('.')[0], full_parameters_string)

    gamma_file_name = file_prefix + '_gamma'
#    gamma_file = open(gamma_file_name, 'w')
    np.savetxt(gamma_file_name, gamma, '%1.5f')

    A_file_name = file_prefix + '_A'
#    A_file = open(A_file_name, 'w')
    np.savetxt(A_file_name, A, fmt='%1.5f')

    B_file_name = file_prefix + '_B'
#    B_file = open(B_file_name, 'w')
    np.savetxt(B_file_name, B, fmt='%1.5f')

    half_B_file_name = file_prefix + '_half_B'
#    gamma_file = open(gamma_file_name, 'w')
    np.savetxt(half_B_file_name, half_B, '%1.5f')

    other_half_B_file_name = file_prefix + '_other_half_B'
#    gamma_file = open(gamma_file_name, 'w')
    np.savetxt(other_half_B_file_name, other_half_B, '%1.5f')

    std_half_B_file_name = file_prefix + '_std_half_B'
#    gamma_file = open(gamma_file_name, 'w')
    np.savetxt(std_half_B_file_name, std_half_B, '%1.5f')

    #open("%s%s_%s_gamma.dat" % (gammas_directory, training_set_file.split('/')[-1].split('.')[0], full_parameters_string), 'w').write(str(gamma).strip('[]').replace('\n', ' '))

    if noise_filtering:
        filtered_gamma, filtered_B, filtered_lamb, P, lamb = get_filtered_gamma_B_lamb_P_and_lamb(
            A, B, half_B, other_half_B, std_half_B, total_phis, num_decoys)
        # gamma_file_name = "%sfiltered_%s_%s_gamma.dat" % (gammas_directory, training_set_file.split('/')[-1].split('.')[0], full_parameters_string)
        # gamma_file = open(gamma_file_name, 'w')
        filtered_gamma_file_name = file_prefix + '_gamma_filtered'
#        filtered_gamma_file = open(filtered_gamma_file_name, 'w')
        np.savetxt(filtered_gamma_file_name, filtered_gamma, fmt='%1.5f')

        filtered_B_file_name = file_prefix + '_B_filtered'
#        filtered_B_file = open(filtered_B_file_name, 'w')
        np.savetxt(filtered_B_file_name, filtered_B, fmt='%1.5f')

        filtered_lamb_file_name = file_prefix + '_lamb_filtered'
#        filtered_lamb_file = open(filtered_lamb_file_name, 'w')
        np.savetxt(filtered_lamb_file_name, filtered_lamb, fmt='%1.5f')

        P_file_name = file_prefix + '_P'
        # print(P)
        # P_file = open(P_file_name, 'wb')
        # np.savetxt(P_file, P, fmt='%1.5f')
        np.savetxt(P_file_name, P, fmt='%1.5f')

        lamb_file_name = file_prefix + '_lamb'
#        lamb_file = open(lamb_file_name, 'w')
        np.savetxt(lamb_file_name, lamb, fmt='%1.5f')

        # open("%sfiltered_%s_%s_gamma.dat" % (gammas_directory, training_set_file.split('/')[-1].split('.')[0], full_parameters_string), 'w').write(str(filtered_gamma).strip('[]').replace('\n', ' '))

    if noise_filtering:
        # return
        return A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb
    else:
        return A, B, gamma



def get_filtered_gamma_B_lamb_P_and_lamb(A, B, half_B, other_half_B, std_half_B, total_phis, num_decoys, noise_iterations=10, relative_error_threshold=0.5, mode=2, setCutoff=None):
    lamb, P = np.linalg.eig(B)
    lamb, P = sort_eigenvalues_and_eigenvectors(lamb, P)
    if not setCutoff:
        cutoff_modes = []
        for i_noise in range(noise_iterations):
            noisy_B = np.zeros((total_phis, total_phis))
            for i in range(total_phis):
                for j in range(i, total_phis):
                    if mode == 1:
                        random_B_ij = np.random.normal(
                            loc=half_B[i][j], scale=std_half_B[i][j] / float(num_decoys))
                    elif mode == 2:
                        random_B_ij = np.random.normal(
                            loc=half_B[i][j], scale=std_half_B[i][j] / float(num_decoys)**0.5)
                    # random_B_ij = np.random.normal(
                    #     loc=half_B[i][j], scale=std_half_B[i][j] / float(num_decoys))
                    # random_B_ij = np.random.normal(
                    #     loc=half_B[i][j], scale=std_half_B[i][j])
                    noisy_B[i][j] = noisy_B[j][i] = random_B_ij - \
                        other_half_B[i][j]

            noisy_lamb, noisy_P = np.linalg.eig(noisy_B)
            noisy_lamb, noisy_P = sort_eigenvalues_and_eigenvectors(
                noisy_lamb, noisy_P)

            try:
                cutoff_mode = np.where(
                    np.abs(lamb - noisy_lamb) / lamb > relative_error_threshold)[0][0]
            except IndexError:
                cutoff_mode = len(lamb)
            cutoff_modes.append(cutoff_mode)

        cutoff_mode = min(cutoff_modes)
    else:
        cutoff_mode = setCutoff
    print(cutoff_mode)

    filtered_lamb = np.copy(lamb)
    filtered_B_inv, filtered_lamb, P = get_filtered_B_inv_lambda_and_P(
        filtered_lamb, cutoff_mode, P)

    filtered_gamma = np.dot(filtered_B_inv, A)
    filtered_B = np.linalg.inv(filtered_B_inv)
    return filtered_gamma, filtered_B, filtered_lamb, P, lamb


def get_filtered_B_inv_lambda_and_P(filtered_lamb, cutoff_mode, P, method='extend_all_after_first_noisy_mode'):
    if method == 'zero_all_after_first_noisy_mode':
        filtered_lamb_inv = 1 / filtered_lamb
        # for "zeroing unreliable eigenvalues"
        filtered_lamb_inv[cutoff_mode:] = 0.0
        filtered_B_inv = np.dot(
            P, np.dot(np.diag(filtered_lamb_inv), np.linalg.inv(P)))
        filtered_lamb = 1 / filtered_lamb_inv
    if method == 'extend_all_after_first_noisy_mode':
        # for "extending lowest reliable eigenvalue"
        filtered_lamb[cutoff_mode:] = filtered_lamb[cutoff_mode - 1]
        filtered_B_inv = np.dot(
            P, np.dot(np.diag(1 / filtered_lamb), np.linalg.inv(P)))

    return filtered_B_inv, filtered_lamb, P


def sort_eigenvalues_and_eigenvectors(eigenvalues, eigenvectors):
    idx = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    return eigenvalues, eigenvectors


def add_virtual_glycine_to_residue(residue):
        # get atom coordinates as vectors
    n = residue['N'].get_vector()
    c = residue['C'].get_vector()
    ca = residue['CA'].get_vector()
    # center at origin
    n = n - ca
    c = c - ca
    # find rotation matrix that rotates n -120 degrees along the ca-c vector
    rot = rotaxis(-np.pi * 120.0 / 180.0, c)
    # apply rotation to ca-n vector
    cb_at_origin = n.left_multiply(rot)
    # put on top of ca atom
    cb = cb_at_origin + ca
    atom = Atom.Atom("CB", cb, 0, 1, " ", " CB ", 0, element="CB")
    residue.add(atom)
    return residue


def add_virtual_glycines_2(structure):
    residues = get_res_list(structure)
    for residue in residues:
        try:
            if residue["CB"] is not None:
                pass
        except KeyError:
            residue = add_virtual_glycine_to_residue(residue)
            # print(residue)
    return structure

def add_virtual_glycines(structure):
    residues = get_res_list(structure)
    for residue in residues:
        if residue.get_resname() == "GLY":
            residue = add_virtual_glycine_to_residue(residue)

    return structure


def save_structure(structure, file_name):
    io = PDBIO()
    io.set_structure(structure)
    io.save(file_name)


def add_virtual_glycines_list(proteins_list_file_name):
    proteins_list = read_column_from_file(proteins_list_file_name, 1)
    error_list_file = open("key_errors.dat", 'a')
    for protein in proteins_list:
        structure = parse_pdb(protein)
        try:
            structure = add_virtual_glycines_2(structure)
        except KeyError:
            error_list_file.write("%s\n" % protein)
            continue
        save_structure(structure, protein + '.pdb')

def evaluate_hamiltonian_native_structures_provided(protein, hamiltonian, native_training_set_file, training_decoy_method, test_decoy_method, use_filtered_gammas=True):
    phi_list = read_phi_list(hamiltonian)
    native_training_set = read_column_from_file(native_training_set_file, 1)

    # read in Hamiltonian;
    # Find out how many total phi_i there are and get full parameter string
    total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string(
        phi_list, native_training_set)

    # read in corresponding gammas
    if use_filtered_gammas:
        gamma_file_name = "%s%s_%s_gamma_filtered" % (
            gammas_directory, native_training_set_file.split('/')[-1].split('.')[0], full_parameters_string)
    else:
        gamma_file_name = "%s%s_%s_gamma" % (gammas_directory, native_training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)

    # Need to filter out the complex number if in the "filtered" mode;
    if use_filtered_gammas:
        gamma = np.loadtxt(gamma_file_name, dtype=complex, converters={
                           0: lambda s: complex(s.decode().replace('+-', '-'))})
    else:
        gamma = np.loadtxt(gamma_file_name)

    # Read in corresponding phis (native);
    phi_native = read_native_phi(protein, phi_list, total_phis)

    # perform dot products to get native energies
    e_native = np.dot(gamma, phi_native)
    return e_native

def evaluate_hamiltonian_decoy_structures_provided(protein, hamiltonian, native_training_set_file, decoy_training_set_file, training_decoy_method, test_decoy_method, use_filtered_gammas=True):
    phi_list = read_phi_list(hamiltonian)
    decoy_training_set = read_column_from_file(decoy_training_set_file, 1)

    # read in Hamiltonian;
    # Find out how many total phi_i there are and get full parameter string
    total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string_decoy_structures_provided(
        phi_list, decoy_training_set)

    # read in corresponding gammas, note there is only one gamma for both the native and decoy sets;
    if use_filtered_gammas:
        gamma_file_name = "%s%s_%s_gamma_filtered" % (
            gammas_directory, native_training_set_file.split('/')[-1].split('.')[0], full_parameters_string)
    else:
        gamma_file_name = "%s%s_%s_gamma" % (gammas_directory, native_training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)

    # Need to filter out the complex number if in the "filtered" mode;
    if use_filtered_gammas:
        gamma = np.loadtxt(gamma_file_name, dtype=complex, converters={
                           0: lambda s: complex(s.decode().replace('+-', '-'))})
    else:
        gamma = np.loadtxt(gamma_file_name)

    # Read in corresponding phis (decoy);
    phi_i_decoy = read_decoy_phi_structures_provided(protein, phi_list, total_phis)

    # perform dot products to get decoy energies
    e_decoy = np.dot(gamma, phi_i_decoy)
    return e_decoy


def evaluate_hamiltonian(protein, hamiltonian, training_set_file, training_decoy_method, test_decoy_method, num_decoys, use_filtered_gammas=True):
    phi_list = read_phi_list(hamiltonian)
    training_set = read_column_from_file(training_set_file, 1)
    # read in Hamiltonian
    # Find out how many total phi_i there are and get full parameter string
    total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string(
        phi_list, training_set)
    # read in corresponding gammas
    if use_filtered_gammas:
        gamma_file_name = "%s%s_%s_gamma_filtered" % (
            gammas_directory, training_set_file.split('/')[-1].split('.')[0], full_parameters_string)
    else:
        gamma_file_name = "%s%s_%s_gamma" % (gammas_directory, training_set_file.split(
            '/')[-1].split('.')[0], full_parameters_string)

    # Need to filter out the complex number if in the "filtered" mode;
    if use_filtered_gammas:
        gamma = np.loadtxt(gamma_file_name, dtype=complex, converters={
                           0: lambda s: complex(s.decode().replace('+-', '-'))})
    else:
        gamma = np.loadtxt(gamma_file_name)

    # read in corresponding phis (native and decoys)
    phi_native = read_native_phi(protein, phi_list, total_phis)
    phi_i_decoy = read_decoy_phis(
        protein, phi_list, total_phis, num_phis, num_decoys, test_decoy_method)
    # perform dot products to get energies (native and decoys)
    e_decoy = np.zeros(num_decoys)
    e_native = np.dot(gamma, phi_native)
    for i_decoy in range(num_decoys):
        e_decoy[i_decoy] = np.dot(gamma, phi_i_decoy[i_decoy])
    e_mg = np.average(e_decoy)
    e_mg_std = np.std(e_decoy)
    # calculate z-score
    z_score = (e_mg - e_native) / e_mg_std
    return z_score, e_native, e_mg, e_mg_std


def evaluate_hamiltonian_wei(protein, hamiltonian, training_set_file, gamma_file_name, test_decoy_method, num_decoys, use_filtered_gammas=True, outputDecoy=False, **kwargs):
    phi_list = read_phi_list(hamiltonian)
    training_set = read_column_from_file(training_set_file, 1)
    # read in Hamiltonian
    # Find out how many total phi_i there are and get full parameter string
    total_phis, full_parameters_string, num_phis = get_total_phis_and_parameter_string(
        phi_list, training_set, **kwargs)
    # print(total_phis, full_parameters_string, num_phis)
    # read in corresponding gammas
    # if use_filtered_gammas:
    #     gamma_file_name = "%s%s_%s_gamma_filtered" % (
    #         gammas_directory, training_set_file.split('/')[-1].split('.')[0], full_parameters_string)
    # else:
    #     pass
    #     gamma_file_name = "%s%s_%s_gamma" % (gammas_directory, training_set_file.split(
    #         '/')[-1].split('.')[0], full_parameters_string)

    # Need to filter out the complex number if in the "filtered" mode;
    # if use_filtered_gammas:
    #     gamma = np.loadtxt(gamma_file_name, dtype=complex, converters={
    #                        0: lambda s: complex(s.decode().replace('+-', '-'))})
    # else:
    #     gamma = np.loadtxt(gamma_file_name)
    # gamma = np.loadtxt(gamma_file_name, dtype=complex, converters={
    #                             0: lambda s: complex(s.decode().replace('+-', '-'))})
    gamma = np.loadtxt(gamma_file_name)
    # read in corresponding phis (native and decoys)
    phi_native = read_native_phi(protein, phi_list, total_phis, **kwargs)
    phi_i_decoy = read_decoy_phis(
        protein, phi_list, total_phis, num_phis, num_decoys, test_decoy_method, **kwargs)
    # perform dot products to get energies (native and decoys)
    e_decoy = np.zeros(num_decoys)
    # print(gamma.shape, phi_native.shape)
    e_native = np.dot(gamma, phi_native)
    for i_decoy in range(num_decoys):
        e_decoy[i_decoy] = np.dot(gamma, phi_i_decoy[i_decoy])
    e_mg = np.average(e_decoy)
    e_mg_std = np.std(e_decoy)
    # calculate z-score
    z_score = (e_mg - e_native) / e_mg_std
    # print(protein, z_score, e_native, e_mg, e_mg_std)
    # if protein == "1r69":
    #     print(phi_i_decoy)
        # print(e_decoy)
    if outputDecoy:
        return e_native, e_decoy
    return z_score, e_native, e_mg, e_mg_std

def validate_hamiltonian_wei(hamiltonian, training_set_file, gamma_file_name, training_decoy_method, num_decoys, test_set_file=None, test_decoy_method=None, use_filtered_gammas=False, **kwargs):
    if test_set_file is None:
        test_set_file = training_set_file
    if test_decoy_method is None:
        test_decoy_method = training_decoy_method
    test_set = read_column_from_file(test_set_file, 1)
    z_scores = []
    e_natives = []
    e_mgs = []
    e_mg_stds = []
    for i, protein in enumerate(test_set):

        z, en, emg, emgstd = evaluate_hamiltonian_wei(
            protein, hamiltonian, training_set_file, gamma_file_name, test_decoy_method, num_decoys, use_filtered_gammas, **kwargs)
        # print(protein, z)
        if np.isnan(z):
            continue
        z_scores.append(z)
        e_natives.append(en)
        e_mgs.append(emg)
        e_mg_stds.append(emgstd)
        if i % 1000 == 0:
            print(i, z)
    data = pd.concat([pd.Series(test_set), pd.Series(z_scores), pd.Series(e_natives), pd.Series(e_mgs), pd.Series(e_mg_stds)], axis=1)
    data.columns = ["Protein", "Z_scores", "E_native", "E_mgs", "Std_mg"]
    data.Z_scores = data.Z_scores.astype(float)
    data.E_native = data.E_native.astype(float)
    return data


# def validate_hamiltonian_decoy_structures_provided(hamiltonian, native_training_set_file, decoy_training_set_file, training_decoy_method, native_test_set_file=None, decoy_test_set_file=None, test_decoy_method=None, use_filtered_gammas=False):
#     if native_test_set_file == None:
#         native_test_set_file = native_training_set_file
#     if decoy_test_set_file == None:
#         decoy_test_set_file = decoy_training_set_file
#     if test_decoy_method == None:
#         test_decoy_method = training_decoy_method
#     native_test_set = read_column_from_file(native_test_set_file, 1)
#     decoy_test_set = read_column_from_file(decoy_test_set_file, 1)
#     z_scores = []
#     e_natives = []
#     e_mgs = []
#     e_mg_stds = []
#     for protein in native_test_set:
#         en = evaluate_hamiltonian_native_structures_provided(
#             protein, hamiltonian, native_test_set_file, training_decoy_method, test_decoy_method, use_filtered_gammas)
#         e_natives.append(en)
#     for protein in decoy_test_set:
#         e_decoy = evaluate_hamiltonian_decoy_structures_provided(
#             protein, hamiltonian, native_test_set_file, decoy_test_set_file, training_decoy_method, test_decoy_method, use_filtered_gammas)
#         e_mgs.append(e_decoy)

#     # Calculate the real averaged <E>mg, <E>n, std(Emg) and zScore;
#     average_e_mg = np.average(e_mgs)
#     std_e_mgs = np.std(e_mgs)
#     average_e_native = np.average(e_natives)
#     z_score_forall = (average_e_mg - average_e_native) / std_e_mgs
#     return z_score_forall, average_e_native, average_e_mg, std_e_mgs, e_mgs, e_natives

# def validate_hamiltonian(hamiltonian, training_set_file, training_decoy_method, num_decoys, test_set_file=None, test_decoy_method=None, use_filtered_gammas=False):
#     if test_set_file == None:
#         test_set_file = training_set_file
#     if test_decoy_method == None:
#         test_decoy_method = training_decoy_method
#     test_set = read_column_from_file(test_set_file, 1)
#     z_scores = []
#     e_natives = []
#     e_mgs = []
#     e_mg_stds = []
#     for protein in test_set:
#         z, en, emg, emgstd = evaluate_hamiltonian(
#             protein, hamiltonian, training_set_file, training_decoy_method, test_decoy_method, num_decoys, use_filtered_gammas)
#         if np.isnan(z):
#             continue
#         z_scores.append(z)
#         e_natives.append(en)
#         e_mgs.append(emg)
#         e_mg_stds.append(emgstd)
#     return z_scores, e_natives, e_mgs, e_mg_stds


###########################################################################
# Conversion function for the AWSEM-PCA project;
# Cause I made a mistake
#
# Conversion formulae: origin_data = coefficient * (x - parameter)
# converted_data = 1/std * (origin_data / coefficient + parameter)

# Written by Xingcheng Lin, 09/28/2018
###########################################################################


def conversion_from_featureScaling_to_scalebyStd(parameter, coefficient, std, data_x):
    parameter = float(parameter)
    coefficient = float(coefficient)
    std = float(std)
    data_x = np.asarray(data_x)

    data_x_converted = 1 / std * (data_x / coefficient + parameter)

    return data_x_converted
