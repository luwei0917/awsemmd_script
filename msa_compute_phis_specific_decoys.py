import os
import sys
import numpy as np
import pandas as pd
from Bio.PDB.Polypeptide import three_to_index
from Bio.PDB.Polypeptide import one_to_three
from Bio.PDB.Polypeptide import three_to_one
import random
import argparse
from functools import partial
import time
# res to index
# index_dic
parameter_folder = "/Users/weilu/Research/database/"
parameter_folder = "/home/wl45/dataset/"

parser = argparse.ArgumentParser(
    description="compute phis by shuffle MSA")

parser.add_argument("pdb", type=str, help="The name of the pdb, and output name")
parser.add_argument("file", type=str, help="location to the pdb file")
parser.add_argument("msa", type=str, help="location to the msa file")
parser.add_argument("to", type=str, help="to the path that stores the results")
parser.add_argument("--topo", type=str, default=None, help="default topo is None")
parser.add_argument("--decoys", type=str, default=None, help="default decoys is None")
parser.add_argument("-m", "--mode", type=int, default=0, help="default mode is 0")
parser.add_argument("-s", "--seed", type=int, default=1, help="random seed")
args = parser.parse_args()

with open('cmd_msa_compute_phis.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')


if os.path.exists(f"{parameter_folder}/gxxxg_index_dic_v5.csv"):
    info = pd.read_csv(f"{parameter_folder}/gxxxg_index_dic_v5.csv", index_col=0)
    interaction_index_dic = {}
    n = 20
    for i, line in info.iterrows():
        # print(i, line["i"])
        interaction_index_dic[f"{line['Direction']}_{line['i']}"] = line["Group"]


def shuffle_string(string):
    list_string = list(string)
    random.shuffle(list_string)
    return ''.join(list_string)

def get_400_based_index(res1_1, res1_2):
    index = three_to_index(res1_1)*20 + three_to_index(res1_2)
    return index
def get_overall_index_v5(index1, index2, direction, interaction_index_dic):
    n = 20
    # plus 1, total parameters. 21*20/2 = 210
    n_shift = 210
    new_index1 = interaction_index_dic[f"{direction}_{index1}"]
    new_index2 = interaction_index_dic[f"{direction}_{index2}"]
    if new_index1 > new_index2:
        new_index1, new_index2 = new_index2, new_index1
    overall_index = ((2*n-(new_index1-1))*(new_index1)/2 + new_index2 - new_index1)
    if direction == "anti":
        overall_index += n_shift
    return int(overall_index)

def get_interaction_index_from_four_residues_v5(res1_1, res1_2, res2_1, res2_2, direction, interaction_index_dic=interaction_index_dic):

    index1 = get_400_based_index(res1_1, res1_2)
    if direction == "parallel":
        index2 = get_400_based_index(res2_1, res2_2)
    elif direction == "anti":
        index2 = get_400_based_index(res2_2, res2_1)
    else:
        print("unknown direction")
        raise
    index = get_overall_index_v5(index1, index2, direction, interaction_index_dic)

    return index


def dis(a, b):
    return ((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)**0.5

def get_side_chain_center_of_mass(res):
    atoms = res.get_atoms()
    total = np.array([0., 0., 0.])
    total_mass = 0
    for atom in atoms:
        if atom.get_name() in ["N", "CA", "C", "O"]:
            continue
        if atom.element == "H":
            continue
        total += atom.mass * atom.get_coord()
        total_mass += atom.mass
        # print(atom.get_name(), atom.get_coord())
    if total_mass == 0:
        x_com = res["CA"].get_coord()
    else:
        x_com = total / total_mass
    return x_com

def get_interaction_distance_com(res1, res2):
    # when Some Amino acids is mutated to GLY.
    # I want it still use the COM of side chain. not the position of CA.
    x1 = get_side_chain_center_of_mass(res1)
    x2 = get_side_chain_center_of_mass(res2)
    return dis(x1, x2)

def get_local_index(residue):
    return residue.get_id()[1]


def get_global_index(residue_list, residue):
    return residue_list.index(residue)


def get_chain(residue):
    return residue.get_parent().get_id()

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

def get_neighbors_within_radius(neighbor_list, residue, radius):
    return neighbor_list.search(get_interaction_atom(residue).get_coord(), radius, level='R')

def get_res_by_globalindex(res_list, index, chain):
    # the res has to be on the same chain as "chain"
    if index < 0:
        return -1
    try:
        res = res_list[index]
    except:
        return -1
    if res.get_parent().get_id() == chain:
        return res
    else:
        return -1
def interaction_well(r, r_min, r_max, kappa):
    return 0.5 * (np.tanh(kappa * (r - r_min)) * np.tanh(kappa * (r_max - r))) + 0.5


from Bio.PDB import PDBParser
from Bio.PDB import NeighborSearch
from Bio.PDB import Selection
from Bio.PDB import PPBuilder
def parse_pdb(pdb_id):
    parser = PDBParser()
    return parser.get_structure(pdb_id, "%s.pdb" % pdb_id)
def is_hetero(residue):
    if residue.id[0] != ' ':
        return True
    else:
        return False
def get_res_list(structure):
    res_list = Selection.unfold_entities(structure, 'R')
    # Get all residues from a structure
    res_list = [residue for residue in res_list if not is_hetero(residue)]
    return res_list
def get_atom_list(structure):
    atom_list = Selection.unfold_entities(structure, 'A')  # A for atoms
    return atom_list
def get_neighbor_list(structure):

    res_list = get_res_list(structure)
    atom_list = [a for a in get_atom_list(
        structure) if not is_hetero(a.get_parent())]
    # print(atom_list)
    neighbor_list = NeighborSearch(atom_list)
    return neighbor_list
def get_sequence_from_structure(structure):
    sequence = ""
    ppb = PPBuilder(radius=10.0)
    for pp in ppb.build_peptides(structure, aa_only=False):
        sequence += '%s\n' % pp.get_sequence()
    return sequence.replace('\n', '')

def get_phi_info_gxxxg_v5_well(res_list, neighbor_list, parameter_list):
    info = []
    min_seq_sep = 10
    r_min = 2.0
    r_max = 6.5
    r_cutoff = 8.5
    kappa = 5
    n_parameters = 420
    info_list = []
    get_distance_between_two_residues = get_interaction_distance_com
    phi_gxxxg_well = np.zeros(n_parameters)
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_cutoff):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)

            for shift_to_res2_2, direction in zip([-4, 4], ["anti", "parallel"]):
                res1_2_globalindex = res1globalindex + 4
                res1_2 = get_res_by_globalindex(res_list, res1_2_globalindex, res1chain)
                # for shift_to_res2_2 in [-4, 4]:

                # consider parallel, and anti-parallel.
                res2_2_globalindex = res2globalindex + shift_to_res2_2
                res2_2 = get_res_by_globalindex(res_list, res2_2_globalindex, res2chain)
                if res1_2 == -1 or res2_2 == -1:
                    continue
                if direction == "parallel":
                    group2index = res2globalindex
                elif direction == "anti":
                    group2index = res2_2_globalindex
                sep = group2index - res1globalindex
                if (res1chain == res2chain and sep >= min_seq_sep) or (res1chain != res2chain and group2index > res1globalindex):
                    rij = get_distance_between_two_residues(res1, res2)

                    rij_2 = get_distance_between_two_residues(res1_2, res2_2)
                    if rij_2 > r_cutoff or rij > r_cutoff:
                        continue

                    res1_name = res1.resname
                    res1_2_name = res1_2.resname
                    res2_name = res2.resname
                    res2_2_name = res2_2.resname
                    # interaction_index = get_interaction_index_from_four_residues_v5(res1.resname, res1_2.resname, res2.resname, res2_2.resname, direction)
                    interaction_index = get_interaction_index_from_four_residues_v5(res1_name, res1_2_name, res2_name, res2_2_name, direction)
                    phi_ = interaction_well(rij, r_min, r_max, kappa) * interaction_well(rij_2, r_min, r_max, kappa)
                    # phi_gxxxg_well[interaction_index] += phi_
                    phi_gxxxg_well[interaction_index] += phi_
                    if phi_ > 1e-5:
                        info.append([round(phi_,4), res1globalindex, res1_2_globalindex, res2globalindex, res2_2_globalindex, direction, res1_name, res1_2_name, res2_name, res2_2_name, interaction_index])
    info = pd.DataFrame(info, columns=["phi", "res1", "res1_2", "res2", "res2_2", "direction", "res1_name", "res1_2_name", "res2_name", "res2_2_name", "interaction_index"])
    return info

def get_overall_index_v6(index1, index2, direction, interaction_index_dic):
    n = 400
    # plus 1, total parameters. 401*400/2 = 80200
    n_shift = 80200
    new_index1 = index1
    new_index2 = index2
    if new_index1 > new_index2:
        new_index1, new_index2 = new_index2, new_index1
    overall_index = ((2*n-(new_index1-1))*(new_index1)/2 + new_index2 - new_index1)
    if direction == "anti":
        overall_index += n_shift
    return int(overall_index)

def get_interaction_index_from_four_residues_v6(res1_1, res1_2, res2_1, res2_2, direction):

    index1 = get_400_based_index(res1_1, res1_2)
    if direction == "parallel":
        index2 = get_400_based_index(res2_1, res2_2)
    elif direction == "anti":
        index2 = get_400_based_index(res2_2, res2_1)
    else:
        print("unknown direction")
        raise
    index = get_overall_index_v6(index1, index2, direction, interaction_index_dic)

    return index

def get_phi_info_gxxxg_v6_well(res_list, neighbor_list, parameter_list, n_parameters=420, get_interaction_index_from_four_residues=None):
    info = []
    min_seq_sep = 10
    r_min = 2.0
    r_max = 6.5
    r_cutoff = 8.5
    kappa = 5
    info_list = []
    get_distance_between_two_residues = get_interaction_distance_com
    phi_gxxxg_well = np.zeros(n_parameters)
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_cutoff):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)

            for shift_to_res2_2, direction in zip([-4, 4], ["anti", "parallel"]):
                res1_2_globalindex = res1globalindex + 4
                res1_2 = get_res_by_globalindex(res_list, res1_2_globalindex, res1chain)
                # for shift_to_res2_2 in [-4, 4]:

                # consider parallel, and anti-parallel.
                res2_2_globalindex = res2globalindex + shift_to_res2_2
                res2_2 = get_res_by_globalindex(res_list, res2_2_globalindex, res2chain)
                if res1_2 == -1 or res2_2 == -1:
                    continue
                if direction == "parallel":
                    group2index = res2globalindex
                elif direction == "anti":
                    group2index = res2_2_globalindex
                sep = group2index - res1globalindex
                if (res1chain == res2chain and sep >= min_seq_sep) or (res1chain != res2chain and group2index > res1globalindex):
                    rij = get_distance_between_two_residues(res1, res2)

                    rij_2 = get_distance_between_two_residues(res1_2, res2_2)
                    if rij_2 > r_cutoff or rij > r_cutoff:
                        continue

                    res1_name = res1.resname
                    res1_2_name = res1_2.resname
                    res2_name = res2.resname
                    res2_2_name = res2_2.resname
                    # interaction_index = get_interaction_index_from_four_residues_v5(res1.resname, res1_2.resname, res2.resname, res2_2.resname, direction)
                    interaction_index = get_interaction_index_from_four_residues(res1_name, res1_2_name, res2_name, res2_2_name, direction)
                    phi_ = interaction_well(rij, r_min, r_max, kappa) * interaction_well(rij_2, r_min, r_max, kappa)
                    # phi_gxxxg_well[interaction_index] += phi_
                    phi_gxxxg_well[interaction_index] += phi_
                    if phi_ > 1e-5:
                        info.append([phi_, res1globalindex, res1_2_globalindex, res2globalindex, res2_2_globalindex, direction, res1_name, res1_2_name, res2_name, res2_2_name, interaction_index])
    info = pd.DataFrame(info, columns=["phi", "res1", "res1_2", "res2", "res2_2", "direction", "res1_name", "res1_2_name", "res2_name", "res2_2_name", "interaction_index"])
    return info

def get_interaction_distance(res1, res2):
    return get_interaction_atom(res1) - get_interaction_atom(res2)

def get_side_chain_atoms(res):
    atoms = res.get_atoms()
    atom_list = []
    if res.resname == "GLY":
        atom_list = list(atoms)
        # atom_list.append(res["CA"])
    for atom in atoms:
        if atom.get_name() in ["N", "CA", "C", "O", "OXT",]:
            continue
        # if atom.element == "H":
        #     continue
        atom_list.append(atom)
    return atom_list

def get_phi_info_gxxxg_v7_well(res_list, neighbor_list, parameter_list, n_parameters=420, get_distance_between_two_residues=None, get_interaction_index_from_four_residues=None):
    info = []
    min_seq_sep = 10
    r_min = 2.0
    r_max = 6.5
    r_cutoff = 8.5
    kappa = 5
    info_list = []
    # get_distance_between_two_residues = get_interaction_distance_com
    phi_gxxxg_well = np.zeros(n_parameters)
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_cutoff):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)

            for shift_to_res2_2, direction in zip([-4, 4], ["anti", "parallel"]):
                res1_2_globalindex = res1globalindex + 4
                res1_2 = get_res_by_globalindex(res_list, res1_2_globalindex, res1chain)
                # for shift_to_res2_2 in [-4, 4]:

                # consider parallel, and anti-parallel.
                res2_2_globalindex = res2globalindex + shift_to_res2_2
                res2_2 = get_res_by_globalindex(res_list, res2_2_globalindex, res2chain)
                if res1_2 == -1 or res2_2 == -1:
                    continue
                if direction == "parallel":
                    group2index = res2globalindex
                elif direction == "anti":
                    group2index = res2_2_globalindex
                sep = group2index - res1globalindex
                if (res1chain == res2chain and sep >= min_seq_sep) or (res1chain != res2chain and group2index > res1globalindex):
                    rij = get_distance_between_two_residues(res1, res2)

                    rij_2 = get_distance_between_two_residues(res1_2, res2_2)
                    if rij_2 > r_cutoff or rij > r_cutoff:
                        continue

                    # # below is new.
                    # # res1_1, and res2_1
                    # in_real_contact = False
                    # atom_list = get_side_chain_atoms(res1)
                    # for atom in atom_list:
                    #     if res2 in neighbor_list.search(atom.get_coord(), 4, level='R'):
                    #         in_real_contact = True
                    # if not in_real_contact:
                    #     # print(pdb, res1.id[1], resName1, res2.id[1], resName2)
                    #     continue
                    # # res1_1, and res2_1
                    # in_real_contact = False
                    # atom_list = get_side_chain_atoms(res1_2)
                    # for atom in atom_list:
                    #     if res2_2 in neighbor_list.search(atom.get_coord(), 4, level='R'):
                    #         in_real_contact = True
                    # if not in_real_contact:
                    #     # print(pdb, res1.id[1], resName1, res2.id[1], resName2)
                    #     continue
                    # # -----
                    res1_name = res1.resname
                    res1_2_name = res1_2.resname
                    res2_name = res2.resname
                    res2_2_name = res2_2.resname
                    # interaction_index = get_interaction_index_from_four_residues_v5(res1.resname, res1_2.resname, res2.resname, res2_2.resname, direction)
                    interaction_index = get_interaction_index_from_four_residues(res1_name, res1_2_name, res2_name, res2_2_name, direction)
                    phi_ = interaction_well(rij, r_min, r_max, kappa) * interaction_well(rij_2, r_min, r_max, kappa)
                    # phi_gxxxg_well[interaction_index] += phi_
                    phi_gxxxg_well[interaction_index] += phi_
                    if phi_ > 1e-5:
                        info.append([phi_, res1globalindex, res1_2_globalindex, res2globalindex, res2_2_globalindex, direction, res1_name, res1_2_name, res2_name, res2_2_name, interaction_index])
    info = pd.DataFrame(info, columns=["phi", "res1", "res1_2", "res2", "res2_2", "direction", "res1_name", "res1_2_name", "res2_name", "res2_2_name", "interaction_index"])
    return info

def encode_four_body_index(res1_1, res1_2, res2_1, res2_2, direction):
    # if anti-parallel, the order by reading from res2_2 to res2_1
    if direction == "anti":
        res2 = three_to_index(res2_2)*20 + three_to_index(res2_1)
        Res2_letter = res2_2 + "_" + res2_1
    elif direction == "parallel":
        res2 = three_to_index(res2_1)*20 + three_to_index(res2_2)
        Res2_letter = res2_1 + "_" + res2_2
    res1 = three_to_index(res1_1)*20 + three_to_index(res1_2)
    Res1_letter = res1_1 + "_" + res1_2
    # if res2 index is smaller than res1, We will swtich the index. to ensure res1 is less than res2.
    if res2 < res1:
        return res2, res1, Res1_letter, Res2_letter
    else:
        return res1, res2, Res1_letter, Res2_letter

def get_interaction_index_dic(fileLocation, cutoff):
    # c = pd.read_csv("/Users/weilu/Research/database/interaction_index_single_chain.csv", index_col=0)
    c = pd.read_csv(fileLocation, index_col=0)
    c_anti = c.query("Direction=='anti'").reset_index(drop=True).reset_index()
    c_parallel = c.query("Direction=='parallel'").reset_index(drop=True).reset_index()

    interaction_index_dic = {}
    shift = cutoff + 1

    for i, line in c_parallel.iterrows():
        direction = line["Direction"]
        res1 = line["Res1"]
        res2 = line["Res2"]
        index = line["index"]
        if index < cutoff:
            interaction_index_dic[f"{direction}_{res1}_{res2}"] = index
        else:
            interaction_index_dic[f"{direction}_{res1}_{res2}"] = cutoff
    for i, line in c_anti.iterrows():
        direction = line["Direction"]
        res1 = line["Res1"]
        res2 = line["Res2"]
        index = line["index"]
        if index < cutoff:
            interaction_index_dic[f"{direction}_{res1}_{res2}"] = index + shift
        else:
            interaction_index_dic[f"{direction}_{res1}_{res2}"] = cutoff + shift

    return interaction_index_dic



def get_interaction_index_from_four_residues_v7(res1_1, res1_2, res2_1, res2_2, direction, interaction_index_dic=None):
    res1, res2, Res1_letter, Res2_letter = encode_four_body_index(res1_1, res1_2, res2_1, res2_2, direction)
    unique_index = f"{direction}_{res1}_{res2}"
    try:
        index = interaction_index_dic[unique_index]
        return index
    except:
        if direction == "anti":
            return 401
        else:
            return 200




six_letter_code_letters = {
    'I':3, 'M':3, 'L':3, 'V':3, 'F':5, 'Y':5, 'W':5, 'G':0, 'P':2, 'C':1, 'A':1, 'S':1, 'T':1, 'N':4, 'H':4, 'Q':4, 'E':4, 'D':4, 'R':4, 'K':4
}


def get_six_letter_based_index(res1_1, res1_2):
    index = six_letter_code_letters[three_to_one(res1_1)]*6 + six_letter_code_letters[three_to_one(res1_2)]
    return index

def get_overall_index_v8(index1, index2, direction, n):
    # n =36
    # plus 1, total parameters. 401*400/2 = 80200
    n_shift = int((n+1)*n/2)
    new_index1 = index1
    new_index2 = index2
    if new_index1 > new_index2:
        new_index1, new_index2 = new_index2, new_index1
    overall_index = ((2*n-(new_index1-1))*(new_index1)/2 + new_index2 - new_index1)
    if direction == "anti":
        overall_index += n_shift
    return int(overall_index)

def get_interaction_index_from_four_residues_v8(res1_1, res1_2, res2_1, res2_2, direction):

    index1 = get_six_letter_based_index(res1_1, res1_2)
    if direction == "parallel":
        index2 = get_six_letter_based_index(res2_1, res2_2)
    elif direction == "anti":
        index2 = get_six_letter_based_index(res2_2, res2_1)
    else:
        print("unknown direction")
        raise
    index = get_overall_index_v8(index1, index2, direction, n=36)

    return index



def prot_water_switchFunc_sigmaWater(rho_i, rho_j, rho_0, kappa):
    return 0.25 * (1 - np.tanh(kappa * (rho_i - rho_0))) * (1 - np.tanh(kappa * (rho_j - rho_0)))


def prot_water_switchFunc_sigmaProt(rho_i, rho_j, rho_0, kappa):
    return 1 - prot_water_switchFunc_sigmaWater(rho_i, rho_j, rho_0, kappa)

def contact_interaction_well(rij, r_min, r_max, kappa, rho_i, rho_j, density_threshold, density_kappa, interactionType):
    if interactionType == "Direct":
        return interaction_well(rij, r_min, r_max, kappa)
    elif interactionType == "HighDensityMediated":
        return     prot_water_switchFunc_sigmaProt(
                    rho_i, rho_j, density_threshold, density_kappa) * interaction_well(rij, r_min, r_max, kappa)
    elif interactionType == "LowDensityMediated":
        return prot_water_switchFunc_sigmaWater(
                    rho_i, rho_j, density_threshold, density_kappa) * interaction_well(rij, r_min, r_max, kappa)
    else:
        print("ERROR:, ", interactionType)

res_type_map = {
    'A': 0,'C': 4,'D': 3,'E': 6,'F': 13,'G': 7,'H': 8,'I': 9,'K': 11,'L': 10,'M': 12,'N': 2,'P': 14,'Q': 5,'R': 1,'S': 15,'T': 16,'V': 19,'W': 17,'Y': 18
}


from Bio.PDB.Polypeptide import three_to_one
def get_direct_contact_interaction_index(res1_name, res2_name, interactionType, burial_i=-1):
    n = 20
    # plus 1, total parameters. 21*20/2 = 210
    n_shift = 210
    index1 = res_type_map[three_to_one(res1_name)]
    if interactionType == "Burial":
        overall_index = 3 * n_shift + burial_i * 20 + index1
        return int(overall_index)
    index2 = res_type_map[three_to_one(res2_name)]
    if index1 > index2:
        index1, index2 = index2, index1
    overall_index = ((2*n-(index1-1))*(index1)/2 + index2 - index1)
    if interactionType == "HighDensityMediated":
        overall_index += n_shift
    if interactionType == "LowDensityMediated":
        overall_index += 2 * n_shift
    return int(overall_index)


def calculate_cb_density(res_list, neighbor_list, min_seq_sep=2, rmin=2.5):
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
                density[res1globalindex] += interaction_well(rij, rmin, 6.5, 5)
    return density

def get_phi_info_contact_well(res_list, neighbor_list, parameter_list,
                               n_parameters=690, get_distance_between_two_residues=None):
    cb_density = calculate_cb_density(res_list, neighbor_list)
    info = []
    min_seq_sep = 10
    r_min_direct = 2.0
    r_max_direct = 6.5
    r_min_mediated = 6.5
    r_max_mediated = 9.5
    r_cutoff = 12.0
    kappa = 5
    density_threshold = 2.6
    density_kappa = 7.0
    rho_table = [[0.0, 3.0], [3.0, 6.0], [6.0, 9.0]]
    burial_kappa= 4.0
    info_list = []
    # get_distance_between_two_residues = get_interaction_distance_com
    phi_contact_well = np.zeros(n_parameters)
    for res1globalindex, res1 in enumerate(res_list):
        res1index = get_local_index(res1)
        res1chain = get_chain(res1)
        rho_i = cb_density[res1globalindex]
        res1_name = res1.resname
        for res2 in get_neighbors_within_radius(neighbor_list, res1, r_cutoff):
            res2index = get_local_index(res2)
            res2chain = get_chain(res2)
            res2globalindex = get_global_index(res_list, res2)
            rho_j = cb_density[res2globalindex]
            res2_name = res2.resname
            sep = res2globalindex - res1globalindex
            if (res1chain == res2chain and sep >= min_seq_sep) or (res1chain != res2chain and res2globalindex > res1globalindex):
                rij = get_distance_between_two_residues(res1, res2)

                for interactionType in ["Direct", "LowDensityMediated", "HighDensityMediated"]:
                    if interactionType == "Direct":
                        r_min, r_max = r_min_direct, r_max_direct
                    else:
                        r_min, r_max = r_min_mediated, r_max_mediated
                    interaction_index = get_direct_contact_interaction_index(res1_name, res2_name, interactionType)
                    phi_ = contact_interaction_well(rij, r_min, r_max, kappa, rho_i, rho_j, density_threshold, density_kappa, interactionType)
                    phi_contact_well[interaction_index] += phi_
                    if phi_ > 1e-5:
                        info.append([phi_, res1globalindex, res2globalindex, interactionType, res1_name, res2_name, interaction_index])
        for i in range(3):
            interactionType = "Burial"
            phi_burial_i = interaction_well(rho_i, rho_table[i][0], rho_table[i][1], burial_kappa)
            interaction_index = get_direct_contact_interaction_index(res1_name, res2_name, interactionType, burial_i=i)
            if phi_burial_i > 1e-5:
                info.append([phi_burial_i, res1globalindex, i, "Burial", res1_name, "NA", interaction_index])
    info = pd.DataFrame(info, columns=["phi", "res1", "res2", "Type", "res1_name", "res2_name", "interaction_index"])
    return info


def get_phis_from_info_and_sequence(info, sequence, n_parameters=420, get_interaction_index_from_four_residues=None, verbose=False):
    phi_gxxxg_well = np.zeros(n_parameters)
    skip_count = 0
    for i, line in info.iterrows():
        direction = line["direction"]

        res1_index = line["res1"]
        res1_2_index = line["res1_2"]
        res2_index = line["res2"]
        res2_2_index = line["res2_2"]
        try:
            res1_name = one_to_three(sequence[res1_index])
            res1_2_name = one_to_three(sequence[res1_2_index])
            res2_name = one_to_three(sequence[res2_index])
            res2_2_name = one_to_three(sequence[res2_2_index])
        except:
            skip_count += 1
            continue
        interaction_index = get_interaction_index_from_four_residues(res1_name, res1_2_name, res2_name, res2_2_name, direction)
        phi = line["phi"]
        phi_gxxxg_well[interaction_index] += phi
    phis_to_return = []
    for i in range(n_parameters):
        phis_to_return.append(phi_gxxxg_well[i])
    if verbose:
        print("total skipped: ", skip_count, "total: ", len(info))
    return phis_to_return


def get_phis_from_info_and_sequence_v2(info, sequence, n_parameters=420, contact_mode=False, get_interaction_index_from_four_residues=None, verbose=False):
    # fix a bug in original v2
    phi_well = np.zeros(n_parameters)
    skip_count = 0
    for i, line in info.iterrows():
        if contact_mode:
            interactionType = line["Type"]
            if interactionType == "Burial":
                res1_index = int(line["res1"])
                try:
                    res1_name = one_to_three(sequence[res1_index])
                except Exception as e:
                    skip_count += 1
                    continue
                burial_i = line["res2"]
                res2_name = None
            else:
                res1_index = line["res1"]
                res2_index = line["res2"]
                try:
                    res1_name = one_to_three(sequence[res1_index])
                    res2_name = one_to_three(sequence[res2_index])
                except:
                    skip_count += 1
                    continue
                burial_i = -1
            interaction_index = get_direct_contact_interaction_index(res1_name, res2_name, interactionType, burial_i=burial_i)
        else:
            direction = line["direction"]

            res1_index = line["res1"]
            res1_2_index = line["res1_2"]
            res2_index = line["res2"]
            res2_2_index = line["res2_2"]
            try:
                res1_name = one_to_three(sequence[res1_index])
                res1_2_name = one_to_three(sequence[res1_2_index])
                res2_name = one_to_three(sequence[res2_index])
                res2_2_name = one_to_three(sequence[res2_2_index])
            except:
                skip_count += 1
                continue
            interaction_index = get_interaction_index_from_four_residues(res1_name, res1_2_name, res2_name, res2_2_name, direction)
        phi = line["phi"]
        phi_well[interaction_index] += phi

    if verbose:
        print("total skipped: ", skip_count, "total: ", len(info))
    return phi_well
def get_phis_from_info_and_sequence_v3(info, sequence, n_parameters=420, shift_index=0, get_interaction_index_from_four_residues=None, verbose=True):
    # fix a bug in original v3.
    phi_well = np.zeros(n_parameters)
    skip_count = 0
    for i, line in info.iterrows():
        interactionTerm = line["interactionTerm"]
        if interactionTerm == "contact":
            interactionType = line["Type"]
            if interactionType == "Burial":
                res1_index = int(line["res1"])
                try:
                    res1_name = one_to_three(sequence[res1_index])
                except Exception as e:
                    skip_count += 1
                    continue
                burial_i = line["res2"]
                res2_name = None
            else:
                res1_index = int(line["res1"])
                res2_index = int(line["res2"])
                try:
                    res1_name = one_to_three(sequence[res1_index])
                    res2_name = one_to_three(sequence[res2_index])
                except Exception as e:
                    # if f"{e}" != "'-'":
                    #     print(f"a{e}a")
                    skip_count += 1
                    continue
                burial_i = -1
            interaction_index = get_direct_contact_interaction_index(res1_name, res2_name, interactionType, burial_i=burial_i)
        elif interactionTerm == "fourBody":
            direction = line["direction"]

            res1_index = int(line["res1"])
            res1_2_index = int(line["res1_2"])
            res2_index = int(line["res2"])
            res2_2_index = int(line["res2_2"])
            try:
                res1_name = one_to_three(sequence[res1_index])
                res1_2_name = one_to_three(sequence[res1_2_index])
                res2_name = one_to_three(sequence[res2_index])
                res2_2_name = one_to_three(sequence[res2_2_index])
            except Exception as e:
                # if f"{e}" != "'-'":
                #     print(f"x{e}x")
                skip_count += 1
                continue
            interaction_index = get_interaction_index_from_four_residues(res1_name, res1_2_name, res2_name, res2_2_name, direction)
            interaction_index += shift_index           # could be changed to a dicitonary, based on the interactionTerm name.
        phi = line["phi"]
        phi_well[interaction_index] += phi


    if verbose:
        print("total skipped: ", skip_count, "total: ", len(info))
    return phi_well

def calculate_A_and_B_single_pdb(average_phi_decoy, phi_native, all_phis):
    A = average_phi_decoy - phi_native
    num_decoys, total_phis = all_phis.shape
    half_B = np.zeros((total_phis, total_phis))
    std_half_B = np.zeros((total_phis, total_phis))
    other_half_B = np.zeros((total_phis, total_phis))

    phis_i = all_phis.reshape(num_decoys, total_phis, 1)
    for j in range(total_phis):
        phis_j = phis_i[:, j].reshape(num_decoys, 1, 1)
        half_B[j] += np.average(phis_i * phis_j, axis=0).reshape(total_phis)
        std_half_B[j] += np.std(phis_i * phis_j, axis=0).reshape(total_phis)

    average_phi = np.average(all_phis, axis=0)
    other_half_B += average_phi.reshape(total_phis, 1) * average_phi.reshape(1, total_phis)
    B = half_B - other_half_B

    return A, B, half_B, other_half_B, std_half_B

def helix_swapping(seq, data, n_swap=1):
    new_seq = list(seq)
    helices_list = data["symbol"].unique()

    # could swap multiple times.
    n_swap_real = np.random.choice(range(1, n_swap + 1))
    for i in range(n_swap_real):
        # randomly choose two helices
        chosen_helices = np.random.choice(helices_list, 2, replace=False)
        # swap their sequence.
        ## get the size of two chosen helices
        h1_size = data.query(f"symbol=='{chosen_helices[0]}'")["count"].values[0]
        h1_start_index = data.query(f"symbol=='{chosen_helices[0]}'")["start_index"].values[0]
        h2_size = data.query(f"symbol=='{chosen_helices[1]}'")["count"].values[0]
        h2_start_index = data.query(f"symbol=='{chosen_helices[1]}'")["start_index"].values[0]
        ## the swap region size in chosen in random.
        swap_size_raw = np.random.randint(10, 31)
        ## swap size must be lower than the helices size.
        swap_size = min(h2_size, min(h1_size, swap_size_raw))

        ## the exact place within the helix is also chosen at random.
        h1_swap_start_shift = np.random.choice(range(h1_size - swap_size + 1))
        h1_swap_start_index = h1_start_index + h1_swap_start_shift
        h1_seq = seq[h1_swap_start_index:h1_swap_start_index+swap_size]

        h2_swap_start_shift = np.random.choice(range(h2_size - swap_size + 1))
        h2_swap_start_index = h2_start_index + h2_swap_start_shift
        h2_seq = seq[h2_swap_start_index:h2_swap_start_index+swap_size]

        ## has a chance the swapping sequence is reversed.
        reverse_chance = 0.2
        if random.random() < reverse_chance:
            h1_seq = h1_seq[::-1]
        if random.random() < reverse_chance:
            h2_seq = h2_seq[::-1]
        new_seq[h1_swap_start_index:h1_swap_start_index+swap_size] = h2_seq
        new_seq[h2_swap_start_index:h2_swap_start_index+swap_size] = h1_seq
    return "".join(new_seq)

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

# fileLocation = "/Users/weilu/Research/server/aug_2020/experimenting_optimization/database/dompdb/5tin_A"
fileLocation = args.file
structure = parse_pdb(fileLocation)
res_list = get_res_list(structure)
neighbor_list = get_neighbor_list(structure)
sequence = get_sequence_from_structure(structure)

pdb = args.pdb
toLocation = args.to

if args.mode == -2:
    # contact term, with four body.
    n_parameters = int(690 + 1332)
    n_msa = 1
    # msa = [msa[0]]
    n_shuffle = 1

    # n_msa = 10
    # # msa = [msa[0]]
    # n_shuffle = 10

    info_ = []
    info = get_phi_info_contact_well(res_list, neighbor_list, [], get_distance_between_two_residues=get_interaction_distance,)
    os.system(f"mkdir -p {toLocation}/info_folder")
    a = info.assign(interactionTerm="contact")
    # info.to_csv(f"{toLocation}/info_folder/{pdb}_contact.csv")
    info_.append(info)

    shift_index = 690
    get_interaction_index_from_four_residues = get_interaction_index_from_four_residues_v8
    get_distance_between_two_residues = get_interaction_distance
    info = get_phi_info_gxxxg_v7_well(res_list, neighbor_list, [], get_distance_between_two_residues=get_distance_between_two_residues, n_parameters=n_parameters, get_interaction_index_from_four_residues=get_interaction_index_from_four_residues)
    # os.system(f"mkdir -p {toLocation}/info_folder")
    b = info.assign(interactionTerm="fourBody")

    get_phis_from_info_and_sequence = partial(get_phis_from_info_and_sequence_v3, shift_index=shift_index)
    info = pd.concat([a,b], sort=False).reset_index(drop=True)
    info.to_csv(f"{toLocation}/info_folder/{pdb}_complete.csv")
    exit()

if args.mode == -1:
    # using 6 letter code.
    n_parameters = 1332  # 36*37
    n_msa = 1
    # msa = [msa[0]]
    n_shuffle = 1
    # print(interaction_index_dic)
    get_interaction_index_from_four_residues = get_interaction_index_from_four_residues_v8
    get_distance_between_two_residues = get_interaction_distance
    info = get_phi_info_gxxxg_v7_well(res_list, neighbor_list, [], get_distance_between_two_residues=get_distance_between_two_residues, n_parameters=n_parameters, get_interaction_index_from_four_residues=get_interaction_index_from_four_residues)
    os.system(f"mkdir -p {toLocation}/info_folder")
    info.to_csv(f"{toLocation}/info_folder/{pdb}.csv")
    exit()

msaFile = args.msa
msa = np.loadtxt(msaFile, dtype=str)

if args.topo:
    # if there is a specified topo file.
    result_topo = str(np.loadtxt(args.topo, dtype=str))

    # gather topo information.
    sym_table = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
    info = []
    for i in range(1, sym_table.index(max(result_topo))+1):
        # print(pdb, i, result_topo.count(sym_table[i]))
        info.append([i, sym_table[i], result_topo.count(sym_table[i]), result_topo.find(sym_table[i])])
    data = pd.DataFrame(info, columns=["i", "symbol", "count", "start_index"])
    # enure all helices size within in 10 to 30.
    topo_data = data.query("count >= 10 and count <= 30").reset_index(drop=True)
# shuffle the string
# each protein, sample 1000 MSA.
# each MSA, shuffle 100 times

if args.mode == 1:
    # instead of shuffling, using helices swapping.
    # set n_max_swap = 5
    # repeat mode 9, but with bug-fixed get_phis_from_info_and_sequence_v2
    # A bit less shuffling.
    # contact term
    n_parameters = 690
    n_msa = 1
    msa = [msa[0]]
    n_decoy = 1000
    withBiased = True
    structure_file_name = args.decoys
    decoy_structures, Qs = read_decoy_structures_andQ(structure_file_name)
    Qs = np.array(Qs)
    available_n_decoy = len(decoy_structures)
    bias = 1 - Qs
    info = get_phi_info_contact_well(res_list, neighbor_list, [], get_distance_between_two_residues=get_interaction_distance,)
    os.system(f"mkdir -p {toLocation}/info_folder")
    info.to_csv(f"{toLocation}/info_folder/{pdb}_contact.csv")
    get_interaction_index_from_four_residues = None
    get_phis_from_info_and_sequence = partial(get_phis_from_info_and_sequence_v2, contact_mode=True)

if args.mode == 2:
    # instead of shuffling, using helices swapping.
    # set n_max_swap = 5
    # repeat mode 9, but with bug-fixed get_phis_from_info_and_sequence_v2
    # A bit less shuffling.
    # contact term
    n_parameters = 690
    n_msa = 100
    # msa = [msa[0]]
    n_decoy = 50
    withBiased = True
    structure_file_name = args.decoys
    decoy_structures, Qs = read_decoy_structures_andQ(structure_file_name)
    Qs = np.array(Qs)
    available_n_decoy = len(decoy_structures)
    bias = 1 - Qs
    info = get_phi_info_contact_well(res_list, neighbor_list, [], get_distance_between_two_residues=get_interaction_distance,)
    os.system(f"mkdir -p {toLocation}/info_folder")
    info.to_csv(f"{toLocation}/info_folder/{pdb}_contact.csv")
    get_interaction_index_from_four_residues = None
    get_phis_from_info_and_sequence = partial(get_phis_from_info_and_sequence_v2, contact_mode=True)

if args.mode == 3:
    # same as mode 1, but n_decoy is only 50.
    # instead of shuffling, using helices swapping.
    # set n_max_swap = 5
    # repeat mode 9, but with bug-fixed get_phis_from_info_and_sequence_v2
    # A bit less shuffling.
    # contact term
    n_parameters = 690
    n_msa = 1
    msa = [msa[0]]
    n_decoy = 50
    withBiased = True
    structure_file_name = args.decoys
    decoy_structures, Qs = read_decoy_structures_andQ(structure_file_name)
    Qs = np.array(Qs)
    available_n_decoy = len(decoy_structures)
    bias = 1 - Qs
    info = get_phi_info_contact_well(res_list, neighbor_list, [], get_distance_between_two_residues=get_interaction_distance,)
    os.system(f"mkdir -p {toLocation}/info_folder")
    info.to_csv(f"{toLocation}/info_folder/{pdb}_contact.csv")
    get_interaction_index_from_four_residues = None
    get_phis_from_info_and_sequence = partial(get_phis_from_info_and_sequence_v2, contact_mode=True)
# if args.mode == 19:
#     # same as mode 10, but with bug-fixed get_phis_from_info_and_sequence_v3
#     # contact term, with four body.
#     n_parameters = int(690 + 1332)
#     n_msa = 50
#     # msa = [msa[0]]
#     n_shuffle = 100

#     # n_msa = 10
#     # # msa = [msa[0]]
#     # n_shuffle = 10
#     n_swap = 3
#     helix_swapping = partial(helix_swapping, n_swap=n_swap)
#     info_ = []
#     info = get_phi_info_contact_well(res_list, neighbor_list, [], get_distance_between_two_residues=get_interaction_distance,)
#     os.system(f"mkdir -p {toLocation}/info_folder")
#     a = info.assign(interactionTerm="contact")
#     # info.to_csv(f"{toLocation}/info_folder/{pdb}_contact.csv")
#     info_.append(info)

#     shift_index = 690
#     get_interaction_index_from_four_residues = get_interaction_index_from_four_residues_v8
#     get_distance_between_two_residues = get_interaction_distance
#     info = get_phi_info_gxxxg_v7_well(res_list, neighbor_list, [], get_distance_between_two_residues=get_distance_between_two_residues, n_parameters=n_parameters, get_interaction_index_from_four_residues=get_interaction_index_from_four_residues)
#     # os.system(f"mkdir -p {toLocation}/info_folder")
#     b = info.assign(interactionTerm="fourBody")

#     get_phis_from_info_and_sequence = partial(get_phis_from_info_and_sequence_v3, shift_index=shift_index)
#     info = pd.concat([a,b], sort=False).reset_index(drop=True)
#     info.to_csv(f"{toLocation}/info_folder/{pdb}_complete.csv")
# exit()
print("n_parameters: ", n_parameters, ", n_msa: ", n_msa, ", n_decoy: ", n_decoy, ", mode: ", args.mode)
# print(info)
print("Compute Phis Starts")
start_time = time.time()
# n_msa = 10
# n_shuffle = 20

total_decoys = n_msa * n_decoy

decoy_phis = np.zeros((total_decoys, n_parameters))
native_phis = np.zeros((n_msa, n_parameters))
random.seed(args.seed)
chosen_msa = random.sample(list(msa), n_msa)
chosen_index = random.sample(range(available_n_decoy), n_decoy)

chosen_decoy_structures = [decoy_structures[i] for i in chosen_index]
chosen_bias = [bias[i] for i in chosen_index]
average_bias = np.average(chosen_bias)

count = 0
for idx, seq in enumerate(chosen_msa):
    native_phis[idx] = get_phis_from_info_and_sequence(info, seq, n_parameters=n_parameters, get_interaction_index_from_four_residues=get_interaction_index_from_four_residues)


for i in range(n_decoy):
    structure_decoy = chosen_decoy_structures[i]
    res_list_decoy = get_res_list(structure_decoy)
    neighbor_list_decoy = get_neighbor_list(structure_decoy)
    sequence_decoy = get_sequence_from_structure(structure_decoy)
    info_decoy = get_phi_info_contact_well(res_list_decoy, neighbor_list_decoy, [], get_distance_between_two_residues=get_interaction_distance,)
    for idx, seq in enumerate(chosen_msa):
        phis = get_phis_from_info_and_sequence(info_decoy, seq, n_parameters=n_parameters, get_interaction_index_from_four_residues=get_interaction_index_from_four_residues)
        decoy_phis[count] = phis
        if withBiased:
            decoy_phis[count] *= chosen_bias[i] / average_bias
        count += 1


# compute A, B, other_half_B, std_half_B
all_phis = decoy_phis
average_phi_decoy = np.average(all_phis, axis=0)
phi_native = np.average(native_phis, axis=0)
A, B, half_B, other_half_B, std_half_B = calculate_A_and_B_single_pdb(average_phi_decoy, phi_native, all_phis)



# toLocation = "/Users/weilu/Research/server/aug_2020/curated_single_chain_optimization/optimization_msa"

# small chance it will save these information.
os.system(f"mkdir -p {toLocation}/decoys")
os.system(f"mkdir -p {toLocation}/phis")
os.system(f"mkdir -p {toLocation}/chosen_msa")
os.system(f"mkdir -p {toLocation}/chosen_bias")
if random.random() < 0.01:
    np.save(f"{toLocation}/phis/{pdb}.npy", decoy_phis)
    np.save(f"{toLocation}/chosen_msa/{pdb}.npy", chosen_msa)
np.save(f"{toLocation}/chosen_bias/{pdb}.npy", chosen_bias)

# free memory
all_phis = None
decoy_phis = None
A_B_dic = {}
A_B_dic["A"] = A
A_B_dic["B"] = B
A_B_dic["half_B"] = half_B
A_B_dic["other_half_B"] = other_half_B
A_B_dic["std_half_B"] = std_half_B
A_B_dic["A_prime"] = average_phi_decoy

os.system(f"mkdir -p {toLocation}/A_B_dic")
np.save(f"{toLocation}/A_B_dic/{pdb}.npy", A_B_dic)

time_taken = time.time() - start_time  # time_taken is in seconds
hours, rest = divmod(time_taken,3600)
minutes, seconds = divmod(rest, 60)
print(f"---{hours} hours {minutes} minutes {seconds} seconds ---")