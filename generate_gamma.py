from pyCodeLib import *
import warnings
import glob
import re
import argparse
import sys
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(
    description="Generate Gammas in optimization.")

# parser.add_argument("proteins", help="The name of the protein list")
parser.add_argument("-p", "--protein_list", type=str, default="protein_list", help="The name of the protein list")
parser.add_argument("-m", "--mode", type=int, default=0)
parser.add_argument("-n", "--n_decoys", type=int, default=1000)
args = parser.parse_args()


if args.mode == 1:
    # from shuffle data.
    # complete_proteins = "iter0.txt"
    complete_proteins = args.protein_list
    n_decoys = args.n_decoys
    A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='shuffle',
                                    withBiased=False, num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, )
if args.mode == 2:
    # complete_proteins = "iter0.txt"
    complete_proteins = "protein_list_complete"
    # n_decoys = 50 * run_n
    n_decoys = args.n_decoys
    A, B, gamma, filtered_B, filtered_gamma, filtered_lamb, P, lamb = calculate_A_B_and_gamma_wl45(complete_proteins, "phi_list.txt", decoy_method='openMM',
                                    withBiased=True, oneMinus=True, decoyBiasName='decoysQ', num_decoys=n_decoys, noise_filtering=True, jackhmmer=False, read=False, mode=0, multiSeq=False, )
