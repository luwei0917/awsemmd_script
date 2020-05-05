from pyCodeLib import *
import warnings
import glob
import re
import argparse
import sys
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(
    description="Generate decoys by shuffling.")

parser.add_argument("proteins", help="The name of the protein list")
parser.add_argument("-m", "--mode", type=int, default=0)
parser.add_argument("-n", "--n_decoys", type=int, default=1000)
args = parser.parse_args()

if args.mode == 1:
    proteins = args.proteins
    n_decoys = args.n_decoys
    generate_decoy_sequences(proteins, methods=['shuffle'], num_decoys=[n_decoys], databaseLocation="../../../")

