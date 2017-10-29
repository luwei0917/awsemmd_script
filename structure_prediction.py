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
from myPersonalFunctions import *
import glob
from small_script.myFunctions import *
from small_script.extract_pdb import *
# Useful codes
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# tr " " "\n"
# sed 1d
# sort -u -k 3
# sed -e 's/+T//'
# awk '$5=-$5' data
# awk '{print $2}' wham.dat |  sed 's/,$//' | sed 1d > qw.dat
mypath = os.environ["PATH"]
os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
my_env = os.environ.copy()

parser = argparse.ArgumentParser(description="This is my playground for current project")
# parser.add_argument("protein", help="the name of protein")
# parser.add_argument("template", help="the name of template file")
parser.add_argument("-t", "--test", help="test ", action="store_true", default=False)
parser.add_argument("-n", "--number", type=int, default=10, help="number of run")
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("-m", "--mode",
                    type=int, default=0)
args = parser.parse_args()


do = os.system
cd = os.chdir

# convert \( 0.0.png 0.2.png 0.4.png  +append \) \( 0.6.png  0.8.png    1.0.png +append \) -background none -append final.png
def test():
    print("don't show me")

if(args.mode == 1):
    location = "/Users/weilu/opt/crystal_structures/membrane_proteins/"
    pdbFileDataBase = "http://opm.phar.umich.edu/pdb/"
    protein_info_list = []
    protein_info_list.append((location, "2xov","A",91,271))
    # protein_info_list.append((location, "1occ","C",71,261))
    # protein_info_list.append((location, "1pv6", "A", 1, 190))
    # protein_info_list.append((location, "2bl2", "A", 12, 156))
    # protein_info_list.append((location, "2bg9", "A", 211, 301))
    # protein_info_list.append((location, "1j4n", "A", 4, 119))
    # protein_info_list.append((location, "1py6", "A", 77, 199))
    for (location, protein, chain, residue_start, residue_end) in protein_info_list:
        do("wget {0}{1} -O ~/opt/crystal_structures/membrane_proteins/original_pdb/{1}".format(pdbFileDataBase, protein+".pdb"))
        extract_pdb(location, protein, chain, residue_start, residue_end)

if(args.mode == 2):
    protein_list = ["1occ", "1pv6", "2bl2", "2bg9", "1j4n", "1py6"]
    for protein in protein_list:
        do("mkdir -p " + protein)
        cd(protein)
        do("mkdir -p " + protein)
        cd(protein)
        do("cp ~/opt/crystal_structures/membrane_proteins/for_simulation/{}.pdb crystal_structure.pdb ".format(protein))
        do("create_project.py {} --frag --membrane > generation.log".format(protein))
        check_and_correct_fragment_memory()
        cd("../..")

if(args.mode ==3):
    protein_list = ["1occ", "1pv6", "2bl2", "2bg9", "1j4n", "1py6"]
    for protein in protein_list:
        structure_prediction_run(protein)
