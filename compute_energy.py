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
from pathlib import Path
from Bio.PDB.PDBParser import PDBParser
from pyCodeLib import *
from compute_energy_helperFunctions import *


parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    automatic copy the template file, \
    run simulation and analysis")

parser.add_argument("protein", help="The name of the protein")
parser.add_argument("-l", "--label", type=str, default="/Users/weilu/opt/ff_contact/1r69/")
parser.add_argument("-f", "--familyFold", action="store_true", default=False)
parser.add_argument("-m", "--membrane", action="store_true", default=False)
args = parser.parse_args()




if(platform.system() == 'Darwin'):
    pre = str(Path.home()) + "/openmmawsem/parameters/"
elif(platform.system() == 'Linux'):
    pre = "/projects/pw8/wl45/openawsem/parameters/"
else:
    print("system unkown")


# pre = "/Users/weilu/Research/server/sep_2019/saved_gammas/no_mix"
# pre = "/Users/weilu/Research/server/sep_2019/saved_gammas/Sep28_trial_5_cutoff500_impose_Aprime_constraint"
if args.membrane:
    gamma_direct, gamma_mediated = read_gamma(f"{pre}/membrane_gamma.dat")
else:
    gamma_direct, gamma_mediated = read_gamma(f"{pre}/gamma.dat")
burial_gamma = np.loadtxt(f"{pre}/burial_gamma.dat").T

gamma_ijm, water_gamma_ijm, protein_gamma_ijm = change_gamma_format(gamma_direct, gamma_mediated)



pdb = (args.protein).split(".")[0]
structure = parse_pdb(pdb)

if False:
    e_pap_aph, e_pap_ap = compute_pap1(structure)
    e_pap_p = compute_pap2(structure)
    print("PAP_APH: ", e_pap_aph, "PAP_AP: ", e_pap_ap, "PAP_P: ", e_pap_p, "PAP_total: ", e_pap_p+e_pap_ap+e_pap_aph)

if False:
    e_membrane = compute_membrane(structure)
    print("Membrane energy", e_membrane)
    # exit()
if False:
    e_debye_huckel = compute_debye_huckel(structure)
    print("Debye Huckel", e_debye_huckel)
# e_chi = compute_chi(structure)
# print("Chi", e_chi)
if False:
    e_mediated = compute_mediated(structure, hasPhosphorylation=True)
    e_direct = compute_direct(structure, hasPhosphorylation=True)
    e_burial = compute_burial(structure, hasPhosphorylation=True)

if True:
    e_mediated = compute_mediated(structure, protein_gamma_ijm, water_gamma_ijm, hasPhosphorylation=False)
    e_direct = compute_direct(structure, gamma_ijm, hasPhosphorylation=False)
    e_burial = compute_burial(structure, burial_gamma, hasPhosphorylation=False)
# print("Mediated, Direct, Mediated+Direct, Burial, Mediated+Direct+Burial")
name_list = "Protein, Mediated, Direct, Mediated+Direct, Burial, Mediated+Direct+Burial".split(",")
out_line = " ".join([f"{i:<20}" for i in name_list])

print(out_line)
# print(e_mediated, e_direct, e_mediated + e_direct, e_burial, e_mediated + e_direct + e_burial)

energy_out_list = [-e_mediated, -e_direct, -(e_mediated+e_direct), -e_burial, -(e_mediated + e_direct + e_burial)]
out_line = " ".join([f"{pdb:<20}"] + ["{0:<20.3f}".format(i) for i in energy_out_list])
print(out_line)

# print(compute_direct_2("8ab.pdb"))
# exit()
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
    energy_out_list = [-e_mediated_ff, -e_direct_ff, -(e_direct_ff+e_mediated_ff), -e_burial_ff, -(e_direct_ff+e_mediated_ff+e_burial_ff)]
    out_line = " ".join([f"{pdb:<20}"] + ["{0:<20.3f}".format(i) for i in energy_out_list])
    print(out_line)


exit()
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