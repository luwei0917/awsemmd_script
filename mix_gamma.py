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
import glob
from time import sleep
import fileinput
import numpy as np
import pandas as pd
from small_script.variable_test import variable_test
from small_script.variable_test2 import variable_test2
import subprocess
from small_script.myFunctions import compute_theta_for_each_helix
from small_script.myFunctions import *
from collections import defaultdict


parser = argparse.ArgumentParser(description="This is my playground for current project")
parser.add_argument("preGamma", help="preGamma")
parser.add_argument("gamma", help="Gamma")
# parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("-a", "--alpha", type=float, default=0.3)
parser.add_argument("-i", "--iter", type=str, default="-1")
parser.add_argument("-o", "--onlyToSimulation", action="store_true", default=False)
parser.add_argument("-s", "--scale", action="store_true", default=False)
args = parser.parse_args()

with open('mix_gamma_cmd.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')

i = args.iter
alpha = args.alpha
pre = "./"
# g = "gammas/proteins_name_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_gamma_filtered"
# Gamma = np.loadtxt(f"/Users/weilu/Research/server/april_2019/optimization_restart_iter{i}/{g}")
# preGamma = np.loadtxt(f"/Users/weilu/Research/server/april_2019/complete_gammas/iter_{pre_i}")
Gamma = np.loadtxt(args.gamma)


if args.onlyToSimulation:
    gamma_for_simulation = pre + f"for_simulation/iteration_{i}_gamma.dat"
    burial_gamma_for_simulation = pre + f"for_simulation/iteration_{i}_burial_gamma.dat"
    gamma_format_convertion_iteration_to_simulation(Gamma, gamma_for_simulation, burial_gamma_for_simulation=burial_gamma_for_simulation)
    exit()

preGamma = np.loadtxt(args.preGamma)

# if alpha == 0.3:
#     mix_gammas_3(pre, Gamma, preGamma, alpha=alpha, scale=args.scale iterGammaName=f"iter_{i}", iteration=f"iter_{i}")
# else:
#     mix_gammas_3(pre, Gamma, preGamma, alpha=alpha, scale=args.scale, iterGammaName=f"iter_{i}_{int(alpha*100)}", iteration=f"iter_{i}")

mix_gammas_3(pre, Gamma, preGamma, alpha=alpha, scale=args.scale, iterGammaName=f"iter_{i}_{int(alpha*100)}", iteration=f"iter_{i}")

os.system(f"mv iteration_iter_{i}_* for_simulation/")
