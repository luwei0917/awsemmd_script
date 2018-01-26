#!/usr/bin/env python3
# @Author: Wei Lu <weilu>
# @Date:   25-Jan-2017
# @Email:  wl45@rice.edu
# @Last modified by:   weilu
# @Last modified time: 25-Jan-2017
# @Copyright: Free


import os
import sys
import argparse
import platform
from datetime import datetime
import imp
from myPersonalFunctions import *
import glob
import numpy as np
import datetime
import pandas as pd
import re
from small_script.myFunctions import readPMF
from small_script.myFunctions import readPMF_2
from small_script.myFunctions import readPMF_basic
# Useful codes
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# tr " " "\n"
# sed 1d
# sort -u -k 3
# sed -e 's/+T//'
parser = argparse.ArgumentParser(description="Analysis code, need run multiple times")

# parser.add_argument("protein", help="the name of protein")
# parser.add_argument("template", help="the name of template file")
parser.add_argument("-t", "--test", help="test ", action="store_true", default=False)
parser.add_argument("-d", "--day", type=str, default="someday")
parser.add_argument("-l", "--label", default="1")
parser.add_argument("--dimension", type=int, default=2)
parser.add_argument("-m", "--mode", type=int, default=0)
parser.add_argument("--force", type=float, default=1.0)
parser.add_argument("-p", "--patch", type=int, default=1)
parser.add_argument("--commons", type=int, default=0)
parser.add_argument("--nsample", type=int, default=2500)
parser.add_argument("--submode", type=int, default=-1)
args = parser.parse_args()


do = os.system
cd = os.chdir

label = args.label
if args.day == "jan22":
    if args.mode ==1:
        all_pmf_list = []
        # simulation_list = ["next_gen_native_based_memb_3_rg_0.2_lipid_0.6_extended",
        #                     "next_gen_native_based_memb_3_rg_0.4_lipid_0.6_extended",
        #                     "next_gen_native_based_memb_3_rg_0.4_lipid_0.6_topology"]
        simulation_list = ["rg_0.3_lipid_0.6_mem_1_350-550"]

        for simulation in simulation_list:
            cd(simulation)
            tmp = readPMF_2(".", mode_list=["1d_dis", "1d_qw", "1d_z"]).assign(submode="0", simulation=simulation)
            all_pmf_list.append(tmp)
            cd("..")
        data = pd.concat(all_pmf_list).reset_index(drop=True)


if(args.test):
    print("hello world")
if args.day == "old":
    if args.mode == -1:
        data = readPMF_basic(".")
    if args.mode == 0:
        data = readPMF(".")
    elif args.mode == 1:
        data = readPMF_2(".")  # read both 1d_dis and 1d_qw
    elif args.mode ==3:
        all_pmf_list = []
        simulation_list = ["memb_3_rg_0.1_lipid_1_extended", "memb_3_rg_0.1_lipid_1_topology"]
        for simulation in simulation_list:
            for mode in range(3):
                cd(f"nov_15_all_freeEnergy_calculation_sample_range_mode_{mode}")
                cd(simulation)
                tmp = readPMF_2(".").assign(submode=mode, simulation=simulation)
                all_pmf_list.append(tmp)
                cd("../..")
        data = pd.concat(all_pmf_list).reset_index(drop=True)
    elif args.mode ==4:
        all_pmf_list = []
        # simulation_list = ["next_gen_native_based_memb_3_rg_0.2_lipid_0.6_extended",
        #                     "next_gen_native_based_memb_3_rg_0.4_lipid_0.6_extended",
        #                     "next_gen_native_based_memb_3_rg_0.4_lipid_0.6_topology"]
        simulation_list = ["no_side_contraint_memb_3_rg_0.4_lipid_0.6_extended_no_energy"]

        cd(f"no_side_contraint_memb_3_rg_0.4_lipid_0.6_extended_0")
        for simulation in simulation_list:
            cd(simulation)
            tmp = readPMF_2(".").assign(submode="0", simulation=simulation)
            all_pmf_list.append(tmp)
            cd("..")
        data = pd.concat(all_pmf_list).reset_index(drop=True)
    # data = readPMF(".")
    elif args.mode ==5:
        all_pmf_list = []
        simulation_list = ["new_next_gen_native_based_memb_3_rg_0.4_lipid_0.6_extended", "new_next_gen_native_based_memb_3_rg_0.4_lipid_0.6_extended_2", "new_next_gen_native_based_memb_3_rg_0.4_lipid_0.6_extendedshort"]
        for simulation in simulation_list:
            for mode in range(3):
                cd(f"nov_18_all_freeEnergy_calculation_sample_range_mode_{mode}")
                cd(simulation)
                tmp = readPMF_2(".").assign(submode=mode, simulation=simulation)
                all_pmf_list.append(tmp)
                cd("../..")
        data = pd.concat(all_pmf_list).reset_index(drop=True)

remove_columns = ['bin']
data = data.drop(remove_columns, axis=1)
fileName = f"/Users/weilu/Research/data/pulling/{datetime.datetime.today().strftime('%d_%h_%H%M%S')}_pmf_{label}.feather"
print(fileName)
data.to_feather(fileName)
