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
from itertools import product
import re
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
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("-l", "--label", default="1")
parser.add_argument("--dimension", type=int, default=2)
parser.add_argument("-m", "--mode", type=int, default=1)
parser.add_argument("--force", type=float, default=1.0)
parser.add_argument("-p", "--patch", type=int, default=1)
parser.add_argument("--commons", type=int, default=0)
parser.add_argument("--nsample", type=int, default=2500)
parser.add_argument("--submode", type=int, default=-1)
args = parser.parse_args()

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir



if(args.test):
    print("hello world")

def expand_grid(dictionary):
    return pd.DataFrame([row for row in product(*dictionary.values())],
                        columns=dictionary.keys())
def readPMF(pre):
    perturbation_table = {0:"original", 1:"p_mem",
                          2:"m_mem", 3:"p_lipid",
                          4:"m_lipid", 5:"p_go",
                          6:"m_go", 7:"p_rg", 8:"m_rg"}
    pmf_list = {
        "perturbation":list(perturbation_table.keys()),
        "force":["0.0", "0.1", "0.2"]
    }
    pmf_list_data = expand_grid(pmf_list)
    all_pmf_list = []
    for index, row in pmf_list_data.iterrows():
        #     print(index)
        #     print("--")
        #     print(row)
        force = row["force"]
        perturbation = row["perturbation"]
        if perturbation == 0:
            location = pre + f"/force_{force}/pmf-*.dat"
            pmf_list = glob.glob(location)
            change = "none"
            upOrDown = "none"
        else:
            location = pre + f"/force_{force}/perturbation-{perturbation}-pmf-*.dat"
            pmf_list = glob.glob(location)
            change = perturbation_table[perturbation].split("_")[-1]
            upOrDown = perturbation_table[perturbation].split("_")[0]
        # print(location)
        name_list = ["f", "df", "e", "s"]
        names = ["bin", "x"] + name_list
        for location in pmf_list:
            # print(location)
            temp = re.findall(r'pmf-(\d+)', location)
            if len(temp) != 1:
                raise ValueError('Not expected to see more than one or none')
            else:
                temp = temp[0]
            data = pd.read_table(location, skiprows=2, sep='\s+', names=names).assign(upOrDown=upOrDown, change=change, force=force, temp=temp, perturbation=perturbation_table[perturbation])
            all_pmf_list.append(data)

    return pd.concat(all_pmf_list).dropna().reset_index()

data = readPMF(".")
remove_columns = ['bin']
data = data.drop(remove_columns, axis=1)
label = args.label
data.to_feather(f"/Users/weilu/Research/data/pulling/{datetime.datetime.today().strftime('%d_%h')}_pmf_{label}.feather")
