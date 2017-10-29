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
import numpy as np
import fileinput
from itertools import product
import pandas as pd
# from run_parameter import *

do = os.system
cd = os.chdir
parser = argparse.ArgumentParser(
    description="modify the template into working input file")

parser.add_argument("--name", default="2xov_eval.in", help="The name of the in file")

args = parser.parse_args()

type_list = {"MY_RG":float,
                "MY_LIPID":float,
                "MY_RUN":int}
modify_list = {"MY_RG":[0.4],
                "MY_LIPID":[2],
                "MY_RUN":list(range(12))}

def change(fileName, from_str, to_str):
    with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
        for line in file:
            tmp = line
            tmp = tmp.replace(from_str, str(to_str))
            print(tmp, end='')

def modify(fileName):
    for key, value in modify_list.items():
        change(fileName, key, value)

def expand_grid(dictionary):
    return pd.DataFrame([row for row in product(*dictionary.values())],
                        columns=dictionary.keys())

# Decide the name of the folders
folder_name_list = []
for key, value_list in modify_list.items():
    # print(key, value_list, len(value_list))
    if len(value_list) > 1:
        folder_name_list.append(key)

data = expand_grid(modify_list)
for index, row in data.iterrows():
    folder = ""
    for name in folder_name_list:
        data = type_list[name](row[name])
        folder += name[3:] + "_{}".format(data)
        # print(type(row[name]))
    do("mkdir -p " + folder)
    cd(folder)
    do("cp ~/opt/2xov_eval/* .")
    for key in modify_list.keys():
        change(args.name, key, type_list[key](row[key]))
        # print(row[key])
    do("sbatch run.slurm")
    cd("..")
    print(folder)

# for key, value_list in modify_list.items():
#     print(key, value_list, len(value_list))
#     if len(value_list) > 1:
#         for value in value_list:
#             my_folder = folder.format(value)
# print(folder)
# do("mkdir -p ")
# modify(args.name)
