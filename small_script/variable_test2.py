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

fileName = "2xov_multi.in"
# if args.rerun == 0:
#     start_from = "read_data data.2xov"
# if args.rerun == 1:
#     start_from = "read_restart restart.extended"


# simulation_steps = 5e7
def expand_grid(dictionary):
    return pd.DataFrame([row for row in product(*dictionary.values())],
                        columns=dictionary.keys())


def variable_test2(k_list=[1],
                      force_ramp_rate_list=[1],
                      memb_k_list=[1],
                      force_list=["ramp"],
                      rg_list=[0.08],
                      pressure_list=[0.1],
                      repeat=1,
                      mode_list=[2],
                      commons=0,
                      temperature_list=[300],
                      start_from_list=["native"],
                      simulation_model_list=["go"]):
    inputs = locals()
    tmp = {}  # some none list one are not important
    # variables that changing determines the name of folder
    folder_name_list = []
    for key,value in inputs.items():
        if isinstance(value, list):
            tmp[key] = value
            if len(value) > 1:
                folder_name_list.append(key)
                # folder_name_template += key.replace("_list", "") + "_{"+key.replace("_list", "")+"}_"
    all_inputs = expand_grid(tmp)

    for index, row in all_inputs.iterrows():
        # print(index, row)
        folder_name_template = ""
        for name in row.index:
            exec(name.replace("_list", "")+"= '"+str(row[name]) + "'")
        # k = row["k_list"]
        # force_ramp_rate = row["force_ramp_rate_list"]
        # memb_k = row["memb_k_list"]
        # force = row["force_list"]
        # rg = row["rg_list"]
        # pressure = row["pressure_list"]
        # mode = row["mode_list"]
        # temperature = row["temperature_list"]
        # simulation_model = row["simulation_model_list"]
        # start_from = row["start_from_list"]

        for name in folder_name_list:
            tmp = row[name]
            folder_name_template += name.replace("_list", "")+f"_{tmp}_"
        # print(folder_name_template)
        if folder_name_template == "":
            folder_name_template = "simulation"
        print(folder_name_template)

        do("mkdir "+folder_name_template)
        do("cp -r 2xov " + folder_name_template + "/")
        cd(folder_name_template + "/2xov")

        if simulation_model == "go":
            fixFile = "fix_backbone_coeff_go.data"
            simulation_steps = 1e7/force_ramp_rate
        if simulation_model == "single":
            fixFile = "fix_backbone_coeff_single.data"
            simulation_steps = 5e6/force_ramp_rate

        if start_from == "native":
            start_from_input = "read_data data.2xov"
        if start_from == "extended":
            start_from_input = "read_restart restart.extended"
        if start_from == "topology":
            start_from_input = "read_restart restart.native_topology"

        with fileinput.FileInput(fixFile, inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace("MY_MEMB_K", str(memb_k)), end='')
        with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
            for line in file:
                tmp = line
                tmp = tmp.replace("MY_K", str(k))
                tmp = tmp.replace("START_FROM", str(start_from_input))
                tmp = tmp.replace("MY_FIX", str(fixFile))
                tmp = tmp.replace("MY_TEMPERATURE", str(temperature))
                tmp = tmp.replace("MY_PRESSURE", str(pressure))
                tmp = tmp.replace("MY_FORCE", str(force))
                tmp = tmp.replace("MY_RG", str(rg))
                tmp = tmp.replace("RATE", str(force_ramp_rate))
                tmp = tmp.replace("SIMULATION_STEPS", str(int(simulation_steps)))
                # tmp = tmp.replace("SIMULATION_STEPS", str(int(simulation_steps/force_ramp_rate)))
                print(tmp, end='')
        cd("..")

        do("run.py -m {} 2xov -n {} --commons {}".format(mode, repeat, commons))
        # do("run.py -m 2 2xov --start extended -n {}".format(repeat))
        cd("..")
