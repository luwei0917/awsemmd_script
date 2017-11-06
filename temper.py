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
import numpy as np
import pandas as pd

def process_data(folder, dis):
    location = folder + "/simulation/dis_{}/0/".format(dis)
    #     print(location)
    all_lipid_list = []
    for i in range(4):
        file = "lipid.{}.dat".format(i)
        lipid = pd.read_csv(location+file).assign(Run = i)
        lipid.columns = lipid.columns.str.strip()
        lipid = lipid[["Steps","Lipid","Run"]]
        all_lipid_list.append(lipid)
    lipid = pd.concat(all_lipid_list)

    all_dis_list = []
    for i in range(4):
        file = "addforce.{}.dat".format(i)
        dis = pd.read_csv(location+file).assign(Run = i)
        dis.columns = dis.columns.str.strip()
        remove_columns = ['AddedForce', 'Dis12', 'Dis34', 'Dis56']
        dis.drop(remove_columns, axis=1,inplace=True)
        all_dis_list.append(dis)
    dis = pd.concat(all_dis_list)

    all_wham_list = []
    for i in range(4):
        file = "wham.{}.dat".format(i)
        wham = pd.read_csv(location+file).assign(Run = i)
        wham.columns = wham.columns.str.strip()
        remove_columns = ['Rg', 'Tc']
        wham = wham.drop(remove_columns, axis=1)
        all_wham_list.append(wham)
    wham = pd.concat(all_wham_list)

    file = "../log.lammps"
    temper = pd.read_table(location+file, skiprows=2, sep=' ')
    temper = temper.melt(id_vars=['Step'], value_vars=['T0', 'T1', 'T2', 'T3'], value_name="Run", var_name="Temp")

    t2 = temper.merge(wham, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]
                        ).sort_values('Step').drop('Steps', axis=1)
    t3 = t2.merge(dis, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]
                     ).sort_values('Step').drop('Steps', axis=1)
    t4 = t3.merge(lipid, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]
                             ).sort_values('Step').drop('Steps', axis=1)
    t4 = t4.assign(TotalE=t4.Energy + t4.Lipid)
    t350 = t4.query('Temp=="T0" & Step > 1e7')
    t350.to_csv(location+"t350.dat", sep=' ', index=False, header=False)
    t400 = t4.query('Temp=="T1" & Step > 1e7')
    t400.to_csv(location+"t400.dat", sep=' ', index=False, header=False)
    t450 = t4.query('Temp=="T2" & Step > 1e7')
    t450.to_csv(location+"t450.dat", sep=' ', index=False, header=False)
    t4.query('Temp=="T3" & Step > 1e7').to_csv(location+"t500.dat", sep=' ', index=False, header=False)

# pre = "/Users/weilu/Research/server/oct_2017/week_of_oct02/freeEnergy"
# folder_list = glob.glob(pre+"/*")
dis_list = np.linspace(30, 130, 51)
folder = "/scratch/wl45/oct_2017/week_of_oct09/no_go_freeEnergy"
for dis in dis_list:
    process_data(folder, dis)
