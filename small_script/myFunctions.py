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
import subprocess
import glob
import re
from small_script.myFunctions_helper import *
import numpy as np
import pandas as pd
import fileinput
from itertools import product
# compute cross Q for every pdb pair in one folder
# parser = argparse.ArgumentParser(description="Compute cross q")
# parser.add_argument("-m", "--mode",
#                     type=int, default=1)

# args = parser.parse_args()
def expand_grid(dictionary):
    return pd.DataFrame([row for row in product(*dictionary.values())],
                        columns=dictionary.keys())

def make_metadata(k=1000.0, temps_list=["450"]):
    cwd = os.getcwd()
    files = glob.glob("../data/*")
    kconstant = k
    with open("metadatafile", "w") as out:
        for oneFile in files:
            tmp = oneFile.split("/")[-1].replace('.dat', '')
            t = tmp.split("_")[1]
            bias = tmp.split("_")[3]
            # print(tmp)
            # if int(float(dis)) > 150:
            #     continue
            if t in temps_list:
                target = "../{} {} {} {}\n".format(oneFile, t, kconstant, bias)
                out.write(target)
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

def readPMF_2(pre):
    mode_list = ["1d_dis", "1d_qw"]
    all_data_list =[]
    for mode in mode_list:
        tmp = readPMF(mode).assign(mode=mode)
        all_data_list.append(tmp)
    return pd.concat(all_data_list).dropna().reset_index()

def shrinkage(n=552, shrink_size=6, max_frame=2000, fileName="dump.lammpstrj"):
    print("Shrinkage: size: {}, max_frame: {}".format(shrink_size, max_frame))
    bashCommand = "wc " + fileName
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    line_number = int(output.decode("utf-8").split()[0])
    print(line_number)
    print(line_number/552)
    # number of atom = 543
    n = 552
    count = 0
    with open("small.lammpstrj", "w") as out:
        with open(fileName, "r") as f:
            for i, line in enumerate(f):
                if (i // n) % shrink_size == 0:
                    if count >= max_frame*n:
                        break
                    count += 1
                    out.write(line)

def compute_theta_for_each_helix(dumpName="../dump.lammpstrj.0"):
    print("This is for 2xov only")
    helices_list = [(94,114), (147,168), (171, 192), (200, 217), (226, 241), (250, 269)]
    atoms_all_frames = read_lammps(dumpName)
    # print(atoms[0])
    # print(len(atoms), len(atoms[0]))
    # helices_angles_all_frames = []
    with open("angles.csv", "w") as out:
        out.write("Frame, Helix, Angle\n")
        for ii, frame in enumerate(atoms_all_frames):
            # helices_angles = []
            for count, (i, j) in enumerate(helices_list):
                # print(i, j)
                i = i-91
                j = j-91
                # end - start
                a = np.array(frame[j]) - np.array(frame[i])
                b = np.array([0, 0, 1])
                angle = a[2]/length(a)  # in form of cos theta
                # helices_angles.append(angle)
                # print(angle)
                out.write("{}, {}, {}\n".format(ii, count+1, angle))
            # helices_angles_all_frames.append(helices_angles)


def structure_prediction_run(protein):
    print(protein)
    protocol_list = ["awsemer", "frag", "er"]
    do = os.system
    cd = os.chdir
    cd(protein)
    # run = "frag"
    for protocol in protocol_list:
        do("rm -r " + protocol)
        do("mkdir -p " + protocol)
        do("cp -r {} {}/".format(protein, protocol))
        cd(protocol)
        cd(protein)
        # do("cp ~/opt/gremlin/protein/{}/gremlin/go_rnativeC* .".format(protein))
        do("cp ~/opt/gremlin/protein/{}/raptor/go_rnativeC* .".format(protein))
        fileName = protein + "_multi.in"
        backboneFile = "fix_backbone_coeff_" + protocol
        with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
            for line in file:
                tmp = line.replace("fix_backbone_coeff_er", backboneFile)
                print(tmp, end='')
        cd("..")
        do("run.py -m 0 -n 20 {}".format(protein))
        cd("..")
    cd("..")
    # do("")


def check_and_correct_fragment_memory():
    with open("tmp.mem", "w") as out:
        with open("fragsLAMW.mem", "r") as f:
            for i in range(4):
                line = next(f)
                out.write(line)
            for line in f:
                gro, _, i, n, _ = line.split()
                delete = False
                # print(gro, i, n)
                # name = gro.split("/")[-1]
                with open(gro, "r") as one:
                    next(one)
                    next(one)
                    all_residues = set()
                    for atom in one:
                        residue, *_ = atom.split()
                        # print(residue)
                        all_residues.add(int(residue))
                    for test in range(int(i), int(i)+int(n)):
                        if test not in all_residues:
                            print("ATTENTION", gro, i, n, "missing:",test)
                            delete = True
                if not delete:
                    out.write(line)
    os.system("mv fragsLAMW.mem fragsLAMW_back")
    os.system("mv tmp.mem fragsLAMW.mem")


def read_temper(n=4, location=".", rerun=-1):
    all_lipid_list = []
    for i in range(n):
        file = "lipid.{}.dat".format(i)
        lipid = pd.read_csv(location+file).assign(Run=i)
        lipid.columns = lipid.columns.str.strip()
        lipid = lipid[["Steps","Lipid","Run"]]
        all_lipid_list.append(lipid)
    lipid = pd.concat(all_lipid_list)

    all_energy_list = []
    for i in range(n):
        file = "energy.{}.dat".format(i)
        energy = pd.read_csv(location+file).assign(Run=i)
        energy.columns = energy.columns.str.strip()
        energy = energy[["Steps", "AMH-Go", "Membrane", "Rg", "Run"]]
        all_energy_list.append(energy)
    energy = pd.concat(all_energy_list)

    all_dis_list = []
    for i in range(n):
        file = "addforce.{}.dat".format(i)
        dis = pd.read_csv(location+file).assign(Run=i)
        dis.columns = dis.columns.str.strip()
        remove_columns = ['AddedForce', 'Dis12', 'Dis34', 'Dis56']
        dis.drop(remove_columns, axis=1,inplace=True)
        all_dis_list.append(dis)
    dis = pd.concat(all_dis_list)

    all_wham_list = []
    for i in range(n):
        file = "wham.{}.dat".format(i)
        wham = pd.read_csv(location+file).assign(Run=i)
        wham.columns = wham.columns.str.strip()
        remove_columns = ['Rg', 'Tc']
        wham = wham.drop(remove_columns, axis=1)
        all_wham_list.append(wham)
    wham = pd.concat(all_wham_list)
    if rerun == -1:
        file = "../log.lammps"
    else:
        file = f"../log{rerun}/log.lammps"
    temper = pd.read_table(location+file, skiprows=2, sep=' ')
    temper = temper.melt(id_vars=['Step'], value_vars=['T' + str(i) for i in range(n)], value_name="Temp", var_name="Run")
    temper["Run"] = temper["Run"].str[1:].astype(int)
    temper["Temp"] = "T" + temper["Temp"].astype(str)
    t2 = temper.merge(wham, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]).sort_values('Step').drop('Steps', axis=1)
    t3 = t2.merge(dis, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]).sort_values('Step').drop('Steps', axis=1)
    t4 = t3.merge(lipid, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]).sort_values('Step').drop('Steps', axis=1)
    t5 = t4.merge(energy, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]).sort_values('Step').drop('Steps', axis=1)
    t6 = t5.assign(TotalE=t5.Energy + t5.Lipid)
    return t6

def process_temper_data(pre, data_folder, folder_list, rerun=-1, n=12, bias="dis"):
    print("process temp data")
    for folder in folder_list:
        simulation_list = glob.glob(pre+folder+f"/simulation/{bias}_*")
        print(pre+folder+f"/simulation/{bias}_*")
        os.system("mkdir -p " + pre+folder+"/data")
        for one_simulation in simulation_list:
            bias_num = one_simulation.split("_")[-1]
            print(bias_num, "!")
            if rerun == -1:
                location = one_simulation + "/0/"
                try:
                    data = read_temper(location=location, n=n)
                    # remove_columns = ['Step', "Run"]
                    # data = data.drop(remove_columns, axis=1)
                    data.reset_index().to_feather(pre+folder+"/data/"+f"{bias}{bias_num}.feather")
                except:
                    print("Unexpected error:", sys.exc_info()[0])
                    print("notrun?", dis)
            else:
                all_data_list = []
                for i in range(rerun):
                    location = one_simulation + f"/{i}/"
                    try:
                        data = read_temper(location=location, n=n, rerun=i)
                        # remove_columns = ['Step', "Run"]
                        # data = data.drop(remove_columns, axis=1)
                        all_data_list.append(data)
                    except:
                        print("Unexpected error:", sys.exc_info()[0])
                        print("notrun?", bias_num)
                try:
                    data = pd.concat(all_data_list)
                    data.reset_index().to_feather(pre+folder+"/data/"+f"{bias}{bias_num}.feather")
                except:
                    print("Unexpected error:", sys.exc_info()[0])
                    print("not data?", bias_num)
    #         temps = list(dic.keys())
        os.system("mv "+pre+folder+"/data "+data_folder+folder)

def move_data(data_folder, freeEnergy_folder, folder, sub_mode_name="", kmem=0.2, klipid=0.1, kgo=0.1, krg=0.2, sample_range_mode=0, bias="dis"):
    print("move data")
    os.system("mkdir -p "+freeEnergy_folder+folder+sub_mode_name+"/data")
    dis_list = glob.glob(data_folder+folder+f"/{bias}*.feather")
    for dis_file in dis_list:
        dis = dis_file.split("/")[-1].replace(bias, '').replace('.feather', '')
        print(dis)
        t6 = pd.read_feather(dis_file)
        remove_columns = ['index']
        t6 = t6.drop(remove_columns, axis=1)
        t6 = t6.assign(TotalE_perturb_mem_p=t6.TotalE + kmem*t6.Membrane)
        t6 = t6.assign(TotalE_perturb_mem_m=t6.TotalE - kmem*t6.Membrane)
        t6 = t6.assign(TotalE_perturb_lipid_p=t6.TotalE + klipid*t6.Lipid)
        t6 = t6.assign(TotalE_perturb_lipid_m=t6.TotalE - klipid*t6.Lipid)
        t6 = t6.assign(TotalE_perturb_go_p=t6.TotalE + kgo*t6["AMH-Go"])
        t6 = t6.assign(TotalE_perturb_go_m=t6.TotalE - kgo*t6["AMH-Go"])
        t6 = t6.assign(TotalE_perturb_rg_p=t6.TotalE + krg*t6.Rg)
        t6 = t6.assign(TotalE_perturb_rg_m=t6.TotalE - krg*t6.Rg)
        # t6["TotalE"] = 0
        dic = {"T0":350, "T1":400, "T2":450, "T3":500, "T4":550, "T5":600, "T6":650, "T7":700, "T8":750, "T9":800, "T10":900, "T11":1000}
        temps = list(dic.values())
        def convert(x):
            return dic[x]
        t6["Temp"] = t6["Temp"].apply(convert)

        for temp in temps:
            if temp > 800:
                continue
            if sample_range_mode == 0:
                tmp = t6.query('Temp=="{}"& Step > 1e7 & Step <= 2e7'.format(temp))
            elif sample_range_mode == 1:
                tmp = t6.query('Temp=="{}"& Step > 2e7 & Step <= 3e7'.format(temp))
            elif sample_range_mode == 2:
                tmp = t6.query('Temp=="{}"& Step > 3e7 & Step <= 4e7'.format(temp))
            tmp.to_csv(freeEnergy_folder+folder+sub_mode_name+"/data/t_{}_{}_{}.dat".format(temp, bias, dis), sep=' ', index=False, header=False)
# def pick_out_and_show():
#     protein_list = ["1occ", "1pv6", "2bl2", "2bg9", "1j4n", "1py6"]
#     for protein in protein_list:
#         frames = [pd.read_csv("{}/awsemer/simulation/{}/0/wham.dat".format(protein, run)).assign(Run=run) for run in range(20)]
#         result = pd.concat(frames)
#         answer = result.iloc[result[' Qw'].argsort()].iloc[-1]
#         print(protein, answer.Steps, answer.Run)
#         os.chdir("{}/awsemer/simulation/{}/0/".format(protein, int(answer.Run)))
#         os.system("show.py --frame {} {} -p".format(int(answer.Steps/4000), protein))
#         os.chdir("../../../../../")
