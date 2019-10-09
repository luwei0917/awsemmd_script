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
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBIO import Select
from Bio.PDB import PDBList
from Bio.PDB.Polypeptide import three_to_one
# from pdbfixer import PDBFixer
# from simtk.openmm.app import PDBFile

# compute cross Q for every pdb pair in one folder
# parser = argparse.ArgumentParser(description="Compute cross q")
# parser.add_argument("-m", "--mode",
#                     type=int, default=1)

# args = parser.parse_args()

def getFromTerminal(CMD):
    return subprocess.Popen(CMD,stdout=subprocess.PIPE,shell=True).communicate()[0].decode()


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [atoi(c) for c in re.split(r'(\d+)', text)]

def get_data(pre, pdb_list, simType="all_simulations", n_rum=30, rerun=1, formatName=True):
    # to get last 50 frame of each run
    _all = []
    for p in pdb_list:
        if formatName:
            name = p.lower()[:4]
        else:
            name = p
        for i in range(n_rum):
            for ii in range(rerun):
                location = pre + f"{simType}/{name}/simulation/{i}/{ii}/wham.dat"
                try:
                    tmp = pd.read_csv(location).tail(50).reset_index()
                    tmp.columns = tmp.columns.str.strip()
                    _all.append(tmp.assign(Run=i, Name=name, Rerun=ii))
                except Exception as e:
                    print(e)
    data = pd.concat(_all)
    data["Run"] = "Run" + data["Run"].astype(str)
    return data

def computeRg(pdb_file, chain="A"):
    # compute Radius of gyration
    # pdb_file = f"/Users/weilu/Research/server/feb_2019/iterative_optimization_new_temp_range/all_simulations/{p}/{p}/crystal_structure.pdb"
    chain_name = chain
    parser = PDBParser()
    structure = parser.get_structure('X', pdb_file)
    chain = list(structure[0][chain_name])
    n = len(chain)
    rg = 0.0
    for i, residue_i in enumerate(chain):
        for j, residue_j in enumerate(chain[i+1:]):
            try:
                r = residue_i["CA"] - residue_j["CA"]
            except:
                print(residue_i, residue_j)
            rg += r**2
    return (rg/(n**2))**0.5

def splitPDB(pre, fileName):
    location = f"{pre}/{fileName}"
    with open(location, "r") as f:
        a = f.readlines()
    print(len(a))
    i = 0
    tmp = ""
    for line in a:
        tmp += line
    #     os.system(f"echo '{line}' >> {pre}frame{i}")
        if line == "END\n":
            with open(f"{pre}/frame{i}.pdb", "w") as out:
                out.write(tmp)
            i += 1
            tmp = ""
        else:
            if line.strip() == "END":
                print("problem")

def read_hydrophobicity_scale(seq, isNew=False):
    seq_dataFrame = pd.DataFrame({"oneLetterCode":list(seq)})
    # HFscales = pd.read_table("~/opt/small_script/Whole_residue_HFscales.txt")
    HFscales = pd.read_csv("~/opt/small_script/Whole_residue_HFscales.txt", sep="\t")
    if not isNew:
        # Octanol Scale
        # new and old difference is at HIS.
        code = {"GLY" : "G", "ALA" : "A", "LEU" : "L", "ILE" : "I",
                "ARG+" : "R", "LYS+" : "K", "MET" : "M", "CYS" : "C",
                "TYR" : "Y", "THR" : "T", "PRO" : "P", "SER" : "S",
                "TRP" : "W", "ASP-" : "D", "GLU-" : "E", "ASN" : "N",
                "GLN" : "Q", "PHE" : "F", "HIS+" : "H", "VAL" : "V",
                "M3L" : "K", "MSE" : "M", "CAS" : "C"}
    else:
        code = {"GLY" : "G", "ALA" : "A", "LEU" : "L", "ILE" : "I",
                "ARG+" : "R", "LYS+" : "K", "MET" : "M", "CYS" : "C",
                "TYR" : "Y", "THR" : "T", "PRO" : "P", "SER" : "S",
                "TRP" : "W", "ASP-" : "D", "GLU-" : "E", "ASN" : "N",
                "GLN" : "Q", "PHE" : "F", "HIS0" : "H", "VAL" : "V",
                "M3L" : "K", "MSE" : "M", "CAS" : "C"}
    HFscales_with_oneLetterCode = HFscales.assign(oneLetterCode=HFscales.AA.str.upper().map(code)).dropna()
    data = seq_dataFrame.merge(HFscales_with_oneLetterCode, on="oneLetterCode", how="left")
    return data

def create_zim(seqFile, isNew=False):
    a = seqFile
    seq = getFromTerminal("cat " + a).rstrip()
    data = read_hydrophobicity_scale(seq, isNew=isNew)
    z = data["DGwoct"].values
    np.savetxt("zim", z, fmt="%.2f")


def expand_grid(dictionary):
    return pd.DataFrame([row for row in product(*dictionary.values())],
                        columns=dictionary.keys())


def duplicate_pdb(From, To, offset_x=0, offset_y=0, offset_z=0, new_chain="B"):
    with open(To, "w") as out:
        with open(From, "r") as f:
            for line in f:
                tmp = list(line)
                if len(tmp) < 26:
                    out.write(line)
                    continue
                atom = line[0:4]
                atomSerialNumber = line[6:11]
                atomName = line[12:16]
                atomResidueName = line[17:20]
                chain = line[21]
                residueNumber = line[22:26]
                # change chain A to B
                # new_chain = "B"

                if atom == "ATOM":
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    # add 40 to the x
                    new_x = x + offset_x
                    new_y = y + offset_y
                    new_z = z + offset_z
                    if new_chain == -1:
                        pass
                    else:
                        tmp[21] = new_chain
                    tmp[30:38] = "{:8.3f}".format(new_x)
                    tmp[38:46] = "{:8.3f}".format(new_y)
                    tmp[46:54] = "{:8.3f}".format(new_z)

                a = "".join(tmp)
                out.write(a)

def fill_chain(From, To, chain="A"):
    # add the missing chain info
    with open(To, "w") as out:
        with open(From, "r") as f:
            for line in f:
                if len(line) < 21:
                    continue
                tmp = list(line)
                if line[0:4] == "ATOM":
                    tmp[21] = "A"
                a = "".join(tmp)
                out.write(a)

def get_two_part_from_prediction(topFile):
    with open(f"{topFile}") as f:
        a = f.readlines()
    assert len(a) == 3
    res_list = []
    first = None
    count = 1
    previousEnd = 0
    # print("g_all = [")
    topo = a[2].strip()
    cutoff = 30 if pdb != "1py6" else 15
    linkerSize = 10 if pdb != "1py6" else 5
    for i, res in enumerate(topo):
        o = "2" if res == "1" else "1"
        if res == "0":
            if len(res_list) > 0:
                # print(f"g{count} =", res_list)
                # print(res_list, ", ")
                count += 1
                last = res_list[-1]
                first = res_list[0] if first is None else first
                span = res_list[0] - previousEnd
                if span > cutoff:
                    # print(f"{pdb} Globular", previousEnd, res_list[0])
                    globular = list(range(previousEnd+linkerSize, res_list[0]-linkerSize))
                previousEnd = last
            res_list = []
        if res == "1":
            res_list.append(i)
    n = len(topo)
    # print(f"{pdb}: size {n}")
    span = n - previousEnd
    if span > cutoff:
        # print(f"{pdb} Globular", previousEnd, n)
        globular = list(range(previousEnd+linkerSize, n-linkerSize))

    membranePart = []
    for i in range(first-5, last+5):
        if i not in globular:
            membranePart.append(i)
    return membranePart, globular
# def get_two_part_from_prediction(probFile):
#     #   probFile = f"TM_pred/{pdb}_PureTM/{pdb}.prob"
#     with open(f"{probFile}") as f:
#         a = f.readlines()
#     res_list = []
#     first = None
#     count = 1
#     previousEnd = 0
#     # print("g_all = [")
#     out = "[\n"
#     for i, line in enumerate(a[3:]):
#         prob = float(line.strip().split()[3])
#         res = "0" if prob < 0.5 else "1"
#         o = "2" if res == "1" else "1"

#         if res == "0":
#             if len(res_list) > 0:
#                 # print(f"g{count} =", res_list)
#                 print(res_list, ", ")
#                 out += f"    {res_list},\n"
#                 count += 1
#                 last = res_list[-1]
#                 first = res_list[0] if first is None else first
#                 span = res_list[0] - previousEnd
#                 if span > 30:
#                     print("Globular", previousEnd, res_list[0])
#                     globular = list(range(previousEnd+10, res_list[0]-10))
#                 previousEnd = last
#             res_list = []
#         if res == "1":
#             res_list.append(i)
#     n = len(a[3:])
#     print(f"size {n}")
#     span = n - previousEnd
#     if span > 30:
#         print(f" Globular", previousEnd, n)
#         globular = list(range(previousEnd+10, n-10))
#     out += "]\n"

#     membranePart = []
#     for i in range(first-5, last+5):
#         if i not in globular:
#             membranePart.append(i)
#     return globular, membranePart

def compute_native_contacts(coords, MAX_OFFSET=4, DISTANCE_CUTOFF=9.5):
    native_coords = np.array(coords)
    a= native_coords[:,np.newaxis]
    dis = np.sqrt(np.sum((a - native_coords)**2, axis=2))

    n = len(dis)
    remove_band = np.eye(n)
    for i in range(1, MAX_OFFSET):
        remove_band += np.eye(n, k=i)
        remove_band += np.eye(n, k=-i)
    dis[remove_band==1] = np.max(dis)
    native_contacts = dis < DISTANCE_CUTOFF
    return native_contacts.astype("int")

def compute_contacts(coords, native_contacts, DISTANCE_CUTOFF=9.5):
    native_coords = np.array(coords)
    a= native_coords[:,np.newaxis]
    dis = np.sqrt(np.sum((a - native_coords)**2, axis=2))
    constacts = dis < DISTANCE_CUTOFF
    constacts = constacts*native_contacts  # remove non native contacts
    return np.sum(constacts, axis=1).astype("float")

def compute_localQ_init(MAX_OFFSET=4, DISTANCE_CUTOFF=9.5):
    from pathlib import Path
    home = str(Path.home())
    struct_id = '2xov'
    filename = os.path.join(home, "opt/pulling/2xov.pdb")
    p = PDBParser(PERMISSIVE=1)
    s = p.get_structure(struct_id, filename)
    chains = s[0].get_list()

    # import pdb file
    native_coords = []
    for chain in chains:
        dis = []
        all_res = []
        for res in chain:
            is_regular_res = res.has_id('CA') and res.has_id('O')
            res_id = res.get_id()[0]
            if (res.get_resname()=='GLY'):
                native_coords.append(res['CA'].get_coord())
            elif (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS') and is_regular_res:
                native_coords.append(res['CB'].get_coord())
            else:
                print('ERROR: irregular residue at %s!' % res)
                exit()
    native_contacts_table = compute_native_contacts(native_coords, MAX_OFFSET, DISTANCE_CUTOFF)

    return native_contacts_table

def compute_localQ(native_contacts_table, pre=".", ii=-1, MAX_OFFSET=4, DISTANCE_CUTOFF=9.5):
    native_contacts = np.sum(native_contacts_table, axis=1).astype("float")
    dump = read_lammps(os.path.join(pre, f"dump.lammpstrj.{ii}"), ca=False)
    localQ_list = []
    for atom in dump:
        contacts = compute_contacts(np.array(atom), native_contacts_table, DISTANCE_CUTOFF=DISTANCE_CUTOFF)
        c = np.divide(contacts, native_contacts, out=np.zeros_like(contacts), where=native_contacts!=0)
        localQ_list.append(c)
    data = pd.DataFrame(localQ_list)
    data.columns = ["Res" + str(i+1) for i in data.columns]
    data.to_csv(os.path.join(pre, f"localQ.{ii}.csv"), index=False)

def readPMF_basic(pre):
    # perturbation_table = {0:"original", 1:"p_mem",
    #                       2:"m_mem", 3:"p_lipid",
    #                       4:"m_lipid", 5:"p_go",
    #                       6:"m_go", 7:"p_rg", 8:"m_rg"}
    perturbation_table = {0:"original", 1:"m_go",
                          2:"p_go", 3:"m_lipid",
                          4:"p_lipid", 5:"m_mem",
                          6:"p_mem", 7:"m_rg", 8:"p_rg"}
    pmf_list = {
        "perturbation":list(perturbation_table.keys())
    }
    pmf_list_data = expand_grid(pmf_list)
    all_pmf_list = []
    for index, row in pmf_list_data.iterrows():
        perturbation = row["perturbation"]
        if perturbation == 0:
            location = pre + f"/pmf-*.dat"
            pmf_list = glob.glob(location)
            change = "none"
            upOrDown = "none"
        else:
            location = pre + f"/perturbation-{perturbation}-pmf-*.dat"
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
            data = pd.read_table(location, skiprows=2, sep='\s+', names=names).assign(upOrDown=upOrDown, change=change, temp=temp, perturbation=perturbation_table[perturbation])
            all_pmf_list.append(data)

    return pd.concat(all_pmf_list).dropna().reset_index()

def make_metadata_3(k=1000.0, temps_list=["450"], i=-1, biasLow=None, biasHigh=None):
    print("make metadata")
    cwd = os.getcwd()
    files = glob.glob(f"../data_{i}/*")
    kconstant = k
    with open("metadatafile", "w") as out:
        for oneFile in sorted(files):
            tmp = oneFile.split("/")[-1].replace('.dat', '')
            t = tmp.split("_")[1]
            bias = tmp.split("_")[3]
            if biasLow:
                if float(bias) < biasLow:
                    continue
            if biasHigh:
                if float(bias) > biasHigh:
                    continue
            # print(tmp)
            # if int(float(dis)) > 150:
            #     continue
            if t in temps_list:
                target = "../{} {} {} {}\n".format(oneFile, t, kconstant, bias)
                out.write(target)


def readPMF(pre, is2d=False, force_list=["0.0", "0.1", "0.2"]):
    # perturbation_table = {0:"original", 1:"p_mem",
    #                       2:"m_mem", 3:"p_lipid",
    #                       4:"m_lipid", 5:"p_go",
    #                       6:"m_go", 7:"p_rg", 8:"m_rg"}
    perturbation_table = {0:"original", 1:"m_go",
                          2:"p_go", 3:"m_lipid",
                          4:"p_lipid", 5:"m_mem",
                          6:"p_mem", 7:"m_rg", 8:"p_rg"}
    pmf_list = {
        "perturbation":list(perturbation_table.keys()),
        "force":force_list
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
        # print(pmf_list)
        name_list = ["f", "df", "e", "s"]
        if is2d:
            names = ["x", "y"] + name_list
        else:
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

def readPMF_2(pre, is2d=0, force_list=["0.0", "0.1", "0.2"]):
    if is2d:
        print("reading 2d pmfs")
    else:
        print("reading 1d dis, qw and z")
    if is2d == 1:
        mode_list = ["2d_qw_dis", "2d_z_dis", "2d_z_qw"]
    elif is2d == 2:
        mode_list = ["quick"]
    else:
        mode_list = ["1d_dis", "1d_qw", "1d_z"]
    all_data_list =[]
    for mode in mode_list:
        tmp = readPMF(mode, is2d, force_list).assign(mode=mode)
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

def compute_theta_for_each_helix(output="angles.csv", dumpName="../dump.lammpstrj.0"):
    print("This is for 2xov only")
    helices_list = [(94,114), (147,168), (171, 192), (200, 217), (226, 241), (250, 269)]
    atoms_all_frames = read_lammps(dumpName)
    # print(atoms[0])
    # print(len(atoms), len(atoms[0]))
    # helices_angles_all_frames = []
    with open(output, "w") as out:
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
                # tmp = line.replace("fix_backbone_coeff_er", backboneFile)
                tmp = line.replace("fix_backbone_coeff_hybrid", backboneFile)
                print(tmp, end='')
        cd("..")
        # do("run.py -m 0 -n 20 {}".format(protein))
        do("run.py -n 20 {}".format(protein))
        cd("..")
    cd("..")
    # do("")


def check_and_correct_fragment_memory_back(fragFile="fragsLAMW.mem"):
    with open("tmp.mem", "w") as out:
        with open(fragFile, "r") as f:
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
                    all_residues = []
                    for atom in one:
                        residue, resType, atomType, *_ = atom.split()
                        # print(residue, resType, atomType)
                        if atomType == "CA":
                            all_residues.append(int(residue))
                    all_residues = np.array(all_residues)
                    for test in range(int(i), int(i)+int(n)):
                        if (test == all_residues).sum() > 1:
                            # In rare case, one res id may have two different possible residues.
                            # on example, pdb 3vpg. chain A, res id 220.
                            # ATOM   1467  N   ARG A 220A      9.151 -20.984  46.737  1.00 31.30           N
                            # ATOM   1468  CA  ARG A 220A      9.120 -19.710  46.027  1.00 31.52           C
                            # ATOM   1469  C   ARG A 220A      9.768 -19.832  44.650  1.00 33.58           C
                            # ATOM   1470  O   ARG A 220A     10.552 -18.973  44.240  1.00 28.91           O
                            # ATOM   1471  CB  ARG A 220A      9.853 -18.641  46.847  1.00 31.58           C
                            # ATOM   1472  CG  ARG A 220A      9.181 -18.295  48.168  1.00 33.55           C
                            # ATOM   1473  CD  ARG A 220A      7.834 -17.651  47.916  1.00 34.70           C
                            # ATOM   1474  NE  ARG A 220A      7.959 -16.526  46.994  1.00 43.05           N
                            # ATOM   1475  CZ  ARG A 220A      6.931 -15.906  46.425  1.00 46.69           C
                            # ATOM   1476  NH1 ARG A 220A      5.691 -16.300  46.683  1.00 39.12           N
                            # ATOM   1477  NH2 ARG A 220A      7.144 -14.898  45.590  1.00 41.15           N
                            # ATOM   1478  N   ALA A 220B      9.429 -20.901  43.936  1.00 33.78           N
                            # ATOM   1479  CA  ALA A 220B      9.979 -21.153  42.608  1.00 32.13           C
                            # ATOM   1480  C   ALA A 220B      9.944 -19.933  41.692  1.00 30.71           C
                            # ATOM   1481  O   ALA A 220B      9.050 -19.088  41.787  1.00 28.56           O
                            # ATOM   1482  CB  ALA A 220B      9.234 -22.310  41.951  1.00 35.20           C
                            print("ATTENTION", gro, i, n, "duplicate:",test)
                            delete = True
                        if test not in all_residues:
                            print("ATTENTION", gro, i, n, "missing:",test)
                            delete = True
                if not delete:
                    out.write(line)
    os.system(f"mv {fragFile} fragsLAMW_back")
    os.system(f"mv tmp.mem {fragFile}")

# def check_and_correct_fragment_memory(fragFile="fragsLAMW.mem"):
#     with open("tmp.mem", "w") as out:
#         with open(fragFile, "r") as f:
#             for i in range(4):
#                 line = next(f)
#                 out.write(line)
#             for line in f:
#                 gro, _, i, n, _ = line.split()
#                 delete = False
#                 # print(gro, i, n)
#                 # name = gro.split("/")[-1]
#                 with open(gro, "r") as one:
#                     next(one)
#                     next(one)
#                     all_residues = set()
#                     for atom in one:
#                         residue, *_ = atom.split()
#                         # print(residue)
#                         all_residues.add(int(residue))
#                     for test in range(int(i), int(i)+int(n)):
#                         if test not in all_residues:
#                             print("ATTENTION", gro, i, n, "missing:",test)
#                             delete = True
#                 if not delete:
#                     out.write(line)
#     os.system(f"mv {fragFile} fragsLAMW_back")
#     os.system(f"mv tmp.mem {fragFile}")

def read_complete_temper_2(n=4, location=".", rerun=-1, qnqc=False, average_z=False, localQ=False, disReal=False, dis_h56=False, goEnergy=False, goEnergy3H=False, goEnergy4H=False):
    all_data_list = []
    for i in range(n):
        file = "lipid.{}.dat".format(i)
        lipid = pd.read_csv(location+file)
        lipid.columns = lipid.columns.str.strip()
        remove_columns = ['Steps']
        lipid = lipid.drop(remove_columns, axis=1)

        file = "rgs.{}.dat".format(i)
        rgs = pd.read_csv(location+file)
        rgs.columns = rgs.columns.str.strip()
        remove_columns = ['Steps']
        rgs = rgs.drop(remove_columns, axis=1)

        file = "energy.{}.dat".format(i)
        energy = pd.read_csv(location+file)
        energy.columns = energy.columns.str.strip()
        energy = energy[["AMH-Go", "Membrane", "Rg"]]
        file = "addforce.{}.dat".format(i)
        dis = pd.read_csv(location+file)
        dis.columns = dis.columns.str.strip()
        remove_columns = ['Steps', 'AddedForce', 'Dis12', 'Dis34', 'Dis56']
        dis.drop(remove_columns, axis=1,inplace=True)


        file = "wham.{}.dat".format(i)
        wham = pd.read_csv(location+file).assign(Run=i)
        wham.columns = wham.columns.str.strip()
        remove_columns = ['Rg', 'Tc']
        wham = wham.drop(remove_columns, axis=1)
        if qnqc:
            qc = pd.read_table(location+f"qc_{i}", names=["qc"])[1:].reset_index(drop=True)
            qn = pd.read_table(location+f"qn_{i}", names=["qn"])[1:].reset_index(drop=True)
            qc2 = pd.read_table(location+f"qc2_{i}", names=["qc2"])[1:].reset_index(drop=True)
            wham = pd.concat([wham, qn, qc, qc2],axis=1)
        # if average_z:
        #     z = pd.read_table(location+f"z_{i}.dat", names=["AverageZ"])[1:].reset_index(drop=True)
        #     wham = pd.concat([wham, z],axis=1)
        if disReal:
            tmp = pd.read_csv(location+f"distance_{i}.dat")[1:].reset_index(drop=True).drop('Steps', axis=1)
            # print(tmp)
            tmp.columns = tmp.columns.str.strip()
            wham = pd.concat([wham, tmp],axis=1)
        if dis_h56:
            tmp = pd.read_csv(location+f"distance_h56_{i}.dat")[1:].reset_index(drop=True).drop('Steps', axis=1)
            tmp1 = pd.read_csv(location+f"distance_h12_{i}.dat")[1:].reset_index(drop=True).drop('Steps', axis=1)
            tmp2 = pd.read_csv(location+f"distance_h34_{i}.dat")[1:].reset_index(drop=True).drop('Steps', axis=1)
            # print(tmp)
            tmp.columns = tmp.columns.str.strip()
            tmp1.columns = tmp1.columns.str.strip()
            tmp2.columns = tmp2.columns.str.strip()
            wham = pd.concat([wham, tmp, tmp1, tmp2],axis=1)
        if average_z:
            z = pd.read_csv(location+f"z_complete_{i}.dat")[1:].reset_index(drop=True)
            z.columns = z.columns.str.strip()
            wham = pd.concat([wham, z],axis=1)
        if localQ:
            all_localQ = pd.read_csv(location+f"localQ.{i}.csv")[1:].reset_index(drop=True)
            wham = pd.concat([wham, all_localQ], axis=1)
        if goEnergy:
            tmp = pd.read_csv(location+f"Go_{i}/goEnergy.dat")[1:].reset_index(drop=True).drop('Steps', axis=1)
            # print(tmp)
            tmp.columns = tmp.columns.str.strip()
            wham = pd.concat([wham, tmp],axis=1)
        if goEnergy3H:
            nEnergy = pd.read_csv(location+f"Go_3helix_{i}/goEnergy.dat")[1:].reset_index(drop=True).drop('Steps', axis=1)
            # print(tmp)
            nEnergy.columns = nEnergy.columns.str.strip()
            wham = pd.concat([wham, nEnergy],axis=1)
        if goEnergy4H:
            nEnergy = pd.read_csv(location+f"Go_4helix_{i}/goEnergy.dat")[1:].reset_index(drop=True).drop('Steps', axis=1)
            # print(tmp)
            nEnergy.columns = nEnergy.columns.str.strip()
            wham = pd.concat([wham, nEnergy],axis=1)
        data = pd.concat([wham, dis, energy, rgs, lipid], axis=1)

        # lipid = lipid[["Steps","Lipid","Run"]]
        all_data_list.append(data)
    data = pd.concat(all_data_list)
    file = f"../log{rerun}/log.lammps"
    temper = pd.read_table(location+file, skiprows=2, sep=' ')
    temper = temper.melt(id_vars=['Step'], value_vars=['T' + str(i) for i in range(n)], value_name="Temp", var_name="Run")
    temper["Run"] = temper["Run"].str[1:].astype(int)
    temper["Temp"] = "T" + temper["Temp"].astype(str)
#     print(temper)
#     print(wham)
    t2 = temper.merge(data, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]).sort_values('Step').drop('Steps', axis=1)
#     print(t2)
    t3 = t2.assign(TotalE=t2.Energy + t2.Lipid)
    return t3.sort_values(["Step", "Run"]).reset_index(drop=True)

def process_complete_temper_data_3(pre, data_folder, folder_list, rerun=-1, end=-1, n=12, bias="dis", qnqc=False, average_z=False, disReal=False, dis_h56=False, localQ=False, goEnergy=False, goEnergy3H=False, goEnergy4H=False, label=""):
    print("process temp data")
    dateAndTime = datetime.today().strftime('%d_%h_%H%M%S')
    for folder in folder_list:
        simulation_list = glob.glob(pre+folder+f"/simulation/{bias}_*")
        print(pre+folder+f"/simulation/{bias}_*")
        os.system("mkdir -p " + pre+folder+"/data")
        # this one only consider rerun >=0, for the case rerun=-1, move log.lammps to log0
        for i in range(rerun, end, -1):
            all_data_list = []
            for one_simulation in simulation_list:
                bias_num = one_simulation.split("_")[-1]
                print(bias_num, "!")

                location = one_simulation + f"/{i}/"
                print(location)
                data = read_complete_temper_2(location=location, n=n, rerun=i, qnqc=qnqc, average_z=average_z, localQ=localQ, disReal=disReal, dis_h56=dis_h56, goEnergy=goEnergy, goEnergy3H=goEnergy3H, goEnergy4H=goEnergy4H)
                print(data.shape)
                # remove_columns = ['Step', "Run"]
                # data = data.drop(remove_columns, axis=1)
                all_data_list.append(data.assign(BiasTo=bias_num))

            data = pd.concat(all_data_list).reset_index(drop=True)
            # if localQ:
            #     print("hi")
            # else:
            #     data.to_csv(os.path.join(pre, folder, f"data/rerun_{i}.csv"))
            # complete_data_list.append(data)
            #         temps = list(dic.keys())
            # complete_data = pd.concat(complete_data_list)
            name = f"rerun_{2*i}_{dateAndTime}.feather"
            data = data.reset_index(drop=True)
            data.query(f'Step > {2*i}e7 & Step <= {2*i+1}e7').reset_index(drop=True).to_feather(pre+folder+"/" + name)
            os.system("cp "+pre+folder+"/" + name + " "+data_folder+label+name)
            name = f"rerun_{2*i+1}_{dateAndTime}.feather"
            data = data.reset_index(drop=True)
            data.query(f'Step > {2*i+1}e7 & Step <= {2*i+2}e7').reset_index(drop=True).to_feather(pre+folder+"/" + name)
            os.system("cp "+pre+folder+"/" + name + " "+data_folder+label+name)



def move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=1, sub_mode_name="", kmem=0.2, klipid=0.1, kgo=0.1, krg=0.2, sample_range_mode=0, biasName="dis", qnqc=False, average_z=0, chosen_mode=0):
    print("move data")
    # dic = {"T_defined":300, "T0":350, "T1":400, "T2":450, "T3":500, "T4":550, "T5":600, "T6":650, "T7":700, "T8":750, "T9":800, "T10":900, "T11":1000}
    if temp_dict_mode == 1:
        dic = {"T0":280, "T1":300, "T2":325, "T3":350, "T4":375, "T5":400, "T6":450, "T7":500, "T8":550, "T9":600, "T10":650, "T11":700}
    if temp_dict_mode == 2:
        dic = {"T0":280, "T1":290, "T2":300, "T3":315, "T4":335, "T5":355, "T6":380, "T7":410, "T8":440, "T9":470, "T10":500, "T11":530}
    if temp_dict_mode == 3:
        dic = {"T0":280, "T1":290, "T2":300, "T3":310, "T4":320, "T5":335, "T6":350, "T7":365, "T8":380, "T9":410, "T10":440, "T11":470}
    if temp_dict_mode == 4:
        dic = {"T0":300, "T1":335, "T2":373, "T3":417, "T4":465, "T5":519, "T6":579, "T7":645, "T8":720, "T9":803, "T10":896, "T11":1000}
    # read in complete.feather
    data_list = []
    for folder in folder_list:
        tmp = pd.read_feather(data_folder + folder +".feather")
        data_list.append(tmp)
    data = pd.concat(data_list)
    os.system("mkdir -p "+freeEnergy_folder+"/"+sub_mode_name+f"/data_{sample_range_mode}")
    for bias, oneBias in data.groupby("BiasTo"):
        for tempSymbol, oneTempAndBias in oneBias.groupby("Temp"):
            temp = dic[tempSymbol]
            if float(temp) > 800:
                continue
            print(f"t_{temp}_{biasName}_{bias}.dat")
            if sample_range_mode == 0:
                queryCmd = 'Step > 0 & Step <= 1e7'
            if sample_range_mode == 1:
                queryCmd = 'Step > 1e7 & Step <= 2e7'
            elif sample_range_mode == 2:
                queryCmd ='Step > 2e7 & Step <= 3e7'
            elif sample_range_mode == 3:
                queryCmd ='Step > 3e7 & Step <= 4e7'
            elif sample_range_mode == 4:
                queryCmd ='Step > 4e7 & Step <= 5e7'
            elif sample_range_mode == 5:
                queryCmd ='Step > 5e7 & Step <= 6e7'
            elif sample_range_mode == 6:
                queryCmd ='Step > 6e7 & Step <= 7e7'
            elif sample_range_mode == 7:
                queryCmd ='Step > 7e7 & Step <= 8e7'
            elif sample_range_mode == -1:
                queryCmd ='Step > 4e7 & Step <= 6e7'
            if sample_range_mode == -2:
                tmp = oneTempAndBias.reset_index(drop=True)
            else:
                tmp = oneTempAndBias.query(queryCmd).reset_index()
            if average_z < 5:
                chosen_list = ["TotalE", "Qw", "Distance"]
            elif average_z == 5:
                chosen_list = ["TotalE", "Qw", "DisReal"]
                chosen_list += ["z_h6"]
            if average_z == 1:
                chosen_list += ["abs_z_average"]
            if average_z == 2 or average_z == 3:
                chosen_list += ["z_h6"]
            if average_z == 3:
                chosen_list += ["DisReal"]
            if average_z == 4:
                tmp["z_h5_and_h6"] = tmp["z_h5"] + tmp["z_h6"]
                chosen_list += ["z_h5_and_h6"]
                chosen_list += ["DisReal"]
            if average_z == 6:
                chosen_list = ["TotalE", "Qw", "DisReal"]
                tmp["z_h5_and_h6"] = tmp["z_h5"] + tmp["z_h6"]
                chosen_list += ["z_h5_and_h6"]
                chosen_list += ["z_h5"]
                chosen_list += ["z_h6"]
                chosen_list += ["Dis_h56"]
            if average_z == 7:
                chosen_list = ["TotalE", "Qw", "DisReal"]
                tmp["z_h56"] = tmp["z_h5"] + tmp["z_h6"]
                tmp["z_h14"] = tmp["z_h1"] + tmp["z_h2"] + tmp["z_h3"] + tmp["z_h4"]
                chosen_list += ["z_h14"]
                chosen_list += ["z_h56"]
                chosen_list += ["z_h5"]
                chosen_list += ["z_h6"]
                chosen_list += ["Dis_h12"]
                chosen_list += ["Dis_h34"]
                chosen_list += ["Dis_h56"]
            if chosen_mode == 0:
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_perturb_mem_p=tmp.TotalE + kmem*tmp.Membrane,
                                        TotalE_perturb_mem_m=tmp.TotalE - kmem*tmp.Membrane,
                                        TotalE_perturb_lipid_p=tmp.TotalE + klipid*tmp.Lipid,
                                        TotalE_perturb_lipid_m=tmp.TotalE - klipid*tmp.Lipid,
                                        TotalE_perturb_go_p=tmp.TotalE + kgo*tmp["AMH-Go"],
                                        TotalE_perturb_go_m=tmp.TotalE - kgo*tmp["AMH-Go"],
                                        TotalE_perturb_rg_p=tmp.TotalE + krg*tmp.Rg,
                                        TotalE_perturb_rg_m=tmp.TotalE - krg*tmp.Rg)
            if chosen_mode == 1:
                chosen_list += ["Res" + str(i+1) for i in range(181)]
                chosen = tmp[chosen_list]
            if chosen_mode == 2:
                chosen_list += ["Res" + str(i+1) for i in range(181)]
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_perturb_go_m=tmp.TotalE - kgo*tmp["AMH-Go"],
                                        TotalE_perturb_go_p=tmp.TotalE + kgo*tmp["AMH-Go"],
                                        TotalE_perturb_lipid_m=tmp.TotalE - klipid*tmp.Lipid,
                                        TotalE_perturb_lipid_p=tmp.TotalE + klipid*tmp.Lipid,
                                        TotalE_perturb_mem_m=tmp.TotalE - kmem*tmp.Membrane,
                                        TotalE_perturb_mem_p=tmp.TotalE + kmem*tmp.Membrane,
                                        TotalE_perturb_rg_m=tmp.TotalE - krg*tmp.Rg,
                                        TotalE_perturb_rg_p=tmp.TotalE + krg*tmp.Rg)
    #         print(tmp.count())
            if chosen_mode == 3:
                chosen_list += ["AMH-Go", "Lipid", "Membrane", "Rg"]
                chosen = tmp[chosen_list]
            if chosen_mode == 4:
                chosen_list += ["Dis_h56"]
                chosen = tmp[chosen_list]
            if chosen_mode == 5:
                chosen_list += ["Dis_h56"]
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_perturb_go_m=tmp.TotalE/10,
                                        TotalE_perturb_go_p=0,
                                        Go=tmp["AMH-Go"])
            if chosen_mode == 6:
                chosen_list += ["Dis_h56"]
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_1=tmp.TotalE + 0.1*tmp.AMH,
                                        TotalE_2=tmp.TotalE + 0.2*tmp.AMH,
                                        TotalE_3=tmp.TotalE + 0.5*tmp.AMH,
                                        TotalE_4=tmp.TotalE + tmp.AMH,
                                        TotalE_5=tmp.AMH)
            if chosen_mode == 7:
                chosen_list += ["Dis_h56"]
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_1=tmp.TotalE + 0.1*tmp.AMH_3H,
                                        TotalE_2=tmp.TotalE + 0.2*tmp.AMH_3H,
                                        TotalE_3=tmp.TotalE + 0.5*tmp.AMH_3H,
                                        TotalE_4=tmp.TotalE + tmp.AMH_3H,
                                        TotalE_5=tmp.TotalE + 0.1*tmp.AMH,
                                        TotalE_6=tmp.TotalE + 0.2*tmp.AMH)
            if chosen_mode == 8:
                # chosen_list += ["Dis_h56"]
                chosen_list += ["z_average"]
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_1=tmp.TotalE + 0.1*tmp.AMH_4H,
                                        TotalE_2=tmp.TotalE + 0.2*tmp.AMH_4H,
                                        TotalE_3=tmp.TotalE + 0.5*tmp.AMH_4H,
                                        TotalE_4=tmp.TotalE + 0.1*tmp.AMH_3H,
                                        TotalE_5=tmp.TotalE + 0.2*tmp.AMH_3H,
                                        TotalE_6=tmp.TotalE + 0.5*tmp.AMH_3H)
            if chosen_mode == 9:
                # chosen_list += ["Dis_h56"]
                chosen_list += ["z_average"]
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_1=tmp.TotalE + 0.1*tmp.AMH_4H,
                                        TotalE_2=tmp.TotalE + 0.2*tmp.AMH_4H,
                                        TotalE_3=tmp.TotalE + 0.5*tmp.AMH_4H)
                chosen = chosen.assign(TotalE_perturb_1go_m=chosen.TotalE_2 - kgo*tmp["AMH-Go"],
                                        TotalE_perturb_1go_p=chosen.TotalE_2 + kgo*tmp["AMH-Go"],
                                        TotalE_perturb_2lipid_m=chosen.TotalE_2 - tmp.Lipid,
                                        TotalE_perturb_2lipid_p=chosen.TotalE_2 + tmp.Lipid,
                                        TotalE_perturb_3mem_m=chosen.TotalE_2 - tmp.Membrane,
                                        TotalE_perturb_3mem_p=chosen.TotalE_2 + tmp.Membrane,
                                        TotalE_perturb_4rg_m=chosen.TotalE_2 - tmp.Rg,
                                        TotalE_perturb_4rg_p=chosen.TotalE_2 + tmp.Rg,
                                        TotalE_perturb_5go=tmp["AMH-Go"],
                                        TotalE_perturb_5lipid=tmp.Lipid,
                                        TotalE_perturb_5mem=tmp.Membrane,
                                        TotalE_perturb_5rg=tmp.Rg)
            if chosen_mode == 10:
                # chosen_list += ["Dis_h56"]
                chosen_list += ["z_average"]
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_1=tmp.TotalE + 0.1*tmp.AMH_4H,
                                        TotalE_2=tmp.TotalE + 0.2*tmp.AMH_4H,
                                        TotalE_3=tmp.TotalE + 0.5*tmp.AMH_4H)
                chosen = chosen.assign(TotalE_perturb_1lipid_m1=chosen.TotalE_2 - 0.1*tmp.Lipid,
                                        TotalE_perturb_1lipid_p1=chosen.TotalE_2 + 0.1*tmp.Lipid,
                                        TotalE_perturb_2lipid_m2=chosen.TotalE_2 - 0.2*tmp.Lipid,
                                        TotalE_perturb_2lipid_p2=chosen.TotalE_2 + 0.2*tmp.Lipid,
                                        TotalE_perturb_3lipid_m3=chosen.TotalE_2 - 0.3*tmp.Lipid,
                                        TotalE_perturb_3lipid_p3=chosen.TotalE_2 + 0.3*tmp.Lipid,
                                        TotalE_perturb_4lipid_m4=chosen.TotalE_2 - 0.5*tmp.Lipid,
                                        TotalE_perturb_4lipid_p4=chosen.TotalE_2 + 0.5*tmp.Lipid,
                                        TotalE_perturb_5go=tmp["AMH-Go"],
                                        TotalE_perturb_5lipid=tmp.Lipid,
                                        TotalE_perturb_5mem=tmp.Membrane,
                                        TotalE_perturb_5rg=tmp.Rg)
            if chosen_mode == 11:
                # chosen_list += ["Dis_h56"]
                chosen_list += ["z_average"]
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_1=tmp.TotalE + 1.1*0.1*tmp.AMH_4H + 0.1*tmp["AMH-Go"],
                                        TotalE_2=tmp.TotalE + 1.1*0.2*tmp.AMH_4H + 0.1*tmp["AMH-Go"],
                                        TotalE_3=tmp.TotalE + 1.1*0.5*tmp.AMH_4H + 0.1*tmp["AMH-Go"])
                chosen = chosen.assign(TotalE_perturb_1lipid_m1=chosen.TotalE_2 - 0.1*tmp.Lipid,
                                        TotalE_perturb_1lipid_p1=chosen.TotalE_2 + 0.1*tmp.Lipid,
                                        TotalE_perturb_2lipid_m2=chosen.TotalE_2 - 0.2*tmp.Lipid,
                                        TotalE_perturb_2lipid_p2=chosen.TotalE_2 + 0.2*tmp.Lipid,
                                        TotalE_perturb_3lipid_m3=chosen.TotalE_2 - 0.1*tmp.Membrane,
                                        TotalE_perturb_3lipid_p3=chosen.TotalE_2 + 0.1*tmp.Membrane,
                                        TotalE_perturb_4lipid_m4=chosen.TotalE_2 - 0.2*tmp.Membrane,
                                        TotalE_perturb_4lipid_p4=chosen.TotalE_2 + 0.2*tmp.Membrane,
                                        TotalE_perturb_5go=tmp["AMH-Go"],
                                        TotalE_perturb_5lipid=tmp.Lipid,
                                        TotalE_perturb_5mem=tmp.Membrane,
                                        TotalE_perturb_5rg=tmp.Rg)
            if chosen_mode == 12:
                chosen = tmp[chosen_list]
                # chosen["z_h56"] = (chosen["z_h5"] + chosen["z_h6"])/2
                chosen = chosen.assign(TotalE_2=tmp.TotalE + 0.2*tmp.AMH_4H,
                                        z_h56=(tmp.z_h5 + tmp.z_h6)/2)
            if chosen_mode == 13:
                chosen_list += ["z_average"]
                chosen = tmp[chosen_list]
                # chosen["z_h56"] = (chosen["z_h5"] + chosen["z_h6"])/2
                force = 0.1
                chosen = chosen.assign(TotalE_2=tmp.TotalE + 0.2*tmp.AMH_4H - (tmp.DisReal - 25.1)*force,
                                        TotalE_3=tmp.TotalE - (tmp.DisReal - 25.1)*force,
                                        TotalE_4=tmp.TotalE + 0.2*tmp.AMH_4H,
                                        TotalE_5=tmp.TotalE + 0.2*tmp.AMH_4H - (tmp.DisReal)*force)
            chosen.to_csv(freeEnergy_folder+"/"+sub_mode_name+f"/data_{sample_range_mode}/t_{temp}_{biasName}_{bias}.dat", sep=' ', index=False, header=False)

    # perturbation_table = {0:"original", 1:"m_go",
    #                       2:"p_go", 3:"m_lipid",
    #                       4:"p_lipid", 5:"m_mem",
    #                       6:"p_mem", 7:"m_rg", 8:"p_rg"}
def compute_average_z(dumpFile, outFile):
    # input dump, output z.dat
    z_list = []
    with open(outFile, "w") as f:
        a = read_lammps(dumpFile)
        for atoms in a:
            b = np.array(atoms)
            z = b.mean(axis=0)[2]
            z_list.append(z)
            f.write(str(z)+"\n")

def compute_average_z_2(dumpFile, outFile):
    # input dump, output z.dat

    helices_list = [(94,114), (147,168), (171, 192), (200, 217), (226, 241), (250, 269)]
    with open(outFile, "w") as f:
        a = read_lammps(dumpFile)
        f.write("z_average, abs_z_average, z_h1, z_h2, z_h3, z_h4, z_h5, z_h6\n")
        for atoms in a:
            b = np.array(atoms)
            z = b.mean(axis=0)[2]
            f.write(str(z)+ ", ")
            z = np.abs(b).mean(axis=0)[2]
            f.write(str(z)+ ", ")
            for count, (i,j) in enumerate(helices_list):
                i = i - 91
                j = j - 91
                z = np.mean(b[i:j], axis=0)[2]
                if count == 5:
                    f.write(str(z))
                else:
                    f.write(str(z)+ ", ")
            f.write("\n")

def read_simulation_2(location=".", i=-1, qnqc=False, average_z=False, localQ=False, disReal=False, **kwargs):
    file = "lipid.dat"
    lipid = pd.read_csv(location+file)
    lipid.columns = lipid.columns.str.strip()
    remove_columns = ['Steps']
    lipid = lipid.drop(remove_columns, axis=1)

    file = "rgs.dat"
    rgs = pd.read_csv(location+file)
    rgs.columns = rgs.columns.str.strip()
    remove_columns = ['Steps']
    rgs = rgs.drop(remove_columns, axis=1)

    file = "energy.dat"
    energy = pd.read_csv(location+file)
    energy.columns = energy.columns.str.strip()
    energy = energy[["AMH-Go", "Membrane", "Rg"]]
    file = "addforce.dat"
    dis = pd.read_csv(location+file)
    dis.columns = dis.columns.str.strip()
    remove_columns = ['Steps', 'AddedForce', 'Dis12', 'Dis34', 'Dis56']
    dis.drop(remove_columns, axis=1,inplace=True)


    file = "wham.dat"
    wham = pd.read_csv(location+file).assign(Run=i)
    wham.columns = wham.columns.str.strip()
    remove_columns = ['Rg', 'Tc']
    wham = wham.drop(remove_columns, axis=1)
    if qnqc:
        qc = pd.read_table(location+f"qc", names=["qc"])[1:].reset_index(drop=True)
        qn = pd.read_table(location+f"qn", names=["qn"])[1:].reset_index(drop=True)
        qc2 = pd.read_table(location+f"qc2", names=["qc2"])[1:].reset_index(drop=True)
        wham = pd.concat([wham, qn, qc, qc2],axis=1)
    # if average_z:
    #     z = pd.read_table(location+f"z_{i}.dat", names=["AverageZ"])[1:].reset_index(drop=True)
    #     wham = pd.concat([wham, z],axis=1)
    if disReal:
        tmp = pd.read_csv(location+f"distance.dat")[1:].reset_index(drop=True).drop('Steps', axis=1)
        # print(tmp)
        tmp.columns = tmp.columns.str.strip()
        wham = pd.concat([wham, tmp],axis=1)
    if average_z:
        z = pd.read_csv(location+f"z_complete.dat")[1:].reset_index(drop=True)
        z.columns = z.columns.str.strip()
        wham = pd.concat([wham, z],axis=1)
    if localQ:
        all_localQ = pd.read_csv(location+f"localQ.csv")[1:].reset_index(drop=True)
        wham = pd.concat([wham, all_localQ], axis=1)

    data = pd.concat([wham, dis, energy, rgs, lipid], axis=1)
    t3 = data.assign(TotalE=data.Energy + data.Lipid)
    return t3.reset_index(drop=True)

def read_folder(location, match="", **kwargs):
    runFolders = os.listdir(location+"/simulation")
    if match == "qbias":
        runFolders = [f for f in runFolders if re.match(r'qbias_[0-9]+', f)]
    else:
        runFolders = [f for f in runFolders if re.match(r'[0-9]+', f)]
    print(runFolders)
    data_list = []
    for run in runFolders:
        tmp = read_simulation_2(location+"/simulation/"+run+"/0/", **kwargs).assign(Run=run)
        data_list.append(tmp)
    return pd.concat(data_list).reset_index(drop=True)

def read_variable_folder(location, match="*_", **kwargs):
    variables = glob.glob(os.path.join(location, match))
    print(variables)
    data_list = []
    for variableFolder in variables:
        tmp = variableFolder.split("/")[-1]
        data_list.append(read_folder(variableFolder, **kwargs).assign(Folder=tmp))
    data = pd.concat(data_list)
    name = f"{datetime.today().strftime('%d_%h_%H%M%S')}.feather"
    data.reset_index(drop=True).to_feather(name)

def getFragPdb(pdbId, i, outFile=None):
    pdb = pdbId + ".pdb"
    if outFile is None:
        outFile = f"{i}_{pdb}"
#     pdb = "1igqB00.pdb"
#     pdbId = pdb.split('.')[0]
    pre = "/Users/weilu/Research/optimization/fragment/"
    database = "/Users/weilu/Research/optimization/fragment/database/dompdb/"
    parser = bio.PDBParser(QUIET=True)
    structure = parser.get_structure("x", os.path.join(database, pdb))
    for model in structure:
        for chain in model:
            all_residues = list(chain)
            io = bio.PDBIO()
            c = bio.Chain.Chain("A")
            c.child_list = all_residues[i:i+9]
#             for ii, res in enumerate(c):
#                 res.id = (' ', ii+1, ' ')
            io.set_structure(c)
            io.save(f'{pre}{outFile}')

# membrane_gamma = np.loadtxt("/Users/weilu/opt/parameters/membrane/gamma.dat")
# membrane_gamma_formated = simulation_to_iteration_gamma(membrane_gamma)
# water_gamma = np.loadtxt("/Users/weilu/opt/parameters/globular_parameters/gamma.dat")
# water_gamma_formated = simulation_to_iteration_gamma(water_gamma)
# ratio = np.std(water_gamma_formated)/np.std(membrane_gamma_formated)
# membrane_gamma_formated_rescaled = membrane_gamma_formated * ratio
def simulation_to_iteration_gamma(membrane_gamma):
    direct = membrane_gamma[:210][:,0]
    protein = membrane_gamma[210:][:,0]
    water = membrane_gamma[210:][:,1]
    return np.concatenate([direct, protein, water])

def gamma_format_convertion_iteration_to_simulation(iteration_gamma, gamma_for_simulation, burial_gamma_for_simulation=None):
    from Bio.PDB.Polypeptide import one_to_three, three_to_one
    res_type_map = {
        'A': 0,
        'C': 4,
        'D': 3,
        'E': 6,
        'F': 13,
        'G': 7,
        'H': 8,
        'I': 9,
        'K': 11,
        'L': 10,
        'M': 12,
        'N': 2,
        'P': 14,
        'Q': 5,
        'R': 1,
        'S': 15,
        'T': 16,
        'V': 19,
        'W': 17,
        'Y': 18
    }
    # res_type_map = gamma_se_map_1_letter = {   'A': 0,  'R': 1,  'N': 2,  'D': 3,  'C': 4,
    #                             'Q': 5,  'E': 6,  'G': 7,  'H': 8,  'I': 9,
    #                             'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14,
    #                             'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19}
    res_type_map_letters = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
                            'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    inverse_res_type_map = dict(list(zip(list(range(20)), res_type_map_letters)))

    # gamma_location = "/Users/weilu/Research/server_backup/jan_2019/optimization/gammas_dec30/cath-dataset-nonredundant-S20Clean_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0_gamma"
    # gamma_for_simulation = "/Users/weilu/Research/server_backup/jan_2019/optimization/iteration_gamma.dat"
    gamma = iteration_gamma
    gamma = -gamma  # caused by tradition.
    # convert gamma to gamma used by simulation
    with open(gamma_for_simulation, "w") as out:
        c = 0
        for i in range(20):
            for j in range(i, 20):
                out.write(f"{gamma[c]:<.5f} {gamma[c]:10.5f}\n")
                c += 1
        out.write("\n")
        for i in range(20):
            for j in range(i, 20):
                # protein, water
                out.write(f"{gamma[c]:<.5f} {gamma[c+210]:10.5f}\n")
                c += 1
    if burial_gamma_for_simulation:
        rhoGamma = pd.DataFrame(gamma[630:].reshape(3,20).T, columns=["rho1", "rho2", "rho3"]).reset_index()
        rhoGamma["oneLetter"] = rhoGamma["index"].apply(lambda x: inverse_res_type_map[x])
        rhoGamma["Residue"] = rhoGamma["index"].apply(lambda x: one_to_three(inverse_res_type_map[x]))
        rhoGamma = rhoGamma[["Residue", "rho1", "rho2", "rho3", "index", "oneLetter"]]
        g = rhoGamma[["rho1", "rho2", "rho3"]].values
        np.savetxt(burial_gamma_for_simulation, g, fmt='%7.4f')


def mix_gammas_3(pre, Gamma, preGamma, alpha=None, scale=True, iterGammaName=None, iteration="5"):
    percent = int(alpha*100)
    if scale:
        scale = np.std(preGamma)/np.std(Gamma)
        iter_gamma = ((1- alpha)*preGamma + alpha*(scale*Gamma)).astype(float)
        iter_gamma *= np.std(preGamma)/np.std(iter_gamma)
    else:
        iter_gamma = ((1- alpha)*preGamma + alpha*(Gamma)).astype(float)

    # pre = "/Users/weilu/Research/server/march_2019/optimization_weighted_by_q_iter1/"
    gamma_for_simulation = pre + f"iteration_{iteration}_gamma_{percent}.dat"
    burial_gamma_for_simulation = pre + f"iteration_{iteration}_burial_gamma_{percent}.dat"
    gamma_format_convertion_iteration_to_simulation(iter_gamma, gamma_for_simulation, burial_gamma_for_simulation=burial_gamma_for_simulation)
    if iterGammaName is not None:
        np.savetxt(pre+iterGammaName, iter_gamma)

# def relocate(fileLocation="frags.mem", toLocation="../fraglib"):
#     # location = "/Users/weilu/Research/server/april_2019/iterative_optimization_new_set_with_frag/all_simulations/1fc2/1fc2"
#     # fileLocation = location + "/frags.mem"
#     # toLocation
#     print(os.getcwd())
#     os.system(f"mkdir -p {toLocation}")
#     a = pd.read_csv(fileLocation, skiprows=4, sep=" ", names=["location", "i", "j", "sep", "w"])
#     b = a["location"].unique()
#     for l in b:
#         out = os.system(f"cp {l} {toLocation}/")
#         if out != 0:
#             print(f"!!Problem!!, {l}")

# def replace(TARGET, FROM, TO):
#     os.system("sed -i.bak 's@{}@{}@g' {}".format(FROM,TO,TARGET))
def replace(TARGET, FROM, TO):
    os.system("sed -i.bak 's@{}@{}@g' {}".format(FROM,TO,TARGET))

def generate_SEQRES(fastaFile):
    # fastaFile = "/Users/weilu/Research/server/april_2019/complete_2xov/P09391.fasta"
    with open(fastaFile) as input_data:
        data = ""
        for line in input_data:
            if(line[0] == ">"):
                print(line)
            elif(line == "\n"):
                pass
            else:
                data += line.strip("\n")
    i = 0
    c = 1
    template = ""
    length = len(data)
    while i < len(data):
        seq = data[i:i+13]
        a = [one_to_three(r) for r in seq]
        ss = " ".join(a)
        template += f"SEQRES {c: >3} A {length: >4}  {ss: <51}          \n"
        i += 13
        c += 1
    return template

def getSeqFromFasta(fastaFile):
    # "/Users/weilu/Research/server/aug_2019/hybrid_protein_simulation/setup/1pv6/1pv6.fasta"
    with open(fastaFile) as f:
        a = f.readlines()
    seq = ""
    for line in a:
        if line[0] == ">":
            continue
        seq += line.strip()
    return seq

def cleanPdb1(pdb_list, chain=None, source=None, toFolder="cleaned_pdbs", formatName=False, verbose=False, removeTwoEndsMissingResidues=True):
    from pdbfixer import PDBFixer
    from simtk.openmm.app import PDBFile
    os.system(f"mkdir -p {toFolder}")
    for pdb_id in pdb_list:
        # print(chain)
        print(pdb_id)
        # pdb = f"{pdb_id.lower()[:4]}"
        # pdbFile = pdb+".pdb"
        if formatName:
            pdb = f"{pdb_id.lower()[:4]}"
        else:
            pdb = pdb_id
        pdbFile = pdb + ".pdb"
        if source is None:
            fromFile = os.path.join("original_pdbs", pdbFile)
        elif source[-4:] == ".pdb":
            fromFile = source
        else:
            fromFile = os.path.join(source, pdbFile)

        # clean pdb
        # os.system("pwd")
        # print(f"--{fromFile}--")
        fixer = PDBFixer(filename=fromFile)
        # remove unwanted chains
        chains = list(fixer.topology.chains())
        print(chains)
        if chain is None:  # None mean deafult is chain A unless specified.
            if len(pdb_id) >= 5:
                Chosen_chain = pdb_id[4]
                # Chosen_chain = pdb_id[4].upper()
            else:
                assert(len(pdb_id) == 4)
                Chosen_chain = "A"
        elif chain == "-1" or chain == -1:
            Chosen_chain = getAllChains(fromFile)
            print(f"Chains: {Chosen_chain}")
        elif chain == "first":
            Chosen_chain = chains[0].id
        else:
            Chosen_chain = chain

        chains_to_remove = [i for i, x in enumerate(chains) if x.id not in Chosen_chain]
        fixer.removeChains(chains_to_remove)

        fixer.findMissingResidues()
        # add missing residues in the middle of a chain, not ones at the start or end of the chain.
        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        if verbose:
            print("missing residues: ", keys)
        if removeTwoEndsMissingResidues:
            for key in list(keys):
                chain_tmp = chains[key[0]]
                if key[1] == 0 or key[1] == len(list(chain_tmp.residues())):
                    del fixer.missingResidues[key]

        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingAtoms()
        try:
            fixer.addMissingAtoms()
        except:
            continue
        fixer.addMissingHydrogens(7.0)
        PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(toFolder, pdbFile), 'w'))

def get_inside_or_not_table(pdb_file):
    parser = PDBParser(PERMISSIVE=1,QUIET=True)
    try:
        structure = parser.get_structure('X', pdb_file)
    except:
        return [0]
    inside_or_not_table = []
    for res in structure.get_residues():
        if res.get_id()[0] != " ":
            continue  # skip
        try:
            res["CA"].get_vector()
        except:
            print(pdb_file, res.get_id())
            return [0]
        inside_or_not_table.append(int(abs(res["CA"].get_vector()[-1]) < 15))
    return inside_or_not_table

def extractTransmembrane(toLocation, location, cutoff=15):
    x = PDBParser().get_structure("x", location)

    class Transmembrane(Select):
        def accept_residue(self, residue):
            if residue.get_id()[0] == ' ' and abs(residue["CA"].get_vector()[-1]) < cutoff:
                return 1
            else:
                return 0

    io = PDBIO()
    io.set_structure(x)
    io.save(toLocation, Transmembrane())

def getSeqFromPDB(location, considerGap=True):
    x = PDBParser().get_structure("x", location)
    seq = ""
    resseqs = []
    preResId = 0
    for res in x.get_residues():
        resId = res.get_id()[1]
        if considerGap and resId != preResId + 1:
            seq += " "
            resseqs.append(-1)
        seq += three_to_one(res.get_resname())
        resseqs.append(res.get_id()[1])
        preResId = resId
    return seq,resseqs


def get_raw_optimization_data(pre, name):
    # pre = "/Users/weilu/Research/server_backup/feb_2019/jan_optimization/gammas/"
    # pp = "cath-dataset-nonredundant-S20Clean_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0"
    # pp = "proteins_name_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0"
    pp = name
    A_name = pp + "_A"
    B_name = pp + "_B"
    B_filtered_name = pp + "_B_filtered"
    P_name = pp + "_P"
    Gamma_name = pp + "_gamma"
    Gamma_filtered_name = pp + "_gamma_filtered"
    Lamb_name = pp + "_lamb"
    Lamb_filtered_name = pp + "_lamb_filtered"

    A_prime_name = pp + "_A_prime"
    A_prime = np.loadtxt(pre+A_prime_name)

    A = np.loadtxt(pre+A_name)
    B = np.loadtxt(pre+B_name)
    B_filtered = np.loadtxt(pre+B_filtered_name, dtype=complex, converters={
                               0: lambda s: complex(s.decode().replace('+-', '-'))})
    Gamma = np.loadtxt(pre+Gamma_name)
    Gamma_filtered = np.loadtxt(pre+Gamma_filtered_name, dtype=complex, converters={
                               0: lambda s: complex(s.decode().replace('+-', '-'))})
    Lamb = np.loadtxt(pre+Lamb_name, dtype=complex, converters={
                               0: lambda s: complex(s.decode().replace('+-', '-'))})
    Lamb_filtered = np.loadtxt(pre+Lamb_filtered_name, dtype=complex, converters={
                               0: lambda s: complex(s.decode().replace('+-', '-'))})

    half_B_name = pp + "_half_B"
    half_B = np.loadtxt(pre+half_B_name)
    other_half_B_name = pp + "_other_half_B"
    other_half_B = np.loadtxt(pre+other_half_B_name)
    std_half_B_name = pp + "_std_half_B"
    std_half_B = np.loadtxt(pre+std_half_B_name)
    return A,B,B_filtered,Gamma,Gamma_filtered,Lamb,Lamb_filtered,half_B,other_half_B,std_half_B,A_prime


def formatBurialGamma(burialGamma):
    # now, positive means favored.
    rhoGamma = pd.DataFrame(-burialGamma.reshape(3,20).T, columns=["rho1", "rho2", "rho3"]).reset_index()
    rhoGamma["oneLetter"] = rhoGamma["index"].apply(lambda x: inverse_res_type_map[x])
    rhoGamma["Residue"] = rhoGamma["index"].apply(lambda x: one_to_three(inverse_res_type_map[x]))
    rhoGamma = rhoGamma[["Residue", "rho1", "rho2", "rho3", "index", "oneLetter"]]
    g = rhoGamma[["rho1", "rho2", "rho3"]].values
    # np.savetxt("/Users/weilu/Research/server/feb_2019/burial_only_gamma.dat", g, fmt='%7.4f')
    # rhoGamma
    rhoGamma["hydrophobicityOrder"] = rhoGamma["oneLetter"].apply(lambda x: hydrophobicity_map[x])
    return rhoGamma.sort_values("hydrophobicityOrder")


def get_MSA_data(a3mFile):
    data = []
    # "/Users/weilu/Research/server/may_2019/family_fold/aligned/1r69.a3m"
    with open(a3mFile, "r") as f:
        for line in f:
            if line[0] == ">":
                continue
            s_new = ""
            seq = line.strip()
            for s in seq:
                if s.islower():
                    continue
                s_new += s
            data.append(s_new)
    return data


def read_gamma(gammaFile):
    data = np.loadtxt(gammaFile)
    gamma_direct = data[:210]
    gamma_mediated = data[210:]
    return gamma_direct, gamma_mediated

def get_gammas(gammaFile="./gamma.dat", memGammaFile="./membrane_gamma_rescaled.dat"):
    gamma_direct, gamma_mediated = read_gamma(gammaFile)
    nwell = 2
    gamma_ijm = np.zeros((nwell, 20, 20))
    water_gamma_ijm = np.zeros((nwell, 20, 20))
    protein_gamma_ijm = np.zeros((nwell, 20, 20))
    m = 0
    count = 0
    for i in range(20):
        for j in range(i, 20):
            gamma_ijm[m][i][j] = gamma_direct[count][0]
            gamma_ijm[m][j][i] = gamma_direct[count][0]
            count += 1
    count = 0
    for i in range(20):
        for j in range(i, 20):
            water_gamma_ijm[m][i][j] = gamma_mediated[count][1]
            water_gamma_ijm[m][j][i] = gamma_mediated[count][1]
            count += 1
    count = 0
    for i in range(20):
        for j in range(i, 20):
            protein_gamma_ijm[m][i][j] = gamma_mediated[count][0]
            protein_gamma_ijm[m][j][i] = gamma_mediated[count][0]
            count += 1
    if memGammaFile:
        mem_gamma_direct, mem_gamma_mediated = read_gamma(memGammaFile)
        m = 1  # membrane environment
        count = 0
        for i in range(20):
            for j in range(i, 20):
                gamma_ijm[m][i][j] = mem_gamma_direct[count][0]
                gamma_ijm[m][j][i] = mem_gamma_direct[count][0]
                count += 1
        count = 0
        for i in range(20):
            for j in range(i, 20):
                water_gamma_ijm[m][i][j] = mem_gamma_mediated[count][1]
                water_gamma_ijm[m][j][i] = mem_gamma_mediated[count][1]
                count += 1
        count = 0
        for i in range(20):
            for j in range(i, 20):
                protein_gamma_ijm[m][i][j] = mem_gamma_mediated[count][0]
                protein_gamma_ijm[m][j][i] = mem_gamma_mediated[count][0]
                count += 1
    return gamma_ijm, water_gamma_ijm, protein_gamma_ijm


def get_ff_dat(data, location=None, gammaLocation=None):
    res_type_map = {
        'A': 0,
        'C': 4,
        'D': 3,
        'E': 6,
        'F': 13,
        'G': 7,
        'H': 8,
        'I': 9,
        'K': 11,
        'L': 10,
        'M': 12,
        'N': 2,
        'P': 14,
        'Q': 5,
        'R': 1,
        'S': 15,
        'T': 16,
        'V': 19,
        'W': 17,
        'Y': 18
    }
    gamma_ijm, water_gamma_ijm, protein_gamma_ijm = get_gammas(f"{gammaLocation}/gamma.dat", memGammaFile=None)
    burial_gamma = np.loadtxt(f"{gammaLocation}/burial_gamma.dat")
    n = len(data[0])
    f_direct = np.zeros((n,n))
    f_water = np.zeros((n,n))
    f_protein = np.zeros((n,n))
    f_burial = np.zeros((n,3))
    for i in range(n):
        for j in range(i+1, n):
            direct = []
            water = []
            protein = []
            for seq in data:
                # seq = data[0]
                if seq[i] == "-" or seq[j] == "-":
                    continue
                if seq[i] == "X" or seq[j] == "X":
                    continue
                if seq[i] == "B" or seq[j] == "B":
                    continue
                if seq[i] == "Z" or seq[j] == "Z":
                    continue
                res1type = res_type_map[seq[i]]
                res2type = res_type_map[seq[j]]
                direct.append(gamma_ijm[0][res1type][res2type])
                water.append(water_gamma_ijm[0][res1type][res2type])
                protein.append(protein_gamma_ijm[0][res1type][res2type])
            f_direct[i][j] += np.average(direct)
            f_water[i][j] += np.average(water)
            f_protein[i][j] += np.average(protein)

            f_direct[j][i] += np.average(direct)
            f_water[j][i] += np.average(water)
            f_protein[j][i] += np.average(protein)

    for i in range(n):
        for j in range(3):
            burial = []
            for seq in data:
                if seq[i] == "-" or seq[i] == "X" or seq[i] == "B" or seq[i] == "Z":
                    continue
                res1type = res_type_map[seq[i]]
                burial.append(burial_gamma[res1type][j])
            f_burial[i][j] += np.average(burial)
    if location:
        np.savetxt(location+"/direct.dat", f_direct)
        np.savetxt(location+"/water.dat", f_water)
        np.savetxt(location+"/protein.dat", f_protein)
        np.savetxt(location+"/burial.dat", f_burial)
    return f_direct, f_water, f_protein, f_burial



def my_reorder(a, first):
    # move first to the top. and keep the rest
    new_order = first.copy()
    for col in a:
        if col not in first:
            new_order.append(col)
    return new_order

def read_pdb(pre, name, run=30, rerun=2):
    all_data = []
    if run == -1:
        run_list = ["native"]
    else:
        run_list = list(range(run))
    for i in run_list:
        if rerun == -1:
            rerun_list = ["rerun"]
        else:
            rerun_list = list(range(rerun))
        for j in rerun_list:
            # pre = "/Users/weilu/Research/server/nov_2018/iterative_optimization_4/all_simulations/"
            location = pre + f"{name}/simulation/{i}/{j}/"
            try:
                wham = pd.read_csv(location+"wham.dat")
            except:
                print(f"PDB: {name}, Run: {i}, Rerun: {j} not exist")
                print(location+"wham.dat")
                continue
            wham.columns = wham.columns.str.strip()
            remove_columns = ['Tc', 'Energy']
            wham = wham.drop(remove_columns, axis=1)
            energy = pd.read_csv(location+"energy.dat")
            energy.columns = energy.columns.str.strip()
            remove_columns = ['Steps', 'Shake', 'Excluded', 'Helix', 'AMH-Go', 'Vec_FM', 'SSB']
            energy = energy.drop(remove_columns, axis=1)
            data = pd.concat([wham, energy], axis=1).assign(Repeat=i, Run=j)
            all_data.append(data)
    data = pd.concat(all_data).reset_index(drop=True)
    data = data.reindex(columns=my_reorder(data.columns, ["Steps", "Qw", "VTotal", "Run", "Repeat"]))
    print(name, len(data))
    return data

def get_complete_data(pre, folder_list, pdb_list, formatName=True, **kwargs):
    complete_all_data = []
    for folder in folder_list:
        # pre = "/Users/weilu/Research/server/april_2019/iterative_optimization_old_set/"
        pre_folder = f"{pre}{folder}/"
        all_data = []
        for p in pdb_list:
            if formatName:
                name = p.lower()[:4]
            else:
                name = p
            tmp = read_pdb(pre_folder, name, **kwargs)
            all_data.append(tmp.assign(Name=name))
        data = pd.concat(all_data)
        complete_all_data.append(data.assign(Folder=folder))
    data = pd.concat(complete_all_data)
    data = data.reindex(columns=my_reorder(data.columns, ["Name", "Folder"]))
    return data


def convertRaptorToInput(pdbID, raptorX_file):
        # pdbID = "2xov_complete_2"
        # read in median distances for pairwise interactions (obtained from analysis of the pdb)
        directory='/Users/weilu/opt/gremlin/'
        distancesCACB=pd.read_csv(directory+'CACBmediandist.dat', delim_whitespace=True, header=None)
        distancesCACA=pd.read_csv(directory+'CACAmediandist.dat', delim_whitespace=True, header=None)
        distancesCBCB=pd.read_csv(directory+'CBCBmediandist.dat', delim_whitespace=True, header=None)
        distancesCACB.columns = ['i', 'j', 'dist']
        distancesCACA.columns = ['i', 'j', 'dist']
        distancesCBCB.columns = ['i', 'j', 'dist']
        # if you want to filter the gremlin data, adjust the parameters below
        filter_threshold=0.5
        column=4
        name='raptorX.5.'


        # make sure that there is a sequence file for the protein and the downloaded gremlin data in the proper directories
        # read and parse raptorX contact prediction data
        # directory = "/Users/weilu/opt/gremlin/protein/" + pdbID + "/"
        # raptorX_file=directory+"raptor." + pdbID + ".dat"
        count=0
        pairs=[]
        seq=''
        with open(raptorX_file) as f:
            for line in f:
                count+=1
                if count<6:
                    continue
                elif line.split()[0]!='END' and len(line.split())==5:
                    pairs.append(line.split())
                elif len(pairs)==0:
                    seq+=line.split()[0]
        n=len(seq)
    #     print(n)
        rnative_matrixCACB=np.ones([n,n])*99
        rnative_matrixCACA=np.ones([n,n])*99
        rnative_matrixCBCB=np.ones([n,n])*99
        for pair in pairs:
            i=int(pair[0])
            j=int(pair[1])
            irestype=seq[i-1]
            jrestype=seq[j-1]
            if float(pair[column]) > filter_threshold:
                if sum((distancesCACB['i']==irestype)&(distancesCACB['j']==jrestype))>0: #check if pair is in correct order
                    well_centerCACB = distancesCACB[(distancesCACB['i']==irestype)&(distancesCACB['j']==jrestype)]['dist'].values[0]
                    well_centerCACA = distancesCACA[(distancesCACA['i']==irestype)&(distancesCACA['j']==jrestype)]['dist'].values[0]
                    well_centerCBCB = distancesCBCB[(distancesCBCB['i']==irestype)&(distancesCBCB['j']==jrestype)]['dist'].values[0]
                else:
                    well_centerCACB = distancesCACB[(distancesCACB['i']==jrestype)&(distancesCACB['j']==irestype)]['dist'].values[0]
                    well_centerCACA = distancesCACA[(distancesCACA['i']==jrestype)&(distancesCACA['j']==irestype)]['dist'].values[0]
                    well_centerCBCB = distancesCBCB[(distancesCBCB['i']==jrestype)&(distancesCBCB['j']==irestype)]['dist'].values[0]

                rnative_matrixCACB[i-1, j-1] = well_centerCACB
                rnative_matrixCACB[j-1, i-1] = well_centerCACB
                rnative_matrixCACA[i-1, j-1] = well_centerCACA
                rnative_matrixCACA[j-1, i-1] = well_centerCACA
                rnative_matrixCBCB[i-1, j-1] = well_centerCBCB
                rnative_matrixCBCB[j-1, i-1] = well_centerCBCB
        import matplotlib.pyplot as plt

        plt.matshow(rnative_matrixCACB)
        # plt.show()
        fig = plt.gcf()
        directory = "/Users/weilu/opt/gremlin/protein/" + pdbID + "/raptor/"
        os.system("mkdir -p " + directory)
        figureDirectory = f"{directory}/contact.png"
        fig.savefig(figureDirectory)
        os.system(f"cp {raptorX_file} {directory}")

        np.savetxt(directory + 'go_rnativeCACB.dat', rnative_matrixCACB, fmt='%10.5f')
        np.savetxt(directory + 'go_rnativeCACA.dat', rnative_matrixCACA, fmt='%10.5f')
        np.savetxt(directory + 'go_rnativeCBCB.dat', rnative_matrixCBCB, fmt='%10.5f')

def get_contactFromDMP(fileLocation, n, threshold=0.2):
    a = np.zeros((n,n))
    c_list = []
    with open(fileLocation, "r") as f:
    #     for i in range(9):
    #         next(f)
        for line in f:
    #         print(line)
            try:
                i,j,_,_,_,p = line.split(" ")
    #             print(i,j,p)
                a[int(i)-1,int(j)-1] = float(p)
                a[int(j)-1,int(i)-1] = float(p)
                if float(p) > threshold:
                    c_list.append([int(i),int(j),float(p)])
            except Exception as e:
                print(e)
                pass
    return a, np.array(c_list)

def convertDMPToInput(pdbID, dmp_file, fasta_file, pre='/Users/weilu/opt/gremlin/'):
        # pdbID = "2xov_complete_2"
        # read in median distances for pairwise interactions (obtained from analysis of the pdb)
        directory=pre
        distancesCACB=pd.read_csv(directory+'CACBmediandist.dat', delim_whitespace=True, header=None)
        distancesCACA=pd.read_csv(directory+'CACAmediandist.dat', delim_whitespace=True, header=None)
        distancesCBCB=pd.read_csv(directory+'CBCBmediandist.dat', delim_whitespace=True, header=None)
        distancesCACB.columns = ['i', 'j', 'dist']
        distancesCACA.columns = ['i', 'j', 'dist']
        distancesCBCB.columns = ['i', 'j', 'dist']
        # if you want to filter the gremlin data, adjust the parameters below
        filter_threshold=0.5
        column=2

        seq = ""
        with open(fasta_file) as f:
            for line in f:
                if line[0] == ">":
                    continue
                seq += line.strip()
        # seq

        n=len(seq)
        _, dmp_pairs = get_contactFromDMP(dmp_file, n=n)

    #     print(n)
        rnative_matrixCACB=np.ones([n,n])*99
        rnative_matrixCACA=np.ones([n,n])*99
        rnative_matrixCBCB=np.ones([n,n])*99
        for pair in dmp_pairs:
            i=int(pair[0])
            j=int(pair[1])
            irestype=seq[i-1]
            jrestype=seq[j-1]
            if float(pair[column]) > filter_threshold:
                if sum((distancesCACB['i']==irestype)&(distancesCACB['j']==jrestype))>0: #check if pair is in correct order
                    well_centerCACB = distancesCACB[(distancesCACB['i']==irestype)&(distancesCACB['j']==jrestype)]['dist'].values[0]
                    well_centerCACA = distancesCACA[(distancesCACA['i']==irestype)&(distancesCACA['j']==jrestype)]['dist'].values[0]
                    well_centerCBCB = distancesCBCB[(distancesCBCB['i']==irestype)&(distancesCBCB['j']==jrestype)]['dist'].values[0]
                else:
                    well_centerCACB = distancesCACB[(distancesCACB['i']==jrestype)&(distancesCACB['j']==irestype)]['dist'].values[0]
                    well_centerCACA = distancesCACA[(distancesCACA['i']==jrestype)&(distancesCACA['j']==irestype)]['dist'].values[0]
                    well_centerCBCB = distancesCBCB[(distancesCBCB['i']==jrestype)&(distancesCBCB['j']==irestype)]['dist'].values[0]

                rnative_matrixCACB[i-1, j-1] = well_centerCACB
                rnative_matrixCACB[j-1, i-1] = well_centerCACB
                rnative_matrixCACA[i-1, j-1] = well_centerCACA
                rnative_matrixCACA[j-1, i-1] = well_centerCACA
                rnative_matrixCBCB[i-1, j-1] = well_centerCBCB
                rnative_matrixCBCB[j-1, i-1] = well_centerCBCB
        import matplotlib.pyplot as plt

        plt.imshow(rnative_matrixCACB, origin=0)
        # plt.show()
        fig = plt.gcf()
        directory = f"{pre}/protein/" + pdbID + "/DMP/"
        os.system("mkdir -p " + directory)
        figureDirectory = f"{directory}/contact.png"
        fig.savefig(figureDirectory)
        os.system(f"cp {dmp_file} {directory}")
        os.system(f"cp {fasta_file} {directory}")

        np.savetxt(directory + 'go_rnativeCACB.dat', rnative_matrixCACB, fmt='%10.5f')
        np.savetxt(directory + 'go_rnativeCACA.dat', rnative_matrixCACA, fmt='%10.5f')
        np.savetxt(directory + 'go_rnativeCBCB.dat', rnative_matrixCBCB, fmt='%10.5f')

def get_PredictedZim(topo, zimFile):
    loc = topo
    with open(loc) as f:
        a = f.readlines()
    assert len(a) % 3 == 0
    chain_count = len(a) // 3
    seq = ""
    for i in range(chain_count):
        seq_i = (a[i*3+2]).strip()
        seq += seq_i
    assert np.alltrue([i in ["0", "1"] for i in seq])

    with open(zimFile, "w") as out:
        for i in seq:
            if i == "0":
                out.write(f"1\n")
            elif i == "1":
                out.write("2\n")
            else:
                raise

def flip_side(side):
#     return side * -1
    if side == "down":
        return "up"
    if side == "up":
        return "down"

def get_PredictedZimSide(topo, zimFile):
    loc = topo
    with open(loc) as f:
        a = f.readlines()
    assert len(a) % 3 == 0
    chain_count = len(a) // 3
    seq = ""
    for i in range(chain_count):
        seq_i = (a[i*3+2]).strip()
        seq += seq_i
    assert np.alltrue([i in ["0", "1"] for i in seq])

    side = "down"
    with open(zimFile, "w") as out:
        for i in seq:
            if i == "0":
                out.write(f"{side}\n")
                inMiddle = False
            elif i == "1":
                out.write("middle\n")
                if not inMiddle:
                    side = flip_side(side)
                    inMiddle = True
            else:
                raise

# def get_inside_or_not_table(pdb_file):
#     parser = PDBParser(PERMISSIVE=1)
#     structure = parser.get_structure('X', pdb_file)
#     inside_or_not_table = []
#     for res in structure.get_residues():
#         inside_or_not_table.append(int(abs(res["CA"].get_vector()[-1]) < 15))
#     return inside_or_not_table

# def getSeqFromPDB(location):
#     x = PDBParser().get_structure("x", location)
#     ppb=PPBuilder()
#     seq = ""
#     for pp in ppb.build_peptides(x):
#         seq += str(pp.get_sequence())
#     return seq

# def relocate(location):
#     # location = "/Users/weilu/Research/server/april_2019/iterative_optimization_new_set_with_frag/all_simulations/1fc2/1fc2"
#     fileLocation = location + "/frags.mem"
#     pre = location + "/../"
#     os.system(f"mkdir -p {pre}/fraglib")
#     with open(fileLocation) as f:
#         next(f)
#         next(f)
#         next(f)
#         next(f)
#         for line in f:
#             out = os.system(f"cp {line.split()[0]} {pre}/fraglib/")
#             if out != 0:
#                 print(f"!!Problem!!, {line.split()[0]}")
# def downloadPdb(pdb_list):
#     os.system("mkdir -p original_pdbs")
#     for pdb_id in pdb_list:
#         pdb = f"{pdb_id.lower()[:4]}"
#         pdbFile = pdb+".pdb"
#         if not os.path.isfile("original_pdbs/"+pdbFile):
#             pdbl = PDBList()
#             name = pdbl.retrieve_pdb_file(pdb, pdir='.', file_format='pdb')
#             os.system(f"mv {name} original_pdbs/{pdbFile}")



# def cleanPdb(pdb_list, chain=None, fromFolder="original_pdbs", toFolder="cleaned_pdbs"):
#     os.system(f"mkdir -p {toFolder}")
#     for pdb_id in pdb_list:
#         pdb = f"{pdb_id.lower()[:4]}"
#         pdbFile = pdb+".pdb"
#         if chain is None:
#             if len(pdb_id) == 5:
#                 Chosen_chain = pdb_id[4].upper()
#             else:
#                 assert(len(pdb_id) == 4)
#                 Chosen_chain = "A"
#         elif chain == "-1" or chain == -1:
#             Chosen_chain = getAllChains(os.path.join(fromFolder, pdbFile))
#         else:
#             Chosen_chain = chain
#         # clean pdb
#         fixer = PDBFixer(filename=os.path.join(fromFolder, pdbFile))
#         # remove unwanted chains
#         chains = list(fixer.topology.chains())
#         chains_to_remove = [i for i, x in enumerate(chains) if x.id not in Chosen_chain]
#         fixer.removeChains(chains_to_remove)

#         fixer.findMissingResidues()
#         # add missing residues in the middle of a chain, not ones at the start or end of the chain.
#         chains = list(fixer.topology.chains())
#         keys = fixer.missingResidues.keys()
#         # print(keys)
#         for key in list(keys):
#             chain = chains[key[0]]
#             if key[1] == 0 or key[1] == len(list(chain.residues())):
#                 del fixer.missingResidues[key]

#         fixer.findNonstandardResidues()
#         fixer.replaceNonstandardResidues()
#         fixer.removeHeterogens(keepWater=False)
#         fixer.findMissingAtoms()
#         fixer.addMissingAtoms()
#         fixer.addMissingHydrogens(7.0)
#         PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(toFolder, pdbFile), 'w'))

# def getAllChains(pdbFile):
#     fixer = PDBFixer(filename=pdbFile)
#     # remove unwanted chains
#     chains = list(fixer.topology.chains())
#     a = ""
#     for i in chains:
#         a += i.id
#     return ''.join(sorted(set(a.upper().replace(" ", ""))))


# def add_chain_to_pymol_pdb(location):
#     # location = "/Users/weilu/Research/server/nov_2018/openMM/random_start/1r69.pdb"
#     with open("tmp", "w") as out:
#         with open(location, "r") as f:
#             for line in f:
#                 info = list(line)
#                 if len(info) > 21:
#                     info[21] = "A"
#                 out.write("".join(info))
#     os.system(f"mv tmp {location}")















# ----------------------------depreciated---------------------------------------



def read_simulation(location):
    file = "lipid.dat"
    lipid = pd.read_csv(location+file)
    lipid.columns = lipid.columns.str.strip()

    file = "energy.dat"
    energy = pd.read_csv(location+file)
    energy.columns = energy.columns.str.strip()
    file = "addforce.dat"
    dis = pd.read_csv(location+file)
    dis.columns = dis.columns.str.strip()
#     remove_columns = ['AddedForce', 'Dis12', 'Dis34', 'Dis56']
    file = "rgs.dat"
    rgs = pd.read_csv(location+file)
    rgs.columns = rgs.columns.str.strip()
    file = "wham.dat"
    wham = pd.read_csv(location+file)
    wham.columns = wham.columns.str.strip()
    remove_columns = ['Rg', 'Tc']
    wham = wham.drop(remove_columns, axis=1)
    data = wham.merge(rgs, how='inner', left_on=["Steps"], right_on=["Steps"]).\
        merge(dis, how='inner', left_on=["Steps"], right_on=["Steps"]).\
        merge(energy, how='inner', left_on=["Steps"], right_on=["Steps"]).\
        merge(lipid, how='inner', left_on=["Steps"], right_on=["Steps"])
    data = data.assign(TotalE=data.Energy + data.Lipid)
    return data


def process_complete_temper_data_2(pre, data_folder, folder_list, rerun=-1, n=12, bias="dis", qnqc=False, average_z=False, localQ=False):
    print("process temp data")
    dateAndTime = datetime.today().strftime('%d_%h_%H%M%S')
    for folder in folder_list:
        simulation_list = glob.glob(pre+folder+f"/simulation/{bias}_*")
        print(pre+folder+f"/simulation/{bias}_*")
        os.system("mkdir -p " + pre+folder+"/data")
        # this one only consider rerun >=0, for the case rerun=-1, move log.lammps to log0
        for i in range(rerun+1):
            all_data_list = []
            for one_simulation in simulation_list:
                bias_num = one_simulation.split("_")[-1]
                print(bias_num, "!")

                location = one_simulation + f"/{i}/"
                print(location)
                data = read_complete_temper_2(location=location, n=n, rerun=i, qnqc=qnqc, average_z=average_z, localQ=localQ)
                print(data.shape)
                # remove_columns = ['Step', "Run"]
                # data = data.drop(remove_columns, axis=1)
                all_data_list.append(data.assign(BiasTo=bias_num))

            data = pd.concat(all_data_list).reset_index(drop=True)
            # if localQ:
            #     print("hi")
            # else:
            #     data.to_csv(os.path.join(pre, folder, f"data/rerun_{i}.csv"))
            # complete_data_list.append(data)
            #         temps = list(dic.keys())
            # complete_data = pd.concat(complete_data_list)
            name = f"rerun_{i}_{dateAndTime}.feather"
            data.reset_index(drop=True).to_feather(pre+folder+"/" + name)
            os.system("cp "+pre+folder+"/" + name + " "+data_folder)

def move_data3(data_folder, freeEnergy_folder, folder, sub_mode_name="", kmem=0.2, klipid=0.1, kgo=0.1, krg=0.2, sample_range_mode=0, biasName="dis", qnqc=False, average_z=0, chosen_mode=0):
    print("move data")
    dic = {"T0":350, "T1":400, "T2":450, "T3":500, "T4":550, "T5":600, "T6":650, "T7":700, "T8":750, "T9":800, "T10":900, "T11":1000}
    # read in complete.feather
    data = pd.read_feather(data_folder + folder +".feather")
    os.system("mkdir -p "+freeEnergy_folder+"/"+sub_mode_name+f"/data_{sample_range_mode}")
    for bias, oneBias in data.groupby("BiasTo"):
        for tempSymbol, oneTempAndBias in oneBias.groupby("Temp"):
            temp = dic[tempSymbol]
            if float(temp) > 800:
                continue
            print(f"t_{temp}_{biasName}_{bias}.dat")
            if sample_range_mode == 0:
                queryCmd = 'Step > 0 & Step <= 1e7'
            if sample_range_mode == 1:
                queryCmd = 'Step > 1e7 & Step <= 2e7'
            elif sample_range_mode == 2:
                queryCmd ='Step > 2e7 & Step <= 3e7'
            elif sample_range_mode == 3:
                queryCmd ='Step > 3e7 & Step <= 4e7'
            elif sample_range_mode == 4:
                queryCmd ='Step > 4e7 & Step <= 5e7'
            elif sample_range_mode == 5:
                queryCmd ='Step > 5e7 & Step <= 6e7'
            elif sample_range_mode == -1:
                queryCmd ='Step > 4e7 & Step <= 6e7'
            tmp = oneTempAndBias.query(queryCmd)
            chosen_list = ["TotalE", "Qw", "Distance"]
            if average_z == 1:
                chosen_list += ["abs_z_average"]
            if average_z == 2:
                chosen_list += ["z_h6"]
            if chosen_mode == 0:
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_perturb_mem_p=tmp.TotalE + kmem*tmp.Membrane,
                                        TotalE_perturb_mem_m=tmp.TotalE - kmem*tmp.Membrane,
                                        TotalE_perturb_lipid_p=tmp.TotalE + klipid*tmp.Lipid,
                                        TotalE_perturb_lipid_m=tmp.TotalE - klipid*tmp.Lipid,
                                        TotalE_perturb_go_p=tmp.TotalE + kgo*tmp["AMH-Go"],
                                        TotalE_perturb_go_m=tmp.TotalE - kgo*tmp["AMH-Go"],
                                        TotalE_perturb_rg_p=tmp.TotalE + krg*tmp.Rg,
                                        TotalE_perturb_rg_m=tmp.TotalE - krg*tmp.Rg)
            if chosen_mode == 1:
                chosen_list += ["Res" + str(i+1) for i in range(181)]
                chosen = tmp[chosen_list]
    #         print(tmp.count())
            chosen.to_csv(freeEnergy_folder+"/"+sub_mode_name+f"/data_{sample_range_mode}/t_{temp}_{biasName}_{bias}.dat", sep=' ', index=False, header=False)


def move_data2(data_folder, freeEnergy_folder, folder, sub_mode_name="", kmem=0.2, klipid=0.1, kgo=0.1, krg=0.2, sample_range_mode=0, biasName="dis", qnqc=False, average_z=0, chosen_mode=0):
    print("move data")
    dic = {"T0":350, "T1":400, "T2":450, "T3":500, "T4":550, "T5":600, "T6":650, "T7":700, "T8":750, "T9":800, "T10":900, "T11":1000}
    # read in complete.feather
    data = pd.read_feather(data_folder + folder +".feather")
    os.system("mkdir -p "+freeEnergy_folder+folder+sub_mode_name+"/data")
    for bias, oneBias in data.groupby("BiasTo"):
        for tempSymbol, oneTempAndBias in oneBias.groupby("Temp"):
            temp = dic[tempSymbol]
            if float(temp) > 800:
                continue
            print(f"t_{temp}_{biasName}_{bias}.dat")
            if sample_range_mode == 0:
                queryCmd = 'Step > 1e7 & Step <= 2e7'
            elif sample_range_mode == 1:
                queryCmd ='Step > 2e7 & Step <= 3e7'
            elif sample_range_mode == 2:
                queryCmd ='Step > 3e7 & Step <= 4e7'
            elif sample_range_mode == 3:
                queryCmd ='Step > 4e7 & Step <= 5e7'
            elif sample_range_mode == 4:
                queryCmd ='Step > 5e7 & Step <= 6e7'
            elif sample_range_mode == -1:
                queryCmd ='Step > 4e7 & Step <= 6e7'
            tmp = oneTempAndBias.query(queryCmd)
            chosen_list = ["TotalE", "Qw", "Distance"]
            if average_z == 1:
                chosen_list += ["abs_z_average"]
            if average_z == 2:
                chosen_list += ["z_h6"]
            if chosen_mode == 0:
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_perturb_mem_p=tmp.TotalE + kmem*tmp.Membrane,
                                        TotalE_perturb_mem_m=tmp.TotalE - kmem*tmp.Membrane,
                                        TotalE_perturb_lipid_p=tmp.TotalE + klipid*tmp.Lipid,
                                        TotalE_perturb_lipid_m=tmp.TotalE - klipid*tmp.Lipid,
                                        TotalE_perturb_go_p=tmp.TotalE + kgo*tmp["AMH-Go"],
                                        TotalE_perturb_go_m=tmp.TotalE - kgo*tmp["AMH-Go"],
                                        TotalE_perturb_rg_p=tmp.TotalE + krg*tmp.Rg,
                                        TotalE_perturb_rg_m=tmp.TotalE - krg*tmp.Rg)
            if chosen_mode == 1:
                chosen_list += ["Res" + str(i+1) for i in range(181)]
                chosen = tmp[chosen_list]
    #         print(tmp.count())
            chosen.to_csv(freeEnergy_folder+folder+sub_mode_name+f"/data/t_{temp}_{biasName}_{bias}.dat", sep=' ', index=False, header=False)
    # chosen

def make_metadata_2(cwd=".", k=1000.0, temps_list=["450"]):
    files = glob.glob("../../data/*")
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
                target = "{} {} {} {}\n".format(oneFile, t, kconstant, bias)
                out.write(target)

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
def read_complete_temper(n=4, location=".", rerun=-1, qnqc=False, average_z=False, localQ=False):
    all_lipid_list = []
    for i in range(n):
        file = "lipid.{}.dat".format(i)
        lipid = pd.read_csv(location+file).assign(Run=i)
        lipid.columns = lipid.columns.str.strip()
        # lipid = lipid[["Steps","Lipid","Run"]]
        all_lipid_list.append(lipid)
    lipid = pd.concat(all_lipid_list)

    all_rgs_list = []
    for i in range(n):
        file = "rgs.{}.dat".format(i)
        rgs = pd.read_csv(location+file).assign(Run=i)
        rgs.columns = rgs.columns.str.strip()
        # lipid = lipid[["Steps","Lipid","Run"]]
        all_rgs_list.append(rgs)
    rgs = pd.concat(all_rgs_list)

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
        if qnqc:
            qc = pd.read_table(location+f"qc_{i}", names=["qc"])[1:].reset_index(drop=True)
            qn = pd.read_table(location+f"qn_{i}", names=["qn"])[1:].reset_index(drop=True)
            qc2 = pd.read_table(location+f"qc2_{i}", names=["qc2"])[1:].reset_index(drop=True)
            wham = pd.concat([wham, qn, qc, qc2],axis=1)
        if average_z:
            z = pd.read_table(location+f"z_{i}.dat", names=["AverageZ"])[1:].reset_index(drop=True)
            wham = pd.concat([wham, z],axis=1)
        if localQ:
            all_localQ = pd.read_csv(location+f"localQ.{i}.csv")[1:].reset_index(drop=True)
            wham = pd.concat([wham, all_localQ], axis=1)
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
    t6 = t5.merge(rgs, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]).sort_values('Step').drop('Steps', axis=1)
    t6 = t6.assign(TotalE=t6.Energy + t6.Lipid)
    return t6
def process_complete_temper_data(pre, data_folder, folder_list, rerun=-1, n=12, bias="dis", qnqc=False, average_z=False, localQ=False):
    print("process temp data")
    for folder in folder_list:
        simulation_list = glob.glob(pre+folder+f"/simulation/{bias}_*")
        print(pre+folder+f"/simulation/{bias}_*")
        os.system("mkdir -p " + pre+folder+"/data")
        complete_data_list = []
        for one_simulation in simulation_list:
            bias_num = one_simulation.split("_")[-1]
            print(bias_num, "!")
            all_data_list = []
            if rerun == -1:
                location = one_simulation + "/0/"
                print(location)
                data = read_complete_temper(location=location, n=n, rerun=rerun, qnqc=qnqc, average_z=average_z, localQ=localQ)
                # remove_columns = ['Step', "Run"]
                # data = data.drop(remove_columns, axis=1)
                all_data_list.append(data)
            else:
                for i in range(rerun+1):
                    location = one_simulation + f"/{i}/"
                    print(location)
                    data = read_complete_temper(location=location, n=n, rerun=i, qnqc=qnqc, average_z=average_z, localQ=localQ)
                    # remove_columns = ['Step', "Run"]
                    # data = data.drop(remove_columns, axis=1)
                    all_data_list.append(data)

            data = pd.concat(all_data_list).assign(BiasTo=bias_num).reset_index(drop=True)
            if localQ:
                print("hi")
            else:
                data.to_csv(os.path.join(pre, folder, f"data/bias_num.csv"))
            complete_data_list.append(data)
    #         temps = list(dic.keys())
        complete_data = pd.concat(complete_data_list)
        name = f"{datetime.today().strftime('%d_%h_%H%M%S')}.feather"
        complete_data.reset_index(drop=True).to_feather(pre+folder+"/" + name)
        os.system("cp "+pre+folder+"/" + name + " "+data_folder+folder+".feather")
def read_temper(n=4, location=".", rerun=-1, qnqc=False):
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
        if qnqc:
            qc = pd.read_table(location+f"qc_{i}", names=["qc"])[1:].reset_index(drop=True)
            qn = pd.read_table(location+f"qn_{i}", names=["qn"])[1:].reset_index(drop=True)
            qc2 = pd.read_table(location+f"qc2_{i}", names=["qc2"])[1:].reset_index(drop=True)
            wham = pd.concat([wham,qn, qc, qc2],axis=1)
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


def process_temper_data(pre, data_folder, folder_list, rerun=-1, n=12, bias="dis", qnqc=False):
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
                    data = read_temper(location=location, n=n, qnqc=qnqc)
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
                        data = read_temper(location=location, n=n, rerun=i, qnqc=qnqc)
                        # remove_columns = ['Step', "Run"]
                        # data = data.drop(remove_columns, axis=1)
                        all_data_list.append(data)
                    except:
                        print("Unexpected error:", sys.exc_info()[0])
                        print("notrun?", bias_num)
                try:
                    data = pd.concat(all_data_list)
                    data.reset_index(drop=True).to_feather(pre+folder+"/data/"+f"{bias}{bias_num}.feather")
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
