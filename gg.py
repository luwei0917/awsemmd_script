#!/usr/bin/env python
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
from myPersonalFunctions import *
import glob
from small_script.myFunctions import *

try:
    OPENAWSEM_LOCATION = os.environ["OPENAWSEM_LOCATION"]
    sys.path.insert(0, OPENAWSEM_LOCATION)
    # print(OPENAWSEM_LOCATION)
except KeyError:
    print("Please set the environment variable name OPENAWSEM_LOCATION.\n Example: export OPENAWSEM_LOCATION='YOUR_OPENAWSEM_LOCATION'")
    exit()

# from openmmawsem import *
from helperFunctions.myFunctions import *

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
parser.add_argument("-f", "--freeEnergy", help="free energy data sort ", action="store_true", default=False)
parser.add_argument("--fix", help="fix ", action="store_true", default=False)
parser.add_argument("--wham", help="wham analysis ", action="store_true", default=False)
parser.add_argument("--wham400", help="wham analysis in temp 400 ", action="store_true", default=False)
parser.add_argument("-p", "--plot", help="plot", action="store_true", default=False)
parser.add_argument("--pull", help="pull ", action="store_true", default=False)
parser.add_argument("--cpull", help="cpull ", action="store_true", default=False)
parser.add_argument("--energy", help="energy ", action="store_true", default=False)
parser.add_argument("--qnqc", help="calculate q of n terminal and q of c terminal ", action="store_true", default=False)
parser.add_argument("-n", "--number", type=int, default=10, help="number of run")
parser.add_argument("-d", "--day", type=str, default="someday")
parser.add_argument("-m", "--mode",type=int, default=0)
args = parser.parse_args()

with open('gg_cmd.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')


def do(cmd):
    subprocess.Popen(cmd, shell=True).wait()
cd = os.chdir

# def pick_structure():
#     with open("show.tcl", "w") as f:
#         structure_index = 1
#         f.write("mol new structure_%s.pdb\n" % structure_index)
#         f.write("mol modcolor 0 [molinfo top] Index\n")
#         f.write("mol modstyle 0 [molinfo top] NewCartoon 0.300000 10.000000 4.100000 0\n")
def pick_structure_generate_show_script(n=2):
    with open("show.pml", "w") as f:
        for structure_index in range(0, n):
            f.write("load structure_%s.pdb\n" % structure_index)
            f.write("cealign structure_0, structure_%s\n" % structure_index)
            f.write("spectrum count, rainbow_rev, structure_%s, byres=1\n" % structure_index)
        # f.write("hide lines, all\n")
        # f.write("show cartoon, all\n")
        # f.write("hide nonbonded, all\n")

dataset = {"old":"1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "),
            "new":"1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "),
            "test":["t089", "t120", "t251", "top7", "1ubq", "t0766", "t0778", "t0782", "t0792", "t0803", "t0815", "t0833", "t0842", "t0844"]}
dataset["combined"] = dataset["old"] + dataset["new"]


if args.day == "may11":
    if args.mode == 1:
        for i in range(1000):
            if i % 100 == 0:
                print(i)
            do(f"python3 ~/opt/compute_energy.py original_1r69_9/frame{i}.pdb  >> data_9_2.info")
    if args.mode == 2:
        for p in dataset["combined"]:
            name = p.lower()[:4]
            a3mFile = f"/Users/weilu/Research/server/may_2019/family_fold/aligned/{name}.a3m"
            pre = f"/Users/weilu/Research/server/may_2019/family_fold/ff_contact/{name}/"
            os.system(f"mkdir -p {pre}")
            data = get_MSA_data(a3mFile)
            print(name, len(data))
            f_direct_2, f_water_2, f_protein_2, f_burial_2 = get_ff_dat(data, location=pre)
if args.day == "may10":
    if args.mode == 1:
        a = glob.glob("cleaned_pdbs/*.pdb")
        do("mkdir -p fasta")
        for line in a:
            pdb = line.split("/")[-1].split(".")[0]
            pdbToFasta(pdb, line, f"fasta/{pdb}.fasta")
            cmd = f"hhblits -i fasta/{pdb}.fasta -d ~/Research/Build/hh-suite/uniclust30_2018_08/uniclust30_2018_08 -o {pdb} -oa3m {pdb}.a3m -n 1"
            do(cmd)
if args.day == "may08":
    if args.mode == 2:
        data = pd.read_csv("protein_info.csv", index_col=0)
        k_rg = 1
        pdb_list = dataset["combined"]
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {name}")
            cd(name)
            do(f"python3 /Users/weilu/openmmawsem/mm_create_project.py {name} --frag --extended")
            rg = data.query(f"Protein == '{name}'")["Rg"].values[0]
            replace(f"params.py", "40", str(rg))
            cd("..")
    if args.mode == 1:
        data = pd.read_csv("protein_info.csv", index_col=0)
        k_rg = 1
        pdb_list = dataset["combined"]
        with open("run.sh", "w") as out:
            for p in pdb_list:

                name = p.lower()[:4]
                cd(name)
                # out.write(f"cd {name}\n")
                for i in range(20):
                    cmd = f"python3 mm_run.py {name} -m 1 --to run_{i} --platform CUDA"
                    os.system(cmd)
                    # out.write(cmd)
                # out.write("cd ..\n")
                cd("..")
if args.day == "may05":
    if args.mode == 1:
        pdb_list = dataset["combined"]
        for pp in pdb_list:
            p = pp.lower()[:4]
            do(f"echo 'minimize       1.0e-4 1.0e-6 100 10000' >>  all_simulations/{p}/{p}/{p}_multi.in")

    if args.mode == 2:
        d = pd.read_csv("seq_info.csv", index_col=0)
        pdb_list = d.query("length < 150 and index % 2 == 0")["protein"].tolist()
        print(pdb_list)
        do("mkdir -p all_simulations")
        # downloadPdb(pdb_list)
        # cleanPdb(pdb_list, chain=None, formatName=False)
        cd("all_simulations")
        for p in pdb_list:
            # name = p.lower()[:4]
            name = p
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            # do("cp ~/opt/crystal_structures/membrane_proteins/for_simulation/{}.pdb crystal_structure.pdb ".format(name))
            do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
            do(f"create_project.py {name} --globular --frag")
            # check_and_correct_fragment_memory("frags.mem")
            relocate(fileLocation="frags.mem", toLocation="../fraglib")
            replace(f"frags.mem", "/Users/weilu/openmmawsem//Gros/", "../../fraglib/")
            protein_length = getFromTerminal("wc ssweight").split()[0]
            print(f"protein: {name}, length: {protein_length}")

            with open("single_frags.mem", "w") as out:
                out.write("[Target]\nquery\n\n[Memories]\n")
                out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")
if args.day == "may03":
    if args.mode == 1:
        # work on pdb fixer of membrane protein
        # pdb_list = ["1uaz", "1m0l", "1vgo", "2ei4"]
        # , '1lv7', '6mlu' bad
        # pdb_list = ['6c70', '5azb', '4r1i', '6bvg', '4zr1', '5o5e', '6eu6', '4pgr']
        # dd = pd.read_csv("/Users/weilu/Research/database/membrane_training_set/chosen_more_data.csv", index_col=0)
        dd = pd.read_csv("/Users/weilu/Research/database/membrane_training_set/chosen_large_data.csv", index_col=0)
        pdb_list = dd.pdbid.tolist()

        rest = []
        need = False
        for pdb in pdb_list:
            if need:
                rest.append(pdb)
            # if pdb == "5xth":
            # if pdb == "4v6m":
            # if pdb == "4v7i":
            if pdb == "6du8":
                need = True
        pdb_list = rest
        # for pdb in pdb_list:
        #     do(f"wget https://opm-assets.storage.googleapis.com/pdb/{pdb}.pdb")
        #     do(f"mv {pdb}.pdb original_pdbs/")

        # for i, pdb in enumerate(pdb_list):
        #     print(i, pdb)
        #     do(f"cp ../Alpha-helical_polytopic/{pdb}.pdb original_pdbs/")
        cleanPdb(pdb_list, chain="first", verbose=True, addMissingResidues=False)

if args.day == "apr30":
    if args.mode == 1:
        # work on pdb fixer of membrane protein
        # pdb_list = ["1uaz", "1m0l", "1vgo", "2ei4"]
        # , '1lv7', '6mlu' bad
        # pdb_list = ['6c70', '5azb', '4r1i', '6bvg', '4zr1', '5o5e', '6eu6', '4pgr']
        pdb_list = ['6aky', '1j4n', '1kpl', '6eu6', '4j05', '4nv6', '2zjs', '6gct', '2c3e', '4b4a', '4llh', '5y79', '5oc0', '3fi1', '4m58', '4a2n', '4uc1', '4qo2', '2q7r', '3h90', '3b4r', '2kdc', '2ksf', '3m73', '2m67', '5t77', '6bvg', '6bbg', '5vkv', '6d9z', '4r1i', '5tcx', '4quv', '4qtn', '5hwy', '3tij', '4av3', '2lop', '2lor', '2lom', '4il3', '3zd0', '5o5e', '4od5', '4o6m', '4pgr', '5jwy', '2mpn', '4q2e', '4x5m', '2mmu', '4rp9', '4zr1', '5a40', '4xu4', '5azb', '5dir', '5wud', '5n6h', '5vre', '5tsa', '5xj5', '6b87', '6cb2', '5zug', '6c70', '6bug', '6bar', '6bhp', '4cad', '6iu3']
        dd = pd.read_csv("/Users/weilu/Research/database/membrane_training_set/chosen.csv", index_col=0)
        pdb_list = dd.pdbid.tolist()
        for pdb in pdb_list:
            do(f"wget https://opm-assets.storage.googleapis.com/pdb/{pdb}.pdb")
            do(f"mv {pdb}.pdb original_pdbs/")
        # for i, pdb in enumerate(pdb_list):
        #     print(i, pdb)
        #     do(f"cp ../Alpha-helical_polytopic/{pdb}.pdb original_pdbs/")
        cleanPdb(pdb_list, chain="first", verbose=True, addMissingResidues=False)

if args.day == "apr18":
    if args.mode == 1:
        # check_and_correct_fragment_memory("frags.mem.bak")
        replace(f"frags.mem", "/Users/weilu/openmmawsem//Gros/", "../../fraglib/")
    if args.mode == 2:
        # do("minimize  1.0e-4 1.0e-6 100 10000")
        d = pd.read_csv("seq_info.csv", index_col=0)
        pdb_list = d.query("length < 150 and index % 2 == 0")["protein"].tolist()
        for p in pdb_list:
            do(f"echo 'minimize       1.0e-4 1.0e-6 100 10000' >>  all_simulations/{p}/{p}/{p}_multi.in")
if args.day == "apr17":
    if args.mode == 1:
        a = glob.glob("/Users/weilu/Research/database/chosen/*.pdb")
        print(len(a))
        for p in a:
            name = p.split("/")[-1]
            do(f"cp {p} original_pdbs/")
    if args.mode == 2:
        # pdb_list = ["2xov"]
        pdb_list = glob.glob("original_pdbs/*.pdb")
        # for p in pdb_list:
        #     add_chain_to_pymol_pdb(p)
        pdb_list = [p.split("/")[-1][:-4] for p in pdb_list]
        do("mkdir -p all_simulations")
        # downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain="first", formatName=False, verbose=False)
    if args.mode == 22:
        # pdb_list = glob.glob("original_pdbs/*.pdb")
        pdb_list = glob.glob("cleaned_pdbs/*.pdb")
        for i, p in enumerate(pdb_list):
            print(i, p.split("/")[-1][:-4], get_PDB_length(p))
    if args.mode == 3:
        # pdb_list = ["2xov"]
        d = pd.read_csv("seq_info.csv", index_col=0)
        pdb_list = d.query("length < 150 and index % 2 == 0")["protein"].tolist()
        print(pdb_list)
        print(len(pdb_list))
        do("mkdir -p all_simulations")
        # downloadPdb(pdb_list)
        # cleanPdb(pdb_list, chain=None, formatName=True)
        cd("all_simulations")
        for p in pdb_list:
            name = p
            # name = p
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            # do("cp ~/opt/crystal_structures/membrane_proteins/for_simulation/{}.pdb crystal_structure.pdb ".format(name))
            do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
            do(f"create_project.py {name} --frag --globular")
            check_and_correct_fragment_memory(fragFile="frags.mem")
            # relocate("./")
            relocate(fileLocation="frags.mem", toLocation="../fraglib")
            # replace(f"frags.mem", "/Users/weilu/opt/script/Gros/", "../fraglib/")
            # replace(f"frags.mem", "../fraglib/", "../../fraglib/")
            # replace(f"frags.mem", "/Users/weilu/openmmawsem/Gros/", "../../fraglib/")
            replace(f"frags.mem", "/Users/weilu/openmmawsem//Gros/", "../../fraglib/")
            # do("cp frags.mem fragsLAMW.mem")
            protein_length = getFromTerminal("wc ssweight").split()[0]
            print(f"protein: {name}, length: {protein_length}")
            with open("single_frags.mem", "w") as out:
                out.write("[Target]\nquery\n\n[Memories]\n")
                out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")

    if args.mode == 4:
        d = pd.read_csv("seq_info.csv", index_col=0)
        pdb_list = d.query("length < 150 and index % 2 == 0")["protein"].tolist()
        cd("all_simulations")
        for p in pdb_list:
            name = p
            cd(f"{name}/{name}")
            do("pwd")
            # replace(f"frags.mem", "/Users/weilu/openmmawsem//Gros/", "../../fraglib/")
            # replace(f"frags.mem", "/Users/weilu/openmmawsem//Gros/", "../../fraglib/")
            # relocate("./"
            relocate(fileLocation="frags.mem.bak", toLocation="../fraglib")
            cd("../..")
            # break
if args.day == "apr14":
    pdb_list = dataset["combined"]
    if args.mode == 1:
        do("mkdir -p all_simulations")
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain=None, formatName=True)
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            # name = p
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            # do("cp ~/opt/crystal_structures/membrane_proteins/for_simulation/{}.pdb crystal_structure.pdb ".format(name))
            do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
            do(f"create_project.py {name} --globular --frag")
            # check_and_correct_fragment_memory("frags.mem")
            relocate(fileLocation="frags.mem", toLocation="../fraglib")
            replace(f"frags.mem", "/Users/weilu/openmmawsem//Gros/", "../../fraglib/")
            protein_length = getFromTerminal("wc ssweight").split()[0]
            print(f"protein: {name}, length: {protein_length}")

            with open("single_frags.mem", "w") as out:
                out.write("[Target]\nquery\n\n[Memories]\n")
                out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")
if args.day == "apr06":
    if args.mode == 1:
        cmd = "calculate_rmsd.py 2xov.pdb 2xov.pdb "
        output = getFromTerminal(cmd)
        print(output)

if args.day == "apr05":
    if args.mode == 1:
        replace(f"frags.mem", "../../../fraglib/", "fraglib/")
if args.day == "apr04":
    if args.mode == 1:
        pdb_list = ["complete"]
        cleanPdb(pdb_list, chain=-1, formatName=False, verbose=True, removeTwoEndsMissingResidues=False)
    if args.mode == 2:
        # pdb_list = ["2xov"]
        pdb_list = dataset["combined"]
        do("mkdir -p all_simulations")
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain=None, formatName=True)
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            # name = p
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            # do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            # do("cp ~/opt/crystal_structures/membrane_proteins/for_simulation/{}.pdb crystal_structure.pdb ".format(name))
            do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
            do(f"create_project.py {name} --frag")
            check_and_correct_fragment_memory(fragFile="frags.mem")
            relocate("./")
            # replace(f"frags.mem", "/Users/weilu/opt/script/Gros/", "../fraglib/")
            # replace(f"frags.mem", "../fraglib/", "../../fraglib/")
            replace(f"frags.mem", "/Users/weilu/opt/script/Gros/", "../../fraglib/")
            do("cp frags.mem fragsLAMW.mem")
            protein_length = getFromTerminal("wc ssweight").split()[0]
            print(f"protein: {name}, length: {protein_length}")
            with open("single_frags.mem", "w") as out:
                out.write("[Target]\nquery\n\n[Memories]\n")
                out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")
if args.day == "apr03":
    if args.mode == 1:
        # downloadPdb(pdb_list)
        pdb_list = ["2bg9", "1j4n", "1py6", "2bl2", "1rhz", "1iwg", "2ic8", "1pv6", "1occ", "1kpl", "2bs2", "1py6", "1u19", "2rh1"]
        pdb_list = ["2xov"]
        pdb_list = ["2bs2", "1py6", "1u19", "2rh1"]
        do("mkdir -p original_pdbs")
        for pdb in pdb_list:
            do(f"wget https://opm-assets.storage.googleapis.com/pdb/{pdb}.pdb")
            do(f"mv {pdb}.pdb original_pdbs/")
        cleanPdb(pdb_list, chain=-1, formatName=False)
    if args.mode == 2:
        # pdb_list = ["2xov"]
        pdb_list = ["1occ", "1pv6", "2bl2", "2bg9", "1j4n", "1py6"]
        do("mkdir -p all_simulations")
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            # do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")

            do("cp ~/opt/crystal_structures/membrane_proteins/for_simulation/{}.pdb crystal_structure.pdb ".format(name))
            do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
            do(f"create_project.py {name} --hybrid --frag")
            check_and_correct_fragment_memory(fragFile="frags.mem")
            relocate("./")

            # replace(f"frags.mem", "/Users/weilu/opt/script/Gros/", "../fraglib/")
            # replace(f"frags.mem", "../fraglib/", "../../fraglib/")
            replace(f"frags.mem", "/Users/weilu/opt/script/Gros/", "../../../fraglib/")
            do("cp frags.mem fragsLAMW.mem")
            protein_length = getFromTerminal("wc ssweight").split()[0]
            print(f"protein: {name}, length: {protein_length}")
            with open("single_frags.mem", "w") as out:
                out.write("[Target]\nquery\n\n[Memories]\n")
                out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")
if args.day == "apr01":
    if args.mode == 1:
        # downloadPdb(pdb_list)
        pdb_list = ["5a63"]
        cleanPdb(pdb_list, chain=-1, formatName=False)
if args.day == "mar31":
    if args.mode == 1:
        iteration = "original"
        percent = 30
        pre = "./"
        iter_gamma = np.loadtxt("/Users/weilu/Research/server/april_2019/complete_gammas/original_gamma")
        gamma_for_simulation = pre + f"original_gamma.dat"
        burial_gamma_for_simulation = pre + f"original_burial_gamma.dat"
        gamma_format_convertion_iteration_to_simulation(iter_gamma, gamma_for_simulation, burial_gamma_for_simulation=burial_gamma_for_simulation)
        do("mv original_*.dat for_simulation/")
    if args.mode == 2:
        pre_i = 1
        i = 2
        pre = "./"
        g = "gammas/proteins_name_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_gamma_filtered"
        Gamma = np.loadtxt(f"/Users/weilu/Research/server/april_2019/optimization_restart_iter{i}/{g}")
        preGamma = np.loadtxt(f"/Users/weilu/Research/server/april_2019/complete_gammas/iter_{pre_i}")
        mix_gammas_3(pre, Gamma, preGamma, alpha=0.3, iterGammaName=f"iter_{i}", iteration=f"iter_{i}")
        do(f"mv iteration_iter_{i}_* for_simulation/")
    if args.mode == 3:
        pre_i = 2
        i = pre_i + 1
        pre = "./"
        g = "gammas/proteins_name_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_gamma_filtered"
        Gamma = np.loadtxt(f"/Users/weilu/Research/server/april_2019/optimization_restart_iter{i}/{g}")
        preGamma = np.loadtxt(f"/Users/weilu/Research/server/april_2019/complete_gammas/iter_{pre_i}")
        mix_gammas_3(pre, Gamma, preGamma, alpha=0.3, iterGammaName=f"iter_{i}", iteration=f"iter_{i}")
        do(f"mv iteration_iter_{i}_* for_simulation/")

    if args.mode == 4:
        pre_i = 2
        i = pre_i + 1
        alpha = 0.9
        pre = "./"
        g = "gammas/proteins_name_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_gamma_filtered"
        Gamma = np.loadtxt(f"/Users/weilu/Research/server/april_2019/optimization_restart_iter{i}/{g}")
        preGamma = np.loadtxt(f"/Users/weilu/Research/server/april_2019/complete_gammas/iter_{pre_i}")
        mix_gammas_3(pre, Gamma, preGamma, alpha=alpha, iterGammaName=f"iter_{i}_{int(alpha*100)}", iteration=f"iter_{i}")
        do(f"mv iteration_iter_{i}_* for_simulation/")
    if args.mode == 5:
        pre_i = 3
        i = pre_i + 1
        alpha = 0.3
        pre = "./"
        g = "gammas/proteins_name_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_gamma_filtered"
        Gamma = np.loadtxt(f"/Users/weilu/Research/server/april_2019/optimization_restart_iter{i}/{g}")
        preGamma = np.loadtxt(f"/Users/weilu/Research/server/april_2019/complete_gammas/iter_{pre_i}")
        if alpha == 0.3:
            mix_gammas_3(pre, Gamma, preGamma, alpha=alpha, iterGammaName=f"iter_{i}", iteration=f"iter_{i}")
        else:
            mix_gammas_3(pre, Gamma, preGamma, alpha=alpha, iterGammaName=f"iter_{i}_{int(alpha*100)}", iteration=f"iter_{i}")
        do(f"mv iteration_iter_{i}_* for_simulation/")

    if args.mode == 6:
        pre_i = 4
        i = pre_i + 1
        alpha = 0.3
        pre = "./"
        g = "gammas/proteins_name_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_gamma_filtered"
        Gamma = np.loadtxt(f"/Users/weilu/Research/server/april_2019/optimization_restart_iter{i}/{g}")
        preGamma = np.loadtxt(f"/Users/weilu/Research/server/april_2019/complete_gammas/iter_{pre_i}")
        if alpha == 0.3:
            mix_gammas_3(pre, Gamma, preGamma, alpha=alpha, iterGammaName=f"iter_{i}", iteration=f"iter_{i}")
        else:
            mix_gammas_3(pre, Gamma, preGamma, alpha=alpha, iterGammaName=f"iter_{i}_{int(alpha*100)}", iteration=f"iter_{i}")
        do(f"mv iteration_iter_{i}_* for_simulation/")

    if args.mode == 7:
        pre_i = 5
        i = pre_i + 1
        alpha = 0.3
        pre = "./"
        g = "gammas/proteins_name_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_gamma_filtered"
        Gamma = np.loadtxt(f"/Users/weilu/Research/server/april_2019/optimization_restart_iter{i}/{g}")
        preGamma = np.loadtxt(f"/Users/weilu/Research/server/april_2019/complete_gammas/iter_{pre_i}")
        if alpha == 0.3:
            mix_gammas_3(pre, Gamma, preGamma, alpha=alpha, iterGammaName=f"iter_{i}", iteration=f"iter_{i}")
        else:
            mix_gammas_3(pre, Gamma, preGamma, alpha=alpha, iterGammaName=f"iter_{i}_{int(alpha*100)}", iteration=f"iter_{i}")
        do(f"mv iteration_iter_{i}_* for_simulation/")

if args.day == "mar30":
    dataset = {"old":("1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "), 40),
                "new":("1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "), 80),
                "test":(["t089", "t120", "t251", "top7", "1ubq", "t0766", "t0778", "t0782", "t0792", "t0803", "t0815", "t0833", "t0842", "t0844"], 40)}

    if args.mode == 1:
        iteration = "0_normalized"
        percent = 30
        pre = "./"
        iter_gamma = np.loadtxt("/Users/weilu/Research/server/april_2019/complete_gammas/iter0_normalized_gamma")
        gamma_for_simulation = pre + f"iteration_{iteration}_gamma.dat"
        burial_gamma_for_simulation = pre + f"iteration_{iteration}_burial_gamma.dat"
        gamma_format_convertion_iteration_to_simulation(iter_gamma, gamma_for_simulation, burial_gamma_for_simulation=burial_gamma_for_simulation)
    if args.mode == 2:
        for tset in ["old", "new"]:
            pdb_list, steps = dataset[tset]
            downloadPdb(pdb_list)
            cleanPdb(pdb_list, chain=None, formatName=True)
            do("mkdir -p all_simulations")
            cd("all_simulations")
            for p in pdb_list:
                name = p.lower()[:4]
                # name = p
                do(f"mkdir -p {name}/{name}")
                cd(f"{name}/{name}")
                do("pwd")
                do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
                do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
                do(f"create_project.py {name} --globular --frag")
                check_and_correct_fragment_memory(fragFile="frags.mem")
                relocate("./")
                replace(f"frags.mem", "/Users/weilu/opt/script/Gros/", "../fraglib/")
                protein_length = getFromTerminal("wc ssweight").split()[0]
                print(f"protein: {name}, length: {protein_length}")
                with open("single_frags.mem", "w") as out:
                    out.write("[Target]\nquery\n\n[Memories]\n")
                    out.write(f"{name}.gro 1 1 {protein_length} 20\n")
                cd("../..")
                # do(f"echo '{name},{protein_length}\n' >> data_info")
            cd("..")
    if args.mode == 3:
        pre_i = 0
        i = 1
        pre = "./"
        g = "gammas/proteins_name_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_gamma_filtered"
        Gamma = np.loadtxt(f"/Users/weilu/Research/server/april_2019/optimization_restart_iter{i}/{g}")
        preGamma = np.loadtxt(f"/Users/weilu/Research/server/april_2019/complete_gammas/iter{pre_i}_normalized_gamma")
        mix_gammas_3(pre, Gamma, preGamma, alpha=0.3, iterGammaName=f"iter_{i}", iteration=f"iter_{i}")
        do("mv iteration_iter_1_* for_simulation/")
if args.day == "mar28":
    dataset = {"old":("1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "), 40),
                "new":("1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "), 80),
                "test":(["t089", "t120", "t251", "top7", "1ubq", "t0766", "t0778", "t0782", "t0792", "t0803", "t0815", "t0833", "t0842", "t0844"], 40)}
    # , "t0846"
    if args.mode == 1:
        pdb_list, steps = dataset["test"]
        do("mkdir -p all_simulations")
        cd("all_simulations")
        for p in pdb_list:
            # name = p.lower()[:4]
            name = p
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
            do(f"create_project.py {name} --globular")
            protein_length = getFromTerminal("wc ssweight").split()[0]
            print(f"protein: {name}, length: {protein_length}")
            with open("frags.mem", "w") as out:
                out.write("[Target]\nquery\n\n[Memories]\n")
                out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")
            do(f"echo '{name},{protein_length}\n' >> data_info")

if args.day == "mar27":
    # pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    dataset = {"old":("1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "), 40),
                "new":("1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "), 80),
                "test":(['1A2J', '1A3H', '1A5Y', '1A8Q', '1AGY', '1AKZ', '1AUZ', '1B1A', '1B31', '1B8X', '1BCO', '1BN6', '1BOH', '1BOI'], 40)}

    pdb_list, steps = dataset["old"]
    pdb_list, steps = dataset["new"]
    pdb_list = ["2lep"]
    if args.mode == 1:
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain=None, formatName=True)
    if args.mode == 2:
        do("mkdir -p all_simulations")
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
            do(f"create_project.py {name} --globular --frag")
            check_and_correct_fragment_memory(fragFile="frags.mem")
            relocate("./")
            # replace(f"frags.mem", "/Users/weilu/opt/script/Gros/", "../fraglib/")
            replace(f"frags.mem", "/Users/weilu/opt/script/Gros/", "../../fraglib/")
            # protein_length = getFromTerminal("wc ssweight").split()[0]
            # print(f"protein: {name}, length: {protein_length}")
            # with open("frags.mem", "w") as out:
            #     out.write("[Target]\nquery\n\n[Memories]\n")
            #     out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            do("cp frags.mem fragsLAMW.mem")
            protein_length = getFromTerminal("wc ssweight").split()[0]
            print(f"protein: {name}, length: {protein_length}")
            with open("single_frags.mem", "w") as out:
                out.write("[Target]\nquery\n\n[Memories]\n")
                out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")
            # break
    if args.mode == 22:
        do("mkdir -p all_simulations")
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
            do(f"create_project.py {name} --globular --frag")
            relocate("./")
            replace(f"frags.mem", "/Users/weilu/opt/script/Gros/", "../../fraglib/")
            # protein_length = getFromTerminal("wc ssweight").split()[0]
            # print(f"protein: {name}, length: {protein_length}")
            # with open("frags.mem", "w") as out:
            #     out.write("[Target]\nquery\n\n[Memories]\n")
            #     out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")

    if args.mode == 3:
        with open("/Users/weilu/Research/database/queriedPDB.dat") as f:
            a = f.readlines()
        a = [i.strip() for i in a]
        pdb_list = a[1:-1]
        # print(a)
        downloadPdb(pdb_list)
        for p in pdb_list:
            if p == "1E4F":
                pa = p + "T"
            elif p == "1FG5":
                pa = p + "N"
            elif p == "1G4W":
                pa = p + "R"
            elif p == "1GKU":
                pa = p + "B"
            elif p == "1J8M":
                pa = p + "F"
            elif p == "1J8Y":
                pa = p + "F"
            elif p == "1JPD":
                pa = p + "X"
            elif p == "1OEM":
                pa = p + "X"
            elif p == "1OEO":
                pa = p + "X"
            elif p == "1W97":
                pa = p + "L"
            elif p == "2PKY":
                pa = p + "X"
            elif p == "2VF4":
                pa = p + "X"
            elif p == "2YXP":
                pa = p + "X"
            else:
                continue
                # pa = p
            pdb_list_a = [pa]
            try:
                cleanPdb(pdb_list_a, chain=None, formatName=True)
            except Exception as e:
                print(f"error, {e}")
    if args.mode == 4:
        pdb_list = ["1E4F", "1FU1", "1XG8"]
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain=-1, formatName=True)
        # cleanPdb(pdb_list, chain=None, formatName=True)
    if args.mode == 5:
        pdb_list = ["t089", "t120", "t251", "top7", "1ubq", "t0766", "t0778", "t0782", "t0792", "t0803", "t0815", "t0833", "t0842", "t0844", "t0846"]
        # downloadPdb(pdb_list)
        # for p in pdb_list:
        #     fill_chain(f"{p}.pdb", f"original_pdbs/{p.lower()}.pdb")
        cleanPdb(pdb_list, chain=-1, formatName=False)
    if args.mode == 6:
        pdb_list = ["1qys", "1ubq"]
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain=-1, formatName=True)
if args.day == "mar26":
    if args.mode == 1:
        pre_i = 4
        pre = "./"
        g = "gammas/proteins_name_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_gamma_filtered"
        Gamma = np.loadtxt(f"/Users/weilu/Research/server/march_2019/optimization_weighted_by_q_iter{pre_i+1}/{g}")
        preGamma = np.loadtxt(f"/Users/weilu/Research/server/march_2019/optimization_weighted_by_q_iter{pre_i}/iter{pre_i}")
        mix_gammas_3(pre, Gamma, preGamma, alpha=0.3, iterGammaName=f"t_{pre_i+1}", iteration=f"test_{pre_i+1}")
    if args.mode == 2:
        pre_i = 4
        i = 6
        pre = "./"
        g = "gammas/proteins_name_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_gamma_filtered"
        Gamma = np.loadtxt(f"/Users/weilu/Research/server/march_2019/optimization_weighted_by_q_iter{i}/{g}")
        preGamma = np.loadtxt(f"/Users/weilu/Research/server/march_2019/optimization_weighted_by_q_iter{pre_i}/iter{pre_i}_thrid_gen")
        mix_gammas_3(pre, Gamma, preGamma, alpha=0.3, iterGammaName=f"t_{i}", iteration=f"test_{i}")
    if args.mode == 3:
        pre_i = 4
        i = 7
        pre = "./"
        g = "gammas/proteins_name_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_gamma_filtered"
        Gamma = np.loadtxt(f"/Users/weilu/Research/server/march_2019/optimization_weighted_by_q_iter{i}/{g}")
        preGamma = np.loadtxt(f"/Users/weilu/Research/server/march_2019/fix_gamma_mix_error/t_6")
        mix_gammas_3(pre, Gamma, preGamma, alpha=0.3, iterGammaName=f"t_{i}_2", iteration=f"test_{i}_2")
    if args.mode == 4:
        pre_i = 4
        i = 7
        pre = "./"
        g = "gammas/proteins_name_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_gamma_filtered"
        Gamma = np.loadtxt(f"/Users/weilu/Research/server/march_2019/optimization_weighted_by_q_iter{i}/{g}")
        mix_gammas_3(pre, Gamma, Gamma, alpha=0.3, iterGammaName=f"t_{i}_3", iteration=f"test_{i}_3")
    if args.mode == 5:
        iteration = "7_normalized"
        percent = 30
        pre = "./"
        iter_gamma = np.loadtxt("/Users/weilu/Research/server/april_2019/fix_gamma_mix_error/t_7_normalized")
        gamma_for_simulation = pre + f"iteration_{iteration}_gamma_{percent}.dat"
        burial_gamma_for_simulation = pre + f"iteration_{iteration}_burial_gamma_{percent}.dat"
        gamma_format_convertion_iteration_to_simulation(iter_gamma, gamma_for_simulation, burial_gamma_for_simulation=burial_gamma_for_simulation)
    if args.mode == 6:
        iteration = "multiSeq"
        percent = 30
        pre = "./"
        iter_gamma = np.loadtxt("/Users/weilu/Research/server/april_2019/fix_gamma_mix_error/multiSeq")
        gamma_for_simulation = pre + f"iteration_{iteration}_gamma_{percent}.dat"
        burial_gamma_for_simulation = pre + f"iteration_{iteration}_burial_gamma_{percent}.dat"
        gamma_format_convertion_iteration_to_simulation(iter_gamma, gamma_for_simulation, burial_gamma_for_simulation=burial_gamma_for_simulation)
if args.day == "mar25":
    pdb_list = pd.read_csv("data_info_2.csv", index_col=0)["FullName"].values
    if args.mode == 1:
        do("mkdir -p all_simulations")
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
            do(f"create_project.py {name} --globular --bias --rgbias")
            protein_length = getFromTerminal("wc ssweight").split()[0]
            print(f"protein: {name}, length: {protein_length}")
            with open("frags.mem", "w") as out:
                out.write("[Target]\nquery\n\n[Memories]\n")
                out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")

if args.day == "mar20":
    if args.mode == 1:
        # pdb_list = ["3GL5A"]
        pdb_list = pd.read_csv("data_info_2.csv", index_col=0)["FullName"].values
        cleanPdb(pdb_list, chain=None)
if args.day == "mar18":
    if args.mode == 1:
        pdb_list = ["3GL5A"]
        cleanPdb(pdb_list, chain=None)

if args.day == "feb27":
    dataset = {"old":("1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", "), 40),
                "new":("1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", "), 80),
                "test":(['1A2J', '1A3H', '1A5Y', '1A8Q', '1AGY', '1AKZ', '1AUZ', '1B1A', '1B31', '1B8X', '1BCO', '1BN6', '1BOH', '1BOI'], 40)}
    pdb_list, steps = dataset["new"]
    if args.mode == 1:
        # pdb_list, steps = dataset["old"]
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain=None)
    if args.mode == 2:
        do("mkdir -p all_simulations")
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
            do(f"create_project.py {name} --globular --bias --rgbias")
            protein_length = getFromTerminal("wc ssweight").split()[0]
            print(f"protein: {name}, length: {protein_length}")
            with open("frags.mem", "w") as out:
                out.write("[Target]\nquery\n\n[Memories]\n")
                out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")

if args.day == "feb25":
    location = "/Users/weilu/Research/database/queriedPDB.dat"
    with open(location) as f:
        n = int(next(f))
        print(n)
        all_lines = f.read().splitlines()
        first_20 = all_lines[:20]
        # print(first_20)
        pdb_list = ['1A2J', '1A3H', '1A5Y', '1A8Q', '1AGY', '1AKZ', '1AUZ', '1B1A', '1B31', '1B8X', '1BCO', '1BN6', '1BOH', '1BOI']
    if args.mode == 1:
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain=None)
    if args.mode == 2:
        do("mkdir -p all_simulations")
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
            do(f"create_project.py {name} --globular --bias --rgbias")
            protein_length = getFromTerminal("wc ssweight").split()[0]
            print(f"protein: {name}, length: {protein_length}")
            with open("frags.mem", "w") as out:
                out.write("[Target]\nquery\n\n[Memories]\n")
                out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")

if args.day == "feb23":
    pdb_list = "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", ")
    # pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    if args.mode == 1:
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain=None)
    if args.mode == 2:
        do("mkdir -p all_simulations")
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
            do(f"create_project.py {name} --globular")
            protein_length = getFromTerminal("wc ssweight").split()[0]
            print(f"protein: {name}, length: {protein_length}")
            with open("frags.mem", "w") as out:
                out.write("[Target]\nquery\n\n[Memories]\n")
                out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")


if args.day == "feb13":
    pdb_list = "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", ")
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    if args.mode == 1:
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain=None)
    if args.mode == 2:
        do("mkdir -p all_simulations")
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
            do(f"create_project.py {name} --globular")
            protein_length = getFromTerminal("wc ssweight").split()[0]
            print(f"protein: {name}, length: {protein_length}")
            with open("frags.mem", "w") as out:
                out.write("[Target]\nquery\n\n[Memories]\n")
                out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")

if args.day == "feb06":
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    if args.mode == 1:
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain=None)
    if args.mode == 2:
        do("mkdir -p all_simulations")
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
            do(f"create_project.py {name} --globular")
            protein_length = getFromTerminal("wc ssweight").split()[0]
            print(f"protein: {name}, length: {protein_length}")
            with open("frags.mem", "w") as out:
                out.write("[Target]\nquery\n\n[Memories]\n")
                out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")


if args.day == "jan14":
    pdb_list = "1FC2C, 1ENH, 2GB1, 2CRO, 1CTF, 4ICB".split(", ")
    if args.mode == 1:
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain=None)
    if args.mode == 2:
        do("mkdir -p all_simulations")
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
            do(f"create_project.py {name} --globular --frag")
            # protein_length = getFromTerminal("wc ssweight").split()[0]
            # print(f"protein: {name}, length: {protein_length}")
            # with open("frags.mem", "w") as out:
            #     out.write("[Target]\nquery\n\n[Memories]\n")
            #     out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")


if args.day == "jan07":
    if args.mode == 1:
        pdb_list = ["5fn2"]
        downloadPdb(pdb_list, membrane_protein=True)
        # downloadPdb(pdb_list, membrane_protein=False)
        print("clean")
        cleanPdb(pdb_list, chain=-1)
    name = "abeta42_1"
    rowN = 1
    columnN = 1
    if args.mode == 2:
        do("rm crystal_structure.pdb")
        to = "tmp1.pdb"

        table = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        count = 0
        for i in range(rowN):
            for j in range(columnN):
                duplicate_pdb("cleaned_pdbs/one_abeta42.pdb", to, offset_x=100, offset_y=-50, offset_z=50.0*j, new_chain="E")
                do(f"cat {to} >> crystal_structure.pdb")
                count += 1
if args.day == "jan06":
    pdb_list = "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", ")
    if args.mode == 1:
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain=None)
    if args.mode == 2:
        do("mkdir -p all_simulations")
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
            do(f"create_project.py {name} --globular")
            protein_length = getFromTerminal("wc ssweight").split()[0]
            print(f"protein: {name}, length: {protein_length}")
            with open("frags.mem", "w") as out:
                out.write("[Target]\nquery\n\n[Memories]\n")
                out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")


if args.day == "jan03":
    if args.mode == 1:
        convert_openMM_to_standard_pdb(fileName="movie.pdb")
if args.day == "dec13":
    if args.mode == 1:
        pdb_list = ["T0958"]
        cleanPdb(pdb_list, chain="-1", fromFile="crystal_structure.pdb", toFolder="cleaned_pdbs")

if args.day == "dec11":
    if args.mode == 1:
        get_frame("movie.pdb", "last_frame.pdb")
        convert_openMM_to_standard_pdb()
        do("~/opt/TMalign/TMalign last_frame.pdb ../crystal_structure.pdb -o result")
        do("cp result_all_atm result_all_atm.pdb")
        do("pymol ~/opt/plot_scripts/tmalign_all.pml")

if args.day == "dec09":
    if args.mode == 1:
        pdb_list = ["5a63"]
        cleanPdb(pdb_list, chain="-1", fromFolder=".", toFolder="cleaned_pdbs")

if args.day == "nov23":
    pdb_list = "5a63".split(", ")
    if args.mode == 1:
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain="ABC")
    if args.mode == 2:
        do("mkdir -p all_simulations")
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            do(f"create_project.py {name} --globular --frag")
            # do(f"create_project.py {name} --globular")
            cd("../..")

if args.day == "nov22":
    pdb_list = "1BRR".split(", ")
    if args.mode == 1:
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain="ABC")
    if args.mode == 2:
        do("mkdir -p all_simulations")
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            do(f"create_project.py {name} --globular --frag")
            # do(f"create_project.py {name} --globular")
            cd("../..")

if args.day == "nov12":
    pdb_list = "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", ")
    if args.mode == 1:
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain=None)
    if args.mode == 2:
        do("mkdir -p all_simulations")
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            do(f"create_project.py {name} --globular --frag")
            cd("../..")


if args.day == "nov09":
    pdb_list = "1R69, 1UTG, 3ICB, 256BA, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", ")
    if args.mode == 1:
        downloadPdb(pdb_list)
        cleanPdb(pdb_list)
    if args.mode == 2:
        do("mkdir -p all_simulations")
        cd("all_simulations")
        for p in pdb_list:
            name = p.lower()[:4]
            do(f"mkdir -p {name}/{name}")
            cd(f"{name}/{name}")
            do("pwd")
            do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")
            do(f"create_project.py {name} --globular")
            cd("../..")


if args.day == "nov02":
    name = "abeta42_2"
    rowN = 2
    columnN = 1
    if args.mode == 1:
        do("rm crystal_structure.pdb")
        to = "tmp1.pdb"

        table = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        count = 0
        for i in range(rowN):
            for j in range(columnN):
                duplicate_pdb("one_abeta42.pdb", to, offset_x=i*40.0, offset_y=j*40.0, offset_z=50.0*j, new_chain=table[count])
                do(f"cat {to} >> crystal_structure.pdb")
                count += 1
    if args.mode == 2:
        do(f"create_project.py {name} --frag --crystal --globular")
        do("cp ~/opt/abeta/ssweight_weihua ssweight")
        # do("cp ~/opt/abeta/zim12 zim")
        do("rm zim")
        for i in range(rowN*columnN):
            do("cat ~/opt/abeta/zim >> zim")
            do("cat ~/opt/abeta/zim >> zimPosition")
        do("cp ~/opt/abeta/seq.gamma .")
        do("cp ~/opt/abeta/fix_backbone_coeff.data .")
        do(f"cp data.crystal data.{name}")
    if args.mode == 3:
        do("python3 ~/opt/small_script/PBC_fixer.py")
        do(f"movie.py {name} -d new.dump")

if args.day == "oct15":
    name = "abeta42_1"
    rowN = 1
    columnN = 1
    if args.mode == 1:
        do("rm crystal_structure.pdb")
        to = "tmp1.pdb"

        table = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        count = 0
        for i in range(rowN):
            for j in range(columnN):
                duplicate_pdb("one_abeta42.pdb", to, offset_x=i*40.0, offset_y=j*40.0, offset_z=50.0*j, new_chain=table[count])
                do(f"cat {to} >> crystal_structure.pdb")
                count += 1
    if args.mode == 2:
        do(f"create_project.py {name} --frag --crystal --globular")
        do("cp ~/opt/abeta/ssweight_weihua ssweight")
        # do("cp ~/opt/abeta/zim12 zim")
        do("rm zim")
        for i in range(rowN*columnN):
            do("cat ~/opt/abeta/zim >> zim")
            do("cat ~/opt/abeta/zim >> zimPosition")
        do("cp ~/opt/abeta/seq.gamma .")
        do("cp ~/opt/abeta/fix_backbone_coeff.data .")
        do(f"cp data.crystal data.{name}")
    if args.mode == 3:
        do("python3 ~/opt/small_script/PBC_fixer.py")
        do(f"movie.py {name} -d new.dump")


if args.day == "sep10":
    name = "abeta42_6"
    rowN = 2
    columnN = 3
    if args.mode == 1:
        do("rm crystal_structure.pdb")
        to = "tmp1.pdb"

        table = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        count = 0
        for i in range(rowN):
            for j in range(columnN):
                duplicate_pdb("one_abeta42.pdb", to, offset_x=i*40.0, offset_y=j*40.0, offset_z=50.0*j-50, new_chain=table[count])
                do(f"cat {to} >> crystal_structure.pdb")
                count += 1
    if args.mode == 2:
        do(f"create_project.py {name} --frag --crystal --globular")
        do("cp ~/opt/abeta/ssweight_weihua ssweight")
        # do("cp ~/opt/abeta/zim12 zim")
        do("rm zim")
        for i in range(rowN*columnN):
            do("cat ~/opt/abeta/zim >> zim")
            do("cat ~/opt/abeta/zim >> zimPosition")
        do("cp ~/opt/abeta/seq.gamma .")
        do("cp ~/opt/abeta/fix_backbone_coeff.data .")
        do(f"cp data.crystal data.{name}")
    if args.mode == 3:
        do("python3 ~/opt/small_script/PBC_fixer.py")
        do(f"movie.py {name} -d new.dump")


if args.day == "aug15":
    cmd_pre = "python2 ~/opt/script/BuildAllAtomsFromLammps.py"
    # location_pre = "/Users/weilu/Research/server/apr_2018/sixth/rg_0.15_lipid_1.0_mem_1_go_0.8/simulation"
    location_pre = "/Users/weilu/Research/server/may_2018/second/simulation"
    # location_pre = "/Users/weilu/Research/server/may_2018/second_long/simulation"
    location_pre = "/Volumes/Wei_backup/GlpG/may_2018_back/second/simulation"
    # location_pre = "/Volumes/Wei_backup/GlpG/may_2018_back/second_long/simulation"
    # cmd = cmd_pre + " " + location + " structure_2 4080 -seq ~/opt/pulling/2xov.seq"

    # tt = pd.read_csv("/Users/weilu/Research/server/jun_2018/low_e_jun01_h56.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Research/server/jun_2018/low_e_jun01_h34.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Research/server/jun_2018/low_e_jun01_h12.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Research/server/jun_2018/low_e_jun01_out.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Research/server/jun_2018/low_e_jun01_pre.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Research/server/jun_2018/low_e_jun01_transition.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Research/server/jun_2018/low_e_jun01_post_transition.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Desktop/56_highest_barrier.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Desktop/56.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Desktop/56_pick.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Desktop/56_pick_2.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Desktop/56_pick_3.csv", index_col=0)  # near native
    # tt = pd.read_csv("/Users/weilu/Desktop/56_pick_4.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Desktop/aug_no_perturb.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Desktop/aug_no_perturb_inter.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Desktop/aug_no_perturb_inter_2.csv", index_col=0)
    tt = pd.read_csv("/Users/weilu/Desktop/native.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Desktop/aug_no_perturb_trans.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Desktop/aug_no_perturb_trans_2.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Desktop/aug_no_perturb_path.csv", index_col=0)
    # rerun = 1
    # sample = tt.sample(5).reset_index(drop=True)
    sample = tt.reset_index(drop=True)
    # sample["Frame"] = ((sample["Step"] - 2e7*rerun)/4000).astype("int")
    sample["rerun"] = (sample["Step"] // 2e7).astype(int)
    sample["Frame"] = ((sample["Step"] % 2e7)/4000).astype("int")
    for index, row in sample.iterrows():
        BiasTo = row["BiasTo"]
        Run = row["Run"]
        Frame = row["Frame"]
        rerun = row["rerun"]
        print(BiasTo, Run, Frame)

        location = location_pre + f"/dis_{BiasTo}/{rerun}/dump.lammpstrj.{int(Run)}"
        cmd = cmd_pre + " " + location + f" structure_{index} {int(Frame)} -seq ~/opt/pulling/2xov.seq"
        print(cmd)
        do(cmd)
    pick_structure_generate_show_script(n=len(sample))

if args.day == "aug10":
    name = "abeta42_2"
    rowN = 2
    columnN = 1
    if args.mode == 1:
        do("rm crystal_structure.pdb")
        to = "tmp1.pdb"

        table = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        count = 0
        for i in range(rowN):
            for j in range(columnN):
                duplicate_pdb("one_abeta42.pdb", to, offset_x=i*40.0, offset_y=j*40.0, offset_z=30.0, new_chain=table[count])
                do(f"cat {to} >> crystal_structure.pdb")
                count += 1
    if args.mode == 2:
        do(f"create_project.py {name} --frag --crystal --globular")
        do("cp ~/opt/abeta/ssweight_weihua ssweight")
        # do("cp ~/opt/abeta/zim12 zim")
        do("rm zim")
        for i in range(rowN*columnN):
            do("cat ~/opt/abeta/zim >> zim")
        do("cp ~/opt/abeta/seq.gamma .")
        do("cp ~/opt/abeta/fix_backbone_coeff.data .")
        do(f"cp data.crystal data.{name}")
    if args.mode == 3:
        do("python3 ~/opt/small_script/PBC_fixer.py")
        do(f"movie.py {name} -d new.dump")


if args.day == "aug04":
    if args.mode == 1:
        do("rm crystal_structure.pdb")
        to = "tmp1.pdb"
        table = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        count = 0
        for i in range(2):
            for j in range(3):
                duplicate_pdb("one_abeta42.pdb", to, offset_x=i*40.0, offset_y=j*40.0, offset_z=30.0, new_chain=table[count])
                do(f"cat {to} >> crystal_structure.pdb")
                count += 1
    if args.mode == 2:
        do("create_project.py abeta42_6 --frag --crystal --globular")
        do("cp ~/opt/abeta/ssweight_weihua ssweight")
        do("cp ~/opt/abeta/zim12 zim")
        do("cp ~/opt/abeta/seq.gamma .")
        do("cp ~/opt/abeta/fix_backbone_coeff.data .")
        do("cp data.crystal data.abeta42_6")
    if args.mode == 3:
        do("python3 ~/opt/small_script/PBC_fixer.py")
        do("movie.py abeta42_6 -d new.dump")

if args.day == "jul31":
    if args.mode == 1:
        do("rm crystal_structure.pdb")
        to = "tmp1.pdb"
        table = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        count = 0
        for i in range(4):
            for j in range(3):
                duplicate_pdb("one_abeta42.pdb", to, offset_x=i*40.0, offset_y=j*40.0, offset_z=30.0, new_chain=table[count])
                do(f"cat {to} >> crystal_structure.pdb")
                count += 1
    if args.mode == 2:
        do("create_project.py abeta42_12 --frag --crystal --globular")
        do("cp ~/opt/abeta/ssweight_weihua ssweight")
        do("cp ~/opt/abeta/zim12 zim")
        do("cp ~/opt/abeta/seq.gamma .")
        do("cp ~/opt/abeta/fix_backbone_coeff.data .")
        do("cp data.crystal data.abeta42_12")
    if args.mode == 3:
        do("python3 ~/opt/small_script/PBC_fixer.py")
        do("movie.py abeta42_12 -d new.dump")
if args.day == "jul20":
    if args.mode == 1:
        do("rm crystal_structure.pdb")
        to = "tmp1.pdb"
        table = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        count = 0
        for i in range(4):
            for j in range(3):
                duplicate_pdb("one_abeta42.pdb", to, offset_x=i*40.0, offset_y=j*40.0, offset_z=30.0, new_chain=table[count])
                do(f"cat {to} >> crystal_structure.pdb")
                count += 1
'''
if args.day == "jul11":
    if args.mode == 1:
        do("rm crystal_structure.pdb")
        to = "tmp1.pdb"
        table = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        count = 0
        for i in range(4):
            for j in range(3):
                duplicate_pdb("one_abeta42.pdb", to, offset_x=i*40.0, offset_y=j*40.0, offset_z=60.0, new_chain=table[count])
                do(f"cat {to} >> crystal_structure.pdb")
                count += 1
if args.day == "jul10":
    if args.mode == 1:
        duplicate_pdb("one_abeta42.pdb", "tmp.pdb", offset_x=40.0, offset_z=60.0, new_chain="B")
        do("cp one_abeta42.pdb crystal_structure.pdb")
        do("cat tmp.pdb >> crystal_structure.pdb")
    if args.mode == 2:
        do("create_project.py abeta42_2 --frag --crystal --globular")
    if args.mode == 3:
        duplicate_pdb("one_abeta42.pdb", "crystal_structure.pdb", offset_x=0.0, offset_z=60.0, new_chain="A")
    if args.mode == 4:
        do("create_project.py abeta42 --frag --crystal --globular")


if args.day == "may28":
    cmd_pre = "python2 ~/opt/script/BuildAllAtomsFromLammps.py"
    # location_pre = "/Users/weilu/Research/server/apr_2018/sixth/rg_0.15_lipid_1.0_mem_1_go_0.8/simulation"
    location_pre = "/Users/weilu/Research/server/may_2018/second/simulation"
    # location_pre = "/Users/weilu/Research/server/may_2018/second_long/simulation"
    # cmd = cmd_pre + " " + location + " structure_2 4080 -seq ~/opt/pulling/2xov.seq"

    tt = pd.read_csv("/Users/weilu/Research/server/jun_2018/low_e_jun01_h56.csv", index_col=0)
    tt = pd.read_csv("/Users/weilu/Research/server/jun_2018/low_e_jun01_h34.csv", index_col=0)
    tt = pd.read_csv("/Users/weilu/Research/server/jun_2018/low_e_jun01_h12.csv", index_col=0)
    tt = pd.read_csv("/Users/weilu/Research/server/jun_2018/low_e_jun01_out.csv", index_col=0)
    tt = pd.read_csv("/Users/weilu/Research/server/jun_2018/low_e_jun01_pre.csv", index_col=0)
    tt = pd.read_csv("/Users/weilu/Research/server/jun_2018/low_e_jun01_transition.csv", index_col=0)
    tt = pd.read_csv("/Users/weilu/Research/server/jun_2018/low_e_jun01_post_transition.csv", index_col=0)
    # rerun = 1
    # sample = tt.sample(5).reset_index(drop=True)
    ample = tt.reset_index(drop=True)
    # sample["Frame"] = ((sample["Step"] - 2e7*rerun)/4000).astype("int")
    sample["rerun"] = (sample["Step"] // 2e7).astype(int)
    sample["Frame"] = ((sample["Step"] % 2e7)/4000).astype("int")
    for index, row in sample.iterrows():
        BiasTo = row["BiasTo"]
        Run = row["Run"]
        Frame = row["Frame"]
        rerun = row["rerun"]
        print(BiasTo, Run, Frame)

        location = location_pre + f"/dis_{BiasTo}/{rerun}/dump.lammpstrj.{int(Run)}"
        cmd = cmd_pre + " " + location + f" structure_{index} {int(Frame)} -seq ~/opt/pulling/2xov.seq"
        print(cmd)
        do(cmd)
    pick_structure_generate_show_script(n=len(sample))
if args.day == "may22":
    cmd_pre = "python2 ~/opt/script/BuildAllAtomsFromLammps.py"
    # location_pre = "/Users/weilu/Research/server/apr_2018/sixth/rg_0.15_lipid_1.0_mem_1_go_0.8/simulation"
    location_pre = "/Users/weilu/Research/server/may_2018/second/simulation"
    # cmd = cmd_pre + " " + location + " structure_2 4080 -seq ~/opt/pulling/2xov.seq"
    # tt = pd.read_csv("/Users/weilu/Research/server/barrier.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Research/server/high_go.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Research/server/rerun3.csv", index_col=0)
    tt = pd.read_csv("/Users/weilu/Research/data/t373_narrow.csv", index_col=0)
    tt = pd.read_csv("/Users/weilu/Research/data/t373_super_narrow.csv", index_col=0)
    tt = pd.read_csv("/Users/weilu/Research/data/selected.csv", index_col=0)
    tt = pd.read_csv("/Users/weilu/Research/data/selected2.csv", index_col=0)
    tt = pd.read_csv("/Users/weilu/Research/data/selected_all.csv", index_col=0)
    tt = pd.read_csv("/Users/weilu/Research/data/constrain_qw.csv", index_col=0)
    tt = pd.read_csv("/Users/weilu/Research/data/constrain_qw_temp.csv", index_col=0)
    # rerun = 1
    sample = tt.sample(5).reset_index(drop=True)
    # sample["Frame"] = ((sample["Step"] - 2e7*rerun)/4000).astype("int")
    sample["rerun"] = (sample["Step"] // 2e7).astype(int)
    sample["Frame"] = ((sample["Step"] % 2e7)/4000).astype("int")
    for index, row in sample.iterrows():
        BiasTo = row["BiasTo"]
        Run = row["Run"]
        Frame = row["Frame"]
        rerun = row["rerun"]
        print(BiasTo, Run, Frame)

        location = location_pre + f"/dis_{BiasTo}/{rerun}/dump.lammpstrj.{int(Run)}"
        cmd = cmd_pre + " " + location + f" structure_{index} {int(Frame)} -seq ~/opt/pulling/2xov.seq"
        print(cmd)
        do(cmd)
    pick_structure_generate_show_script(n=len(sample))
if args.day == "may13":
    protein_list = ["T0949", "T0950"]
    for protein in protein_list:
        for i in range(1,4):
            do(f"cp -r {protein}/post-processing/model.{i}/lowTstructure/ ~/Dropbox/IAAWSEM/casp13/selection/{protein}_model_{i}")
if args.day == "apr20":
    cmd_pre = "python2 ~/opt/script/BuildAllAtomsFromLammps.py"
    location_pre = "/Users/weilu/Research/server/apr_2018/sixth/rg_0.15_lipid_1.0_mem_1_go_0.8/simulation"
    # cmd = cmd_pre + " " + location + " structure_2 4080 -seq ~/opt/pulling/2xov.seq"
    # tt = pd.read_csv("/Users/weilu/Research/server/barrier.csv", index_col=0)
    # tt = pd.read_csv("/Users/weilu/Research/server/high_go.csv", index_col=0)
    tt = pd.read_csv("/Users/weilu/Research/server/rerun3.csv", index_col=0)
    rerun = 1
    sample = tt.sample(5).reset_index(drop=True)
    sample["Frame"] = ((sample["Step"] - 2e7*rerun)/4000).astype("int")
    for index, row in sample.iterrows():
        BiasTo = row["BiasTo"]
        Run = row["Run"]
        Frame = row["Frame"]
        print(BiasTo, Run, Frame)

        location = location_pre + f"/dis_{BiasTo}/{rerun}/dump.lammpstrj.{int(Run)}"
        cmd = cmd_pre + " " + location + f" structure_{index} {int(Frame)} -seq ~/opt/pulling/2xov.seq"
        print(cmd)
        do(cmd)
    pick_structure_generate_show_script(n=len(sample))

if args.day == "mar26":
    cmd_pre = "python2 ~/opt/script/BuildAllAtomsFromLammps.py"
    location_pre = "/Users/weilu/Research/server/mar_2018/eighth/force_0.03_rg_0.15_lipid_1.0_mem_1_go_0.8/simulation"
    # cmd = cmd_pre + " " + location + " structure_2 4080 -seq ~/opt/pulling/2xov.seq"
    tt = pd.read_csv("/Users/weilu/Research/server/mar_2018/05_week/pick_structure/highE.csv", index_col=0)
    rerun = 0
    sample = tt.sample(5).reset_index(drop=True)
    sample["Frame"] = ((sample["Step"] - 2e7*rerun)/4000).astype("int")
    for index, row in sample.iterrows():
        BiasTo = row["BiasTo"]
        Run = row["Run"]
        Frame = row["Frame"]
        print(BiasTo, Run, Frame)

        location = location_pre + f"/dis_{BiasTo}/{rerun}/dump.lammpstrj.{int(Run)}"
        cmd = cmd_pre + " " + location + f" structure_{index} {int(Frame)} -seq ~/opt/pulling/2xov.seq"
        print(cmd)
        do(cmd)
    pick_structure_generate_show_script(n=len(sample))
if args.day == "mar04":
    if args.mode == 1:
        dic = {"T0784":"3u6g", "T0792":"3rcoa", "T0815":"3fh1", "T0833":"4r8o", "T251":"1j0f", "T0803":"4pe2"}
        protein = "T0803"
        model = dic[protein] + "_part_enhanced"
        do(f"sed '/TER/d' {model}.pdb > {model}_remove_ter.pdb")
        do(f"python2 ~/opt/small_script/CalcQValueFromTwoPdb_2.py {protein}.pdb {model}_remove_ter.pdb | tail -n 1 > qw")
        do(f"~/opt/TMalign/TMscore {model}_remove_ter.pdb {protein}.pdb | grep 'Structure1' -A8 > RMSD_and_GDT")
if args.day == "feb15":
    if args.mode == 1:
        with open("metadata", "w") as f:
            dis_list = glob.glob("/Users/weilu/Research/server/feb_2018/week_of_feb05/rg_0.1_lipid_1.0_mem_1/simulation/dis_*")
            for dis in dis_list:
                for idx in range(12):
                    fileName = dis + f"/1/data_{idx}.dat {idx}\n"
                    f.write(fileName)
    if args.mode == 2:
        dis_list = glob.glob("dis_*")
        for dis in dis_list:
            cd(dis)
            cd("1")
            for i in range(12):
                do("awk '{print $2}'" + f" wham.{i}.dat |  sed 's/,$//' > qw.{i}.dat")
                do(f"paste qw.{i}.dat z_{i}.dat > data_{i}.dat")
            cd("../..")
    if args.mode == 3:
        with open("metadata", "w") as f:
            dis_list = glob.glob("/Users/weilu/Research/server/feb_2018/week_of_feb05/rg_0.1_lipid_1.0_mem_1/simulation/dis_*")
            dis_list = ["/Users/weilu/Research/server/feb_2018/week_of_feb05/rg_0.1_lipid_1.0_mem_1/simulation/dis_60.0"]
            for dis in dis_list:
                for idx in range(12):
                    fileName = dis + f"/1/data_{idx}.dat {idx}\n"
                    f.write(fileName)
if args.day == "feb14":
    if(args.mode == 1):
        print("add pdb")
        for i in range(0, 20):
            # print(os.getcwd())
            cd(str(i))
            do("cp 2xov.pdb 0/")
            cd("..")
    if(args.mode == 2):
        print("Extract qw and distance info.")
        for i in range(20):
            cd(str(i))
            cd("0")
            do("awk '{print $2}' wham.dat |  sed 's/,$//' | sed 1d > qw.dat")
            do("awk '{print $2}' addforce.dat |  sed 's/,$//' | sed 1d > distance.dat")
            cd("../..")

    # do("python2 ~/opt/small_script/CalcLocalDistanceStats.py 2xov directory_list out")
    if(args.mode == 3):
        print("create directory_list")
        with open("directory_list", "w") as f:
            for i in range(0, 20):
                # print(os.getcwd())
                location = os.getcwd() + "/../"
                f.write(location+str(i)+"/0\n")

    if(args.mode == 4):
        print("create directory_list")
        for i in range(0, 20):
            with open(str(i), "w") as f:
                # print(os.getcwd())
                location = os.getcwd() + "/../"
                f.write(location+str(i)+"/0\n")
            do("python2 ~/opt/small_script/CalcLocalDistanceStats.py 2xov {0} out_{0}".format(i))

if args.day == "jan25":
    if args.mode == 1:
        compute_average_z(f"dump.lammpstrj", f"z.dat")
if args.day == "jan19":
    if args.mode == 1:
        read_data_Pulling("pressure_0.5_rg_0.2_mem_1")
    if args.mode == 2:
        read_data_Pulling("pressure_1")

if args.day == "jan17":
    if args.mode == 1:
        read_data_Pulling("pressure_0.5_rg_0.2_mem_1")
    if args.mode == 2:
        read_data_Pulling("pressure_0.5")

if args.day == "jan12":
    if args.mode == 1:
        read_data_Pulling()

if args.day == "jan05":
    if args.mode == 1:
        read_data_Pulling()
if args.day == "nov21":
    if args.mode == 1:
        simulation_list = ["memb_3_rg_0.1_lipid_1_extended"]
        for simulation in simulation_list:
            for mode in range(3):
                cd(f"nov_15_all_freeEnergy_calculation_sample_range_mode_{mode}")
                cd(simulation)
                do(f"read_pmf.py -m 1 -l {simulation}_{mode}")
if args.day == "nov15":
    if args.mode == 1:
        simulation_list = ["memb_3_rg_0.1_lipid_1_extended"]
        for simulation in simulation_list:
            for mode in range(3):
                cd(f"nov_15_all_freeEnergy_calculation_sample_range_mode_{mode}")
                cd(simulation)
                do(f"read_pmf.py -m 1 -l {simulation}_{mode}")
'''
if args.day == "old":
    #  convert \( 0.0.png 0.2.png 0.4.png  +append \) \( 0.6.png  0.8.png    1.0.png +append \) -background none -append final.png
    if(args.test):
        check_and_correct_fragment_memory()

    if(args.mode == 10):
        source_img = "final.untitled.*.jpg"
        # do('convert ' + source_img + ' -pointsize 14 -draw "fill black text 1,11 \'some text\' " ' + "test.gif")
        for i in range(1,1301, 1):
            do("convert final.2xov.{:05d}.jpg -pointsize 50 -annotate +800+100 'Current force: {:5.2f}pN' -gravity West {:05d}.jpg".format(i, 69.7*(i*4e-4), i))
    if(args.mode == 9):
        print("create directory_list")
        with open("directory_list", "w") as f:
            force_list = [0.3, 0.35, 0.4]
            for force in force_list:
                for i in range(0, 20):
                    # print(os.getcwd())
                    location = os.getcwd() + "/../force_{}_/simulation/".format(force)
                    f.write(location+str(i)+"/0\n")

    if(args.mode == 8):
        print("mode: {}".format(args.mode))
        compute_theta_for_each_helix()
    if(args.mode == 7):
        # n_list = [10, 16, 20, 21, 31, 33, 4, 45, 53, 55, 56, 6, 60, 7, 74, 75, 80, 90, 92]
        n_list = [32, 41, 47, 48, 57, 58, 6, 63, 69, 72, 82, 96]
        n_list = [0,1,10,11,12,13,14,15,16,17,18,19,2]
        for i,n in enumerate(n_list):
            cd("{}/0".format(n))
            do("python3 ~/opt/small_script/last_n_frame.py -n 1 2xov.")
            cd("../..")
        cd("..")
        cd("end_frame_that_folded")
        for i,n in enumerate(n_list):
            do("cp ../simulation/{0}/0/frames/*.pdb {0}.pdb".format(n))
    if(args.mode == 6):
        n_list = [16, 10, 13, 16,15, 3, 8, 11, 6, 7, 20, 13, 18, 9,6, 19,4,14, 11,14]
        print(len(n_list))
        for i,n in enumerate(n_list):
            do("cp run{0}/frames/{1}.pdb selection/run{0}.pdb".format(i+1,n-1))

    if(args.mode == 5):
        n = 20
        for i in range(1, n+1):
            cd("run{}".format(i))
            do("python3 ~/opt/small_script/last_n_frame.py -n 20 T0815. -m 2")
            cd("frames")
            do("python3 ~/opt/small_script/cross_q.py -m 3")
            do("cp ~/opt/small_script/heatmap_script.m .")
            cd("../..")
    if(args.mode == 1):
        print("create frustration_censored_contacts.dat")
        with open("frustration_censored_contacts.dat", "w") as f:
            for i in range(1,182):
                # f.write("65 {}\n".format(i))
                f.write("116 {}\n".format(i))

    if(args.mode == 2):
        print("Extract qw and distance info.")
        for i in range(20):
            cd(str(i))
            cd("0")
            do("awk '{print $2}' wham.dat |  sed 's/,$//' | sed 1d > qw.dat")
            do("awk '{print $2}' addforce.dat |  sed 's/,$//' | sed 1d > distance.dat")
            cd("../..")

    # do("python2 ~/opt/small_script/CalcLocalDistanceStats.py 2xov directory_list out")
    if(args.mode == 3):
        print("create directory_list")
        with open("directory_list", "w") as f:
            for i in range(0, 20):
                # print(os.getcwd())
                location = os.getcwd() + "/../"
                f.write(location+str(i)+"/0\n")

    if(args.mode == 4):
        print("create directory_list")
        for i in range(0, 20):
            with open(str(i), "w") as f:
                # print(os.getcwd())
                location = os.getcwd() + "/../"
                f.write(location+str(i)+"/0\n")
            do("python2 ~/opt/small_script/CalcLocalDistanceStats.py 2xov {0} out_{0}".format(i))


# if(args.test):
#     if(args.mode == 3):
#         do("cp ../2xov.pdb .")
#         do("python2 ~/opt/script/CalcLocalQTrajectory.py 2xov dump.lammpstrj localQ_trajectory")
#     if(args.mode == 0):
#         run = 4
#         cd(str(run))
#         cd("0")
#         do("movie.py 2xov")
#         do("/Applications/VMD\ 1.9.3.app/Contents/MacOS/startup.command -e 2xov_movie_bicelle.tcl")
#     if(args.mode == 1):
#         folder = "memb_0_force_ramp_rg_0"
#         cd(folder)
#         cd("0")
#         do("tail -n 4 addforce.dat")
#         do("tail wham.dat")
#         do("movie.py 2xov")
#         do("/Applications/VMD\ 1.9.3.app/Contents/MacOS/startup.command -e 2xov_movie_bicelle.tcl")
#     if(args.mode == 2):
#         files = glob.glob("memb_*")
#         # print(files)
#         for folder in files:
#             print(folder)
#             cd(folder)
#             cd("0")
#             do("movie.py 2xov")
#             do("/Applications/VMD\ 1.9.3.app/Contents/MacOS/startup.command -e 2xov_movie_bicelle.tcl")
#             cd("../..")
#     # for i in range(0, 20):
#     #     cd(str(i))
#     #     do("python3 ~/opt/aawsem_show.py --casp -m 2 T0782.")
#     #     cd("..")
    def test():
        print("don't show me")
        # force_list = [round(i*0.1,2) for i in range(20)]
        # for force in force_list:
        #     cd("{}".format(force))
        #     do("plotcontour.py pmf-300.dat")
        #     do("cp test.png ../result/{}.png".format(force))
        #     cd("..")

        with open("metadata", "w") as f:
            temp_list = [325, 350]
            i_range = range(40)
            for temp in temp_list:
                for i in i_range:
                    f.write("/Users/weilu/Research/server/project/freeEnergy_2xov/qValue_v4/simulation/{}/{}/halfdata\n".format(temp, i))

        # n = args.number
        # folder_list = sorted(glob.glob("wham*"))
        # for folder in folder_list:
        #     print(folder)
        #     os.system("tail -n+2 {}/cv-300-400-10.dat | sort -r -k 2 | head -n1".format(folder))

                # for i in range(2):
                #     address = folder + "/simulation/" + str(i)
                #     f.write(address+"  \n")

def calQnQc():
    qn_start = 0
    qn_end = 80
    qc_start = 80
    qc_end = 181
    qc2_start = 130
    qc2_end = 181
    os.system("python2 ~/opt/CalcQnQc.py 2xov.pdb dump.lammpstrj {} 0.15 {} {}".format("qn", qn_start, qn_end))
    os.system("python2 ~/opt/CalcQnQc.py 2xov.pdb dump.lammpstrj {} 0.15 {} {}".format("qc", qc_start, qc_end))
    os.system("python2 ~/opt/CalcQnQc.py 2xov.pdb dump.lammpstrj {} 0.15 {} {}".format("qc2", qc2_start, qc2_end))
    size1 = file_len("qn")
    size2 = file_len("qc")
    size3 = file_len("qc2")
    # os.system("paste qn qc > qnqc")
    # if(size1 < 400 or size2 < 400 or size3 < 400):
    #     raise ValueError('file length too small')
    # os.system("head -n 4000 qn > qn_all")
    # os.system("head -n 4000 qc > qc_all")
    # os.system("head -n 4000 qc2 > qc2_all")
    # os.system("tail -n 2000 qn_all > qn_half")
    # os.system("tail -n 2000 qc_all > qc_half")
    # os.system("tail -n 2000 qc2_all > qc2_half")
    # os.system("paste qn_half qc_half qc2_half ")
if(args.qnqc):
    calQnQc()


def calQo():
    os.system("python2 ~/opt/script/CalcQValue_multi.py 2lhc dump.lammpstrj qo 1")


def cpull():
    os.system("cp ~/opt/pulling/cfx.gp .")
    os.system("gnuplot cfx.gp")
    os.system("open cf_extension.*")
if(args.cpull):
    cpull()


def pull():
    os.system("cp ~/opt/pulling/fx.gp .")
    os.system("gnuplot fx.gp")
    os.system("open f_extension.pdf")
if(args.pull):
    pull()


def energy():
    in_file_name = "qw.dat"
    os.system("cp ~/opt/gagb/temp.gp .")
    out_file = "test.pdf"
    out = "'out_file_name=\"{}\"'".format(out_file)
    in_file = "gb_sequence/wham.dat_2"
    gp_in = "'in_file_name=\"{}\"'".format(in_file)
    # in_file_2 = "0/wham.dat"
    # gp_in_2 = "'in_file_name_2=\"{}\"'".format(in_file_2)
    os.system("gnuplot -e {} -e {} -e {} temp.gp ".format(gp_in, gp_in_2, out))
    os.system("open test.pdf ")
# def wham_analysis():
#     os.system("mkdir -p wham")
#     os.system("rm wham/all.dat")
#     os.system("rm wham/*_total")
#     os.chdir("300")
#     for i in range(20):
#         # os.system("cat {}/halfdata.dat >> ../wham/all.dat".format(i))
#         # os.system("tail -n+3 rerun_"+str(i)+"/wham.dat | awk '{print $3}' | tail -n 5000 >> ../wham/p_total")
#         os.system("tail -n+2 " + str(i) + "/wham.dat | awk '{print $3}' | tail -n 5000 >> ../wham/p_total")
#         os.system("tail -n+2 " + str(i) + "/wham.dat | awk '{print $5}' | tail -n 5000 >> ../wham/e_total")
#         os.system("tail -n+2 " + str(i) + "/wham.dat | awk '{print $2}' | tail -n 5000 >> ../wham/qw_total")
#
#         os.chdir(str(i))
#         os.system("cp ~/opt/gagb/2lhc_part.pdb .")
#         os.system("awk '{if ((!((($1>0 && $1<25) || ($1>159 && $1<200) ) && $3>-10)  ) ) print }' dump.lammpstrj > data_test")
#         os.system("python2 ~/opt/script/CalcQValue.py 2lhc_part.pdb data_test test")
#         os.system("tail -n +2 test > qw_ga.dat")
#         os.chdir("..")
#         os.system("tail -n 5000 " + str(i) + "/qw_ga.dat >> ../wham/qw_ga_total")
#         # os.system("tail -n+3 rerun_{}/wham.dat | awk '{print $3}' >> ../wham/qo_total".format(str(i)))
#     os.chdir("../wham")
#     # os.system("awk '{print $2}' all.dat > Qw_total")
#     # os.system("awk '{print $3}' all.dat > Qgb_total")
#     # os.system("awk '{print $1}' all.dat > qwa_total")
#     # os.system("awk '{print $4}' all.dat > e_total")
#     # os.system("awk '{print $6}' all.dat > p_total")
#     os.system("cp ~/opt/wham_analysis/*.m .")
#     os.chdir("..")
#     os.system("~/opt/script/wham/fused_calc_cv.sc wham/ 2lhd 20 300 250 350 10 50 200 0.05 1")
# protein_name = args.template.split('_', 1)[-1].strip('/')
# protein_name = args.protein.strip('/')
    # name = "ga_2m"
# exec(open("config.py").read())
# n = number_of_run
# steps = simulation_steps
# # protein_name = protein_name
#
# import experiment_analysis
# os.system("cp ~/opt/gagb/energy2.gp .")
# os.system("gnuplot energy2.gp ")
# os.system("open energy2.pdf ")


def wham_analysis():
    os.system("mkdir -p wham")
    os.system("rm wham/all.dat")
    os.system("rm wham/*_total")
    os.chdir("350")
    for i in range(50):

        # os.system("cat {}/halfdata.dat >> ../wham/all.dat".format(i))
        # os.system("tail -n+3 rerun_"+str(i)+"/wham.dat | awk '{print $3}' | tail -n 5000 >> ../wham/p_total")
        os.system("tail -n+2 " + str(i) + "/wham.dat | awk '{print $4}' | tail -n 2000 >> ../wham/p_total")
        os.system("tail -n+2 " + str(i) + "/wham.dat | awk '{print $6}' | tail -n 2000 >> ../wham/e_total")
        os.system("tail -n+2 " + str(i) + "/wham.dat | awk '{print $2}' | tail -n 2000 >> ../wham/qw_total")
        os.system("tail -n 2000 " + str(i) + "/energy.log | awk '{print $18-$13}' >> ../wham/e_total")
        os.chdir(str(i))
        os.system("cp ~/opt/gagb/2lhc_part.pdb .")
        os.system("awk '{if ((!((($1>0 && $1<25) || ($1>159 && $1<200) ) && $3>-10)  ) ) print }' dump.lammpstrj > data_test")
        os.system("python2 ~/opt/script/CalcQValue.py 2lhc_part.pdb data_test test")
        os.system("tail -n +2 test > qw_ga.dat")
        os.system("python2 ~/opt/script/CalcQValue_multi.py 2LHC dump.lammpstrj ga_qo.dat 1")

        os.chdir("..")

        os.system("tail -n 2000 " + str(i) + "/ga_qo.dat >> ../wham/qo_ga_total")
        os.system("tail -n 2000 " + str(i) + "/qw_ga.dat >> ../wham/qw_ga_total")
        # os.system("tail -n+3 rerun_{}/wham.dat | awk '{print $3}' >> ../wham/qo_total".format(str(i)))
    os.chdir("../wham")
    # os.system("awk '{print $2}' all.dat > Qw_total")
    # os.system("awk '{print $3}' all.dat > Qgb_total")
    # os.system("awk '{print $1}' all.dat > qwa_total")
    # os.system("awk '{print $4}' all.dat > e_total")
    # os.system("awk '{print $6}' all.dat > p_total")
    os.system("cp ~/opt/wham_analysis/*.m .")
    os.chdir("..")
    os.system("~/opt/script/wham/fused_calc_cv.sc wham/ 2lhd 50 350 300 400 10 60 100 0.02 1")


def wham_analysis400():
    os.system("mkdir -p wham400")
    os.system("rm wham400/all.dat")
    os.system("rm wham400/*_total")
    os.chdir("400")
    for i in range(18):
        os.system("cat {}/halfdata.dat >> ../wham400/all.dat".format(i))
        os.system("tail -n+3 rerun_"+str(i)+"/wham.dat | awk '{print $3}' | tail -n 5000 >> ../wham400/p_total")

    os.chdir("../wham400")
    os.system("awk '{print $2}' all.dat > Qw_total")
    os.system("awk '{print $3}' all.dat > Qgb_total")
    os.system("awk '{print $1}' all.dat > qwa_total")
    os.system("awk '{print $4}' all.dat > e_total")
    os.system("cp ~/opt/wham_analysis/*.m .")
    os.chdir("..")
    os.system("~/opt/script/wham/fused_calc_cv.sc wham400/ 2lhd 18 400 350 450 10 50 200 0.05 0.9")


# def free_energy_analysis():
#     temp_list = [350]
#     n = 50
#     for temp in temp_list:
#         os.chdir(str(temp))
#         for i in range(n):
#             os.chdir(str(i))
#             os.system("awk '{print $6}' wham.dat | tail -n +2 > e.dat")
#             # os.system("awk '{print $3}' wham.dat | tail -n +2 > p2.dat")
#             os.system("cp ~/opt/gagb/*.pdb .")
#             # os.system("awk '{if ((!((($1>0 && $1<25) || ($1>159 && $1<200) ) && $3>-10)  ) ) print }' dump.lammpstrj > data_test")
#             # os.system("python2 ~/opt/script/CalcQValue.py 2lhc_part.pdb data_test test")
#             # os.system("tail -n +2 test > q_ga_part.dat")
#             # os.system("python2 ~/opt/script/CalcQValue.py 2lhc dump.lammpstrj q_ga_included.dat")
#             # os.system("tail -n +2 q_ga_included.dat > q_ga.dat")
#             # os.system("python2 ~/opt/script/CalcQValue.py 2lhd dump.lammpstrj q_gb_included.dat")
#             # os.system("tail -n +2 q_gb_included.dat > q_gb.dat")
#             os.system("cp ~/opt/gagb/nativecoords_g* .")
#             os.system("mv nativecoords_ga.dat nativecoords.dat")
#             os.system("python2 ~/opt/script/CalcQValue_multi.py 2lhc dump.lammpstrj qo_ga 1")
#             os.system("tail -n +2 qo_ga > qo_ga.dat")
#             os.system("mv nativecoords_gb.dat nativecoords.dat")
#             os.system("python2 ~/opt/script/CalcQValue_multi.py 2lhd dump.lammpstrj qo_gb 1")
#             os.system("tail -n +2 qo_gb > qo_gb.dat")
#             os.system("paste qo_ga.dat  qo_gb.dat e.dat > data.dat")
#             os.system("tail -n 2000 data.dat > halfdata.dat")
#             os.chdir("..")
#         os.chdir("..")
def free_energy_analysis():
    temp_list = [350]
    n = 40

    for temp in temp_list:
        os.chdir(str(temp))
        os.system("mkdir -p wham")
        os.system("rm wham/all.dat")
        os.system("mkdir -p wham_half")
        os.system("rm wham_half/all.dat")
        for i in range(n):
            os.chdir(str(i))
            os.system("tail -n +2 wham.dat > data.dat")
            os.system("tail -n 4000 data.dat > halfdata.dat")
            # os.system("cat data.dat >> ../wham/all.dat")
            # os.system("cat halfdata.dat >> ../wham_half/all.dat")
            # os.system("awk '{print $6}' wham.dat | tail -n +2 > e.dat")
            # # os.system("awk '{print $3}' wham.dat | tail -n +2 > p2.dat")
            # os.system("cp ~/opt/gagb/*.pdb .")
            # # os.system("awk '{if ((!((($1>0 && $1<25) || ($1>159 && $1<200) ) && $3>-10)  ) ) print }' dump.lammpstrj > data_test")
            # # os.system("python2 ~/opt/script/CalcQValue.py 2lhc_part.pdb data_test test")
            # # os.system("tail -n +2 test > q_ga_part.dat")
            # # os.system("python2 ~/opt/script/CalcQValue.py 2lhc dump.lammpstrj q_ga_included.dat")
            # # os.system("tail -n +2 q_ga_included.dat > q_ga.dat")
            # # os.system("python2 ~/opt/script/CalcQValue.py 2lhd dump.lammpstrj q_gb_included.dat")
            # # os.system("tail -n +2 q_gb_included.dat > q_gb.dat")
            # os.system("cp ~/opt/gagb/nativecoords_g* .")
            # os.system("mv nativecoords_ga.dat nativecoords.dat")
            # os.system("python2 ~/opt/script/CalcQValue_multi.py 2lhc dump.lammpstrj qo_ga 1")
            # os.system("tail -n +2 qo_ga > qo_ga.dat")
            # os.system("mv nativecoords_gb.dat nativecoords.dat")
            # os.system("python2 ~/opt/script/CalcQValue_multi.py 2lhd dump.lammpstrj qo_gb 1")
            # os.system("tail -n +2 qo_gb > qo_gb.dat")
            # os.system("paste qo_ga.dat  qo_gb.dat e.dat > data.dat")
            # os.system("tail -n 2000 data.dat > halfdata.dat")
            os.chdir("..")

        # os.chdir("wham")
        # os.system("awk '{print $2}' all.dat > qo_ga_total")
        # os.system("awk '{print $3}' all.dat > qo_gb_total")
        # os.system("awk '{print $3}' all.dat > p_total")
        # os.system("awk '{print $7}' all.dat > e_total")
        # os.system("cp ~/opt/wham_analysis/*.m .")
        # os.chdir("..")
        # os.system("~/opt/script/wham/fused_calc_cv.sc wham 2lhc 40 350 300 400 10 20 56 0.1 0.9")
        #
        # os.chdir("wham_half")
        # os.system("awk '{print $2}' all.dat > qo_ga_total")
        # os.system("awk '{print $3}' all.dat > qo_gb_total")
        # os.system("awk '{print $3}' all.dat > p_total")
        # os.system("awk '{print $7}' all.dat > e_total")
        # os.system("cp ~/opt/wham_analysis/*.m .")
        # os.chdir("..")
        # os.system("~/opt/script/wham/fused_calc_cv.sc wham_half 2lhc 40 350 300 400 10 20 56 0.1 0.9")

        os.chdir("..")

if(args.freeEnergy):
    free_energy_analysis()


def plot():
    print("Plotting")
    # os.system("python2 ~/opt/script/CalcQValue.py 2lhc dump2.lammpstrj qw.dat")
    os.system("cp ~/opt/temp.gp .")
    out_file = "test.pdf"
    out = "'out_file_name=\"{}\"'".format(out_file)
    in_file = "gb_sequence/wham.dat_2"
    gp_in = "'in_file_name=\"{}\"'".format(in_file)
    in_file_2 = "0/wham.dat"
    gp_in_2 = "'in_file_name_2=\"{}\"'".format(in_file_2)
    os.system("gnuplot -e {} -e {} -e {} temp.gp ".format(gp_in, gp_in_2, out))
    os.system("open test.pdf ")


def fix():
    n = 20
    # for i range(n):
    #     os.system("tail -n 2000 energy.log | awk '{print $18-$13}' >> ../wham")
    # os.chdir("analysis")
    # os.system("rm highest_q_gb")
    # os.system("rm highest_q")
    # for i in range(n):
    #     os.chdir(str(i))
    #     # os.system("tail -n 5000 wham.dat > halfdata.dat")
    #     os.system("cp ~/opt/gagb/2lhc.pdb .")
    #     os.system("cp ~/opt/gagb/2lhd.pdb .")
    #     os.system("python2 ~/opt/script/CalcQValue.py 2lhc.pdb dump.lammpstrj ga")
    #     os.system("tail -n 1000 ga | sort | tail -n 1 > ga_highest")
    #     os.system("cat ga_highest >> ../highest_q")
    #
    #     os.system("python2 ~/opt/script/CalcQValue.py 2lhd dump.lammpstrj gb")
    #     os.system("tail -n 1000 gb | sort | tail -n 1 > gb_highest")
    #     os.system("cat gb_highest >> ../highest_q_gb")
    #     os.chdir("..")


def rerun():
    n = 20
    for i in range(n):
        os.system("cp -r {0} rerun_{0}".format(str(i)))
        os.chdir("rerun_"+str(i))
        os.system("cp ~/opt/gagb/rerun.slurm .")
        os.system("cp ~/opt/gagb/2lhc_rerun.in .")
        os.system("sbatch rerun.slurm")
        os.chdir("..")


def rerun_wham_analysis():
    os.system("mkdir -p wham")
    os.system("rm wham/all.dat")
    os.system("rm wham/*_total")
    os.chdir("300")
    for i in range(20):
        # os.system("cat {}/halfdata.dat >> ../wham/all.dat".format(i))
        # os.system("tail -n+3 rerun_"+str(i)+"/wham.dat | awk '{print $3}' | tail -n 5000 >> ../wham/p_total")
        os.system("tail -n+2 rerun_" + str(i) + "/wham.dat | awk '{print $3}' | tail -n 2000 >> ../wham/p_total")
        os.system("tail -n+2 rerun_" + str(i) + "/wham.dat | awk '{print $6}' | tail -n 2000 >> ../wham/e_total")
        os.system("tail -n+2 rerun_" + str(i) + "/wham.dat | awk '{print $2}' | tail -n 2000 >> ../wham/qw_total")

        # os.chdir(str(i))
        # os.system("cp ~/opt/gagb/2lhc_part.pdb .")
        # os.system("awk '{if ((!((($1>0 && $1<25) || ($1>159 && $1<200) ) && $3>-10)  ) ) print }' dump.lammpstrj > data_test")
        # os.system("python2 ~/opt/script/CalcQValue.py 2lhc_part.pdb data_test test")
        # os.system("tail -n +2 test > qw_ga.dat")
        # os.chdir("..")
        os.system("tail -n 2000 " + str(i) + "/qw_ga.dat >> ../wham/qw_ga_total")
        # os.system("tail -n+3 rerun_{}/wham.dat | awk '{print $3}' >> ../wham/qo_total".format(str(i)))
    os.chdir("../wham")
    # os.system("awk '{print $2}' all.dat > Qw_total")
    # os.system("awk '{print $3}' all.dat > Qgb_total")
    # os.system("awk '{print $1}' all.dat > qwa_total")
    # os.system("awk '{print $4}' all.dat > e_total")
    # os.system("awk '{print $6}' all.dat > p_total")
    os.system("cp ~/opt/wham_analysis/*.m .")
    os.chdir("..")
    os.system("~/opt/script/wham/fused_calc_cv.sc wham/ 2lhd 20 300 250 350 10 50 200 0.05 1")
#rerun()
# rerun_wham_analysis()


#rerun()

if(args.wham):
    wham_analysis()
if(args.wham400):
    wham_analysis400()

if(args.fix):
    fix()
if(args.plot):
    plot()
## -------------Pulling--------
# os.system("cp ~/opt/small_script/springForce.plt .")
# os.system("cp ~/opt/small_script/springForce_smooth.plt .")
# os.system("gnuplot springForce.plt")
# os.system("gnuplot springForce_smooth.plt")
# os.system("open springForce.pdf")
# os.system("open springForce_smooth.pdf")
# SpringConstant_list = [0.0001, 0.00001, 0.000001, 0.0000001]
# for SpringConstant in SpringConstant_list:
#     name = "spring"+str(SpringConstant)
#     os.system("mkdir "+name)
#     os.chdir(name)
#     os.system("cp -r ../2xov/ .")
#     os.system("cp ../variables.dat .")
#     os.chdir("2xov")
#     os.system(  # replace SIMULATION_STEPS with specific steps
#         "sed -i.bak 's/SpringForce/'" +
#         str(SpringConstant) +
#         "'/g' "+protein_name+".in")
#     os.chdir("..")
#     os.system("run.py 2xov/ -s 2")
#     os.chdir("..")
# # -----------------GAGB------------------------------
# # number_of_run_list = [2,  4, 8, 16]
# number_of_run_list = [5, 8, 32]
# for n in number_of_run_list:
#     name = "ga_"+str(n)+"m"
#     # os.system("mkdir "+name)
#     os.system("cp -r 2lhd.pdb "+name)
#
#     # os.system("cp -r 2lhc variables.dat "+name)
#     os.chdir(name)
#     for i in range(20):
#         os.chdir("analysis/"+str(i))
#         os.system("cp ../../2lhd.pdb .")
#         os.system("python2 ~/opt/script/CalcQValue.py 2lhd.pdb dump.lammpstrj q_gb.dat")
#         os.system("python2 ~/opt/script/CalcQValue.py 2lhc.pdb dump.lammpstrj q_ga.dat")
#         os.system("cp ~/opt/small_script/qw_gagb.plt .")
#         os.system("gnuplot qw_gagb.plt")
#         os.system("mv qw_gagb.pdf ../../results/qw_gagb_{0}.pdf".format(str(i)))
#         os.chdir("../..")
#     os.chdir("..")
#
# for n in number_of_run_list:
#     name = "ga_"+str(n)+"m"
#     # os.system("mkdir "+name)
#     os.system("cp -r 2lhd.pdb "+name)
#
#     # os.system("cp -r 2lhc variables.dat "+name)
#     os.chdir(name)
#     for i in range(20):
#         os.chdir("analysis/"+str(i))
#         os.system("paste q_ga.dat q_gb.dat > q_gagb.dat")
#         os.system("cp ~/opt/small_script/qw_ga-gb.plt .")
#         os.system("gnuplot qw_ga-gb.plt")
#         os.system("mv qw_ga-gb.pdf ../../results/qw_ga-gb_{0}.pdf".format(str(i)))
#         os.chdir("../..")
#     os.system("cp ~/opt/small_script/qw_ga_all.plt .")
#     os.system("gnuplot qw_ga_all.plt")
#     os.system("cp ~/opt/small_script/qw_gb_all.plt .")
#     os.system("gnuplot qw_gb_all.plt")
#     os.system("cp ~/opt/small_script/qw_diff_all.plt .")
#     os.system("gnuplot qw_diff_all.plt")
#     os.chdir("..")
# # -----------------GAGB------------------------------

# # -----------------GAGB------------------------------
# # number_of_run_list = [2,  4, 8, 16]
# # number_of_run_list = [5, 8, 32]
# number_of_run_list = [32]
# for n in number_of_run_list:
#     name = "gb_"+str(n)+"m"
#     # os.system("mkdir "+name)
#     os.system("cp -r 2lhc.pdb "+name)
#     os.system("cp -r 2lhd.pdb "+name)
#     # os.system("cp -r 2lhc variables.dat "+name)
#     os.chdir(name)
#     for i in range(20):
#         os.chdir("analysis/"+str(i))
#         os.system("cp ../../2lhc.pdb .")
#         os.system("cp ../../2lhd.pdb .")
#         os.system("python2 ~/opt/script/CalcQValue.py 2lhd.pdb dump.lammpstrj q_gb.dat")
#         os.system("python2 ~/opt/script/CalcQValue.py 2lhc.pdb dump.lammpstrj q_ga.dat")
#         os.system("cp ~/opt/small_script/qw_gagb.plt .")
#         os.system("gnuplot qw_gagb.plt")
#         os.system("mv qw_gagb.pdf ../../results/qw_gagb_{0}.pdf".format(str(i)))
#         os.chdir("../..")
#     os.chdir("..")
#
# for n in number_of_run_list:
#     name = "gb_"+str(n)+"m"
#     # os.system("mkdir "+name)
#     os.system("cp -r 2lhd.pdb "+name)
#
#     # os.system("cp -r 2lhc variables.dat "+name)
#     os.chdir(name)
#     for i in range(20):
#         os.chdir("analysis/"+str(i))
#         os.system("paste q_ga.dat q_gb.dat > q_gagb.dat")
#         os.system("cp ~/opt/small_script/qw_ga-gb.plt .")
#         os.system("gnuplot qw_ga-gb.plt")
#         os.system("mv qw_ga-gb.pdf ../../results/qw_ga-gb_{0}.pdf".format(str(i)))
#         os.chdir("../..")
#     os.system("cp ~/opt/small_script/qw_ga_all.plt .")
#     os.system("gnuplot qw_ga_all.plt")
#     os.system("cp ~/opt/small_script/qw_gb_all.plt .")
#     os.system("gnuplot qw_gb_all.plt")
#     os.system("cp ~/opt/small_script/qw_diff_all.plt .")
#     os.system("gnuplot qw_diff_all.plt")
#     os.chdir("..")
# # -----------------GAGB------------------------------


# simulation_steps = 4 * 10**6
# warm_up_steps = 10 * 10**5
#
# seed(datetime.now())
# n= 20
# vmd = "/Applications/VMD\ 1.9.2.app/Contents/MacOS/startup.command"
#
# os.system("BuildAllAtomsFromLammps.py dump.lammpstrj movie")
# os.system("cp ~/opt/plot_scripts/2xov_movie.tcl .")
# os.system(vmd+" -e 2xov_movie.tcl ")
# os.system("mkdir -p MyResults")
# for i in range(n):
#     print(i)
#     os.chdir("analysis/"+str(i))
#     os.system("cp ~/opt/plot_scripts/2xov_movie_screenshot.tcl .")
#     os.system(vmd+" -e 2xov_movie_screenshot.tcl")
#     os.system("cp frame1000.tga ../../MyResults/frame"+str(i)+"_1000.tga")
#     #os.system("cp frame450.tga ../Results/frame"+folder_name+"_450.tga")
#     # os.system("movie.py "+protein_name)
#     os.chdir("../..")
#     # analysis

# folder_name = ""
# result_folder = "WeiLu_Aug_07"

# protein_list = ['T089', 'T120', 'T251', 'TOP7', '1UBQ']
# sublist = ['']
# # sublist = ['_ha', '_he']
# # sublist = ['_lp', '_he_lp']
# # folder_list = []
# for protein in protein_list:
#     for sub in sublist:
#         folder_name = protein+sub
#         os.chdir(folder_name)
#         os.chdir("best_2nd")
#         os.system("pymol ~/opt/plot_scripts/align.pml > matrix.dat")
#         os.system("head -n 70 matrix.dat | tail -n 20 > cealign_matrix.dat")
#         # for i in range(19, -1, -1):
#         #     os.system("mv {}.pdb {}.pdb".format(i, i+1))
#         os.chdir("../..")
# os.chdir(protein)
# os.chdir("best_1st")
# os.system("python3 ~/opt/small_script/cross_q.py")
# os.chdir("..")
# os.chdir("best_2nd")
# os.system("python3 ~/opt/small_script/cross_q.py")
# os.chdir("..")
# os.chdir("..")
# n = 3
# for i in range(n):
#     # simulation set up
#     folder_name = str(i)
#     os.system("mkdir -p "+folder_name)
#     os.system("cp -r "+args.protein+"* "+folder_name)
#     os.chdir(folder_name)
#     os.system("cp ../../helix_less/simulation/"+str(i)+"/restart.4000000 .")
#     os.system(  # replace SIMULATION_STEPS with specific steps
#         "sed -i.bak 's/WARM_UP_STEPS/'" +
#         str(warm_up_steps) +
#         "'/g' "+protein_name+".in")
#     os.system(  # replace RANDOM with a radnom number
#             "sed -i.bak 's/RANDOM/'" +
#             str(randint(1, 10**6)) +
#             "'/g' "+protein_name+".in")
#     os.system(  # replace SIMULATION_STEPS with specific steps
#             "sed -i.bak 's/SIMULATION_STEPS/'" +
#             str(simulation_steps) +
#             "'/g' "+protein_name+".in")
# # if(platform.system() == 'Darwin'):
# #     os.system("/Users/weilu/Documents/lammps-9Oct12_modified/src/lmp_serial \
# #     < "+protein_name+".in")
#     if(platform.system() == 'Darwin'):
#         os.system("/Users/weilu/Documents/lammps-9Oct12_modified/src/lmp_serial \
#         < "+protein_name+".in")
#     elif(platform.system() == 'Linux'):
#         os.system("cp ~/opt/run.slurm .")
#         os.system(  # replace PROTEIN with pdb name
#                 "sed -i.bak 's/PROTEIN/'" +
#                 protein_name +
#                 "'/g' run.slurm")
#         os.system("sbatch run.slurm")
#     else:
#         print("system unkown")
#     os.chdir("..")
# exit(1)

# w_helix_list = [0.1, 0.5, 1, 1.5]
# m_helix_list = [0.1, 0.5, 1, 1.5]
#
# for i in range(len(w_helix_list)):
#     w = w_helix_list[i]
#     for j in range(len(m_helix_list)):
#
#         # m = m_helix_list[j]
#         folder_name = str(i)+"_"+str(j)
#         # os.system("cd "folder_name)
#         os.chdir(folder_name)
#         # os.system("analysis.py 2xov/")
#         # os.system("echo "+folder_name+" >> ../all")
#         os.system("sort -k 3 analysis/list_of_max_q > ../data/"+folder_name)
#         os.chdir("..")
#         # os.system("mkdir "+folder_name)
#         # os.chdir(folder_name)
#         # os.system("cp -r ../2xov .")
#         # os.chdir("2xov")
#         # os.system(
#         #         "sed -i.bak 's/W_HELIX/'" +
#         #         str(w) +
#         #         "'/g' fix_backbone_coeff.data")
#         # os.system(
#         #         "sed -i.bak 's/M_HELIX/'" +
#         #         str(m) +
#         #         "'/g' fix_backbone_coeff.data")
#         # os.chdir("..")
#         # os.system("run.py 2xov/ -n 5")

# os.system("cp ~/opt/gg.py this_gg.py")
# for i in range(5):
#     os.system("mkdir "+str(i))
#     os.chdir(str(i))
#     os.system("cp -r ../2xov/ .")
#     os.system("cp ../../2xov_strong_single_memory_600to500/simulation/"+str(i)+"/restart.2000000 2xov/")
#     os.system("run.py -s 4 -n 2 2xov/")
#     os.chdir("..")

# # rama_list = [6, 8, 16]
# # rama_list = [4]
# melt_t_list = [400, 500, 600]
# for variable in melt_t_list:
#     folder_name = str(variable)
#     os.system("mkdir "+folder_name)
#     os.chdir(folder_name)
#     os.system("cp -r ../1qjp .")
#     os.chdir("1qjp")
#     os.system(
#             "sed -i.bak 's/MELTT/'" +
#             str(variable) +
#             "'/g' 1qjp.in")
#     os.chdir("..")
#     # os.system("pwd")
#     os.system("run.py 1qjp/ -n 5 -s 5")
#     os.chdir("..")
# os.system("cp ~/opt/gg.py this_gg.py")
#
# exec (open("config.py").read())
# n = number_of_run
# steps = simulation_steps
#
# protein_name = args.protein.strip('/')
#
# temp = 400
# folder_name = "{}_t{}_q100_test11".format(protein_name, str(temp))
# print("all going to "+folder_name)
# os.system("mkdir -p "+folder_name)
# os.system("rm -f "+folder_name + "/*")
# command = 'cat simulation/{}/%d/wham11 \
# >> {}/all_wham.dat'.format(temp, folder_name)
# # cal rmsd
# os.chdir("simulation/"+str(temp))
# for i in range(n):
#     os.chdir(str(i))
#     os.system("awk '{print>\"file1\"(NR>(n/2)?2:1)}' n=\"$(wc -l <file1)\" file1")
#     os.system("cat file11 >> ../../../"+folder_name+"/rmsd_total")
#     # os.system("sed 1d wham.dat > wham1d.dat")
#     os.system("awk '{print>\"wham1\"(NR>(n/2)?2:1)}' n=\"$(wc -l <wham1)\" wham1")
#     os.chdir("..")
# os.chdir("../..")
# for i in range(n):
#     cmd = command % i
#     os.system(cmd)
# os.chdir(folder_name)
# os.system("awk '{print $2}' all_wham.dat > Qw_total")
# os.system("awk '{print $3}' all_wham.dat > rg_total")
# os.system("awk '{print $4}' all_wham.dat > p_total")
# os.system("awk '{print $5}' all_wham.dat > tc_total")
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# os.system("cp ~/opt/wham_analysis/*.m .")
# os.chdir("..")
# os.system("~/opt/script/wham/fused_calc_cv.sc {} top7 50 400 350 450 5 50 100 0 0.98".format(folder_name))
#
#
# folder_name = "{}_t{}_q100_test12".format(protein_name, str(temp))
# print("all going to "+folder_name)
# os.system("mkdir -p "+folder_name)
# os.system("rm -f "+folder_name + "/*")
# command = 'cat simulation/{}/%d/wham12 \
# >> {}/all_wham.dat'.format(temp, folder_name)
# # cal rmsd
# os.chdir("simulation/"+str(temp))
# for i in range(n):
#     os.chdir(str(i))
#     os.system("cat file12 >> ../../../"+folder_name+"/rmsd_total")
#     os.chdir("..")
# os.chdir("../..")
# for i in range(n):
#     cmd = command % i
#     os.system(cmd)
# os.chdir(folder_name)
# os.system("awk '{print $2}' all_wham.dat > Qw_total")
# os.system("awk '{print $3}' all_wham.dat > rg_total")
# os.system("awk '{print $4}' all_wham.dat > p_total")
# os.system("awk '{print $5}' all_wham.dat > tc_total")
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# os.system("cp ~/opt/wham_analysis/*.m .")
# os.chdir("..")
#
#
#
# os.system("~/opt/script/wham/fused_calc_cv.sc {} top7 50 400 350 450 5 50 100 0 0.98".format(folder_name))

#
# result_folder = "WeiLu_Aug_07"
# os.system("mkdir -p "+result_folder)
# protein_list = ['T089', 'T120', 'T251', 'top7', '1UBQ']
# # sublist = ['_ha', '_he']
# sublist = ['_lp', '_he_lp']
# folder_list = []
# for protein in protein_list:
#     for sub in sublist:
#         folder_list += [protein+sub]
# print(folder_list)
# # exit(1)
# # awk '{print>'file'(NR>(n/2)?2:1)}' n='$(wc -l <test)' test
# for folder in folder_list:
#     print(folder)
#     os.chdir(folder)
#     exec (open("config.py").read())
#     n = number_of_run
#     steps = simulation_steps
#     os.system("mkdir -p ../{}/".format(result_folder)+folder+"/best_q")
#     os.system("sort analysis/list_of_max_q > ../{}/q_".format(result_folder)+folder+".dat")
#     for i in range(n):
#         # move
#         os.chdir("analysis/"+str(i))
#         os.system("cp chosen.pdb ../../../{}/".format(result_folder) + folder+"/best_q/"+str(i)+".pdb")
#         os.chdir("../..")
#     os.chdir("..")


# result_folder = "WeiLu_Aug_07"
# os.system("mkdir -p "+result_folder)
# protein_list = ['T089', 'T120', 'T251', 'top7', '1UBQ']
# # sublist = ['_ha', '_he']
# sublist = ['_lp', '_he_lp']
# folder_list = []
# for protein in protein_list:
#     for sub in sublist:
#         folder_list += [protein+sub]
# print(folder_list)
# # exit(1)
#
# for folder in folder_list:
#     print(folder)
#     os.chdir(folder)
#     exec (open("config.py").read())
#     n = number_of_run
#     steps = simulation_steps
#     os.system("mkdir -p ../{}/".format(result_folder)+folder+"/best_q")
#     os.system("sort analysis/list_of_max_q > ../{}/q_".format(result_folder)+folder+".dat")
#     for i in range(n):
#         # move
#         os.chdir("analysis/"+str(i))
#         os.system("cp chosen.pdb ../../../{}/".format(result_folder) + folder+"/best_q/"+str(i)+".pdb")
#         os.chdir("../..")
#     os.chdir("..")
