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
parser.add_argument("-l", "--label", type=str, default="label")
args = parser.parse_args()

with open('gg_cmd.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')


def do(cmd, get=False, show=True):
    if get:
        out = subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True).communicate()[0].decode()
        if show:
            print(out, end="")
        return out
    else:
        return subprocess.Popen(cmd, shell=True).wait()
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
dataset["may13"] = ['1r69', '3icb', '256b', '4cpv', '2mhr', '1mba', '2fha', '1fc2', '1enh', '2gb1', '2cro', '1ctf', '4icb']
dataset["membrane"] = ["2bg9", "1j4n", "1py6_SD", "2bl2", "1rhz", "1iwg", "2ic8", "1pv6", "1occ", "1kpl", "2bs2", "1py6", "1u19"]
dataset["hybrid"] = ["2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91", "4nv6", "4p79", "5dsg", "6g7o", "6a93", "2jo1", "1py6", "1pv6", "1u19"]
dataset["optimization_cath"] = ['1a75A00', '1bekA01', '1bqbA02', '1cpcB00', '1cscA02', '1cy5A00', '1dv5A00', '1e8yA05', '1evyA02', '1in4A03', '1l1fA03', '1vq8P01', '1xmkA00', '1zcaA02', '2grhA00', '2ii2A04', '2q6fB03', '2wh6A00', '3g0vA00', '3geuA00', '3h99A03', '3hrdD02', '3ju5A01', '3p1wA03', '4cxfA01', '4i2aA01', '4i4tB03', '4i6uB00', '5kn9A02']
dataset["optimization_v2"] = ['1e0m', '1w4e', '1e0g', '2wqg', '1jo8', '1fex', '2l6r', '1c8c', '1g6p', '1mjc', '2jmc', '1hdn', '1st7', '1n88', '1d6o', '2ga5', '1j5u', '3o4d']

def create_project_for_pdb_list(pdb_list, frag=False):
    do("mkdir -p setup")
    cd("setup")
    for pdb in pdb_list:
        do(f"mkdir -p {pdb}")
        cd(pdb)
        if frag:
            addFrag = "--frag"
        do(f"mm_create_project.py ../../cleaned_pdbs/{pdb}.pdb --extended --membrane {addFrag}")
        cd("..")

def get_two_part_from_eye_seperation(pdb, data):
    row = data.query(f"Protein == '{pdb}'")
    assert len(row) == 1
    row = row.iloc[0]
    glob_start, glob_end = row["Range"].split("-")
    length = row["Length"]
    glob_start = int(glob_start)
    glob_end = int(glob_end)
    print(pdb, glob_start, glob_end)
    GlobularPart = list(range(glob_start, glob_end+1))
    MembranePart = []
    for i in range(1, length+1):
        if i not in GlobularPart:
            MembranePart.append(i)
    return GlobularPart, MembranePart


def mix_frag(globular, mem, GlobularPart, MembranePart):
    out = """[Target]
query

[Memories]
"""
    n = len(mem)//20
    for i in range(n):
        if i in MembranePart:
            out += "".join(mem[i*20:(i+1)*20])
        else:
            out += "".join(globular[i*20:(i+1)*20])

    with open("HA_combined.mem", "w") as o:
        o.write(out)


def get_aligned_info(p1, p2):
    # do(f"~/opt/TMalign/TMalign {p1} {p2} > tm_log")
    cmd = f"~/opt/TMalign/TMalign {p1} {p2} | grep -A 1 'Aligned'"
    output = getFromTerminal(cmd)
    # print(output)
    line = output.split("\n")[0]
    line2 = output.split("\n")[1]
    aligned_length,rmsd,seqid = line.split(",")
    aligned_length = int(aligned_length.split("=")[1])
    rmsd = float(rmsd.split("=")[1])
    tmscore = float(line2.split(",")[0].split("=")[-1].strip().split(" ")[0])
    seqid = float(seqid.split("=")[-1])
    # print("aligned_length, rmsd, tmscore, seqid")
    # print(aligned_length, rmsd, tmscore, seqid)
    return aligned_length, rmsd, tmscore, seqid

def convert_frag_to_cbd(fragFile):
    # a = "frags.mem"
    a = fragFile
    with open(a) as f:
        b = f.readlines()
    pdb_list = [c.split()[0].split("/")[-1].split(".")[0] for c in b[4:]]

    # pdb_list = ["2ix5a"]
    pdb_list = list(set(pdb_list))
    print(len(pdb_list))
    for pdb in pdb_list:
        print(pdb)
        upperPDB = pdb[:4].upper()
        pre = f"/Users/weilu/openmmawsem/PDBs/"
        toPre = "frag_in_cbd"
        do(f"mkdir -p {toPre}/pdbs")
        do(f"mkdir -p {toPre}/frags")
        fromFile = f"{pre}/{upperPDB}.pdb"
        toFile = f"{toPre}/pdbs/cbd_{upperPDB}.pdb"
        convert_all_atom_pdb_to_cbd_representation(fromFile, toFile)
        # cmd = f"python ~/openmmawsem/helperFunctions/Pdb2Gro.py {toPre}/pdbs/cbd_{upperPDB}.pdb {toPre}/{pdb}.gro"
        pdbFile = toFile
        groFile = f"{toPre}/frags/{pdb}.gro"
        chainID = pdb[-1]
        Pdb2Gro(pdbFile, groFile, chainID.upper())
        do(f"cp {fragFile} cbd_frags.mem")
        replace(f"cbd_frags.mem", "./fraglib/", f"{toPre}/frags/")

if args.day == "sep20":
    if args.mode == 1:
        pdb_list = ["2bg9"]
        # he frag memeory
        for pdb in pdb_list:
            frag_folder = f"frag_memory/{pdb}"
            do(f"mkdir -p {frag_folder}")
            cd(frag_folder)
            do(f"cp ../../setups/{pdb}/{pdb}.fasta .")
            do(f"python ~/openmmawsem/helperFunctions/MultCha_prepFrags_index.py ~/openmmawsem/database/cullpdb_pc80_res3.0_R1.0_d160504_chains29712 {pdb}.fasta 20 1 9 > logfile")

            sys.path.insert(0, "/Users/weilu/openmmawsem")
            # name = "1iwg"
            import openmmawsem
            import helperFunctions.myFunctions
            helperFunctions.myFunctions.check_and_correct_fragment_memory("frags.mem")
            helperFunctions.myFunctions.relocate(fileLocation="frags.mem", toLocation="fraglib")
            # print(f"{__location__}//Gros/")
            helperFunctions.myFunctions.replace(f"frags.mem", f"{__location__}//Gros/", f"../../frag_memory/{pdb}/fraglib/")
            do(f"cp frags.mem ../../setups/{pdb}/he.mem")
            cd("../..")

if args.day == "nov16":
    if args.mode == 2:
        pdb_list = ['5ncq', '1uaz', '5lwe', '6e9o', '5i6x', '3vvo', '4tph', '6m20', '6d32']
        new_target_set = pd.read_csv("/Users/weilu/Research/data/openMM/membrane_protein_target_set_run1_submode1_11-15.csv").reset_index(drop=True)
        y = "Q"
        # d = pd.concat([data, previous_data])
        # d = pd.concat([data, mixed_iter2, iter1_v1])
        d = pd.concat([new_target_set])
        d = d.query("Steps > 1000").reset_index(drop=True)
        t = d.groupby(["Protein", "Folder"])[y].idxmax().reset_index()
        max_Q_data = d.iloc[t[y].to_list()].reset_index(drop=True)
        # sub_data = max_Q_data
        max_Q_data["Protein_sorted"] = pd.Categorical(max_Q_data["Protein"], pdb_list)
        print(max_Q_data)
        for i, line in max_Q_data.iterrows():
            pdb = line["Protein"]
            run = line["Run"]
            folder = line["Folder"]
            # do(f"scp -r wl45@nots.rice.edu:/scratch/wl45/nov_2020/membrane_protein_target_set/{folder}/{pdb}/{run}/movie.dcd {pdb}_{folder}_{run}.movie.dcd")
            # do(f"scp -r wl45@nots.rice.edu:/scratch/wl45/nov_2020/membrane_protein_target_set/setups/{pdb}/minimization-openmmawsem.pdb {pdb}.pdb")
            do(f"scp -r wl45@nots.rice.edu:/scratch/wl45/nov_2020/membrane_protein_target_set/{folder}/{pdb}/{run}/lastFrame.pdb last_frame_{pdb}_{folder}_{run}.pdb")
    if args.mode == 1:
        pdb_list = ["2bg9", "1j4n", "1py6_SD", "2bl2", "1iwg", "2ic8", "1pv6", "1occ", "1kpl", "2bs2", "1py6", "1u19"]
        cmd = " ".join([f"setups/{pdb}/{pdb}-cleaned.pdb" for pdb in pdb_list])
        do(f"pymol {cmd}")
if args.day == "nov14":
    if args.mode == 1:
        pdb_list = manual_chosen = ['1uaz', '5lwe', '6e9o', '5i6x', '3vvo', '4tph', '6m20', '6d32']
        for pdb in pdb_list:
            a = open(f"/Users/weilu/Research/server/nov_2020/membrane_protein_target_set/test/qucik_start_native/{pdb}/native.pdb").readlines()
            keep = False
            info = []
            for line in a:
                if len(line) < 5:
                    continue
                if keep:
                    info.append(line)
                if line[:5] == "MODEL":
                    # print(line)
                    # print(line.split())
                    if line.split()[1] == "2":
                        keep = True
            with open(f"/Users/weilu/Research/server/nov_2020/membrane_protein_target_set/setups/{pdb}/minimization-openmmawsem.pdb", "w") as out:
                for line in info:
                    out.write(line)

if args.day == "nov12":
    if args.mode == 1:
        # zim to side
        # PredictedZim to PredictedZimSide
        # 3 to down, 2 to middle, 1 to up
        pdb_list = ["5ncq"]
        pdb_list = manual_chosen = ['1uaz', '5lwe', '6e9o', '5i6x', '3vvo', '4tph', '6m20', '6d32']
        for pdb in pdb_list:
            a = np.loadtxt(f"/Users/weilu/Research/server/nov_2020/membrane_protein_target_set/setups/{pdb}/zimPosition", dtype=int)
            with open(f"/Users/weilu/Research/server/nov_2020/membrane_protein_target_set/setups/{pdb}/zimPositionSide", "w") as out:
                for res in a:
                    if res == 1:
                        out.write("up\n")
                    if res == 2:
                        out.write("middle\n")
                    if res == 3:
                        out.write("down\n")


if args.day == "nov11":
    if args.mode == 3:
        pdb_list = manual_chosen = ['1uaz', '5lwe', '6e9o', '5i6x', '3vvo', '4tph', '6m20', '6d32']
        for pdb in pdb_list:
            cmd = f"mkdir -p setups/{pdb}"
            do(cmd)
            cd(f"setups/{pdb}")
            cmd = f"python ~/openmmawsem/mm_create_project.py ~/Dropbox/gxxxg/membrane_protein_target_set_curation/cleaned_manual_adjust/{pdb}.pdb --hybrid"
            do(cmd)
            cd("../..")
    if args.mode == 2:
        # pdb_list = ["5ncq"]
        pdb_list = manual_chosen = ['1uaz', '5lwe', '6e9o', '5i6x', '3vvo', '4tph', '6m20', '6d32']
        source = "/Users/weilu/Dropbox/gxxxg/membrane_protein_target_set_curation/manual_adjust/"
        toFolder = "/Users/weilu/Dropbox/gxxxg/membrane_protein_target_set_curation/cleaned_manual_adjust/"
        cleanPdb(pdb_list, source=source, addMissingResidues=False, toFolder=toFolder, chain=-1, formatName=False, removeTwoEndsMissingResidues=True)
    if args.mode == 1:
        chosen_list = pd.read_csv("/Users/weilu/Dropbox/gxxxg/chosen_one_from_each_family", index_col=0)
        pdb_list = chosen_list["pdb"].to_list()
        for pdb in pdb_list:
            do(f"cp ../cleaned_membrane_part/{pdb}.pdb {pdb}_cleaned.pdb")
            do(f"cp ../original/{pdb}.pdb {pdb}_original.pdb")
            do(f"pymol {pdb}_cleaned.pdb {pdb}_original.pdb")
if args.day == "nov10":
    if args.mode == 2:
        data = pd.read_csv("/Users/weilu/Dropbox/gxxxg/curated_membrane_protein_target_set_curation.csv", index_col=0)
        os.system("mkdir -p /Users/weilu/Dropbox/gxxxg/membrane_protein_target_set_curation/membrane_part")
        pdb_list = data["pdbid"].to_list()
        pdb_list = [pdb[2:6] for pdb in pdb_list]
        source = "/Users/weilu/Dropbox/gxxxg/membrane_protein_target_set_curation/membrane_part/"
        toFolder = "/Users/weilu/Dropbox/gxxxg/membrane_protein_target_set_curation/cleaned_membrane_part/"
        cleanPdb(pdb_list, source=source, addMissingResidues=False, toFolder=toFolder, chain=-1, formatName=False, removeTwoEndsMissingResidues=True)
    if args.mode == 1:
        data = pd.read_csv("/Users/weilu/Dropbox/gxxxg/curated_membrane_protein_target_set_curation.csv", index_col=0)
        os.system("mkdir -p /Users/weilu/Dropbox/gxxxg/membrane_protein_target_set_curation/membrane_part")
        for pdb in data["pdbid"]:
            # print(pdb)
            pdbName = pdb[2:6]
            do(f"python ~/opt/small_script/extract_membrane_part.py /Users/weilu/Dropbox/gxxxg/membrane_protein_target_set_curation/original/{pdbName}.pdb /Users/weilu/Dropbox/gxxxg/membrane_protein_target_set_curation/membrane_part/{pdbName}.pdb -c 17")

if args.day == "jul06":
    if args.mode == 1:
        # get cbd openAWSEM input.
        # pre = "./"
        pdb_list = dataset["membrane"]
        for pdb in pdb_list:
            pre = f"setups/{pdb}/"
            original_openAWSEM_input = f"{pre}/{pdb}-openmmawsem.pdb"
            new_openAWSEM_input = f"{pre}/cbd-openmmawsem.pdb"
            all_atom_pdb_file = f"{pre}/crystal_structure-cleaned.pdb"
            replace_CB_coord_with_CBD_for_openAWSEM_input(original_openAWSEM_input, new_openAWSEM_input, all_atom_pdb_file)

            # get single memory in cbd format
            # pre = "./"
            fromFile = f"{pre}/crystal_structure.pdb"
            toFile = f"{pre}/cbd_{pdb}.pdb"
            convert_all_atom_pdb_to_cbd_representation(fromFile, toFile)
            cmd = f"python ~/openmmawsem/helperFunctions/Pdb2Gro.py {pre}/cbd_{pdb}.pdb {pre}/cbd_{pdb}.gro"
            do(cmd)
            do(f"cp {pre}/single_frags.mem {pre}/cbd_single_frags.mem")
            # replace(f"{pre}/cbd_single_frags.mem", pdb, f"cbd_{pdb}")
            replace(f"{pre}/cbd_single_frags.mem", f"{pdb}_A", f"cbd_{pdb}")

if args.day == "jul01":
    # pdb_list = ["2jo1"]
    pdb_list = dataset["membrane"]
    if args.mode == 1:
        # downloadPdb(pdb_list)
        # pdb_list = ["2bg9", "1j4n", "1py6", "2bl2", "1rhz", "1iwg", "2ic8", "1pv6", "1occ", "1kpl", "2bs2", "1py6", "1u19", "2rh1"]
        # pdb_list = ["2xov"]
        # pdb_list = ["2bs2", "1py6", "1u19", "2rh1"]
        do("mkdir -p original_pdbs")
        for pdb in pdb_list:
            do(f"wget https://opm-assets.storage.googleapis.com/pdb/{pdb}.pdb")
            do(f"mv {pdb}.pdb original_pdbs/")
        # cleanPdb(pdb_list, source="original_pdbs", chain=-1, formatName=False)
    if args.mode == 2:
        do("mkdir -p extracted")
        from small_script.extract_pdb import *
        protein_info_list = []
        protein_info_list.append(("1py6_SD", "A",77, 199))
        protein_info_list.append(("1py6", "A",5, 231))
        protein_info_list.append(("2ic8", "A",91,272))
        protein_info_list.append(("1occ", "C",71,261))
        protein_info_list.append(("1pv6", "A", 1, 190))
        protein_info_list.append(("1j4n", "A",4,119))
        protein_info_list.append(("2bs2", "C", 21,237))
        protein_info_list.append(("2bl2", "A", 12, 156))
        protein_info_list.append(("2bg9", "A", 211, 301))
        protein_info_list.append(("1iwg", "A", 330, 497))
        protein_info_list.append(("1rhz", "A", 23, 188))
        protein_info_list.append(("1kpl", "A", 31, 233))
        protein_info_list.append(("1u19", "A", 33, 310))
        for (protein, chain, residue_start, residue_end) in protein_info_list:
            # do(f"wget {0}{1} -O ~/opt/crystal_structures/membrane_proteins/original_pdb/{1}".format(pdbFileDataBase, protein+".pdb"))
            # do(f"wget {0}{1} -O {1}".format(pdbFileDataBase, protein+".pdb"))
            extract_pdb("./", protein, chain, residue_start, residue_end)
        # cleanPdb(pdb_list, chain=-1, formatName=False)
        do("rm tmp.pdb")
    if args.mode == 3:
        cleanPdb(pdb_list, source="extracted", chain=-1, formatName=False)
    if args.mode == 4:
        for pdb in pdb_list:
            # do(f"getSize.py extracted/{pdb}.pdb")
            do(f"getSize.py cleaned_pdbs/{pdb}.pdb")
    if args.mode == 5:
        for pdb in pdb_list:
            do(f"mkdir -p setups/{pdb}")
            cd(f"setups/{pdb}")
            do(f"~/openmmawsem/mm_create_project.py ../../extracted/{pdb}.pdb --hybrid")
            cd("../..")
    if args.mode == 6:
        for pdb in pdb_list:
            # do(f"mkdir -p setups/{pdb}")
            # cd(f"setups/{pdb}")
            do(f"python mm_run.py setups/{pdb}/{pdb} --to native/{pdb} -s 1e5 -r 2000 --tempStart 300 ")
            # cd("../..")
if args.day == "jun21":
    if args.mode == 1:
        # fromPdb = "1mnn.pdb"
        fromPdb = "1mnn-openmmawsem.pdb"
        toPdb = "shifted_1mnn.pdb"
        rotation_and_translation(fromPdb, toPdb, rotation_axis=(0,0,0), degree=0, translation=(100,0,0))

if args.day == "jun20":
    pdb_list = ['4y60', '5ke8', '1a1j', '5lxu', '1skn', '6a2h']
    if args.mode == 1:
        for pdb in pdb_list:
            folder = f"setups/{pdb}"
            do(f"mkdir -p {folder}")
            cd(folder)
            do(f"python /Users/weilu/openmmawsem/mm_create_project.py ../../../protein_DNA/cleaned_pdbs/{pdb}.pdb --extended --frag")
            cd("../..")
            # do(f"cp gamma_noCysCys.dat {folder}/")

if args.day == "jun12":
    pdb = "1su4"
    if args.mode == 1:
        # get cbd openAWSEM input.
        pre = "./"
        original_openAWSEM_input = f"{pre}/{pdb}-openmmawsem.pdb"
        new_openAWSEM_input = f"{pre}/cbd-openmmawsem.pdb"
        all_atom_pdb_file = f"{pre}/crystal_structure-cleaned.pdb"
        replace_CB_coord_with_CBD_for_openAWSEM_input(original_openAWSEM_input, new_openAWSEM_input, all_atom_pdb_file)
    if args.mode == 2:
        # get single memory in cbd format
        pre = "./"
        fromFile = f"{pre}/crystal_structure.pdb"
        toFile = f"{pre}/cbd_{pdb}.pdb"
        convert_all_atom_pdb_to_cbd_representation(fromFile, toFile)
        cmd = f"python ~/openmmawsem/helperFunctions/Pdb2Gro.py {pre}/cbd_{pdb}.pdb {pre}/cbd_{pdb}.gro"
        do(cmd)
        do(f"cp {pre}/single_frags.mem {pre}/cbd_single_frags.mem")
        replace(f"{pre}/cbd_single_frags.mem", pdb, f"cbd_{pdb}")

if args.day == "may16":
    # disulfide setups and annealing.
    pdb_list = ["1ppb", "1fs3", "1bpi", "1hn4", "1lmm", "1tcg"]
    pdb_list = ["1ppb_H"]
    pdb_list = ["1hn4_A"]
    if args.mode == 1:
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, source="original_pdbs/", addMissingResidues=True, toFolder="cleaned_pdbs/", chain=-1, formatName=False, removeTwoEndsMissingResidues=True)
    if args.mode == 2:
        for pdb in pdb_list:
            folder = f"setups/{pdb}"
            do(f"mkdir -p {folder}")
            cd(folder)
            do(f"mm_create_project.py ../../cleaned_pdbs/{pdb}.pdb --extended")
            cd("../..")
            do(f"cp gamma_noCysCys.dat {folder}/")
    if args.mode == 3:
        # ho frag memeory
        for pdb in pdb_list:
            frag_folder = f"frag_memory/{pdb}"
            do(f"mkdir -p {frag_folder}")
            cd(frag_folder)
            do(f"cp ../../setups/{pdb}/{pdb}.fasta .")
            do(f"python ~/openmmawsem/helperFunctions/MultCha_prepFrags_index.py ~/openmmawsem/database/cullpdb_pc80_res3.0_R1.0_d160504_chains29712 {pdb}.fasta 20 2 9 > logfile")

            sys.path.insert(0, "/Users/weilu/openmmawsem")
            # name = "1iwg"
            import openmmawsem
            import helperFunctions.myFunctions
            helperFunctions.myFunctions.check_and_correct_fragment_memory("frags.mem")

            do(f"cp frags.mem ../../setups/{pdb}/ho.mem")
            cd("../..")
    if args.mode == 4:
        for pdb in pdb_list:
            do(f"python mm_run.py setups/{pdb}/{pdb} --to native/{pdb}_ho -f forces_setup_gamma_noCysCys.py --subMode 31 --reportFrequency 10 -s 20")
    if args.mode == 5:
        pdb = "1lmm"
        do(f"python mm_run.py setups/{pdb}/{pdb} --to test/{pdb}_ho -f forces_setup_gamma_noCysCys.py --subMode 31 --reportFrequency 10 -s 20")
if args.day == "nov29":
    pdb = "1hn4"
    pdb = "1lmm"
    pdb = "1tcg"
    pdb = args.label
    if args.mode == 4:
        for pdb in ["1tcg", "1lmm"]:
            do(f"gg.py -d nov29 -m 1 -l {pdb}")
            do(f"gg.py -d nov29 -m 2 -l {pdb}")
            do(f"gg.py -d nov29 -m 3 -l {pdb}")
    if args.mode == 1:
        folder = f"setups/{pdb}"
        do(f"mkdir -p {folder}")
        cd(folder)
        do(f"mm_create_project.py ../../original_pdbs/{pdb}.pdb --extended")
        cd("../..")
        do(f"cp gamma_noCysCys.dat {folder}/")
    if args.mode == 2:
        # ho frag memeory
        frag_folder = f"frag_memory/{pdb}"
        do(f"mkdir -p {frag_folder}")
        cd(frag_folder)
        do(f"cp ../../setups/{pdb}/{pdb}.fasta .")
        do(f"python ~/openmmawsem/helperFunctions/MultCha_prepFrags_index.py ~/openmmawsem/database/cullpdb_pc80_res3.0_R1.0_d160504_chains29712 {pdb}.fasta 20 2 9 > logfile")

        sys.path.insert(0, "/Users/weilu/openmmawsem")
        # name = "1iwg"
        import openmmawsem
        import helperFunctions.myFunctions
        helperFunctions.myFunctions.check_and_correct_fragment_memory("frags.mem")

        do(f"cp frags.mem ../../setups/{pdb}/ho.mem")
        cd("../..")
    if args.mode == 3:
        do(f"python mm_run.py setups/{pdb}/{pdb} --to native/{pdb}_single -f forces_setup_gamma_noCysCys.py --subMode 1 --reportFrequency 10 -s 20")
        do(f"python mm_run.py setups/{pdb}/{pdb} --to native/{pdb}_ho -f forces_setup_gamma_noCysCys.py --subMode 31 --reportFrequency 10 -s 20")
    # if args.mode == 2:
    #     do(f"cp frags.mem ../../setups/{pdb}/ho.mem")
    # if args.mode == 3:
    #     sys.path.insert(0, "/Users/weilu/openmmawsem")
    #     # name = "1iwg"
    #     import openmmawsem
    #     import helperFunctions.myFunctions
    #     helperFunctions.myFunctions.check_and_correct_fragment_memory("frags.mem")


if args.day == "may12":
    if args.mode == 1:
        parser = PDBParser(QUIET=True)
        fileName = "/Users/weilu/Research/server/may_week2_2020/cornichon/source/6ud8.pdb"
        s = parser.get_structure("X", fileName)


        class ResSelect(Select):
            def accept_residue(self, residue):
                try:
                    if residue["CA"].get_coord()[2] < 30:
                        return 1
                    else:
                        return 0
                except:
                    if residue.id[0] == "H_ZK1":
                        return 0
                    else:
                        return 1

        io = PDBIO()
        io.set_structure(s)
        io.save('/Users/weilu/Research/server/may_week2_2020/cornichon/source/cutoff_z30.pdb', ResSelect())
    if args.mode == 2:
        do("pdb_selchain -A,B,C,D,E cutoff_z30.pdb > chainABCDE.pdb")
        # do("pdb_selchain -E 20200120_171355_2866_6ud8.pdb > chainE.pdb")
        # do("pdb_selchain -A,B,C,D 20200120_171355_2866_6ud8.pdb > chain_ABCD.pdb")
    if args.mode == 3:
        # pdb_list = ["chainE"]
        pdb_list = ["20200120_171355_2866_6ud8"]
        cleanPdb(pdb_list, source="./", addMissingResidues=True, toFolder="complete_chainE/", chain="E", formatName=False, removeTwoEndsMissingResidues=False, verbose=True)
    if args.mode == 4:
        # pdb_list = ["chainE"]
        pdb_list = ["combine_chain_ABCD_and_complete_chainE"]
        cleanPdb(pdb_list, source="./", removeHeterogens=False, addMissingResidues=False, toFolder="cleaned/", chain="ABCDE", formatName=False, removeTwoEndsMissingResidues=False, verbose=True)
        # then manual combine those two

if args.day == "may08":

    pdb = "6T7A"
    chosen = ["5KL6", "6CEU", "6TEM", "6L49", "6S01", "6KE9"]
    pdb_list = chosen
    if args.mode == 11:
        downloadPdb(pdb_list)
    if args.mode == 1:
        #
        for p in pdb_list:
            print(p)
            pdb = p.lower()
            # do(f"download.py {pdb}")
            do(f"mkdir -p setups/{pdb}")
            cd(f"setups/{pdb}")
            # cmd = f"python setup_DNA.py original_pdbs/{pdb}.pdb -o cleaned/{pdb}_clean.pdb"
            cmd = f"mm_create_project.py ../../source/original_pdbs/{pdb}.pdb"
            do(cmd)
            cd("../..")
            # cmd = f"grep 'S   D' cleaned/{pdb}_clean.pdb | wc"
            # print("DNA")
            # line = do(cmd, get=True, show=True).strip()
            # do(f"echo '{line}' > data/{pdb}_DNA.dat")

            # cmd = f"grep 'CA' cleaned/{pdb}_clean.pdb | wc"
            # print("Protein")
            # line = do(cmd, get=True, show=True).strip()
            # do(f"echo '{line}' > data/{pdb}_Protein.dat")
    if args.mode == 2:
        for p in pdb_list:
            pdb = p.lower()
            do(f"cp cleaned/{pdb}_clean.pdb chosen/")

if args.day == "may07":
    if args.mode == 1:
        pdb_list = ["1akr", "1opd", "1ptf", "1tig", "1tmy", "2acy", "5nul"]
        for pdb in pdb_list:
            do(f"python mm_evaluate_native.py setups/{pdb}/{pdb} --to native_frag/{pdb} --subMode 0")
    if args.mode == 2:
        import sys
        sys.path.insert(0, "/Users/weilu/openmmawsem/helperFunctions/")
        from Pdb2GroLib import *
        pdb_list = ["1akr", "1opd", "1ptf", "1tig", "1tmy", "2acy", "5nul"]
        for pdb in pdb_list:
            cd(f"setups/{pdb}")
            fragFile = "frags.mem"
            convert_frag_to_cbd(fragFile)
            cd("../..")
    if args.mode == 3:
        pdb_list = ["1akr", "1opd", "1ptf", "1tig", "1tmy", "2acy", "5nul"]
        for pdb in pdb_list:
            do(f"python mm_evaluate_native.py setups/{pdb}/{pdb} --to native_frag_cbd/{pdb} --subMode 14 -f forces_setup_cbd.py")
if args.day == "may06":
    pdb_list = ['5ITH', '6D1T', '5NNX', '4X9J', '5KL5', '5KL2', '5KL7', '5KL6', '5KL3', '5D5W', '5ZK1', '5D8K', '6OGK', '6NCE', '5GZB', '6AKP', '6AKO', '4XZF', '5L0M', '5BT2', '6JVY', '6A8R', '5D5U', '5K5J', '5K5I', '5VC9', '5K5H', '6F1K', '6TCE', '4USG', '5MEY', '5EG0', '6E94', '6E93', '6DFY', '5EGO', '5Z2T', '6BLW', '5HOD', '4XRM', '6TBZ', '5BNG', '5T01', '5W9Q', '6U81', '6DFB', '5VMU', '5VMW', '5VMV', '5VMY', '5VMX', '5VMZ', '6DF8', '6DF5', '6DF9', '6V8U', '6DFC', '6DFA', '5VA7', '6E8C', '5VA0', '6BQU', '6Q6R', '6C1T', '5UC6', '5EGB', '5ZFW', '5ZFZ', '5ZFY', '6MDX', '6MG3', '6MG2', '5Z6Z', '6CEV', '6CCG', '6CC8', '6CEU', '5LUX', '5TD5', '5CC0', '6CNQ', '4WK8', '6FBR', '6BUX', '4CN2', '4CN5', '6FBQ', '5EYO', '6LFF', '5BUA', '6L6Q', '6L6L', '6ES2', '5LTY', '6ES3', '5CHZ', '4OFA', '4OFH', '4OFE', '5KEG', '6IIR', '6OGJ', '6C1Y', '5EMC', '5EMP', '5EMQ', '4Y0F', '5KL4', '6NCM', '6O3T', '5YBD', '6FQP', '6CNP', '5D8F', '6IIT', '6IIS', '6EL8', '5YJ3', '5BNH', '4XEG', '5CYS', '5JXY', '4YJ0', '4Z7B', '4Z7Z', '5E69', '5E6B', '5E6A', '4Z47', '5E6D', '5E6C', '4Z3A', '4YO2', '5DUI', '5ZYV', '6KI6', '6OD4', '5I50', '6T78', '5T2W', '6U15', '5HF7', '6U16', '6U17', '5D5V', '5ZU1', '5ODG', '5OD6', '6B0O', '6B0P', '6OEB', '6OEA', '6OE7', '6B0Q', '5FF8', '4XRS', '4ZBN', '6FZS', '6OD5', '6B0R']
    pdb_list = ['5KL6', '6CEU', '4OFA', '5E6A', '5ZOE', '5V09', '5D9Y', '5KFL', '4ROD', '6P94', '6P09', '6BKG', '5XFQ', '5ZKJ', '5SWW', '5ZQF', '5NWA', '6EDB', '5YX2', '6KE9', '5W2C', '5W1C', '5ITU', '5E5A', '5XM0', '5Y0D', '6SE6', '5IIO', '6C0W', '6T7C', '5U1C', '6O96', '6NJ9', '6NOG', '6S01', '6PWX', '5N96', '5N8U', '6ASW', '6T90', '6R25', '5GSE', '6CIL', '6ERG', '6UPK', '4U7D', '3JBX', '6UPL', '6R8Z', '3JBW', '6CG0', '6R92']

    pdb = "6T7A"
    chosen = ["5KL6", "6CEU", "6TEM", "6L49", "6S01", "6KE9"]
    pdb_list = chosen
    if args.mode == 11:
        downloadPdb(pdb_list)
    if args.mode == 1:
        #
        for p in pdb_list:
            print(p)
            pdb = p.lower()
            # do(f"download.py {pdb}")
            cmd = f"python setup_DNA.py original_pdbs/{pdb}.pdb -o cleaned/{pdb}_clean.pdb"
            do(cmd)

            cmd = f"grep 'S   D' cleaned/{pdb}_clean.pdb | wc"
            print("DNA")
            line = do(cmd, get=True, show=True).strip()
            do(f"echo '{line}' > data/{pdb}_DNA.dat")

            cmd = f"grep 'CA' cleaned/{pdb}_clean.pdb | wc"
            print("Protein")
            line = do(cmd, get=True, show=True).strip()
            do(f"echo '{line}' > data/{pdb}_Protein.dat")
    if args.mode == 2:
        for p in pdb_list:
            pdb = p.lower()
            do(f"cp cleaned/{pdb}_clean.pdb chosen/")
if args.day == "may22":
    if args.mode == 1:
        pdb_list = ["1akr", "1opd", "1ptf", "1tig", "1tmy", "2acy", "5nul"]
        for pdb in pdb_list:
            name = pdb
            cd(f"setups/{pdb}")
            protein_length = getFromTerminal("wc ssweight").split()[0]
            print(f"protein: {name}, length: {protein_length}")

            new_line = f"cbd_{name}.gro 1 1 {protein_length} 1"
            do("cp cbd_frags.mem cbd_frags_include_native.mem")
            do(f"echo {new_line} >> cbd_frags_include_native.mem")
            cd("../..")

if args.day == "may05":
    # cornichon
    if args.mode == 1:
        pdb_list = ["6ud8"]
        cleanPdb(pdb_list, source="./", addMissingResidues=True, toFolder="cleaned_pdbs/", chain="E", formatName=False, removeTwoEndsMissingResidues=True)
    if args.mode == 2:
        do(f"mm_create_project.py ../source/chainE_no_filling.pdb --extended --frag --membrane")
    if args.mode == 3:
        # get cbd openAWSEM input.
        pre = "./"
        pdb = "chainE_no_filling"
        original_openAWSEM_input = f"{pre}/{pdb}-openmmawsem.pdb"
        new_openAWSEM_input = f"{pre}/cbd-openmmawsem.pdb"
        all_atom_pdb_file = f"{pre}/crystal_structure-cleaned.pdb"
        replace_CB_coord_with_CBD_for_openAWSEM_input(original_openAWSEM_input, new_openAWSEM_input, all_atom_pdb_file)
    if args.mode == 4:
        # get single memory in cbd format
        pre = "./"
        pdb = "chainE_no_filling"
        fromFile = f"{pre}/crystal_structure.pdb"
        toFile = f"{pre}/cbd_{pdb}.pdb"
        convert_all_atom_pdb_to_cbd_representation(fromFile, toFile)
        cmd = f"python ~/openmmawsem/helperFunctions/Pdb2Gro.py {pre}/cbd_{pdb}.pdb {pre}/cbd_{pdb}.gro"
        do(cmd)
        do(f"cp {pre}/single_frags.mem {pre}/cbd_single_frags.mem")
        replace(f"{pre}/cbd_single_frags.mem", pdb, f"cbd_{pdb}")
    if args.mode == 5:
        # convert all pdb in fras.mem
        a = "frags.mem"
        with open(a) as f:
            b = f.readlines()
        pdb_list = [c.split()[0].split("/")[-1].split(".")[0] for c in b[4:]]
        import sys
        sys.path.insert(0, "/Users/weilu/openmmawsem/helperFunctions/")
        from Pdb2GroLib import *
        # pdb_list = ["2ix5a"]
        pdb_list = list(set(pdb_list))
        print(len(pdb_list))
        for pdb in pdb_list:
            print(pdb)
            upperPDB = pdb[:4].upper()
            pre = f"/Users/weilu/openmmawsem/PDBs/"
            toPre = "frag_in_cbd"
            do(f"mkdir -p {toPre}/pdbs")
            do(f"mkdir -p {toPre}/frags")
            fromFile = f"{pre}/{upperPDB}.pdb"
            toFile = f"{toPre}/pdbs/cbd_{upperPDB}.pdb"
            convert_all_atom_pdb_to_cbd_representation(fromFile, toFile)
            # cmd = f"python ~/openmmawsem/helperFunctions/Pdb2Gro.py {toPre}/pdbs/cbd_{upperPDB}.pdb {toPre}/{pdb}.gro"
            pdbFile = toFile
            groFile = f"{toPre}/frags/{pdb}.gro"
            chainID = pdb[-1]
            Pdb2Gro(pdbFile, groFile, chainID.upper())
    if args.mode == 6:
        pre = "./"
        toPre = "frag_in_cbd"
        do(f"cp {pre}/frags.mem {pre}/cbd_frags.mem")
        replace(f"{pre}/cbd_frags.mem", "./fraglib/", f"{toPre}/frags/")

if args.day == "may02":
    if args.mode == 1:
        pdb_list = ["1r69", "1akr", "1opd", "1ptf", "1tig", "1tmy", "2acy", "5nul"]
        for pdb in pdb_list:
            pre = f"setups/{pdb}"
            fromFile = f"{pre}/crystal_structure.pdb"
            toFile = f"{pre}/cbd_{pdb}.pdb"
            convert_all_atom_pdb_to_cbd_representation(fromFile, toFile)
            cmd = f"python ~/openmmawsem/helperFunctions/Pdb2Gro.py {pre}/cbd_{pdb}.pdb {pre}/cbd_{pdb}.gro"
            do(cmd)
            do(f"cp {pre}/single_frags.mem {pre}/cbd_single_frags.mem")
            replace(f"{pre}/cbd_single_frags.mem", pdb, f"cbd_{pdb}")
    if args.mode == 2:
        # pre = "/Users/weilu/Research/server/feb_2020/casp13_targets/setups/T0953s2-D1"
        pdb_list = ["1r69", "1akr", "1opd", "1ptf", "1tig", "1tmy", "2acy", "5nul"]
        for pdb in pdb_list:
            pre = f"setups/{pdb}"
            original_openAWSEM_input = f"{pre}/{pdb}-openmmawsem.pdb"
            new_openAWSEM_input = f"{pre}/cbd-openmmawsem.pdb"
            all_atom_pdb_file = f"{pre}/crystal_structure-cleaned.pdb"
            replace_CB_coord_with_CBD_for_openAWSEM_input(original_openAWSEM_input, new_openAWSEM_input, all_atom_pdb_file)
if args.day == "may01":
    if args.mode == 1:
        pdb_list = ["1akr", "1opd", "1ptf", "1tig", "1tmy", "2acy", "5nul"]
        for pdb in pdb_list:
            folder = f"setups/{pdb}"
            do(f"mkdir -p {folder}")
            cd(folder)
            do(f"mm_create_project.py ../../source/{pdb}.pdb --extended --frag")
            cd("../..")
    if args.mode == 2:
        pdb_list = ["1akr", "1opd", "1ptf", "1tig", "1tmy", "2acy", "5nul"]
        for pdb in pdb_list:
            do(f"python mm_evaluate_native.py setups/{pdb}/{pdb} --to native/{pdb}")

if args.day == "apr29":
    if args.mode == 1:
        # pre = "/Users/weilu/Research/server/feb_2020/casp13_targets/setups/T0953s2-D1"
        pre = "./chain_E"
        pdb = "chainE"
        original_openAWSEM_input = f"{pre}/{pdb}-openmmawsem.pdb"
        new_openAWSEM_input = f"{pre}/cbd-openmmawsem.pdb"
        all_atom_pdb_file = f"{pre}/crystal_structure-cleaned.pdb"
        replace_CB_coord_with_CBD_for_openAWSEM_input(original_openAWSEM_input, new_openAWSEM_input, all_atom_pdb_file)
    if args.mode == 2:
        a = "/Users/weilu/Research/server/apr_2020/cornichon_cbd/chain_E/frags.mem"
        with open(a) as f:
            b = f.readlines()
        pdb_list = [c.split()[0].split("/")[-1].split(".")[0] for c in b[4:]]
        import sys
        sys.path.insert(0, "/Users/weilu/openmmawsem/helperFunctions/")
        from Pdb2GroLib import *
        # pdb_list = ["2ix5a"]
        pdb_list = list(set(pdb_list))
        print(len(pdb_list))
        for pdb in pdb_list:
            print(pdb)
            upperPDB = pdb[:4].upper()
            pre = f"/Users/weilu/openmmawsem/PDBs/"
            toPre = "/Users/weilu/Research/server/apr_2020/cornichon_cbd/frag_in_cbd"
            fromFile = f"{pre}/{upperPDB}.pdb"
            toFile = f"{toPre}/pdbs/cbd_{upperPDB}.pdb"
            convert_all_atom_pdb_to_cbd_representation(fromFile, toFile)
            # cmd = f"python ~/openmmawsem/helperFunctions/Pdb2Gro.py {toPre}/pdbs/cbd_{upperPDB}.pdb {toPre}/{pdb}.gro"
            pdbFile = toFile
            groFile = f"{toPre}/frags/{pdb}.gro"
            chainID = pdb[-1]
            Pdb2Gro(pdbFile, groFile, chainID.upper())


if args.day == "apr24":
    if args.mode == 1:
        relocate(fileLocation="frags.mem", toLocation="../fraglib")
        replace(f"frags.mem", "/Users/weilu/openmmawsem//Gros/", "../../fraglib/")
        protein_length = getFromTerminal("wc ssweight").split()[0]
        # print(f"protein: {name}, length: {protein_length}")
    if args.mode == 2:
        sys.path.insert(0, "/Users/weilu/openmmawsem")
        # name = "1iwg"
        import openmmawsem
        import helperFunctions.myFunctions
        for i in range(7):
            print(i)
            pdb = f"movie_frame_{i}.pdb"
            input_pdb_filename, cleaned_pdb_filename = openmmawsem.prepare_pdb(pdb, "T", keepIds=False)
            openmmawsem.ensure_atom_order(input_pdb_filename)
    if args.mode == 3:
        # do("echo 'ENDMODEL\nMODEL  {i+1}\n' >> openAWSEM_formatted_movie.pdb")
        do("rm openAWSEM_formatted_movie.pdb")
        for i in range(7):
            pdb = f"movie_frame_{i}-openmmawsem.pdb "
            do(f"cat {pdb} >> openAWSEM_formatted_movie.pdb")
            # do(f"sed '$d' {pdb} >> openAWSEM_formatted_movie.pdb")
            # do("echo 'ENDMODEL\nMODEL  {i+1}\n' >> openAWSEM_formatted_movie.pdb")
if args.day == "apr16":
    if args.mode == 1:
        # pdb = "6e67A"
        for pdb in dataset["hybrid"]:
            part = "globular"
            # part = "membrane"
            with open(f"plot_script/show_{pdb}_{part}.pml", "w") as out:
                # out.write(f"load native_{part}.pdb\n")
                a = f"{pdb}_best_{part}"
                out.write(f"load ../best_Q_structures/{a}.pdb\n")
                b = f"{pdb}_native_{part}"
                out.write(f"load ../native_structures/{b}.pdb\n")
                out.write(f'cealign {a}, {b}\ncmd.spectrum("count",selection="({a})&*/CA")\ncmd.spectrum("count",selection="({b})&*/CA")\n')
                out.write("orient\n")
                out.write("set ray_opaque_background, on\n")
                # out.write("ray 1000, 1000\n")
                # out.write(f"save {pdb}_{part}.png\n")
                # out.write("exit\n")
    if args.mode == 2:
        part = "globular"
        for pdb in dataset["hybrid"]:
            do(f"pymol show_{pdb}_{part}.pml")
    if args.mode == 3:
        # pdb = "6e67A"
        for pdb in dataset["hybrid"]:
            part = "membrane"
            # part = "membrane"
            with open(f"plot_script/show_{pdb}_{part}.pml", "w") as out:
                # out.write(f"load native_{part}.pdb\n")
                a = f"{pdb}_best_{part}"
                out.write(f"load ../best_Q_membrane_structures/{a}.pdb\n")
                b = f"{pdb}_native_{part}"
                out.write(f"load ../native_structures/{b}.pdb\n")
                out.write(f'cealign {a}, {b}\ncmd.spectrum("count",selection="({a})&*/CA")\ncmd.spectrum("count",selection="({b})&*/CA")\n')
                out.write("orient\n")
                out.write("set ray_opaque_background, on\n")
                # out.write("ray 1000, 1000\n")
                # out.write(f"save {pdb}_{part}.png\n")
                # out.write("exit\n")
    if args.mode == 4:
        part = "membrane"
        for pdb in dataset["hybrid"]:
            do(f"pymol show_{pdb}_{part}.pml")

        fromFile = f"{pre}/crystal_structure.pdb"
        toFile = f"{pre}/cbd_{pdb}.pdb"
        convert_all_atom_pdb_to_cbd_representation(fromFile, toFile)
        cmd = f"python ~/openmmawsem/helperFunctions/Pdb2Gro.py {pre}/cbd_{pdb}.pdb {pre}/cbd_{pdb}.gro"
        do(cmd)
        do(f"cp {pre}/single_frags.mem {pre}/cbd_single_frags.mem")
        replace(f"{pre}/cbd_single_frags.mem", pdb, f"cbd_{pdb}")

if args.day == "mar06":
    data = pd.read_csv("training_set.csv")
    specific_decoys = data.query("Length < 150 and Length > 70").reset_index(drop=True)
    pdb_list = specific_decoys["Protein"].to_list()
    pdb_list = [a.lower() for a in pdb_list]
    for pdb in pdb_list:
        pre = f"/Users/weilu/Research/server/mar_2020/mass_iterative_run/setups/{pdb}"
        fromFile = f"{pre}/crystal_structure.pdb"
        toFile = f"{pre}/cbd_{pdb}.pdb"
        convert_all_atom_pdb_to_cbd_representation(fromFile, toFile)
        cmd = f"python ~/openmmawsem/helperFunctions/Pdb2Gro.py {pre}/cbd_{pdb}.pdb {pre}/cbd_{pdb}.gro"
        do(cmd)
        do(f"cp {pre}/single_frags.mem {pre}/cbd_single_frags.mem")
        replace(f"{pre}/cbd_single_frags.mem", pdb, f"cbd_{pdb}")

if args.day == "feb25":
    if args.mode == 1:
        # pre = "/Users/weilu/Research/server/feb_2020/casp13_targets/setups/T0953s2-D1"
        pre = "."
        pdb = "1r69"
        original_openAWSEM_input = f"{pre}/{pdb}-openmmawsem.pdb"
        new_openAWSEM_input = f"{pre}/cbd-openmmawsem.pdb"
        all_atom_pdb_file = f"{pre}/crystal_structure-cleaned.pdb"
        replace_CB_coord_with_CBD_for_openAWSEM_input(original_openAWSEM_input, new_openAWSEM_input, all_atom_pdb_file)
if args.day == "feb20":
    if args.mode == 1:
        convert_all_atom_pdb_to_cbd_representation("crystal_structure.pdb", "cbd.pdb")
    if args.mode == 2:
        cmd = "python ~/opt/compute_energy.py lastFrame.pdb -g ../../../setups/1r69/iter0_gamma.dat -b ../../../setups/1r69/iter0_burial_gamma.dat"
        do(cmd)
    if args.mode == 3:
        pre = "/Users/weilu/Research/server/feb_2020/casp13_targets/setups/T0953s2-D1"
        pdb = "T0953s2-D1"
        original_openAWSEM_input = f"{pre}/{pdb}-openmmawsem.pdb"
        new_openAWSEM_input = f"{pre}/cbd-openmmawsem.pdb"
        all_atom_pdb_file = f"{pre}/crystal_structure-cleaned.pdb"
        replace_CB_coord_with_CBD_for_openAWSEM_input(original_openAWSEM_input, new_openAWSEM_input, all_atom_pdb_file)
    if args.mode == 4:
        cmd = "python mm_run.py setups/256b/cbd-openmmawsem.pdb --to native/256_cbd_submode_7 -s 1e3 --reportFrequency 100 -f forces_setup.py --fromOpenMMPDB --fasta crystal_structure.fasta --subMode 7"
        do("cmd")
if args.day == "feb17":
    pdb_list = ["T0951-D1", "T0953s2-D1", "T0955-D1", "T0957s1-D1", "T0957s1-D2", "T0958-D1", "T0960-D5", "T0963-D3", "T0968s1-D1", "T1008-D1"]
    if args.mode == 1:
        cleanPdb(pdb_list, chain="first")
    if args.mode == 2:
        for pdb in pdb_list:
            folder = f"setups/{pdb}"
            do(f"mkdir -p {folder}")
            cd(folder)
            do(f"mm_create_project.py ../../cleaned_pdbs/{pdb}.pdb --extended --frag")
            cd("../..")

if args.day == "feb16":
    pdb_list = dataset["may13"]
    if args.mode == 1:
        for pdb in pdb_list:
            folder = f"setups/{pdb}"
            do(f"mkdir -p {folder}")
            cd(folder)
            do(f"mm_create_project.py ../../cleaned_pdbs/previous_paper_pdbs/{pdb}.pdb --extended --frag")
            cd("../..")
    # pdb_list = "1R69, 1UZC, 1UTG, 3ICB, 1BG8, 1N2X, 256B, 4CPV, 1CCR, 2MHR, 1MBA, 2FHA".split(", ")
    # pdb_list = ['1r69', '3icb', '256b', '4cpv', '2mhr', '1mba', '2fha', '1fc2', '1enh', '2gb1', '2cro', '1ctf', '4icb']
    pdb_list = []
    pdb_list += ["1uzc", "1ccr", "1jwe"]
    # pdb_list += ["T0172_2"]
    if args.mode == 2:
        # pdb_list, steps = dataset["old"]
        downloadPdb(pdb_list)
        cleanPdb(pdb_list, chain=None)
    pdb_list = ["T0172_2"]
    if args.mode == 3:
        cleanPdb(pdb_list, chain="first")
    pdb_list = ["1uzc", "1ccr", "1jwe", "T0172_2"]
    if args.mode == 4:
        for pdb in pdb_list:
            folder = f"setups/{pdb}"
            do(f"mkdir -p {folder}")
            cd(folder)
            do(f"mm_create_project.py ../../cleaned_pdbs/{pdb}.pdb --extended --frag")
            cd("../..")
if args.day == "feb12":
    if args.mode == 1:
        pdb_list = ['1l0oC00']
        cleanPdb(pdb_list, source="my_CATH_database", toFolder="cleaned_pdbs", chain=-1, formatName=False)

if args.day == "feb06":
    pdb_list = ["6f45"]
    if args.mode == 1:
        for pdb in pdb_list:
            do(f"download.py {pdb}")
            do(f"mv {pdb} original_pdbs")
    if args.mode == 2:
        cleanPdb(pdb_list, source="original_pdbs", toFolder="cleaned_pdbs", chain="ABC", formatName=False)
    if args.mode == 3:
        for pdb in pdb_list:
            do(f"mkdir -p setups/{pdb}")
            cd(f"setups/{pdb}")
            do(f"mm_create_project.py ../../cleaned_pdbs/{pdb}.pdb --extended --frag")
            cd("../..")
    if args.mode == 4:
        for pdb in pdb_list:
            do(f"python mm_run.py setups/{pdb}/{pdb} --to native/{pdb} -s 2e5 --reportFrequency 2000 -f forces_setup.py --tempStart 300 --tempEnd 300")

if args.day == "feb03":
    fileName = "/Users/weilu/Research/optimization/chang_database/training_set.txt"
    pdb_list = []
    with open(fileName) as f:
        for line in f:
            pdbs = line.split()
            pdb_list += pdbs
    data = pd.read_csv("/Users/weilu/Research/optimization/chang_database/training_set.csv")
    specific_decoys = data.query("Length < 150 and Length > 70").reset_index(drop=True)
    pdb_list = specific_decoys["Protein"].to_list()
    pdb_list = [a.lower() for a in pdb_list]
    if args.mode == 1:
        do("mkdir -p training_set")
        cd("training_set")
        for pdb in pdb_list:
            do(f"download.py {pdb}")
    if args.mode == 2:
        cleanPdb(pdb_list, source="training_set/", toFolder="cleaned_pdbs/", chain="A", formatName=False, removeTwoEndsMissingResidues=True)
    if args.mode == 3:
        data = pd.read_csv("/Users/weilu/Research/optimization/chang_database/training_set.csv")
        specific_decoys = data.query("Length < 150 and Length > 70").reset_index(drop=True)
        pdb_list = specific_decoys["Protein"].to_list()
        pdb_list = [a.lower() for a in pdb_list]
        for pdb in pdb_list:
            do(f"mkdir -p setups/{pdb}")
            cd(f"setups/{pdb}")
            do(f"mm_create_project.py ../../cleaned_pdbs/{pdb}.pdb --extended")
            cd("../..")

if args.day == "jan20":
    if args.mode == 1:
        parser = PDBParser(QUIET=True)
        fileName = "/Users/weilu/Research/server/jan_2020/include_small_molecular/cleaned_pdbs/chain_ABCD.pdb"
        s = parser.get_structure("X", fileName)


        class ResSelect(Select):
            def accept_residue(self, residue):
                try:
                    if residue["CA"].get_coord()[2] < 32:
                        return 1
                    else:
                        return 0
                except:
                    if residue.id[0] == "H_ZK1":
                        return 0
                    else:
                        return 1

        io = PDBIO()
        io.set_structure(s)
        io.save('/Users/weilu/Research/server/jan_2020/include_small_molecular/cleaned_pdbs/ABCD_cutout_LBD.pdb', ResSelect())
    if args.mode == 2:
        do("pdb_selchain -E 20200120_171355_2866_6ud8.pdb > chainE.pdb")
        do("pdb_selchain -A,B,C,D 20200120_171355_2866_6ud8.pdb > chain_ABCD.pdb")
    if args.mode == 3:
        # pdb_list = ["chainE"]
        pdb_list = ["20200120_171355_2866_6ud8"]
        cleanPdb(pdb_list, source="./", addMissingResidues=True, toFolder="complete_chainE/", chain="E", formatName=False, removeTwoEndsMissingResidues=False, verbose=True)
    if args.mode == 4:
        # pdb_list = ["chainE"]
        pdb_list = ["combine_chain_ABCD_and_complete_chainE"]
        cleanPdb(pdb_list, source="./", removeHeterogens=False, addMissingResidues=False, toFolder="cleaned/", chain="ABCDE", formatName=False, removeTwoEndsMissingResidues=False, verbose=True)
        # then manual combine those two
if args.day == "jan18":
    pdb_list = dataset["optimization_v2"]
    if args.mode == 1:
        for pdb in pdb_list:
            do(f"python mm_run.py setups/{pdb}/{pdb} --to test/{pdb} -s 1e3 -f forces_setup_single_and_frag.py --reportFrequency 200 --subMode 0")
if args.day == "jan03":
    pdb_list = dataset["optimization_v2"]
    # cleanPdb(pdb_list, source="original_pdbs/first_test_set", toFolder="cleaned_pdbs/first_test_set", chain="A", formatName=False)
    if args.mode == 1:
        for pdb in pdb_list:
            do(f"mkdir -p setups/{pdb}")
            cd(f"setups/{pdb}")
            do(f"mm_create_project.py ../../cleaned_pdbs/first_test_set/{pdb}.pdb --extended --frag")
            cd("../..")
    if args.mode == 2:
        for pdb in pdb_list:
            do(f"python mm_run.py setups/{pdb}/{pdb} --to native/{pdb} -s 1e4 --reportFrequency 1000")
'''

########################################2019################################################
if args.day == "dec27":
    if args.mode == 1:
        pdb_list = ["6ud8_ABCD_noFill"]
        cleanPdb(pdb_list, source="original_pdbs/", addMissingResidues=False, toFolder="cleaned_pdbs/", chain="ABCD", formatName=False, removeTwoEndsMissingResidues=True)
    if args.mode == 2:
        pdb_list = ["6ud8_F_complete_not_embeded"]
        cleanPdb(pdb_list, source="original_pdbs/", addMissingResidues=True, toFolder="cleaned_pdbs/", chain="F", formatName=False, removeTwoEndsMissingResidues=False)
        do("pdb_rplchain -A:F 6ud8_F_complete_not_embeded.pdb > 6ud8_F_complete_not_embeded_asF.pdb")
    if args.mode == 3:
        pdb_list = ["my_ABCDF"]
        cleanPdb(pdb_list, source="original_pdbs/", addMissingResidues=False, toFolder="cleaned_pdbs/", chain="ABCDF", formatName=False, removeTwoEndsMissingResidues=True)
if args.day == "dec26":
    if args.mode == 1:
        # data = pd.read_csv("/Users/weilu/Research/server/dec_2019/optimization_database/for_optimization_dataset.csv", index_col=0)
        # pdb_list = data["Protein"].to_list()
        pdb_list = dataset["optimization_cath"]
        for pdb in pdb_list:
            do(f"mkdir -p {pdb}")
            cd(pdb)
            do(f"mm_create_project.py /Volumes/Wei_backup/cath_dataset/dompdb/{pdb}.pdb --extended")
            cd("..")
    if args.mode == 2:
        data = pd.read_csv("/Users/weilu/Research/server/dec_2019/optimization_database/for_optimization_dataset.csv", index_col=0)
        pdb_list = data["Protein"].to_list()
        for pdb in pdb_list:
            do(f"cp /Volumes/Wei_backup/cath_dataset/dompdb/{pdb}.pdb  .")
if args.day == "dec23":
    pdb_list = ["cut_out_LBD"]
    if args.mode == 1:
        cleanPdb(pdb_list, source="original_pdbs/", toFolder="cleaned_pdbs/", chain="ABCD", formatName=False, removeTwoEndsMissingResidues=True)

if args.day == "dec16":
    # pdb_list = ["6ud8_F_complete"]
    pdb_list = ["6ud8_ABCD"]
    if args.mode == 1:
        pdb_list = ["6ud8_ABCD"]
        cleanPdb(pdb_list, source="original_pdbs/", toFolder="cleaned_pdbs/", chain="ABCD", formatName=False, removeTwoEndsMissingResidues=True)
        # cleanPdb(pdb_list, source="original_pdbs/", toFolder="cleaned_pdbs/", chain="F", formatName=False, removeTwoEndsMissingResidues=False)
    if args.mode == 2:
        print("convert predictedZim")
        # name = "6ud8_F"
        name = "6ud8_F_complete"
        topo_name = f"TM_pred/{name}_topo"
        get_PredictedZim(topo_name, f"setups/{name}/PredictedZim")
        get_PredictedZimSide(topo_name, f"setups/{name}/PredictedZimSide")
    if args.mode == 5:
        name = "6ud8_F"
        name = "6ud8_F_complete"
        topo_name = f"TM_pred/crystal_topo_side"
        get_PredictedZim(topo_name, f"setups/{name}/PredictedZim_4TH")
        get_PredictedZimSide_v2(topo_name, f"setups/{name}/PredictedZimSide_4TH")
    if args.mode == 3:
        # pdb = "6ud8_F"
        pdb = "6ud8_F_complete"
        folder = "with_topology"
        do(f"mkdir -p {folder}")
        toFile = f"{folder}/forces_setup_{pdb}.py"
        fromFile = "forces_setup.py"
        do(f"cp {fromFile} {toFile}")
        pre = "."
        loc = f"{pre}/TM_pred/{pdb}_topo"
        with open(loc) as f:
            a = f.readlines()
        assert len(a) % 3 == 0
        chain_count = len(a) // 3
        seq = ""
        for i in range(chain_count):
            seq_i = (a[i*3+2]).strip()
            seq += seq_i
        assert np.alltrue([i in ["0", "1"] for i in seq])

        res_list = []
        first = None
        count = 1
        previousEnd = 0
        # print("g_all = [")
        out = "[\n"
        for i, res in enumerate(seq):
            if res == "0":
                if len(res_list) > 0:
                    # print(f"g{count} =", res_list)
                    print(res_list, ", ")
                    out += f"    {res_list},\n"
                    count += 1
                    last = res_list[-1]
                    first = res_list[0] if first is None else first
                    span = res_list[0] - previousEnd
                    previousEnd = last
                res_list = []
            if res == "1":
                res_list.append(i)
        out += "]\n"
        with fileinput.FileInput(toFile, inplace=True) as file:
            for line in file:
                tmp = line.replace("GALL", out).replace("FIRST", str(first)).replace("LAST", str(last))
                print(tmp, end='')
if args.day == "dec12":
    d = pd.read_csv("/Users/weilu/Research/server/dec_2019/iterative_optimization/original_pdbs/first_test_set.csv", index_col=0)
    pdb_list = d.PDB.str.lower().to_list()
    if args.mode == 1:
        cleanPdb(pdb_list, source="original_pdbs/first_test_set", toFolder="cleaned_pdbs/first_test_set", chain="A", formatName=False)
    if args.mode == 2:
        for pdb in pdb_list:
            do(f"mkdir -p setups/{pdb}")
            cd(f"setups/{pdb}")
            do(f"mm_create_project.py ../../cleaned_pdbs/first_test_set/{pdb}.pdb --extended")
            cd("../..")
    if args.mode == 3:
        for pdb in pdb_list:
            do(f"python mm_run.py setups/{pdb}/{pdb} --to native/{pdb} -s 1e4 --reportFrequency 1000")
if args.day == "dec11":
    if args.mode == 1:
        for pdb in dataset["may13"]:
            do(f"download.py {pdb}")
    if args.mode == 2:
        d = pd.read_csv("/Users/weilu/Research/server/dec_2019/iterative_optimization/original_pdbs/first_test_set.csv", index_col=0)
        for pdb in d.PDB:
            do(f"download.py {pdb}")
    if args.mode == 3:
        pdb_list = dataset["may13"]
        cleanPdb(pdb_list, source="original_pdbs/previous_paper_pdbs", toFolder="cleaned_pdbs/previous_paper_pdbs", chain="first", formatName=False)
    if args.mode == 4:
        pdb_list = dataset["may13"]
        for pdb in pdb_list:
            pdbToFasta(pdb, f"cleaned_pdbs/previous_paper_pdbs/{pdb}.pdb", f"cleaned_pdbs/previous_paper_pdbs/fasta/{pdb}.fasta", chains="A")

if args.day == "dec10":
    if args.mode == 1:
        # get unfinished job.
        # pre = "/Users/weilu/Research/server/oct_2019/membrane_optimization"
        cmd = f"grep 'CANCELLED' outs/slurm-*"
        has_cancelled = getFromTerminal(cmd)
        aa = list(set([b.split(":")[0] for b in has_cancelled.split("\n")]))
        print(aa)
        skip_pdbs = []
        for a in aa:
            cmd = f"grep -A1 'phi_pairwise_contact_well' {a}"
            c = getFromTerminal(cmd)
            if c == '':
                continue
            pdb = eval(c.split("\n")[1].strip())[0]
            skip_pdbs.append(pdb)
            # print(pdb)
        print(len(skip_pdbs))

if args.day == "dec02":
    if args.mode == 1:
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
if args.day == "nov29":
    pdb = "1hn4"
    pdb = "1lmm"
    pdb = "1tcg"
    pdb = args.label
    if args.mode == 4:
        for pdb in ["1tcg", "1lmm"]:
            do(f"gg.py -d nov29 -m 1 -l {pdb}")
            do(f"gg.py -d nov29 -m 2 -l {pdb}")
            do(f"gg.py -d nov29 -m 3 -l {pdb}")
    if args.mode == 1:
        folder = f"setups/{pdb}"
        do(f"mkdir -p {folder}")
        cd(folder)
        do(f"mm_create_project.py ../../original_pdbs/{pdb}.pdb --extended")
        cd("../..")
        do(f"cp gamma_noCysCys.dat {folder}/")
    if args.mode == 2:
        # ho frag memeory
        frag_folder = f"frag_memory/{pdb}"
        do(f"mkdir -p {frag_folder}")
        cd(frag_folder)
        do(f"cp ../../setups/{pdb}/{pdb}.fasta .")
        do(f"python ~/openmmawsem/helperFunctions/MultCha_prepFrags_index.py ~/openmmawsem/database/cullpdb_pc80_res3.0_R1.0_d160504_chains29712 {pdb}.fasta 20 2 9 > logfile")

        sys.path.insert(0, "/Users/weilu/openmmawsem")
        # name = "1iwg"
        import openmmawsem
        import helperFunctions.myFunctions
        helperFunctions.myFunctions.check_and_correct_fragment_memory("frags.mem")

        do(f"cp frags.mem ../../setups/{pdb}/ho.mem")
        cd("../..")
    if args.mode == 3:
        do(f"python mm_run.py setups/{pdb}/{pdb} --to native/{pdb}_single -f forces_setup_gamma_noCysCys.py --subMode 1 --reportFrequency 10 -s 20")
        do(f"python mm_run.py setups/{pdb}/{pdb} --to native/{pdb}_ho -f forces_setup_gamma_noCysCys.py --subMode 31 --reportFrequency 10 -s 20")
    # if args.mode == 2:
    #     do(f"cp frags.mem ../../setups/{pdb}/ho.mem")
    # if args.mode == 3:
    #     sys.path.insert(0, "/Users/weilu/openmmawsem")
    #     # name = "1iwg"
    #     import openmmawsem
    #     import helperFunctions.myFunctions
    #     helperFunctions.myFunctions.check_and_correct_fragment_memory("frags.mem")

if args.day == "nov26":
    if args.mode == 1:
        # ho frag memeory
        pdb = "1bpi"
        do("python ~/openmmawsem/helperFunctions/MultCha_prepFrags_index.py ~/openmmawsem/database/cullpdb_pc80_res3.0_R1.0_d160504_chains29712 {pdb}.fasta 20 2 9 > logfile")
    if args.mode == 2:
        do("cp frags.mem ../../setups/1bpi/ho.mem")
if args.day == "nov25":
    if args.mode == 1:
        sys.path.insert(0, "/Users/weilu/openmmawsem")
        # name = "1iwg"
        import openmmawsem
        import helperFunctions.myFunctions
        helperFunctions.myFunctions.check_and_correct_fragment_memory("frags.mem")

if args.day == "nov20":
    if args.mode == 1:
        # test new contact potential
        name = "6j14"
        # a3mFile = f"/Users/weilu/Research/server/may_2019/family_fold/{name}.a3m"
        pre = f"ff_contact_2/"
        os.system(f"mkdir -p {pre}")
        # data = get_MSA_data(a3mFile)
        data = ['QLQQSGAELVRPGASVKLSCKALGDTFTDYEIHWVKQTPVHGLEWIGVIHPGSGGTVYNQKFKGKATLTADKYSSTAYMELSSLTSEDSAVYYCTREGMNTDWYFDVWGAGTTVTVSDILMTQDELSLPVSLGDQASISCRSSQTIVHTNGNTYLEWYLQKPGQSPKLLIYKVSNRFSGVPDRFSGSGSGTYFTLKISRLEAEDLGVYYCFQGSHVPYTFGGGTKLEMKNPPTFSPALLVVTEGDNATFTCSFSNTSESFVLNWYRMSPSNQTDKLAAFPEDRSQPGQDCRFRVTQLPNGRDFHMSVVRARRNDSGTYLCGAISLAPKAQIKESLRAELRVTER']
        print(name, len(data))
        f_direct_2, f_water_2, f_protein_2, f_burial_2 = get_index_based_gamma(data, location=pre, gammaLocation="/Users/weilu/openmmawsem/parameters/", hasPhosphorylation=True)

if args.day == "nov13":
    pdb_list = ["1igd", "2sni", "1snb", "3il8", "1ubi", "1pht", "1poh", "1tig", "2acy", "1frd", "1opc", "1rds", "3chy", "5nul"]
    if args.mode == 1:
        for pdb in pdb_list:
            do(f"download.py {pdb}")
    if args.mode == 2:
        for pdb in pdb_list:
            folder = f"setups/{pdb}"
            do(f"mkdir -p {folder}")
            cd(folder)
            do(f"mm_create_project.py ../../original_pdbs/{pdb}.pdb --extended --frag")
            cd("../..")

if args.day == "nov10":
    if args.mode == 1:
        # do("download.py 1fs3")
        pdb = "1fs3"
        folder = f"setups/{pdb}"
        do(f"mkdir -p {folder}")
        cd(folder)
        do(f"mm_create_project.py ../../{pdb}.pdb --extended --frag")
        # do(f"mm_create_project.py ../../original/{pdb}.fasta")
        cd("../..")
if args.day == "nov09":
    # benchmark_list = ["6jdq", "6d3r"]
    # benchmark_list = ["6jdq", "5uak"]
    benchmark_list = ["5uak"]
    if args.mode == 1:
        for pdb in benchmark_list:
            folder = f"setups/{pdb}"
            do(f"mkdir -p {folder}")
            cd(folder)
            do(f"mm_create_project.py ../../original/{pdb}.pdb --verbose")
            # do(f"mm_create_project.py ../../original/{pdb}.fasta")
            cd("../..")

if args.day == "nov04":
    benchmark_list = ["fmt", "6n7n", "2is1"]
    if args.mode == 1:
        for pdb in benchmark_list:
            folder = f"setups/{pdb}"
            do(f"mkdir -p {folder}")
            cd(folder)
            do(f"mm_create_project.py ../../original/{pdb}.pdb")
            # do(f"mm_create_project.py ../../original/{pdb}.fasta")
            cd("../..")
if args.day == "oct30":
    fromPdb = "1a91.pdb"
    toPdb = "1a91_rotation_3.pdb"
    rotation_and_translation(fromPdb, toPdb, rotation_axis=(1,0,0), degree=0, translation=(0,0,10))

if args.day == "oct25":
    # benchmark_list = ["sos", "6g57", "5m2o", "unk2", "pex22", "6q64", "6btc", "1tt8"]
    # benchmark_list = ['6n7n_monomer_proteindna']   # include DNA
    benchmark_list = ['p53_tetramer_proteinonly', '6n7n_hexamer_proteinonly', '5xjy', '5uar', '5dqq']
    # benchmark_list = ['6n7n_hexamer_proteinonly', '5xjy', '5uar', '5dqq']
    if args.mode == 1:
        # a = generate_SEQRES("/Users/weilu/Research/server/oct_2019/benchmark/20191025_benchmark/1tt8.fasta")
        # print(a)
        for pdb in benchmark_list:
            a = generate_SEQRES(f"/Users/weilu/Research/server/oct_2019/benchmark/20191025_benchmark/original/{pdb}.fasta")
            with open(pdb+".pdb", "w") as o:
                o.write(a)
    if args.mode == 2:
        # a = generate_SEQRES("/Users/weilu/Research/server/oct_2019/benchmark/20191025_benchmark/1tt8.fasta")
        # print(a)
        for pdb in benchmark_list:
            folder = f"setups/{pdb}"
            do(f"mkdir -p {folder}")
            cd(folder)
            do(f"mm_create_project.py ../../original/{pdb}.pdb --extended")
            # do(f"mm_create_project.py ../../original/{pdb}.fasta")
            cd("../..")
    if args.mode == 3:
        for pdb in benchmark_list:
            # do(f"cp setups/{pdb}/cleaned_pdbs/{pdb}.pdb cleaned_pdbs/")
            # do(f"cp setups/{pdb}/crystal_structure-cleaned.pdb cleaned_pdbs/{pdb}.pdb")
            do(f"cp setups/{pdb}/{pdb}.fasta cleaned_pdbs/")
    if args.mode == 4:
        for pdb in benchmark_list:
            do(f"python mm_evaluate_native.py setups/{pdb}/{pdb} --to native/{pdb}")
if args.day == "oct19":
    if args.mode == 1:
        sys.path.insert(0, "/Users/weilu/openmmawsem")
        # name = "1iwg"
        import openmmawsem
        import helperFunctions.myFunctions
        for name in dataset["membrane"]:
            cd(f"setups/{name}")
            do(f"python3 /Users/weilu/openmmawsem/helperFunctions/fasta2pdb.py extended -f {name}.fasta")
            helperFunctions.myFunctions.add_chain_to_pymol_pdb("extended.pdb")  # only work for one chain only now
            input_pdb_filename, cleaned_pdb_filename = openmmawsem.prepare_pdb("extended.pdb", "A")
            openmmawsem.ensure_atom_order(input_pdb_filename)
            cd("../..")
    if args.mode == 2:
        for pdb in dataset["membrane"]:
            do(f"cp ../original_pdbs/{pdb}.pdb {pdb}_complete.pdb")
    if args.mode == 3:
        for pdb in dataset["membrane"]:
            print(pdb)
            do(f"python ~/opt/small_script/extract_membrane_part.py {pdb}_complete.pdb {pdb}_membrane.pdb")
    if args.mode == 4:
        sys.path.insert(0, "/Users/weilu/openmmawsem")
        # name = "1iwg"
        import openmmawsem
        import helperFunctions.myFunctions
        for name in dataset["membrane"]:
            print(name)
            cd(f"setups/{name}")
            do("grep -E 'CB|CA  GLY' crystal_structure-cleaned.pdb > cbs.data")
            do("""awk '{if($9>15) print "1"; else if($9<-15) print "3"; else print "2"}'  cbs.data  > zimPosition""")
            helperFunctions.myFunctions.create_zim(f"crystal_structure.fasta", tableLocation=f"/Users/weilu/openmmawsem/helperFunctions")
            cd("../..")

if args.day == "oct16":
    pdb_list = ["4nv6", "4p79", "5dsg", "6g7o", "6a93", "5xpd", "3kp9", "4a2n", "5d91", "2jo1"]
    pdb_list += ["1py6", "1pv6", "1u19"]
    pdb_list += ["2xov_complete", "6e67A"]
    pdb_list = ["6e67", "5xpd", "3kp9", "4a2n", "5d91", "4nv6", "4p79", "5dsg", "6g7o", "6a93", "1py6", "1pv6", "1u19"]
    if args.mode == 1:
        do("mkdir -p original_pdbs")
        for pdb in pdb_list:
            do(f"wget https://opm-assets.storage.googleapis.com/pdb/{pdb}.pdb")
            do(f"mv {pdb}.pdb original_pdbs/")
        cleanPdb(pdb_list, source="original_pdbs", chain=-1, formatName=False)
    if args.mode == 2:
        create_project_for_pdb_list(pdb_list, frag=True)

if args.day == "oct09":
    if args.mode == 1:
        # convert Porter5 format
        # pdb = "cannabinoid_receptor"
        # name = "serotonin_1A_receptor"
        # name = "cannabinoid_receptor"
        name = "beta_2_adrenergic_receptor"
        from_secondary = f"secondary/{name}/{name}_PROP/{name}.ss3"
        to_ssweight = f"{name}/setup/ssweight"
        print("convert ssweight")

        data = pd.read_csv(from_secondary, comment="#", names=["i", "Res", "ss3", "Helix", "Sheet", "Coil"], sep="\s+")
        # print(data)
        with open(to_ssweight, "w") as out:
            for i, line in data.iterrows():
                if line["ss3"] == "H":
                    out.write("1.0 0.0\n")
                if line["ss3"] == "E":
                    out.write("0.0 1.0\n")
                if line["ss3"] == "C":
                    out.write("0.0 0.0\n")

        fasta_file = f"../original_fasta_files/{name}.fasta"
        DMP_file = f"DMP/{name}/{name}.deepmetapsicov.con"
        convertDMPToInput(name, DMP_file, fasta_file)
        do(f"cp ~/opt/gremlin/protein/{name}/DMP/go_rnativeC* {name}/setup/")
        # name = "serotonin_1A_receptor"
        print("convert predictedZim")
        topo_name = f"TM_pred/{name}/{name}_topo"
        get_PredictedZim(topo_name, f"{name}/setup/PredictedZim")
        get_PredictedZimSide(topo_name, f"{name}/setup/PredictedZimSide")

if args.day == "oct08":
    if args.mode == 1:
        name = "1su4.pdb"
        print(name , computeRg(name))
    if args.mode == 2:
        do("python3 -m jpredapi submit --mode=single --format=fasta --file=cannabinoid_receptor.fasta ")
        do("python3 -m jpredapi status --results=jpred_sspred/results --extract --jobid=jp_05vaP_r")
        # jp_05vaP_r
    if args.mode == 3:
        for pdb in dataset["membrane"]:
            do(f"mkdir -p setups/{pdb}")
            cd(f"setups/{pdb}")
            do(f"mm_create_project.py ../../extracted/{pdb}.pdb")
            cd("../..")
    if args.mode == 4:
        # pdb = "6e67B"
        # pdb = args.label
        # pdb = "cannabinoid_receptor"
        brain_damage = 1
        do("mkdir -p frag_database")
        cd("frag_database")
        for pdb in dataset["membrane"]:
            do(f"mkdir {pdb}_HE")
            cd(f"{pdb}_HE")

            # brain_damage = 0
            # do("mkdir -p frag_database")
            # cd("frag_database")
            # do(f"mkdir {pdb}_HA")
            # cd(f"{pdb}_HA")

            # do(f"cp ../../setup/{pdb}/{pdb}.fasta .")
            # do(f"cp ../../../original_fasta_files/{pdb}.fasta .")
            do(f"cp ../../setups/{pdb}/{pdb}.fasta .")
            do("mkdir globularPart")
            cd("globularPart")
            do(f"python ~/openmmawsem/helperFunctions/MultCha_prepFrags_index.py ~/openmmawsem/database/cullpdb_pc80_res3.0_R1.0_d160504_chains29712 ../{pdb}.fasta 20 {brain_damage} 9 > logfile")
            check_and_correct_fragment_memory(fragFile="frags.mem")
            cd("..")

            do("mkdir membranePart")
            cd("membranePart")
            do(f"python ~/openmmawsem/helperFunctions/MultCha_prepFrags_index.py ~/openmmawsem/database/membrane_database/cullpdb_pc25_res3.0_R0.3_d190618_chains497 ../{pdb}.fasta 20 {brain_damage} 9 > logfile")
            check_and_correct_fragment_memory(fragFile="frags.mem")
            # do(f"python ~/openmmawsem/helperFunctions/MultCha_prepFrags_index.py ~/Research/optimization/fragment/self_culled/cullpdb_pc25_res3.0_R0.3_d190618_chains497 ../{pdb}.fasta 20 {brain_damage} 9 > logfile")
            cd("..")
            cd("..")
    if args.mode == 5:
        check_and_correct_fragment_memory(fragFile="frags.mem")
        # relocate(fileLocation="combined.mem", toLocation="fraglib")
        # replace(f"combined.mem", f"/Users/weilu/openmmawsem//Gros/", "./fraglib/")
    if args.mode == 6:
        # oct08
        for pdb in dataset["membrane"]:
            cmd = f"cp frag_database/{pdb}_HE/globularPart/frags.mem setups/{pdb}/globular.mem"
            do(cmd)
            cmd = f"cp frag_database/{pdb}_HE/membranePart/frags.mem setups/{pdb}/membrane.mem"
            do(cmd)
    if args.mode == 7:
        # oct08
        for pdb in dataset["membrane"]:
            print(pdb)
            cmd = f"python mm_evaluate_native.py setups/{pdb}/{pdb} --to native_contact_energy/{pdb} --subMode 0"
            do(cmd)
            # cmd = f"python mm_evaluate_native.py setups/{pdb}/{pdb} --to native_energy/{pdb}_globular --subMode 1"
            # do(cmd)
            # cmd = f"python mm_evaluate_native.py setups/{pdb}/{pdb} --to native_energy/{pdb}_membrane --subMode 2"
            # do(cmd)
if args.day == "oct05":
    # if args.day == "jun26":
    if args.mode == 1:
        # pdb = "6e67B"
        # pdb = args.label
        pdb = "cannabinoid_receptor"
        brain_damage = 1
        do("mkdir -p frag_database")
        cd("frag_database")
        do(f"mkdir {pdb}_HE")
        cd(f"{pdb}_HE")
        # brain_damage = 0
        # do("mkdir -p frag_database")
        # cd("frag_database")
        # do(f"mkdir {pdb}_HA")
        # cd(f"{pdb}_HA")

        # do(f"cp ../../setup/{pdb}/{pdb}.fasta .")
        do(f"cp ../../../original_fasta_files/{pdb}.fasta .")
        do("mkdir globularPart")
        cd("globularPart")
        do(f"python ~/openmmawsem/helperFunctions/MultCha_prepFrags_index.py ~/openmmawsem/database/cullpdb_pc80_res3.0_R1.0_d160504_chains29712 ../{pdb}.fasta 20 {brain_damage} 9 > logfile")
        cd("..")

        do("mkdir membranePart")
        cd("membranePart")
        do(f"python ~/openmmawsem/helperFunctions/MultCha_prepFrags_index.py ~/openmmawsem/database/membrane_database/cullpdb_pc25_res3.0_R0.3_d190618_chains497 ../{pdb}.fasta 20 {brain_damage} 9 > logfile")
        # do(f"python ~/openmmawsem/helperFunctions/MultCha_prepFrags_index.py ~/Research/optimization/fragment/self_culled/cullpdb_pc25_res3.0_R0.3_d190618_chains497 ../{pdb}.fasta 20 {brain_damage} 9 > logfile")
        cd("..")

    if args.mode == 2:
        check_and_correct_fragment_memory(fragFile="combined.mem")
        relocate(fileLocation="combined.mem", toLocation="fraglib")
        replace(f"combined.mem", f"/Users/weilu/openmmawsem//Gros/", "./fraglib/")
    if args.mode == 3:
        relocate(fileLocation="combined.mem", toLocation="fraglib")
    if args.mode == 4:
        replace(f"combined.mem", f"/Users/weilu/openmmawsem//Gros/", "./fraglib/")
    if args.mode == 5:
        probFile= "/Users/weilu/Research/server/jun_2019/simluation_hybrid/TM_pred/6e67A_PureTM/6e67A.prob"
        GlobularPart, MembranePart = get_two_part_from_prediction(probFile)
    if args.mode == 6:
        # crystal = "../5u09_clean.pdb"
        crystal = "4iar_clean.pdb"
        aligned_length, rmsd, tmscore, seqid = get_aligned_info(crystal, "0/lastFrame.pdb")
        print(aligned_length, rmsd, tmscore, seqid)
    if args.mode == 7:
        do("~/opt/TMalign/TMalign awsem.pdb {}.pdb -o result".format(protein_name))
        do("cp result_all_atm result_all_atm.pdb")
        do("pymol ~/opt/plot_scripts/tmalign_all.pml")
    if args.mode == 8:
        # name = "serotonin_1A_receptor"
        name = "cannabinoid_receptor"
        cryslal_look_up_data = {"serotonin_1A_receptor":"/Users/weilu/opt/gremlin/protein/serotonin_1A_receptor/DMP/4iar_clean.pdb",
                                "cannabinoid_receptor":"/Users/weilu/opt/gremlin/protein/cannabinoid_receptor/DMP/5tgz_clean.pdb"}
        template =cryslal_look_up_data[name]
        tm_info = open("tm_info.dat", "w")
        tm_info.write("i, aligned_length, rmsd, tmscore, seqid\n")
        for i in range(20):
            print(i)
            cd(f"{i}")
            aligned_length, rmsd, tmscore, seqid = get_aligned_info(template, "lastFrame.pdb")
            print(aligned_length, rmsd, tmscore, seqid)
            # template = "../../5u09_clean.pdb"
            target = "lastFrame.pdb"
            do(f"~/opt/TMalign/TMalign {template} {target} -o result > tm_log")
            do("cp result_all_atm result_all_atm.pdb")
            do("cp ~/opt/plot_scripts/tmalign_all.pml .")
            cd("..")
            tm_info.write(f"{i}, {aligned_length}, {rmsd}, {tmscore}, {seqid}\n")
            # do("pymol ~/opt/plot_scripts/tmalign_all.pml")
        tm_info.close()
    if args.mode == 9:
        #  gg.py -d oct05 -m 9
        import matplotlib.pyplot as plt
        # name = "cannabinoid_receptor"
        # n = 472
        # fileLocation = f"/Users/weilu/Research/server/oct_2019/draw_contact_for_DMP/{name}.deepmetapsicov.con"
        name = "serotonin_1A_receptor"
        n = 422
        # fileLocation = f"/Users/weilu/Research/server/oct_2019/GPCRs_reorder/simulation_setups/DMP/{name}/{name}.deepmetapsicov.con"
        fileLocation = f"/Users/weilu/opt/gremlin/protein/{name}/DMP/{name}.deepmetapsicov.con"

        a, _ = get_contactFromDMP(fileLocation, n)

        # pdbFile = "/Users/weilu/Research/server/oct_2019/draw_contact_for_DMP/5xr8_clean.pdb"
        # pdbFile = "../5u09_clean.pdb"
        cryslal_look_up_data = {"serotonin_1A_receptor":"/Users/weilu/opt/gremlin/protein/serotonin_1A_receptor/DMP/4iar_clean.pdb",
                                "cannabinoid_receptor":"/Users/weilu/opt/gremlin/protein/cannabinoid_receptor/DMP/5tgz_clean.pdb"}
        crystal =cryslal_look_up_data[name]
        data = getContactMapFromPDB(crystal, n)
        DMP_cutoff = 0.5
        # DMP_cutoff = 0.1
        # DMP_cutoff = 0.05
        do("mkdir -p figures")
        t_s = a.astype(float)
        plt.imshow(t_s, origin="bottom", cmap="Greys")
        plt.colorbar()
        plt.savefig(f"figures/DMP.png", dpi=300)
        plt.figure()
        for i in range(20):
            print(i)
            # pdbFile = "/Users/weilu/Research/server/oct_2019/draw_contact_for_DMP/lastFrame.pdb"
            pdbFile = f"{i}/lastFrame.pdb"
            lastFrame = getContactMapFromPDB(pdbFile, n)

            data_s = data.astype(float)
            t_s = lastFrame.astype(float)
            combined = data_s + t_s * 2
            upper = combined * np.tri(n, k=-1)

            data_s = data.astype(float)
            t_s = (a>DMP_cutoff).astype(float)
            combined = data_s + t_s * 2
            lower = combined * (1 - np.tri(n, k=0))
            combined = upper + lower
            from matplotlib import colors
            # red is in crystal but not in predicted.
            # blue is in predicted but not in crystal
            # black is in both
            cmap = colors.ListedColormap(['white', '#e6194B', '#4363d8', '#bfef45'])
            # cmap = colors.ListedColormap(['white', 'red', 'blue', '#bfef45'])

            bounds=[-1,0.1, 1.1, 2.1, 3.1]
            norm = colors.BoundaryNorm(bounds, cmap.N)
            plt.imshow(combined, origin="bottom", cmap=cmap, norm=norm)
            aligned_length, rmsd, tmscore, seqid = get_aligned_info(crystal, f"{i}/lastFrame.pdb")
            tm_line = f"{i}, Aligned: {aligned_length}, RMSD: {rmsd}, TMscore: {tmscore}, {seqid}\n"
            plt.title(tm_line+"Red is only in crystal, Blue is only in Predicted(upper), DMP in lower part")

            # plt.savefig("contact_5u09.png")
            plt.savefig(f"figures/contact_{i}_cutoff_{DMP_cutoff}.png", dpi=300)
        # plt.show()
    if args.mode == 11:
        for i in range(0, 20):
            print(i)
            cd(f"{i}")
            do("movie.py -m 5 123")
            cd("..")
    if args.mode == 12:
        do("cp -r serotonin_1A_receptor ~/Research/server/oct_2019/GPCRs_reorder/simulation_setups/")

if args.day == "oct02":
    if args.mode == 1:
        # using Porter5 to generate ss3 file.
        name = "serotonin_1A_receptor"

        # fastaFile = "beta_2_adrenergic_receptor.fasta"
        do("mkdir -p porter5")
        cd("porter5")
        sourceFasta = f"../../original_fasta_files/{name}.fasta"
        do(f"cp {sourceFasta} .")
        fastaFile = f"{name}.fasta"
        do(f"python ~/Documents/Porter5/Porter5.py -i {fastaFile} --cpu 2")

    if args.mode == 2:
        # setup simulation folder
        do("mkdir -p setup")
        cd("setup")
        fastaFile = "beta_2_adrenergic_receptor.fasta"
        do(f"mm_create_project.py ../{fastaFile} --hybrid")

if args.day == "sep28":
    if args.mode == 1:
        # convert Porter5 format
        # pdb = "cannabinoid_receptor"
        # pdb = "beta_2_adrenergic_receptor"
        # name = "serotonin_1A_receptor"
        name = "cannabinoid_receptor"
        from_secondary = f"secondary/{name}/{name}_PROP/{name}.ss3"
        to_ssweight = f"{name}/setup/ssweight"

        data = pd.read_csv(from_secondary, comment="#", names=["i", "Res", "ss3", "Helix", "Sheet", "Coil"], sep="\s+")
        # print(data)
        with open(to_ssweight, "w") as out:
            for i, line in data.iterrows():
                if line["ss3"] == "H":
                    out.write("1.0 0.0\n")
                if line["ss3"] == "E":
                    out.write("0.0 1.0\n")
                if line["ss3"] == "C":
                    out.write("0.0 0.0\n")
        cd("..")
        # cd("porter5")
        # data = pd.read_csv(from_secondary, sep="\t")
        # data = data.reset_index(drop=True).reset_index()
        # def chosen(data):
        #     top = max(data["Helix"], data["Sheet"], data["Coil"])
        #     if data["Helix"] == top:
        #         return 0
        #     elif data["Sheet"] == top:
        #         return 1
        #     else:
        #         return 2
        # data["chosen"] = data.apply(chosen, axis=1)
        # with open("ssweight", "w") as out:
        #     for i, line in data.iterrows():
        #         if line["chosen"] == 0:
        #             out.write("1 0\n")
        #         if line["chosen"] == 1:
        #             out.write("0 1\n")
        #         if line["chosen"] == 2:
        #             out.write("0 0\n")



    if args.mode == 2:
        # convert DMP result
        # pdbID = args.label
        name = "serotonin_1A_receptor"
        # fasta_file = f"{name}.fasta"
        fasta_file = f"../original_fasta_files/{name}.fasta"
        DMP_file = f"DMP/{name}/{name}.deepmetapsicov.con"
        convertDMPToInput(name, DMP_file, fasta_file)
        do(f"cp ~/opt/gremlin/protein/{name}/DMP/go_rnativeC* {name}/setup/")

if args.day == "sep20":
    if args.mode == 1:
        iteration = "trial_1"
        pre = "./trial_1/"
        do(f"mkdir -p {pre}")
        iter_gamma = np.loadtxt("/Users/weilu/Research/server/sep_2019/peptide_optimization/saved_gammas/mixed_original_and_cutoff100_impose_Aprime_constraint")
        gamma_for_simulation = pre + f"gamma.dat"
        burial_gamma_for_simulation = pre + f"burial_gamma.dat"
        gamma_format_convertion_iteration_to_simulation(iter_gamma, gamma_for_simulation, burial_gamma_for_simulation=burial_gamma_for_simulation)
        # do("mv original_*.dat for_simulation/")

if args.day == "aug20":
    if args.mode == 1:
        get_PredictedZim("4rws_topo", "PredictedZim")
if args.day == "aug19":
    if args.mode == 1:
        for pdb in dataset["hybrid"]:
            # do(f"pymol show_{pdb}_globular.pml")
            do(f"pymol show_{pdb}_membrane.pml")
    if args.mode == 2:
        infoLocation = "/Users/weilu/Research/database/hybrid_prediction_database/length_info.csv"
        info = pd.read_csv(infoLocation, index_col=0)
        # get_two_part_from_eye
        part_info = pd.read_csv("/Users/weilu/Research/database/hybrid_prediction_database/part_info.csv", names=["Protein", "Range"])
        part_info = part_info.merge(info, on="Protein")
        # pdb = "1pv6"
        folder = "set_topology_force"
        do(f"mkdir -p {folder}")
        for pdb in dataset["hybrid"]:
            GlobularPart, MembranePart = get_two_part_from_eye_seperation(pdb, part_info)
            toFile = f"{folder}/forces_setup_{pdb}.py"
            fromFile = "forces_setup.py"
            do(f"cp {fromFile} {toFile}")
            # pre = "/Users/weilu/Research/server/aug_2019/second_hybrid_protein_simulation/"
            pre = "."
            loc = f"{pre}/TM_pred/{pdb}_topo"
            with open(loc) as f:
                a = f.readlines()
            assert len(a) % 3 == 0
            chain_count = len(a) // 3
            seq = ""
            for i in range(chain_count):
                seq_i = (a[i*3+2]).strip()
                seq += seq_i
            assert np.alltrue([i in ["0", "1"] for i in seq])

            res_list = []
            first = None
            count = 1
            previousEnd = 0
            # print("g_all = [")
            out = "[\n"
            for i, res in enumerate(seq):
                if res == "0":
                    if len(res_list) > 0:
                        # print(f"g{count} =", res_list)
                        print(res_list, ", ")
                        out += f"    {res_list},\n"
                        count += 1
                        last = res_list[-1]
                        first = res_list[0] if first is None else first
                        span = res_list[0] - previousEnd
                        previousEnd = last
                    res_list = []
                if res == "1":
                    res_list.append(i)
            out += "]\n"
            with fileinput.FileInput(toFile, inplace=True) as file:
                for line in file:
                    tmp = line.replace("GALL", out).replace("FIRST", str(first)).replace("LAST", str(last))
                    tmp = tmp.replace("RESMEMB", f"{MembranePart}")
                    tmp = tmp.replace("RESGLOBULAR", f"{GlobularPart}")
                    print(tmp, end='')
if args.day == "aug16":
    if args.mode == 1:
        name = args.label
        check_and_correct_fragment_memory(fragFile=f"{name}.mem")
        relocate(fileLocation=f"{name}.mem", toLocation=f"fraglib_{name}")
        replace(f"{name}.mem", f"/Users/weilu/openmmawsem//Gros/", f"./fraglib_{name}/")

if args.day == "aug15":
    if args.mode == 1:
        pdb_list = ["1fs3"]
        for pdb in pdb_list:
            pdbID = pdb
            raptorX_file = f"contactmap.txt"
            convertRaptorToInput(pdbID, raptorX_file)
            print(pdb)
            do(f"wc ~/opt/gremlin/protein/{pdb}/raptor/go_rnativeCACA.dat")
            # do(f"wc ~/Research/server/jul_2019/hybrid_protein_simulation/setup/{pdb}/ssweight")
            print("--")
    if args.mode == 2:
        check_and_correct_fragment_memory(fragFile="new_frags_HA.mem")
        relocate(fileLocation="new_frags_HA.mem", toLocation="fraglib_new_HA")
        replace(f"new_frags_HA.mem", f"/Users/weilu/openmmawsem//Gros/", "./fraglib_new_HA/")

if args.day == "aug14":
    if args.mode == 1:
        # pdb = "6e67A"
        for pdb in dataset["hybrid"]:
            # part = "globular"
            part = "membrane"
            with open(f"plot_script/show_{pdb}_{part}.pml", "w") as out:
                # out.write(f"load native_{part}.pdb\n")
                a = f"{pdb}_best_{part}"
                out.write(f"load ../best_Q_structures/{a}.pdb\n")
                b = f"{pdb}_native_{part}"
                out.write(f"load ../native_structures/{b}.pdb\n")
                out.write(f'cealign {a}, {b}\ncmd.spectrum("count",selection="({a})&*/CA")\ncmd.spectrum("count",selection="({b})&*/CA")\n')
                out.write("orient\n")
    if args.mode == 2:
        part = "globular"
        for pdb in dataset["hybrid"]:
            do(f"pymol show_{pdb}_{part}.pml")

if args.day == "aug12":
    if args.mode == 1:
        check_and_correct_fragment_memory(fragFile="frags_HA.mem")
        relocate(fileLocation="frags_HA.mem", toLocation="fraglib_HA")
        replace(f"frags_HA.mem", f"/Users/weilu/openmmawsem//Gros/", "./fraglib_HA/")

if args.day == "aug11":
    if args.mode == 1:
        infoLocation = "/Users/weilu/Research/database/hybrid_prediction_database/length_info.csv"
        info = pd.read_csv(infoLocation, index_col=0)
        # get_two_part_from_eye
        part_info = pd.read_csv("/Users/weilu/Research/database/hybrid_prediction_database/part_info.csv", names=["Protein", "Range"])
        part_info = part_info.merge(info, on="Protein")
        # pdb = "1pv6"
        do("mkdir -p forces_recompute")
        for pdb in dataset["hybrid"]:
            GlobularPart, MembranePart = get_two_part_from_eye_seperation(pdb, part_info)
            forceLocation = "forces_recompute"
            do(f"cp forces/forces_setup_{pdb}.py {forceLocation}/forces_setup_{pdb}.py")
            with fileinput.FileInput(f"{forceLocation}/forces_setup_{pdb}.py", inplace=True) as file:
                for line in file:
                    fromLine = 'MembranePart ='
                    toLine = f'MembranePart ={MembranePart}\n#MembranePart ='
                    tmp = line.replace(fromLine, toLine)
                    fromLine = 'GlobularPart ='
                    toLine = f'GlobularPart ={GlobularPart}\n#GlobularPart ='
                    tmp = tmp.replace(fromLine, toLine)
                    print(tmp, end='')


if args.day == "aug10":
    def myRead(fileLocation, skip=0, strip=False):
        # read fragments.
        with open(fileLocation) as f:
            a = f.readlines()
        a = a[skip:]
        if strip:
            a = [i.strip() for i in a]
        return a
    if args.mode == 1:
        for pdb in dataset["hybrid"]:
            do(f"gg.py -d jun26 -m 1 -l {pdb}")
    if args.mode == 4:
        infoLocation = "/Users/weilu/Research/database/hybrid_prediction_database/length_info.csv"
        info = pd.read_csv(infoLocation, index_col=0)
        # get_two_part_from_eye
        part_info = pd.read_csv("/Users/weilu/Research/database/hybrid_prediction_database/part_info.csv", names=["Protein", "Range"])
        part_info = part_info.merge(info, on="Protein")

        for pdb in dataset["hybrid"]:
            cd(pdb+"_HA")
            GlobularPart, MembranePart = get_two_part_from_eye_seperation(pdb, part_info)
            # print(pdb)
            # print(GlobularPart)
            # print(MembranePart)
            fileLocation = "membranePart/frags.mem"
            mem = myRead(fileLocation, skip=4)
            fileLocation = "globularPart/frags.mem"
            globular = myRead(fileLocation, skip=4)
            mix_frag(globular, mem, GlobularPart, MembranePart)
            check_and_correct_fragment_memory(fragFile="HA_combined.mem")
            relocate(fileLocation="HA_combined.mem", toLocation="HA_frags")
            replace(f"HA_combined.mem", f"/Users/weilu/openmmawsem//Gros/", "./HA_frags/")
            cd("..")

#     if args.mode == 2:
#         # pdb = "6e67A"
#         pdb = args.label
#         fileLocation = "membranePart/frags.mem"
#         mem = myRead(fileLocation, skip=4)
#         fileLocation = "globularPart/frags.mem"
#         globular = myRead(fileLocation, skip=4)

#         probFile= f"/Users/weilu/Research/server/jun_2019/simluation_hybrid/TM_pred/{pdb}_PureTM/{pdb}.prob"
#         GlobularPart, MembranePart = get_two_part_from_prediction(probFile)

#         out = """[Target]
# query

# [Memories]
# """
#         n = len(mem)//20
#         for i in range(n):
#             if i in MembranePart:
#                 out += "".join(mem[i*20:(i+1)*20])
#             else:
#                 out += "".join(globular[i*20:(i+1)*20])

#         with open("HA_combined.mem", "w") as o:
#             o.write(out)

#     if args.mode == 3:
#         # pdb = "6e67A"
#         pdb = args.label
#         check_and_correct_fragment_memory(fragFile="HA_combined.mem")
#         relocate(fileLocation="HA_combined.mem", toLocation="HA_frags")
#         replace(f"HA_combined.mem", f"/Users/weilu/openmmawsem//Gros/", "./HA_frags/")
    # if args.mode == 11:
    #     cmd = "pymol "
    #     for pdb in pdb_list:
    #         cmd += f"setup/{pdb}/{pdb}.pdb "
    #     do(cmd)

if args.day == "jun26":
    if args.mode == 1:
        # pdb = "6e67B"
        pdb = args.label
        brain_damage = 1
        do("mkdir -p frag_database")
        cd("frag_database")
        do(f"mkdir {pdb}_HE")
        cd(f"{pdb}_HE")
        # brain_damage = 0
        # do("mkdir -p frag_database")
        # cd("frag_database")
        # do(f"mkdir {pdb}_HA")
        # cd(f"{pdb}_HA")

        do(f"cp ../../setup/{pdb}/{pdb}.fasta .")
        do("mkdir globularPart")
        cd("globularPart")
        do(f"python ~/openmmawsem/helperFunctions/MultCha_prepFrags_index.py ~/openmmawsem/database/cullpdb_pc80_res3.0_R1.0_d160504_chains29712 ../{pdb}.fasta 20 {brain_damage} 9 > logfile")
        cd("..")

        do("mkdir membranePart")
        cd("membranePart")
        do(f"python ~/openmmawsem/helperFunctions/MultCha_prepFrags_index.py ~/Research/optimization/fragment/self_culled/cullpdb_pc25_res3.0_R0.3_d190618_chains497 ../{pdb}.fasta 20 {brain_damage} 9 > logfile")
        cd("..")

    if args.mode == 2:
        check_and_correct_fragment_memory(fragFile="combined.mem")
        relocate(fileLocation="combined.mem", toLocation="fraglib")
        replace(f"combined.mem", f"/Users/weilu/openmmawsem//Gros/", "./fraglib/")
    if args.mode == 3:
        relocate(fileLocation="combined.mem", toLocation="fraglib")
    if args.mode == 4:
        replace(f"combined.mem", f"/Users/weilu/openmmawsem//Gros/", "./fraglib/")
    if args.mode == 5:
        probFile= "/Users/weilu/Research/server/jun_2019/simluation_hybrid/TM_pred/6e67A_PureTM/6e67A.prob"
        GlobularPart, MembranePart = get_two_part_from_prediction(probFile)

if args.day == "aug01":
    if args.mode == 1:
        pdb = "lysozyme"
        line = "lysozyme.pdb"
        pdbToFasta(pdb, line, f"{pdb}.fasta")
    if args.mode == 2:
        pdb_list = ["4rws"]
        for pdb in pdb_list:
            pdbID = pdb
            raptorX_file = f"{pdbID}.txt"
            convertRaptorToInput(pdbID, raptorX_file)
            print(pdb)
            do(f"wc ~/opt/gremlin/protein/{pdb}/raptor/go_rnativeCACA.dat")
            # do(f"wc ~/Research/server/jul_2019/hybrid_protein_simulation/setup/{pdb}/ssweight")
            print("--")
if args.day == "jul27":
    if args.mode == 1:
        # pdb_list = ["5d91", "2xov_complete", "3kp9", "4a2n"]
        pdb_list = ["2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91"]
        pdb_list = ["4nv6", "4p79", "5dsg", "6g7o", "6a93", "2jo1", "1py6", "1pv6", "1u19"]
        for pdb in pdb_list:
            do(f"cp ~/opt/gremlin/protein/{pdb}/raptor/go_rnativeC* {pdb}/")

if args.day == "jul22":
    # two chain assembly.
    pdb_list = ["6iu3"]
    if args.mode == 1:
        do("mkdir -p original_pdbs")
        for pdb in pdb_list:
            do(f"wget https://opm-assets.storage.googleapis.com/pdb/{pdb}.pdb")
            do(f"mv {pdb}.pdb original_pdbs/")
        cleanPdb(pdb_list, source="original_pdbs", chain=-1, formatName=False)
        create_project_for_pdb_list(pdb_list, frag=True)

if args.day == "jul21":
    pdb_list = ["4nv6", "4p79", "5dsg", "6g7o", "6a93", "5xpd", "3kp9", "4a2n", "5d91", "2jo1"]
    pdb_list += ["1py6", "1pv6", "1u19"]
    pdb_list += ["2xov_complete", "6e67A"]
    if args.mode == 1:
        do("mkdir -p original_pdbs")
        for pdb in pdb_list:
            do(f"wget https://opm-assets.storage.googleapis.com/pdb/{pdb}.pdb")
            do(f"mv {pdb}.pdb original_pdbs/")
        cleanPdb(pdb_list, source="original_pdbs", chain=-1, formatName=False)
    if args.mode == 2:
        create_project_for_pdb_list(pdb_list, frag=True)
if args.day == "jul16":
    from Bio.PDB import *

    class ExtractResidues(Select):

        def __init__(self, ResidueIndexGroup, resList):
            super(ExtractResidues, self).__init__()
            self.ResidueIndexGroup = ResidueIndexGroup
            self.resList = resList

        def accept_residue(self, residue):
            if self.resList.index(residue) in self.ResidueIndexGroup:
                return True
            else:
                return False

    def extractResidues(structure, toName, ResidueIndexGroup):
        resList = list(structure.get_residues())
        io = PDBIO()
        io.set_structure(structure)
        io.save(toName, ExtractResidues(ResidueIndexGroup, resList))

    if args.mode == 1:
        # pdb = "2xov_complete"
        pdb = args.label
        # probFile= f"/Users/weilu/Research/server/jun_2019/simluation_hybrid/TM_pred/{pdb}_PureTM/{pdb}.prob"
        probFile= f"/Users/weilu/Research/server/jul_2019/group1_hybrid_simulation/TM_pred/{pdb}_PureTM/{pdb}.prob"
        GlobularPart, MembranePart = get_two_part_from_prediction(probFile)
        if pdb == "2xov_complete":
            GlobularPart = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62]

        # print(GlobularPart, MembranePart)
        # fileLocation = "lastFrame.pdb"


        maxQList = glob.glob("max_Q*.pdb")
        for maxQ in maxQList:
            fileLocation = maxQ
            parser = PDBParser()
            structure = parser.get_structure('X', fileLocation)
            extractResidues(structure, f"{maxQ}_globular.pdb", GlobularPart)
            extractResidues(structure, f"{maxQ}_membrane.pdb", MembranePart)

        fileLocation = "crystal_structure.pdb"
        parser = PDBParser()
        structure = parser.get_structure('X', fileLocation)
        extractResidues(structure, "native_globular.pdb", GlobularPart)
        extractResidues(structure, "native_membrane.pdb", MembranePart)

        fileLocation = "lastFrame.pdb"
        parser = PDBParser()
        structure = parser.get_structure('X', fileLocation)
        extractResidues(structure, "lastFrame_globular.pdb", GlobularPart)
        extractResidues(structure, "lastFrame_membrane.pdb", MembranePart)

    if args.mode == 2:
        location = "info.dat"
        t = pd.read_csv(location, sep="\s+")
        # skip first two.
        t = t.query("Steps > 1").reset_index(drop=True)
        frame = t["Q_wat"].idxmax()

        location = "movie.pdb"
        with open(location) as f:
            a = f.readlines()
        n = len(a)
        # get the position of every model title
        model_title_index_list = []
        for i in range(n):
            if len(a[i]) >= 5 and a[i][:5] == "MODEL":
                model_title_index = i
                model_title_index_list.append(model_title_index)
        model_title_index_list.append(n)
        check_array = np.diff(model_title_index_list)
        if not np.allclose(check_array, check_array[0]):
            print("!!!! Someting is wrong  !!!!")
            print(check_array)
        else:
            size = check_array[0]
        with open(f"max_Q_wat_frame_{frame}.pdb", "w") as out:
            out.write("".join(a[size*frame:size*(frame+1)]))

    if args.mode == 3:
        # pdb = "6e67A"
        pdb = args.label
        part = "globular"
        with open(f"show_{pdb}_{part}.pml", "w") as out:
            out.write(f"load native_{part}.pdb\n")
            out.write(f"load lastFrame_{part}.pdb\n")
            out.write(f'cealign native_{part}, lastFrame_{part}\ncmd.spectrum("count",selection="(lastFrame_{part})&*/CA")\ncmd.spectrum("count",selection="(native_{part})&*/CA")\n')
            out.write("orient\n")

if args.day == "jul03":
    if args.mode == 1:
        pdb_list_A = ["4nv6", "4p79", "5dsg", "6g7o", "6a93"]
        pdb_list_B = ["4zyo", "5n6m"]
        pdb_list_C = ["4rws", "4xt3", "5uiw"]  # two chains.
        pdb_list_D = ["5uig"]
        pdb_list_E = ["2rh1"]
        pdb_list = pdb_list_A + pdb_list_B + pdb_list_C + pdb_list_D + pdb_list_E
        pdb_list = ["5uiw"]
        pdb_list = pdb_list_C
        do("mkdir -p original_pdbs")
        for pdb in pdb_list:
            do(f"wget https://opm-assets.storage.googleapis.com/pdb/{pdb}.pdb")
            do(f"mv {pdb}.pdb original_pdbs/")
        cleanPdb(pdb_list, source="original_pdbs", chain=-1, formatName=False)
    if args.mode == 2:
        # two chain, in connection.
        pdb_list_C = ["4rws", "4xt3", "5uiw"]
        pdb_list = ["4rws"]
        create_project_for_pdb_list(pdb_list)

if args.day == "jul02":
    pdb_list = ["2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91", "2jo1"]
    if args.mode == 1:
        # pdb = "5d91"
        # part = "membrane"
        # part = "globular"
        # pdb_list = ["5d91", "2xov_complete", "3kp9", "4a2n"]
        # pdb_list = ["2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91"]
        pdb_list = ["2xov_complete", "5xpd", "3kp9", "4a2n", "5d91"]
        part_list = ["globular", "membrane"]
        part_list = ["globular"]
        do("mkdir -p best_pml")
        for pdb in pdb_list:
            for part in part_list:
                if part == "globular":
                    loc = "wat"
                else:
                    loc = "mem"
                with open(f"best_pml/show_{part}_{pdb}.pml", "w") as out:
                    out.write(f"load {pdb}/0/native_{part}.pdb\n")
                    out.write(f"load Q_{loc}_max/{pdb}_best_{part}.pdb\n")
                    out.write(f'cealign native_{part}, {pdb}_best_{part}\ncmd.spectrum("count",selection="({pdb}_best_{part})&*/CA")\ncmd.spectrum("count",selection="(native_{part})&*/CA")\n')
                    out.write("orient\n")
                    out.write("ray 1000, 1000\n")
                    out.write(f"save snapshot/{part}_{pdb}.png")
                    # out.write("exit")
                cmd = f"pymol best_pml/show_{part}_{pdb}.pml"
                out = getFromTerminal(cmd)
                with open(f"log/{part}_{pdb}", "w") as o:
                    o.write(out)
                # print(out)

if args.day == "jul01":
    pdb_list = ["2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91", "2jo1"]
    if args.mode == 6:
        # pdb = "5d91"
        # part = "membrane"
        # part = "globular"
        # pdb_list = ["5d91", "2xov_complete", "3kp9", "4a2n"]
        pdb_list = ["2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91"]
        part_list = ["membrane", "globular"]
        do("mkdir -p pml")
        for pdb in pdb_list:
            for part in part_list:
                with open(f"pml/show_{pdb}_{part}.pml", "w") as out:
                    out.write(f"load {pdb}/0/native_{part}.pdb\n")
                    out.write(f"load {pdb}/0/lastFrame_{part}.pdb\n")
                    out.write(f'cealign native_{part}, lastFrame_{part}\ncmd.spectrum("count",selection="(lastFrame_{part})&*/CA")\ncmd.spectrum("count",selection="(native_{part})&*/CA")\n')
                    out.write("orient\n")
                # do("pymol 5d91/0/lastFrame_membrane.pdb 5d91/0/native_membrane.pdb")
    if args.mode == 5:
        # pdb_list = ["5d91", "2xov_complete", "3kp9", "4a2n"]
        pdb_list = ["2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91"]
        for pdb in pdb_list:
            cd(pdb+"/0")
            do("movie.py -m 4 123")
            do(f"gg.py -d jun29 -m 2 -l {pdb}")
            do(f"gg.py -d jun29 -m 1 -l {pdb}")
            cd("../..")
    if args.mode == 4:
        for pdb in pdb_list:
            with fileinput.FileInput(f"forces_setup_{pdb}.py", inplace=True) as file:
                for line in file:
                    tmp = line.replace("# er_term(oa)", "er_term(oa)")
                    print(tmp, end='')

    if args.mode == 3:
        # pdb_list = ["5d91", "2xov_complete", "3kp9", "4a2n"]
        pdb_list = ["2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91"]
        for pdb in pdb_list:
            do(f"cp ~/opt/gremlin/protein/{pdb}/raptor/go_rnativeC* {pdb}/")
    if args.mode == 1:
        for pdb in pdb_list:
            print(pdb)
            do(f"cat {pdb}/{pdb}.fasta")
    if args.mode == 2:

        # pdbID=args.label
        # raptorX_file = args.label
        pdb_list = ["4nv6", "4p79", "5dsg", "6g7o", "6a93", "2jo1", "1py6", "1pv6", "1u19"]
        for pdb in pdb_list:
            pdbID = pdb
            raptorX_file = f"{pdbID}/contactmap.txt"
            convertRaptorToInput(pdbID, raptorX_file)
            print(pdb)
            do(f"wc ~/opt/gremlin/protein/{pdb}/raptor/go_rnativeCACA.dat")
            do(f"wc ~/Research/server/jul_2019/hybrid_protein_simulation/setup/{pdb}/ssweight")
            print("--")
if args.day == "jun30":
    if args.mode == 1:
        pdb_list = ["lastFrame"]
        cleanPdb(pdb_list, source="./", chain=-1, formatName=False)

if args.day == "jun28":
    # time.sleep(60*60*2)
    pdb_list = ["2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91", "2jo1"]

if args.day == "jun27":
    pdb_list = ["5xpd", "3kp9", "4a2n", "5d91", "5uiw", "6akg", "4xt3"]
    if args.mode == 1:
        # pdb_list = ['1xrd', '5uiw', '3kp9', '3e9j', '2n7r', '5tcx', '4aw6', '5d91', '4zyo', '5ktf', '5mm0', '6bms']
        # pdb_list = ['2n7r', '1jo5', '5ktf', '2moz', '3zd0', '2lor', '2ksr', '3wkv', '4p79', '4a2n', '5tcx', '4b4a', '3wo7', '3tx3', '3ddl', '5jwy', '3kp9', '5xpd', '4zr1', '6gci', '6a2j', '4jr9', '4il3', '6ids', '5n6m', '6bug', '4dji']
        # a = pd.read_csv("/Users/weilu/Research/database/hybrid_prediction_database/picked2.csv", index_col=0)
        # pdb_list = a["Protein"].to_list()
        print(pdb_list)
        do("mkdir -p original_pdbs")
        for pdb in pdb_list:
            do(f"wget https://opm-assets.storage.googleapis.com/pdb/{pdb}.pdb")
            do(f"mv {pdb}.pdb original_pdbs/")
        cleanPdb(pdb_list, source="original_pdbs", chain=-1, formatName=False)
    if args.mode == 2:
        do("mkdir -p setup")
        cd("setup")
        for pdb in pdb_list:
            do(f"mkdir -p {pdb}")
            cd(pdb)
            do(f"mm_create_project.py ../../cleaned_pdbs/{pdb}.pdb --extended --membrane")
            cd("..")

# if args.day == "jun26":
#     if args.mode == 1:
#         # pdb = "6e67B"
#         pdb = args.label
#         # do("mkdir frag_database")
#         # cd("frag_database")
#         do(f"mkdir {pdb}_HA")
#         cd(f"{pdb}_HA")
#         do(f"cp ../../setup/{pdb}/{pdb}.fasta .")
#         do("mkdir globularPart")
#         cd("globularPart")
#         do(f"python ~/openmmawsem/helperFunctions/MultCha_prepFrags_index.py ~/openmmawsem/database/cullpdb_pc80_res3.0_R1.0_d160504_chains29712 ../{pdb}.fasta 20 0 9 > logfile")
#         cd("..")

#         do("mkdir membranePart")
#         cd("membranePart")
#         do(f"python ~/openmmawsem/helperFunctions/MultCha_prepFrags_index.py ~/Research/optimization/fragment/self_culled/cullpdb_pc25_res3.0_R0.3_d190618_chains497 ../{pdb}.fasta 20 0 9 > logfile")
#         cd("..")

#     if args.mode == 2:
#         check_and_correct_fragment_memory(fragFile="combined.mem")
#         relocate(fileLocation="combined.mem", toLocation="fraglib")
#         replace(f"combined.mem", f"/Users/weilu/openmmawsem//Gros/", "./fraglib/")
#     if args.mode == 3:
#         relocate(fileLocation="combined.mem", toLocation="fraglib")
#     if args.mode == 4:
#         replace(f"combined.mem", f"/Users/weilu/openmmawsem//Gros/", "./fraglib/")
#     if args.mode == 5:
#         probFile= "/Users/weilu/Research/server/jun_2019/simluation_hybrid/TM_pred/6e67A_PureTM/6e67A.prob"
#         GlobularPart, MembranePart = get_two_part_from_prediction(probFile)

if args.day == "jun23":
    if args.mode == 1:
        pdb_list = ["2mja", "2lep"]
        downloadPdb(pdb_list)
        # do("cp zimPosition back_zimPosition")

if args.day == "jun17":
    if args.mode == 3:
        name = "lastframe_m1.pdb"
        print(name , computeRg(name))
        name = "lastframe_im1.pdb"
        print(name , computeRg(name))
        name = "chainB.pdb"
        print(name , computeRg(name, chain="B"))
    if args.mode == 1:
        pdb_list = ["6e67A"]
        cleanPdb(pdb_list, source="original_pdbs", chain="A", formatName=False)
        pdb_list = ["6e67B"]
        cleanPdb(pdb_list, source="original_pdbs", chain="B", formatName=False)
    if args.mode == 2:
        # pdb_list = ["2xov"]
        # pdb_list = ["1occ", "1pv6", "2bl2", "2bg9", "1j4n", "1py6"]
        # pdb_list = ["2bg9"]
        pdb_list = ["6e67A", "6e67B"]
        do("mkdir -p setup")
        cd("setup")
        for pdb in pdb_list:
            do(f"mkdir -p {pdb}")
            cd(pdb)
            do(f"mm_create_project.py ../../cleaned_pdbs/{pdb}.pdb --extended --membrane --frag")
            cd("..")

if args.day == "jun15":
    if args.mode == 1:
        do("rm crystal_structure.pdb")
        to = "tmp1.pdb"

        duplicate_pdb("2xov_complete.pdb", to, offset_x=0, offset_y=0, offset_z=0, new_chain="A")
        do(f"cat {to} >> crystal_structure.pdb")
        duplicate_pdb("2xov_complete.pdb", to, offset_x=50, offset_y=0, offset_z=0, new_chain="B")
        do(f"cat {to} >> crystal_structure.pdb")

if args.day == "jun13":
    # pdb_list = ["2jo1"]
    pdb_list = dataset["membrane"]
    if args.mode == 1:
        # downloadPdb(pdb_list)
        # pdb_list = ["2bg9", "1j4n", "1py6", "2bl2", "1rhz", "1iwg", "2ic8", "1pv6", "1occ", "1kpl", "2bs2", "1py6", "1u19", "2rh1"]
        # pdb_list = ["2xov"]
        # pdb_list = ["2bs2", "1py6", "1u19", "2rh1"]
        do("mkdir -p original_pdbs")
        for pdb in pdb_list:
            do(f"wget https://opm-assets.storage.googleapis.com/pdb/{pdb}.pdb")
            do(f"mv {pdb}.pdb original_pdbs/")
        cleanPdb(pdb_list, source="original_pdbs", chain=-1, formatName=False)
    if args.mode == 2:
        do("mkdir -p setup")
        cd("setup")

        for pdb in pdb_list:
            do(f"mkdir -p {pdb}")
            cd(pdb)
            do(f"mm_create_project.py ../../cleaned_pdbs/{pdb}.pdb --extended --membrane --frag")
            cd("..")
if args.day == "jun08":
    if args.mode == 1:
        # a = glob.glob("cleaned_pdbs/*.pdb")
        a = glob.glob("dompdb/*.pdb")
        do("mkdir -p fasta")
        for line in a:
            print(line)
            pdb = line.split("/")[-1].split(".")[0]
            pdbToFasta(pdb, line, f"fasta/{pdb}.fasta")
            cmd = f"hhblits -i fasta/{pdb}.fasta -d ~/Research/Build/hh-suite/uniclust30_2018_08/uniclust30_2018_08 -o {pdb} -oa3m {pdb}.a3m -n 2"
            do(cmd)
            break

if args.day == "jun06":
    pdb_list = ["2bg9", "1j4n", "1py6_SD", "2bl2", "1rhz", "1iwg", "2ic8", "1pv6", "1occ", "1kpl", "2bs2", "1py6", "1u19"]
    if args.mode == 1:
        # downloadPdb(pdb_list)
        pdb_list = ["2bg9", "1j4n", "1py6", "2bl2", "1rhz", "1iwg", "2ic8", "1pv6", "1occ", "1kpl", "2bs2", "1py6", "1u19", "2rh1"]
        # pdb_list = ["2xov"]
        # pdb_list = ["2bs2", "1py6", "1u19", "2rh1"]
        do("mkdir -p original_pdbs")
        for pdb in pdb_list:
            do(f"wget https://opm-assets.storage.googleapis.com/pdb/{pdb}.pdb")
            do(f"mv {pdb}.pdb original_pdbs/")
    if args.mode == 2:
        from small_script.extract_pdb import *
        protein_info_list = []
        protein_info_list.append(("1py6_SD", "A",77, 199))
        protein_info_list.append(("1py6", "A",5, 231))
        protein_info_list.append(("2ic8", "A",91,272))
        protein_info_list.append(("1occ", "C",71,261))
        protein_info_list.append(("1pv6", "A", 1, 190))
        protein_info_list.append(("1j4n", "A",4,119))
        protein_info_list.append(("2bs2", "C", 21,237))
        protein_info_list.append(("2bl2", "A", 12, 156))
        protein_info_list.append(("2bg9", "A", 211, 301))
        protein_info_list.append(("1iwg", "A", 330, 497))
        protein_info_list.append(("1rhz", "A", 23, 188))
        protein_info_list.append(("1kpl", "A", 31, 233))
        protein_info_list.append(("1u19", "A", 33, 310))
        for (protein, chain, residue_start, residue_end) in protein_info_list:
            # do(f"wget {0}{1} -O ~/opt/crystal_structures/membrane_proteins/original_pdb/{1}".format(pdbFileDataBase, protein+".pdb"))
            # do(f"wget {0}{1} -O {1}".format(pdbFileDataBase, protein+".pdb"))
            extract_pdb("./", protein, chain, residue_start, residue_end)
        # cleanPdb(pdb_list, chain=-1, formatName=False)
        do("rm tmp.pdb")
    if args.mode == 3:
        cleanPdb(pdb_list, source="extracted", chain=-1, formatName=False)
    if args.mode == 4:
        for pdb in pdb_list:
            # do(f"getSize.py extracted/{pdb}.pdb")
            do(f"getSize.py cleaned_pdbs/{pdb}.pdb")
    if args.mode == 5:
        # pdb_list = ["2xov"]
        # pdb_list = ["1occ", "1pv6", "2bl2", "2bg9", "1j4n", "1py6"]
        # pdb_list = ["2bg9"]
        do("mkdir -p all_simulations")
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

            do(f"create_project.py {name} --membrane")

            # do(f"create_project.py {name} --membrane --frag")
            # # check_and_correct_fragment_memory("frags.mem")
            # relocate(fileLocation="frags.mem", toLocation="../fraglib")
            # replace(f"frags.mem", "/Users/weilu/openmmawsem//Gros/", "../../fraglib/")
            protein_length = getFromTerminal("wc ssweight").split()[0]
            print(f"protein: {name}, length: {protein_length}")

            with open("single_frags.mem", "w") as out:
                out.write("[Target]\nquery\n\n[Memories]\n")
                out.write(f"{name}.gro 1 1 {protein_length} 20\n")
            cd("../..")

        # for p in pdb_list:
        #     name = p.lower()[:4]
        #     do(f"mkdir -p {name}/{name}")
        #     cd(f"{name}/{name}")
        #     do("pwd")
        #     # do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")

        #     do("cp ~/opt/crystal_structures/membrane_proteins/for_simulation/{}.pdb crystal_structure.pdb ".format(name))
        #     do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
        #     do(f"create_project.py {name} --hybrid --frag")
        #     check_and_correct_fragment_memory(fragFile="frags.mem")
        #     relocate("./")

        #     # replace(f"frags.mem", "/Users/weilu/opt/script/Gros/", "../fraglib/")
        #     # replace(f"frags.mem", "../fraglib/", "../../fraglib/")
        #     replace(f"frags.mem", "/Users/weilu/opt/script/Gros/", "../../../fraglib/")
        #     do("cp frags.mem fragsLAMW.mem")
        #     protein_length = getFromTerminal("wc ssweight").split()[0]
        #     print(f"protein: {name}, length: {protein_length}")
        #     with open("single_frags.mem", "w") as out:
        #         out.write("[Target]\nquery\n\n[Memories]\n")
        #         out.write(f"{name}.gro 1 1 {protein_length} 20\n")
        #     cd("../..")
if args.day == "jun04":
    if args.mode == 1:
        for i in range(4):
            for j in range(i,4):
                do(f"show_gamma.py correct_ni_{i}_nj_{j} -o i{i}j{j}")


if args.day == "may25":
    # pdb_list = dataset["may13"]
    pdb_list = dataset["membrane"]
    if args.mode == 1:
        do("mkdir -p setup")
        cd("setup")

        for pdb in pdb_list:
            do(f"mkdir -p {pdb}")
            cd(pdb)
            do(f"mm_create_project.py ../../cleaned_pdbs/{pdb}.pdb --extended --membrane")
            cd("..")

    if args.mode == 2:
        pdb_list = dataset["may13"]
        for pdb in pdb_list:
            do(f"mkdir -p {pdb}")
            do(f"python mm_run.py setup/{pdb}/extended --to testRun/{pdb} -m 1 -s 200000")
            do(f"python mm_analysis.py setup/{pdb}/extended -t testRun/{pdb}/movie.dcd")
    if args.mode == 3:
        do(f"python mm_run.py setup/1r69/1r69 --to testRun/1r69_native -s 10")
        do(f"python mm_analysis.py setup/1r69/1r69 -t testRun/1r69_native/native.pdb -o testRun/1r69_native/native.dat")

if args.day == "may24":
    if args.mode == 1:
        a = pd.read_csv("/Users/weilu/Research/server/may_2019/multiSeq_iteration/chosen_may24.csv", index_col=0)
        pdb_list = a["FullName"].tolist()
        print(pdb_list)
        cleanPdb(pdb_list, chain="first", source="database/dompdb", verbose=True, addMissingResidues=False)

if args.day == "may20":
    if args.mode == 1:
        # test new contact potential
        name = "1r69"
        a3mFile = f"/Users/weilu/Research/server/may_2019/family_fold/{name}.a3m"
        pre = f"/Users/weilu/Research/server/may_2019/family_fold/ff_contact_2/{name}/"
        os.system(f"mkdir -p {pre}")
        data = get_MSA_data(a3mFile)
        print(name, len(data))
        f_direct_2, f_water_2, f_protein_2, f_burial_2 = get_ff_dat(data, location=pre, gammaLocation="/Users/weilu/Research/server/may_2019/openMM_2/1r69/")
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
# if args.day == "may10":
#     if args.mode == 1:
#         a = glob.glob("cleaned_pdbs/*.pdb")
#         do("mkdir -p fasta")
#         for line in a:
#             pdb = line.split("/")[-1].split(".")[0]
#             pdbToFasta(pdb, line, f"fasta/{pdb}.fasta")
#             cmd = f"hhblits -i fasta/{pdb}.fasta -d ~/Research/Build/hh-suite/uniclust30_2018_08/uniclust30_2018_08 -o {pdb} -oa3m {pdb}.a3m -n 1"
#             do(cmd)
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
# if args.day == "apr03":
#     if args.mode == 1:
#         # downloadPdb(pdb_list)
#         pdb_list = ["2bg9", "1j4n", "1py6", "2bl2", "1rhz", "1iwg", "2ic8", "1pv6", "1occ", "1kpl", "2bs2", "1py6", "1u19", "2rh1"]
#         pdb_list = ["2xov"]
#         pdb_list = ["2bs2", "1py6", "1u19", "2rh1"]
#         do("mkdir -p original_pdbs")
#         for pdb in pdb_list:
#             do(f"wget https://opm-assets.storage.googleapis.com/pdb/{pdb}.pdb")
#             do(f"mv {pdb}.pdb original_pdbs/")
#         cleanPdb(pdb_list, chain=-1, formatName=False)
#     if args.mode == 2:
#         # pdb_list = ["2xov"]
#         pdb_list = ["1occ", "1pv6", "2bl2", "2bg9", "1j4n", "1py6"]
#         do("mkdir -p all_simulations")
#         cd("all_simulations")
#         for p in pdb_list:
#             name = p.lower()[:4]
#             do(f"mkdir -p {name}/{name}")
#             cd(f"{name}/{name}")
#             do("pwd")
#             # do(f"cp ../../../cleaned_pdbs/{name}.pdb crystal_structure.pdb")

#             do("cp ~/opt/crystal_structures/membrane_proteins/for_simulation/{}.pdb crystal_structure.pdb ".format(name))
#             do(f"python2 ~/opt/script/Pdb2Gro.py crystal_structure.pdb {name}.gro")
#             do(f"create_project.py {name} --hybrid --frag")
#             check_and_correct_fragment_memory(fragFile="frags.mem")
#             relocate("./")

#             # replace(f"frags.mem", "/Users/weilu/opt/script/Gros/", "../fraglib/")
#             # replace(f"frags.mem", "../fraglib/", "../../fraglib/")
#             replace(f"frags.mem", "/Users/weilu/opt/script/Gros/", "../../../fraglib/")
#             do("cp frags.mem fragsLAMW.mem")
#             protein_length = getFromTerminal("wc ssweight").split()[0]
#             print(f"protein: {name}, length: {protein_length}")
#             with open("single_frags.mem", "w") as out:
#                 out.write("[Target]\nquery\n\n[Memories]\n")
#                 out.write(f"{name}.gro 1 1 {protein_length} 20\n")
#             cd("../..")
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

'''

'''
########----------------------2018-------------------------#######
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
