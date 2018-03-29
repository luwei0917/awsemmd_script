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
from small_script.variable_test import variable_test
import subprocess
from small_script.myFunctions import compute_theta_for_each_helix
from small_script.myFunctions import *
from small_script.myFunctions import make_metadata
# Useful codes
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# tr " " "\n"
# sed 1d
# sort -u -k 3
# sed -e 's/+T//'
# import re
# numbers = re.compile(r'(\d+)')
# def numericalSort(value):
#     parts = numbers.split(value)
#     parts[1::2] = map(int, parts[1::2])
#     return parts
# mypath = os.environ["PATH"]
# os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
# my_env = os.environ.copy()

parser = argparse.ArgumentParser(description="This is my playground for current project")
parser.add_argument("-r", "--run", help="test mode",
                    action="store_true")
parser.add_argument("-s", "--see", help="test mode",
                    action="store_true")
# parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("-m", "--mode", type=int, default=0)
parser.add_argument("-d", "--day", type=str, default="someday")
parser.add_argument("-t", "--test", action="store_true", default=False)
args = parser.parse_args()

if args.test:
    do = print
else:
    do = os.system
cd = os.chdir

base_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=01:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun {}\n'''

def replace(TARGET, FROM, TO):
    do("sed -i.bak 's/{}/{}/g' {}".format(FROM,TO,TARGET))

def getFromTerminal(CMD):
    return subprocess.Popen(CMD,stdout=subprocess.PIPE,shell=True).communicate()[0].decode()

def continueRunConvertion(n=12, rerun=0):
    rerun_plus_one = rerun + 1
    do(f"cp 2xov_0.in 2xov_{rerun_plus_one}.in")
    fileName = f"2xov_{rerun_plus_one}.in"
    replace(fileName, "variable r world "+ " ".join(str(i) for i in list(range(n))), "")
    replace(fileName, "# read_restart restart.25000000", "variable r world "+ " ".join(str(i) for i in list(range(n))))
    initial_steps = 20000000 * rerun_plus_one
    replace(fileName, "read_restart restart.extended", f"read_restart restart.$r.{initial_steps}")
    replace(fileName, "read_restart restart.native_topology", f"read_restart restart.$r.{initial_steps}")
    replace(fileName, "0\/", f"{rerun_plus_one}\/")
    cmd = 'tail -n 1 log.lammps | cut -d" " -f2-'
    line = getFromTerminal(cmd).rstrip()
    replace(fileName, "reset_timestep	0", "variable w world " + line)
    replace(fileName, "fix     xbias all colvars colvars.x output x", "fix     xbias all colvars colvars.x output x.$r")
    cmd = f'grep "temper" 2xov_{rerun_plus_one}.in'
    line = getFromTerminal(cmd).rstrip()
    replace(fileName, line, line + " $w")

def scancel_jobs_in_folder(folder):
    cd(folder)
    cmd = "find -name 'slurm-*' | rev | awk -F'[-.]' '{print $2}' | rev"
    lines = getFromTerminal(cmd).splitlines()
    for line in lines:
        print(line)
        do("scancel " + line)
    cd("..")

quick_slurm = '''#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:30:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun python3 ~/opt/davinci_gg.py -d mar13 -m 5
'''

localQ_slurm = '''#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:30:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun python3 ~/opt/davinci_gg.py -d mar13 -m 3
'''


if args.day == "mar13":
    if args.mode == 7:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 1
            i_plus_one = i +1
            # do(f"mkdir -p log{i}")
            # do(f"mv log.* log{i}/")
            # do(f"cp log{i}/log.lammps .")
            # do(f"cp x.* log{i}/")
            continueRunConvertion(n=12, rerun=i)
            do(f"mkdir {i_plus_one}")
            do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
    if args.mode == 6:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"

        freeEnergy_folder = f"third_combined_expectedDistance_freeEnergy/"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"second_rerun_{i}_12_Mar_211030" for i in range(2,4)]
        folder_list = [f"third_rerun_{i}_14_Mar_015209" for i in range(2,4)]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, sample_range_mode=-2, sub_mode_name=temp_mode, average_z=2, chosen_mode=0)


        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=-2)
                    nsample = len(folder_list)*2500
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 0 --nsample {nsample} --submode 2")
                    cd("..")
                cd("..")
        cd("..")
    if args.mode ==5:
        for i in range(12):
            compute_average_z_2(f"dump.lammpstrj.{i}", f"z_complete_{i}.dat")
    if args.mode == 4:
        print("compute localQ")
        # print(native_contacts_table)
        # cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # sim_list = ["0"]
        sim_list = ["0", "1"]
        for sim in sim_list:
            for folder in simulation_list:
                cd(folder)
                cd(sim)
                print(folder)
                with open("localQ.slurm", "w") as f:
                    f.write(localQ_slurm)
                    # f.write(localQ_slurm.replace("ctbp-common", "commons"))
                do("sbatch localQ.slurm")
                cd("../..")
    if args.mode == 3:
        native_contacts_table = compute_localQ_init()
        for i in range(12):
            compute_localQ(native_contacts_table, pre=".", ii=i)
    if args.mode == 2:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        sim_list = ["0", "1"]
        for sim in sim_list:
            for folder in simulation_list:
                cd(folder)
                cd(sim)
                print(folder)
                with open("computeZ.slurm", "w") as f:
                    f.write(quick_slurm)
                    # f.write(quick_slurm.replace("ctbp-common", "commons"))
                do("sbatch computeZ.slurm")
                cd("../..")
    if args.mode == 1:
        pre = "/scratch/wl45/"
        data_folder = "/scratch/wl45/all_data_folder/"
        folder_list = ["rg_0.1_lipid_1.0_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data_3(pre, data_folder, folder_list, rerun=1, average_z=True, localQ=True, label="third_")
if args.day == "mar09":
    if args.mode == 1:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 0
            i_plus_one = i +1
            do(f"mkdir -p log{i}")
            do(f"mv log.* log{i}/")
            do(f"cp log{i}/log.lammps .")
            do(f"cp x.* log{i}/")
            continueRunConvertion(n=12, rerun=i)
            do(f"mkdir {i_plus_one}")

            do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
if args.day == "mar05":
    if args.mode == 1:
        print("cp data files.")
        protein_list = ["T0766", "1mba", "T0784", "T0792", "T0803", "T0815", "T0833", "T0251"]
        for protein in protein_list:
            # cmd = f"cp -r /work/cms16/xl23/shared/IAAWSEM/MC_DATA_28Feb2018/{protein}/AWSEM_energy/AWSEM_energy.log {protein}_awsem.log"
            cmd = f"cp /work/cms16/xl23/shared/IAAWSEM/AWSEM_HO_Results/protein_pool/02282018/{protein.lower()}/iter/post-processing/noCST/qw/pca/lowTstructure/Qw.short.out {protein}_qw.txt"
            do(cmd)
            cmd = f"cp /work/cms16/xl23/shared/IAAWSEM/AWSEM_HO_Results/protein_pool/02282018/{protein.lower()}/iter/post-processing/noCST/qw/pca/lowTstructure/rmsd-angstrom.short.xvg {protein}_rmsd.txt"
            do(cmd)
if args.day == "mar03":
    if args.mode == 1:
        print("cp data files.")
        protein_list = ["T0766", "1mba", "T0784", "T0792", "T0803", "T0815", "T0833", "T0251"]
        for protein in protein_list:
            # cmd = f"cp -r /work/cms16/xl23/shared/IAAWSEM/MC_DATA_28Feb2018/{protein}/AWSEM_energy/AWSEM_energy.log {protein}_awsem.log"
            cmd = f"cp /work/cms16/xl23/shared/IAAWSEM/AWSEM_HO_Results/protein_pool/02282018/{protein.lower()}/iter/post-processing/noCST/qw/pca/lowTstructure/rwplusScore.short.txt {protein}_rw.txt"

            do(cmd)
    if args.mode == 2:
        a = pd.read_csv("/scratch/wl45/structure_selector_mar03/old_best_by_prediction.csv")
        a = a.assign(index=a.Step-1)
        for name, data in a.groupby("Name"):
            if name == "T0251":
                nn = "T251"
            else:
                nn = name
            print(name)
            do(f"mkdir {name}")
        #     print(data["index"])
            for i in data["index"]:
                do(f"cp /work/cms16/xl23/shared/IAAWSEM/MC_DATA_24Aug/{nn.upper()}/lowTstructure/lowTstructure{i}.pdb {name}/chosen_{i}.pdb")
    if args.mode == 3:
        a = pd.read_csv("/scratch/wl45/structure_selector_mar03/best_by_prediction_correction.csv")
        # a = a.assign(index=a.Step-1)
        for name, data in a.groupby("Name"):
            if name == "T0251":
                nn = "T251"
            else:
                nn = name
            print(name)
            do(f"mkdir {name}")
        #     print(data["index"])
            for i in data["index"]:
                do(f"cp /work/cms16/xl23/shared/IAAWSEM/MC_DATA_28Feb2018/{nn}/lowTstructure/lowTstructure{i}.pdb {name}/pca_chosen_{i}.pdb")
    if args.mode == 4:
        a = pd.read_csv("/scratch/wl45/structure_selector_mar03/best_by_prediction_based_on_new.csv")
        # a = a.assign(index=a.Step-1)
        for name, data in a.groupby("Name"):
            if name == "1MBA":
                nn = "1mba"
            else:
                nn = name
            print(name)
            do(f"mkdir {name}_new")
        #     print(data["index"])
            for i in data["index"]:
                do(f"cp /work/cms16/xl23/shared/IAAWSEM/MC_DATA_28Feb2018/{nn}/lowTstructure/lowTstructure{i}.pdb {name}_new/pca_chosen_{i}.pdb")
if args.day == "mar01":
    if args.mode == 1:
        cd("simulation")
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        sim_list = ["0", "1"]
        for sim in sim_list:
            for folder in simulation_list:
                cd(folder)
                cd(sim)
                print(folder)
                with open("computeZ.slurm", "w") as f:
                    f.write(quick_slurm)
                    # f.write(quick_slurm.replace("ctbp-common", "commons"))
                do("sbatch computeZ.slurm")
                cd("../..")
    if args.mode ==2:
        for i in range(12):
            compute_average_z(f"dump.lammpstrj.{i}", f"z_{i}.dat")
    if args.mode == 3:
        pre = "/scratch/wl45/"
        data_folder = "/scratch/wl45//all_data_folder/"
        folder_list = ["rg_0.2_lipid_1.0_mem_1"]
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_complete_temper_data(pre, data_folder, folder_list, rerun=1, average_z=True)
    if args.mode == 4:
        simulation_list = glob.glob("dis_*")
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 1
            i_plus_one = i +1
            do(f"mkdir -p log{i}")
            do(f"mv log.* log{i}/")
            do(f"cp log{i}/log.lammps .")
            do(f"cp x.* log{i}/")

            # continueRunConvertion(n=12, rerun=i)
            # do(f"mkdir {i_plus_one}")

            # run_slurm = base_run_slurm.format(i_plus_one)
            # with open(f"run_{i_plus_one}.slurm", "w") as r:
            #     r.write(run_slurm)
            # do(f"sbatch run_{i_plus_one}.slurm")

            # do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            # do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
    if args.mode == 5:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10", "1d_z":"12", "2d_z_qw":"13", "2d_z_dis":"14"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(3):
            freeEnergy_folder = f"freeEnergy_rg_0.1_lipid_1.0_mem_1_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["rg_0.2_lipid_1.0_mem_1"]
            # submode_list = ["_no_energy"]
            # submode_list = ["", "only_500"]
            # submode_list = ["350", "400", "450", "500", "550"]
            temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    move_data2(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=temp_mode, average_z=True)

            cd(freeEnergy_folder)
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    cd(folder+temp_mode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        make_metadata(temps_list=temp_list,k=0.02)
                        do("pulling_analysis.py -m {} --commons 0 --nsample 2500 --submode 5".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
if args.day == "dec03":
    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(3):
            freeEnergy_folder = f"dec02_no_side_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["no_side_contraint_memb_3_rg_0.4_lipid_0.6_extended"]
            # submode_list = ["_no_energy"]
            # submode_list = ["", "only_500"]
            # submode_list = ["350", "400", "450", "500", "550"]
            temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    move_data(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=temp_mode)

            cd(freeEnergy_folder)
            for temp_mode, temp_list in temp_dic.items():
                for folder in folder_list:
                    cd(folder+temp_mode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        for temp in temp_list:
                            make_metadata(temps_list=temp_list,k=0.02)
                            do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
if args.day == "dec02":
    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(3):
            freeEnergy_folder = f"dec02_no_side_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["no_side_contraint_memb_3_rg_0.4_lipid_0.6_extended"]
            # submode_list = ["_no_energy"]
            # submode_list = ["", "only_500"]
            submode_list = ["350", "400", "450", "500", "550"]
            for submode in submode_list:
                for folder in folder_list:
                    move_data(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=submode)

            cd(freeEnergy_folder)
            for submode in submode_list:
                for folder in folder_list:
                    cd(folder+submode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        for temp in temp_list:
                            make_metadata(temps_list=[submode],k=0.02)
                            do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
if args.day == "nov29":
    if args.mode == 3:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(3):
            freeEnergy_folder = f"nov29_no_side_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["no_side_contraint_memb_3_rg_0.4_lipid_0.6_extended"]
            # submode_list = ["_no_energy"]
            submode_list = ["", "only_500"]
            for submode in submode_list:
                for folder in folder_list:
                    move_data(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=submode)

            cd(freeEnergy_folder)
            for submode in submode_list:
                for folder in folder_list:
                    cd(folder+submode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        for temp in temp_list:
                            if submode == "only_500":
                                make_metadata(temps_list=[500],k=0.02)
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                            if submode == "_no_energy":
                                do("make_metadata.py -m 18 -k 0.05")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                            if submode == "":
                                do("make_metadata.py -m 18 -k 0.02")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                            if submode == "short":
                                do("make_metadata.py -m 21")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 3".format(mode))
                            elif submode == "low_t_":
                                do("make_metadata.py -m 20")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 4".format(mode))
                        cd("..")
                    cd("..")
            cd("..")

    if args.mode == 2:
        make_metadata(temps_list=[400])
    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"15", "1d_dis":"14", "1d_qw":"13"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(3):
            freeEnergy_folder = f"nov29_q_bias_temper_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["q_bias_temper_new"]
            # submode_list = ["_no_energy"]
            submode_list = [""]
            for submode in submode_list:
                for folder in folder_list:
                    move_data(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=submode, bias="qbias")

            cd(freeEnergy_folder)
            for submode in submode_list:
                for folder in folder_list:
                    cd(folder+submode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        for temp in temp_list:
                            if submode == "_no_energy":
                                do("make_metadata.py -m 18 -k 0.02")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                            if submode == "":
                                do("make_metadata.py -m 18 -k 1000")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                            if submode == "short":
                                do("make_metadata.py -m 21")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 3".format(mode))
                            elif submode == "low_t_":
                                do("make_metadata.py -m 20")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 4".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
if args.day == "nov28":
    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(1):
            freeEnergy_folder = f"nov28_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["bias_0.05_memb_3_rg_0.4_lipid_0.6_extended"]
            # submode_list = ["_no_energy"]
            submode_list = [""]
            for submode in submode_list:
                for folder in folder_list:
                    move_data(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=submode)

            cd(freeEnergy_folder)
            for submode in submode_list:
                for folder in folder_list:
                    cd(folder+submode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        for temp in temp_list:
                            if submode == "_no_energy":
                                do("make_metadata.py -m 18 -k 0.05")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                            if submode == "":
                                do("make_metadata.py -m 18 -k 0.05")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                            if submode == "short":
                                do("make_metadata.py -m 21")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 3".format(mode))
                            elif submode == "low_t_":
                                do("make_metadata.py -m 20")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 4".format(mode))
                        cd("..")
                    cd("..")
            cd("..")

if args.day == "nov27":
    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"15", "1d_dis":"14", "1d_qw":"13"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(1):
            freeEnergy_folder = f"q_bias_temper_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["q_bias_temper_new"]
            # submode_list = ["_no_energy"]
            submode_list = [""]
            for submode in submode_list:
                for folder in folder_list:
                    move_data(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=submode, bias="qbias")

            cd(freeEnergy_folder)
            for submode in submode_list:
                for folder in folder_list:
                    cd(folder+submode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        for temp in temp_list:
                            if submode == "_no_energy":
                                do("make_metadata.py -m 18 -k 0.02")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                            if submode == "":
                                do("make_metadata.py -m 18 -k 1000")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                            if submode == "short":
                                do("make_metadata.py -m 21")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 3".format(mode))
                            elif submode == "low_t_":
                                do("make_metadata.py -m 20")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 4".format(mode))
                        cd("..")
                    cd("..")
            cd("..")

if args.day == "nov21":
    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(1):
            freeEnergy_folder = f"no_side_contraint_memb_3_rg_0.4_lipid_0.6_extended_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["no_side_contraint_memb_3_rg_0.4_lipid_0.6_extended"]
            # submode_list = ["_no_energy"]
            submode_list = [""]
            for submode in submode_list:
                for folder in folder_list:
                    move_data(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=submode)

            cd(freeEnergy_folder)
            for submode in submode_list:
                for folder in folder_list:
                    cd(folder+submode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        for temp in temp_list:
                            if submode == "_no_energy":
                                do("make_metadata.py -m 18 -k 0.02")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                            if submode == "":
                                do("make_metadata.py -m 18 -k 0.02")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                            if submode == "short":
                                do("make_metadata.py -m 21")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 3".format(mode))
                            elif submode == "low_t_":
                                do("make_metadata.py -m 20")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 4".format(mode))
                        cd("..")
                    cd("..")
            cd("..")
if args.day == "nov19":
    if args.mode == 2:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(3):
            freeEnergy_folder = f"nov_18_all_freeEnergy_calculation_sample_range_mode_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["new_next_gen_native_based_memb_3_rg_0.4_lipid_0.6_extended"]
            submode_list = ["_no_energy"]
            for submode in submode_list:
                for folder in folder_list:
                    move_data(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=submode)


    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(3):
            freeEnergy_folder = f"nov_18_all_freeEnergy_calculation_sample_range_mode_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["new_next_gen_native_based_memb_3_rg_0.4_lipid_0.6_extended"]
            submode_list = ["short"]
            for submode in submode_list:
                for folder in folder_list:
                    move_data(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode, sub_mode_name=submode)

            cd(freeEnergy_folder)
            for submode in submode_list:
                for folder in folder_list:
                    cd(folder+submode)
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        for temp in temp_list:
                            if submode == "":
                                do("make_metadata.py -m 18 -k 0.02")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                            if submode == "short":
                                do("make_metadata.py -m 21")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 3".format(mode))
                            elif submode == "low_t_":
                                do("make_metadata.py -m 20")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 4".format(mode))
                        cd("..")
                    cd("..")
            cd("..")

if args.day == "nov18":
    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(3):
            freeEnergy_folder = f"nov_18_all_freeEnergy_calculation_sample_range_mode_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["new_next_gen_native_based_memb_3_rg_0.4_lipid_0.6_extended"]
            for folder in folder_list:
                move_data(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode)
            submode_list = [""]
            cd(freeEnergy_folder)
            for folder in folder_list:
                cd(folder)
                for submode in submode_list:
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = submode+bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        for temp in temp_list:
                            if submode == "":
                                do("make_metadata.py -m 18 -k 0.02")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                            if submode == "short":
                                do("make_metadata.py -m 19")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 3".format(mode))
                            elif submode == "low_t_":
                                do("make_metadata.py -m 20")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 4".format(mode))
                        cd("..")
                cd("..")
            cd("..")
    if args.mode == 2:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(3):
            freeEnergy_folder = f"nov_18_all_freeEnergy_calculation_sample_range_mode_{sample_range_mode}_2/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["new_next_gen_native_based_memb_3_rg_0.4_lipid_0.6_extended"]
            for folder in folder_list:
                move_data(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode)
            submode_list = [""]
            cd(freeEnergy_folder)
            for folder in folder_list:
                cd(folder)
                for submode in submode_list:
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = submode+bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        for temp in temp_list:
                            if submode == "":
                                do("make_metadata.py -m 18 -k 0.02")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                            if submode == "short":
                                do("make_metadata.py -m 19")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 3".format(mode))
                            elif submode == "low_t_":
                                do("make_metadata.py -m 20")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 4".format(mode))
                        cd("..")
                cd("..")
            cd("..")
if args.day == "nov15":
    if args.mode == 2:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        freeEnergy_folder = "nov_15_all_freeEnergy_calculation/"
        # do("mkdir " + freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "next_gen_native_based_memb_3_rg_0.2_lipid_0.6_extended",
        #                 "next_gen_native_based_memb_3_rg_0.4_lipid_0.6_extended",
        #                 "next_gen_native_based_memb_3_rg_0.4_lipid_0.6_topology",
        #                 "stronger_bias_for_expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_topology",
                        "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended",
                        "next_gen_native_based_memb_3_rg_0.2_lipid_0.6_extended",
                        "next_gen_native_based_memb_3_rg_0.4_lipid_0.6_extended",
                        "next_gen_native_based_memb_3_rg_0.4_lipid_0.6_topology",
                        "stronger_bias_for_expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        for folder in folder_list:
            move_data(data_folder, freeEnergy_folder, folder, sample_range_mode=0)
        submode_list = [""]
        cd(freeEnergy_folder)
        for folder in folder_list:
            cd(folder)
            for submode in submode_list:
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = submode+bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    for temp in temp_list:
                        if submode == "":
                            if folder == "stronger_bias_for_expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended":
                                do("make_metadata.py -m 18 -k 0.05")
                            else:
                                do("make_metadata.py -m 18 -k 0.02")
                            do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                    cd("..")
            cd("..")

    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        for sample_range_mode in range(3):
            freeEnergy_folder = f"nov_15_all_freeEnergy_calculation_sample_range_mode_{sample_range_mode}/"
            # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
            folder_list = ["memb_3_rg_0.1_lipid_1_topology"]
            for folder in folder_list:
                move_data(data_folder, freeEnergy_folder, folder, sample_range_mode=sample_range_mode)
            submode_list = [""]
            cd(freeEnergy_folder)
            for folder in folder_list:
                cd(folder)
                for submode in submode_list:
                    for bias, mode in bias_list.items():
                        # name = "low_t_" + bias
                        name = submode+bias
                        print(name)
                        do("rm -r "+name)
                        do("mkdir -p " + name)
                        cd(name)
                        for temp in temp_list:
                            if submode == "":
                                do("make_metadata.py -m 18 -k 0.02")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                            if submode == "short":
                                do("make_metadata.py -m 19")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 3".format(mode))
                            elif submode == "low_t_":
                                do("make_metadata.py -m 20")
                                do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 4".format(mode))
                        cd("..")
                cd("..")
            cd("..")
if args.day == "nov14":
    if args.mode == 2:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        freeEnergy_folder = "nov_14_all_freeEnergy_calculation/"
        folder_list = ["stronger_bias_for_expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        for folder in folder_list:
            move_data(data_folder, freeEnergy_folder, folder)
        submode_list = [""]
        cd(freeEnergy_folder)
        for folder in folder_list:
            cd(folder)
            for submode in submode_list:
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = submode+bias
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    for temp in temp_list:
                        if submode == "":
                            do("make_metadata.py -m 18 -k 0.05")
                            do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                        if submode == "short":
                            do("make_metadata.py -m 19")
                            do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 3".format(mode))
                        elif submode == "low_t_":
                            do("make_metadata.py -m 20")
                            do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 4".format(mode))
                    cd("..")
            cd("..")
    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        freeEnergy_folder = "all_freeEnergy_calculation_nov14/"
        folder_list = ["next_gen_native_based_memb_3_rg_0.4_lipid_0.6_topology"]
        for folder in folder_list:
            move_data(data_folder, freeEnergy_folder, folder)
        submode_list = [""]
        cd(freeEnergy_folder)
        for folder in folder_list:
            cd(folder)
            for submode in submode_list:
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = submode+bias
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    for temp in temp_list:
                        if submode == "":
                            do("make_metadata.py -m 18")
                            do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                        if submode == "short":
                            do("make_metadata.py -m 19")
                            do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 3".format(mode))
                        elif submode == "low_t_":
                            do("make_metadata.py -m 20")
                            do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 4".format(mode))
                    cd("..")
            cd("..")


if args.day == "nov12":
    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        freeEnergy_folder = "all_freeEnergy_calculation_nov11/"
        folder_list = ["next_gen_native_based_memb_3_rg_0.4_lipid_0.6_extended"]
        # for folder in folder_list:
        #     move_data(data_folder, freeEnergy_folder, folder)
        submode_list = ["low_t_"]
        cd(freeEnergy_folder)
        for folder in folder_list:
            cd(folder)
            for submode in submode_list:
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = submode+bias
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    for temp in temp_list:
                        if submode == "":
                            do("make_metadata.py -m 18")
                            do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                        if submode == "short":
                            do("make_metadata.py -m 19")
                            do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 3".format(mode))
                        elif submode == "low_t_":
                            do("make_metadata.py -m 20")
                            do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 4".format(mode))
                    cd("..")
            cd("..")

if args.day == "nov11":
    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        freeEnergy_folder = "all_freeEnergy_calculation_nov11/"
        folder_list = ["next_gen_native_based_memb_3_rg_0.2_lipid_0.6_extended"]
        for folder in folder_list:
            move_data(data_folder, freeEnergy_folder, folder)
        submode_list = [""]
        cd(freeEnergy_folder)
        for folder in folder_list:
            cd(folder)
            for submode in submode_list:
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = submode+bias
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    for temp in temp_list:
                        if submode == "":
                            do("make_metadata.py -m 18")
                            do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                        # elif submode == "low_t_":
                        #     do("make_metadata.py -m 16")
                        #     do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 3".format(mode))
                    cd("..")
            cd("..")
    if args.mode == 2:
        temp_list = ["all"]
        # bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        bias_list = {"1d_dis":"9"}
        data_folder = "all_data_folder/"
        freeEnergy_folder = "all_freeEnergy_calculation_nov11_2/"
        folder_list = ["next_gen_native_based_memb_3_rg_0.2_lipid_0.6_extended"]
        for folder in folder_list:
            move_data(data_folder, freeEnergy_folder, folder, krg=1)
        submode_list = [""]
        cd(freeEnergy_folder)
        for folder in folder_list:
            cd(folder)
            for submode in submode_list:
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = submode+bias
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    for temp in temp_list:
                        if submode == "":
                            do("make_metadata.py -m 19")
                            do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 3".format(mode))
                        # elif submode == "low_t_":
                        #     do("make_metadata.py -m 16")
                        #     do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 3".format(mode))
                    cd("..")
            cd("..")
if args.day == "nov09":
    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        freeEnergy_folder = "all_freeEnergy_calculation_nov09/"
        folder_list = ["expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        # for folder in folder_list:
        #     move_data(data_folder, freeEnergy_folder, folder, krg=5, klipid=0.8, kgo=0.5)
        submode_list = [""]
        cd(freeEnergy_folder)
        for folder in folder_list:
            cd(folder)
            for submode in submode_list:
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = submode+bias
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    for temp in temp_list:
                        if submode == "":
                            do("make_metadata.py -m 17")
                            do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                        elif submode == "low_t_":
                            do("make_metadata.py -m 16")
                            do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 3".format(mode))
                    cd("..")
            cd("..")
if args.day == "nov08":
    if args.mode == 2:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        # bias_list = {"2d_qw_dis":"11"}
        # bias_list = {"1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        freeEnergy_folder = "all_freeEnergy_calculation/"
        # folder = "rgWidth_memb_3_rg_0.1_lipid_1_extended"
        folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_topology",
                        "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        folder_list = ["expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
                        "rgWidth_memb_3_rg_0.1_lipid_1_topology",
                        "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        # for folder in folder_list:
        #     move_data(data_folder, freeEnergy_folder, folder)
        submode_list = ["", "low_t_"]
        cd("all_freeEnergy_calculation")
        for folder in folder_list:
            cd(folder)
            for submode in submode_list:
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = submode+bias
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    for temp in temp_list:
                        if submode == "":
                            do("make_metadata.py -m 17")
                            do("pulling_analysis.py -m {} --commons 0 --nsample 2500 --submode 2".format(mode))
                        elif submode == "low_t_":
                            do("make_metadata.py -m 16")
                            do("pulling_analysis.py -m {} --commons 0 --nsample 2500 --submode 3".format(mode))
                    cd("..")
            cd("..")
    if args.mode == 1:
        pre = "/scratch/wl45/nov_2017/06nov/"
        data_folder = "/scratch/wl45/nov_2017/06nov/all_data_folder/"
        folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
                        "rgWidth_memb_3_rg_0.1_lipid_1_topology",
                        "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        process_temper_data(pre, data_folder, folder_list)
if args.day == "nov07":
    if args.mode == 5:
        scancel_jobs_in_folder(".")
    if args.mode == 4:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        # bias_list = {"2d_qw_dis":"11"}
        # bias_list = {"1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        freeEnergy_folder = "all_freeEnergy_calculation/"
        # folder = "rgWidth_memb_3_rg_0.1_lipid_1_extended"
        folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_topology",
                        "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        folder_list = ["expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
                        "rgWidth_memb_3_rg_0.1_lipid_1_topology",
                        "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]
        # for folder in folder_list:
        #     move_data(data_folder, freeEnergy_folder, folder)
        submode_list = ["", "low_t_"]
        cd("all_freeEnergy_calculation")
        for folder in folder_list:
            cd(folder)
            for submode in submode_list:
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = submode+bias
                    do("mkdir -p " + name)
                    cd(name)
                    for temp in temp_list:
                        # do("make_metadata.py -m 16")
                        do("pulling_analysis.py -m {} --commons 0 --nsample 2500 --submode 2".format(mode))
                    cd("..")
            cd("..")
    if args.mode == 3:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        # bias_list = {"2d_qw_dis":"11"}
        # bias_list = {"1d_dis":"9", "1d_qw":"10"}
        data_folder = "all_data_folder/"
        freeEnergy_folder = "all_freeEnergy_calculation/"
        # folder = "rgWidth_memb_3_rg_0.1_lipid_1_extended"
        folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_topology",
                        "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]

        for folder in folder_list:
            move_data(data_folder, freeEnergy_folder, folder)
        cd("all_freeEnergy_calculation")
        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    do("make_metadata.py -m 15")
                    do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 2".format(mode))
                cd("..")
            cd("..")
    if args.mode == 2:
        data_folder = "all_data_folder/"
        freeEnergy_folder = "all_freeEnergy_calculation/"
        folder = "rgWidth_memb_3_rg_0.1_lipid_1_extended"
        move_data(data_folder, freeEnergy_folder, folder)

    if args.mode == 1:
        pre = "/scratch/wl45/nov_2017/06nov/"
        data_folder = "/scratch/wl45/nov_2017/06nov/all_data_folder/"
        folder_list = [
            "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"
        ]
        process_temper_data(pre, data_folder, folder_list)

if args.day == "nov05":
    if args.mode == 1:
        temp_list = ["all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}
        # bias_list = {"2d_qw_dis":"11"}
        # bias_list = {"1d_dis":"9", "1d_qw":"10"}
        for bias, mode in bias_list.items():
            do("mkdir -p " + bias)
            cd(bias)
            for temp in temp_list:
                do("make_metadata.py -m 12")
                do("mkdir t_" + temp)
                do("mv metadatafile t_" + temp)
                cd("t_" + temp)
                do("pulling_analysis.py -m {} --commons 0 --nsample 2500 --submode 1".format(mode))
                cd("..")
            cd("..")
        # cmd = "find -name 'slurm-*' | rev | awk -F'[-.]' '{print $2}' | rev"
        # lines = getFromTerminal(cmd).splitlines()
        # for line in lines:
        #     print(line)
        #     do("scancel " + line)

if args.day == "oct31":
    if args.mode == 4:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        folder_list = ["memb_3_rg_0.1_lipid_1_extended", "memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended", "rgWidth_memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "500"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    if temp != "all":
                        do("make_metadata.py -m 10 --submode 3 -t " + temp)
                    else:
                        do("make_metadata.py -m 9 --submode 3")
                    do("mkdir t_" + temp)
                    do("mv metadatafile t_" + temp)
                    cd("t_" + temp)
                    do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 1".format(mode))
                    cd("..")
                cd("..")
            cd("..")
    if args.mode == 3:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        # folder_list = ["memb_2_rg_0.1_lipid_1_extended", "memb_2_rg_0.1_lipid_1_topology"]
        folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended", "rgWidth_memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "500"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    if temp != "all":
                        do("make_metadata.py -m 10 --submode 3 -t " + temp)
                    else:
                        do("make_metadata.py -m 9 --submode 3")
                    do("mkdir t_" + temp)
                    do("mv metadatafile t_" + temp)
                    cd("t_" + temp)
                    do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 1".format(mode))
                    cd("..")
                cd("..")
            cd("..")
    if args.mode == 2:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        # folder_list = ["memb_2_rg_0.1_lipid_1_extended", "memb_2_rg_0.1_lipid_1_topology"]
        folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended", "rgWidth_memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "500"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("rm -r " + bias)
            cd("..")
    if args.mode == 1:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        # folder_list = ["memb_2_rg_0.1_lipid_1_extended", "memb_2_rg_0.1_lipid_1_topology"]
        folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended", "rgWidth_memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "500"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    if temp != "all":
                        do("make_metadata.py -m 10 --submode 1 -t " + temp)
                    else:
                        do("make_metadata.py -m 9 --submode 1")
                    do("mkdir t_" + temp)
                    do("mv metadatafile t_" + temp)
                    cd("t_" + temp)
                    do("pulling_analysis.py -m {} --commons 0 --nsample 2500 --submode 1".format(mode))
                    cd("..")
                cd("..")
            cd("..")

if args.day == "oct28":
    if args.mode == 1:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        # folder_list = ["memb_2_rg_0.1_lipid_1_extended", "memb_2_rg_0.1_lipid_1_topology"]
        folder_list = ["memb_3_rg_0.1_lipid_1_extended", "memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "500"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    if temp != "all":
                        do("make_metadata.py -m 10 --submode 1 -t " + temp)
                    else:
                        do("make_metadata.py -m 9 --submode 1")
                    do("mkdir t_" + temp)
                    do("mv metadatafile t_" + temp)
                    cd("t_" + temp)
                    do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 1".format(mode))
                    cd("..")
                cd("..")
            cd("..")

if args.day == "oct25":
    if args.mode == 1:
        compute_theta_for_each_helix()
        print("Done")
    if args.mode == 2:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        # folder_list = ["memb_2_rg_0.1_lipid_1_extended", "memb_2_rg_0.1_lipid_1_topology"]
        folder_list = ["memb_3_rg_0.1_lipid_1_extended", "memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "500"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    if temp != "all":
                        do("make_metadata.py -m 10 -t " + temp)
                    else:
                        do("make_metadata.py -m 9")
                    do("mkdir t_" + temp)
                    do("mv metadatafile t_" + temp)
                    cd("t_" + temp)
                    do("pulling_analysis.py -m {} --commons 1 --nsample 2500 --submode 1".format(mode))
                    cd("..")
                cd("..")
            cd("..")
if args.day == "oct24":
    if args.mode == 2:
        cmd = 'tail -n 1 log.lammps | cut -d" " -f2-'
        line = getFromTerminal(cmd).rstrip()
        print(line)
    if args.mode == 1:
        dis_list = glob.glob("dis_*")
        print(dis_list)
        for dis in dis_list:
            cd(dis)
            do("cp log0/log.lammps .")
            continueRunConvertion()
            # do("mkdir log0")
            # do("mv log.* log0/")
            # do("mkdir 1")

            do("sed 's/2xov_0/2xov_1/g' run_0.slurm > run_1.slurm")
            do("sbatch run_1.slurm")
            cd("..")
        # continueRunConvertion()

if args.day == "oct21":
    if args.mode == 1:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        # folder_list = ["memb_2_rg_0.1_lipid_1_extended", "memb_2_rg_0.1_lipid_1_topology"]
        folder_list = ["memb_3_rg_0.1_lipid_1_extended", "memb_3_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "450"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    if temp != "all":
                        do("make_metadata.py -m 10 -t " + temp)
                    else:
                        do("make_metadata.py -m 9")
                    do("mkdir t_" + temp)
                    do("mv metadatafile t_" + temp)
                    cd("t_" + temp)
                    do("pulling_analysis.py -m {} --commons 0 --nsample 2500".format(mode))
                    cd("..")
                cd("..")
            cd("..")
    if args.mode == 2:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        folder_list = ["memb_2_rg_0.1_lipid_1_extended", "memb_2_rg_0.1_lipid_1_topology"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["important_all"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    do("make_metadata.py -m 9")
                    do("rm t_t_" +temp)
                    do("mkdir -p t_t_" + temp)
                    do("mv metadatafile t_t_" + temp)
                    cd("t_t_" + temp)
                    do("pulling_analysis.py -m {} --commons 0 --nsample 2500".format(mode))
                    cd("..")
                cd("..")
            cd("..")
if args.day == "oct17":
    if args.mode == 1:
        print("can it fold")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 1, 2]
        force_ramp_rate_list=[1]
        temperature_list=[350, 500]
        memb_k_list = [1, 2, 4]
        rg_list = [0, 0.1, 0.4]
        force_list = [0.0]
        repeat = 2
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=1)
if args.day == "oct16":
    if args.mode == 1:
        # folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        folder_list = ["more_higher_temp"]
        # folder_list = ["more_higher_temp_topology"]
        temp_list = ["all", "350", "400", "450", "500"]
        bias_list = {"2d_qw_dis":"11", "1d_dis":"9", "1d_qw":"10"}

        for folder in folder_list:
            cd(folder)
            for bias, mode in bias_list.items():
                do("mkdir -p " + bias)
                cd(bias)
                for temp in temp_list:
                    if temp != "all":
                        do("make_metadata.py -m 10 -t " + temp)
                    else:
                        do("make_metadata.py -m 9")
                    do("mkdir t" + temp)
                    do("mv metadatafile t" + temp)
                    cd("t" + temp)
                    do("pulling_analysis.py -m {} --commons 1 --nsample 4000".format(mode))
                    cd("..")
                cd("..")
            cd("..")
if args.day == "oct13":
    if args.mode == 2:
        print("strong membrane fore ramp")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [1, 2]
        force_ramp_rate_list=[1]
        temperature_list=[350]
        memb_k_list = [1, 2, 4, 8, 16]
        rg_list = [0, 0.4]
        force_list = [0.0]
        repeat = 10
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=0)
    if args.mode == 1:
        folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        temp_list = ["all", "350", "400", "450", "500"]
        for folder in folder_list:
            cd(folder)

            do("mkdir -p 2d_qw_dis")
            cd("2d_qw_dis")
            for temp in temp_list:
                if temp != "all":
                    do("make_metadata.py -m 10 -t " + temp)
                else:
                    do("make_metadata.py -m 9")
                do("mkdir t" + temp)
                do("mv metadatafile t" + temp)
                cd("t" + temp)
                do("pulling_analysis.py -m 11 --commons 1")
                cd("..")
            cd("..")
            cd("..")
if args.day == "oct12":
    if args.mode == 1:
        temp_list = ["400", "450", "500", "all"]
        for temp in temp_list:
            if temp != "all":
                do("make_metadata.py -m 10 -t " + temp)
            else:
                do("make_metadata.py -m 9")
            do("mkdir t" + temp)
            do("mv metadatafile t" + temp)
            cd("t" + temp)
            do("pulling_analysis.py -m 9")
            cd("..")
    if args.mode == 2:
        folder_list = ["rg_0.4_lipid_2_extended", "rg_0.4_lipid_2_topology"]
        temp_list = ["350", "400", "450", "500", "all"]
        for folder in folder_list:
            cd(folder)
            # do("mkdir -p 1d_dis")
            # cd("1d_dis")
            # for temp in temp_list:
            #     if temp != "all":
            #         do("make_metadata.py -m 10 -t " + temp)
            #     else:
            #         do("make_metadata.py -m 9")
            #     do("mkdir t" + temp)
            #     do("mv metadatafile t" + temp)
            #     cd("t" + temp)
            #     do("pulling_analysis.py -m 9")
            #     cd("..")
            # cd("..")
            do("mkdir -p 1d_qw")
            cd("1d_qw")
            for temp in temp_list:
                if temp != "all":
                    do("make_metadata.py -m 10 -t " + temp)
                else:
                    do("make_metadata.py -m 9")
                do("mkdir t" + temp)
                do("mv metadatafile t" + temp)
                cd("t" + temp)
                do("pulling_analysis.py -m 10")
                cd("..")
            cd("..")
            cd("..")
    if args.mode == 3:
        print("high temp refold")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [2]
        force_ramp_rate_list=[1000]
        temperature_list=[100]
        memb_k_list = [1, 2, 4, 8, 16, 32]
        rg_list = [0.4]
        force_list = [0.0]
        repeat = 1
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat,
                        commons=1)

if args.day == "oct10":
    if args.mode == 5:
        do("mkdir combined")
        cd("combined")
        do("make_metadata.py -m 9")
        do("pulling_analysis.py -m 9 --commons 1")
    if args.mode == 4:
        temp_list = [350, 400, 450, 500]
        for temp in temp_list:
            do("mkdir t" + str(temp))
            cd("t" + str(temp))
            do("make_metadata.py -m 10 -t {}".format(temp))
            do("pulling_analysis.py -m 9 --commons 1")
            cd("..")
    if args.mode == 3:
        do("make_metadata.py -m 10 -t 400")
        do("pulling_analysis.py -m 9")
    if args.mode == 2:
        print("high temp refold")
        start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [2]
        force_ramp_rate_list=[1]
        temperature_list=[500, 550, 600]
        memb_k_list = [1]
        rg_list = [0.4]
        force_list = [0.0]
        repeat = 50
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat)
    if args.mode == 1:
        fold_list = ["rg_0.4_lipid_2_temp_300_extended", "rg_0.4_lipid_2_temp_300_topology",
                        "rg_0.4_lipid_2_temp_400_extended", "rg_0.4_lipid_2_temp_400_topology"]
        for folder in fold_list:
            cd(folder)
            cd("simulation")
            do("pulling_prepare.py -m 8 --submode 1")
            cd("..")
            do("mkdir freeEnergy_2")
            cd("freeEnergy_2")
            do("make_metadata.py -m 3 -k 0.02 -t 300 --submode 1")
            do("pulling_analysis.py -m 6")
            cd("../..")

if args.day == "oct05":
    if args.mode == 5:
        print("no go, constant force, refold")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [2]
        force_ramp_rate_list=[0.5]
        temperature_list=[400]
        memb_k_list = [0, 1, 2, 4]
        rg_list = [0.4]
        force_list = [0.0]
        repeat = 10
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat)
        # # start_from_list=["native", "extended", "topology"]
        # # start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        # mode_list = [3]  # lipid mediated interaction
        # # pressure_list = [0, 0.1, 1.0]
        # pressure_list = [0, 2]
        # force_ramp_rate_list=[0.5]
        # temperature_list=[400]
        # memb_k_list = [0, 1, 2, 4]
        # rg_list = [0, 0.4]
        # force_list = [0.0]
        # repeat = 10
        # variable_test(temperature_list=temperature_list,
        #                 start_from_list=start_from_list,
        #                 rg_list=rg_list,
        #                 memb_k_list=memb_k_list,
        #                 mode_list=mode_list,
        #                 pressure_list=pressure_list,
        #                 force_ramp_rate_list=force_ramp_rate_list,
        #                 force_list=force_list,
        #                 repeat=repeat)

    if args.mode == 4:
        print("membrane effect")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 0.1, 1, 2]
        force_ramp_rate_list=[100]
        temperature_list=[200, 300, 400, 500]
        memb_k_list = [1, 2, 4, 8, 16]
        rg_list = [0, 0.08, 0.4]
        force_list =["ramp"]
        repeat = 2
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat)
    if args.mode == 1:
        fold_list = ["rg_0.4_lipid_2_temp_300_extended", "rg_0.4_lipid_2_temp_300_topology",
                        "rg_0.4_lipid_2_temp_400_extended", "rg_0.4_lipid_2_temp_400_topology"]
        for folder in fold_list:
            cd(folder)
            cd("simulation")
            do("pulling_prepare.py -m 8 --submode 1")
            cd("..")
            do("mkdir freeEnergy_2")
            cd("freeEnergy_2")
            do("make_metadata.py -m 3 -k 0.02 -t 300 --submode 1")
            do("pulling_analysis.py -m 6")
            cd("../..")
    if args.mode == 2:
        print("2d wham")
        # cd("simulation")
        # do("pulling_prepare.py -m 8")
        # cd("..")
        fold_list = ["rg_0.4_lipid_2_temp_300_extended", "rg_0.4_lipid_2_temp_300_topology",
                        "rg_0.4_lipid_2_temp_400_extended", "rg_0.4_lipid_2_temp_400_topology"]
        for folder in fold_list:
            cd(folder)
            name = "wham2d"
            do("mkdir {}".format(name))
            cd(name)
            do("make_metadata.py -m 3 -k 0.02 -t 300 --submode 1")
            do("pulling_analysis.py -m 5")
            cd("../..")
    if args.mode == 3:
        print("qw wham")
        # cd("simulation")
        # do("pulling_prepare.py -m 8")
        # cd("..")
        fold_list = ["rg_0.4_lipid_2_temp_300_extended", "rg_0.4_lipid_2_temp_300_topology",
                        "rg_0.4_lipid_2_temp_400_extended", "rg_0.4_lipid_2_temp_400_topology"]
        for folder in fold_list:
            cd(folder)
            name = "qw_2"
            do("mkdir {}".format(name))
            cd(name)
            do("make_metadata.py -m 3 -k 0.02 -t 300 --submode 1")
            do("pulling_analysis.py -m 8")
            cd("../..")
if args.day == "oct04":
    if args.mode == 1:
        # cd("simulation")
        # do("pulling_prepare.py -m 8")
        # cd("..")
        do("mv freeEnergy old_freeEnergy")
        do("mkdir freeEnergy")
        cd("freeEnergy")
        do("make_metadata.py -m 3 -k 0.02 -t 300")
        do("pulling_analysis.py -m 6")
    if args.mode == 2:
        # cd("simulation")
        # do("pulling_prepare.py -m 8")
        # cd("..")
        fold_list = ["rg_0.4_lipid_2_temp_300_extended", "rg_0.4_lipid_2_temp_300_topology",
                        "rg_0.4_lipid_2_temp_400_extended", "rg_0.4_lipid_2_temp_400_topology"]
        for folder in fold_list:
            cd(folder)
            name = "qw"
            do("mkdir {}".format(name))
            cd(name)
            do("make_metadata.py -m 3 -k 0.02 -t 300")
            do("pulling_analysis.py -m 8")
            cd("../..")
if args.day == "oct03":
    if args.mode == 1:
        print("p3 force ramp")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0]
        force_ramp_rate_list=[0.5]
        temperature_list=[200, 300, 400]
        memb_k_list = [1]
        rg_list = [0, 0.08, 0.4]
        force_list =["ramp"]
        repeat = 5
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat)
    if args.mode == 2:
        print("p3 folding temperature")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 2]
        force_ramp_rate_list=[0.5]
        temperature_list=["ramp"]
        memb_k_list = [1, 2, 4]
        rg_list = [0, 0.08, 0.4]
        force_list =["ramp"]
        repeat = 2
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat)
if args.day == "sep26":
    if args.mode == 1:
        print("constant force")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0]
        force_ramp_rate_list=[0.5]
        temperature_list=[200, 300, 400]
        memb_k_list = [1]
        rg_list = [0, 0.08, 0.4]
        force_list =[0.1, 0.2, 0.25, 0.3, 0.4]
        repeat = 5
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat)
if args.day == "sep25":
    if args.mode == 5:
        print("constant force")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        # start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [2]
        force_ramp_rate_list=[0.5]
        temperature_list=[300]
        memb_k_list = [1]
        rg_list = [0.4]
        force_list =[0.6, 0.7, 0.75, 0.8, 0.85]
        repeat = 20
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat)
    if args.mode == 4:
        print("start with basic")
        # start_from_list=["native", "extended", "topology"]
        # start_from_list=["native"]
        start_from_list=["extended", "topology"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [2]
        force_ramp_rate_list=[0.5]
        temperature_list=[300]
        memb_k_list = [1]
        rg_list = [0.4]
        force_list =[0.0, 0.1, 0.2, 0.3]
        repeat = 50
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        force_list=force_list,
                        repeat=repeat)

    if args.mode == 3:
        print("start with basic")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0]
        force_ramp_rate_list=[1]
        temperature_list=[200, 300]
        memb_k_list = [1]
        rg_list = [0.08]
        repeat = 5
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)
    if args.mode == 2:
        print("Force ramp with high rg")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 1, 2, 3]
        force_ramp_rate_list=[1]
        temperature_list=[300]
        memb_k_list = [1, 2]
        rg_list = [0, 0.2, 0.4, 0.8, 1.6]
        repeat = 3
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)
    if args.mode == 1:
        print("Force ramp without go")
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]  # lipid mediated interaction
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 1, 2, 3]
        force_ramp_rate_list=[1]
        temperature_list=[300]
        memb_k_list = [0, 1, 2, 4]
        rg_list = [0, 0.08, 0.2, 0.4, 1, 2]
        repeat = 3
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)
if args.day == "sep24":
    if args.mode == 1:
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 1, 2, 3]
        force_ramp_rate_list=[1]
        temperature_list=[200, 300]
        memb_k_list = [1, 2, 4]
        rg_list = [0, 0.02, 0.08, 0.2]
        repeat = 3
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        rg_list=rg_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)

if args.day == "sep23":
    if args.mode == 1:
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 2, 4, 6, 8]
        force_ramp_rate_list=[1]
        temperature_list=[200, 300, 400]
        memb_k_list = [1, 2, 4, 8, 16]
        repeat = 5
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)

if args.day == "sep20":
    if args.mode == 3:
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 0.1, 1, 2, 4]
        force_ramp_rate_list=[1]
        temperature_list=[300]
        memb_k_list = [1, 2, 4, 8, 16]
        repeat = 5
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)
    if args.mode == 1:
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0]
        force_ramp_rate_list=[10]
        temperature_list=["ramp"]
        memb_k_list = [1, 2, 4, 8, 16]
        repeat = 2
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        memb_k_list=memb_k_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)
    if args.mode == 2:
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0]
        force_ramp_rate_list=[10]
        temperature_list=["ramp"]
        repeat = 5
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)
if args.day == "sep17":
    if args.mode == 1:
        # start_from_list=["native", "extended", "topology"]
        start_from_list=["native"]
        mode_list = [3]
        # pressure_list = [0, 0.1, 1.0]
        pressure_list = [0, 0.1, 1, 2, 4]
        force_ramp_rate_list=[1, 10]
        temperature_list=[230, 300]
        repeat = 5
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        mode_list=mode_list,
                        pressure_list=pressure_list,
                        force_ramp_rate_list=force_ramp_rate_list,
                        repeat=repeat)
if args.day == "sep11":
    if args.mode == 1:
        # folding tempearture
        pressure_list = [0.1]
        rg_list = [0.08]
        force_list =["ramp"]
        temperature_list=[230, 300]
        memb_k_list=[1]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go"]
        force_ramp_rate_list=[1]
        variable_test(temperature_list=temperature_list,
                        memb_k_list=memb_k_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=100,force_ramp_rate_list=force_ramp_rate_list)

if args.day == "sep10":
    if args.mode == 1:
        # folding tempearture
        pressure_list = [0.4, 0.8]
        rg_list = [0.4]
        force_list =[0.3, 0.35, 0.4]
        temperature_list=[230, 300]
        memb_k_list=[1, 4]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go"]
        force_ramp_rate_list=[2]
        variable_test(temperature_list=temperature_list,
                        memb_k_list=memb_k_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=10,force_ramp_rate_list=force_ramp_rate_list)
    if args.mode == 2:
        # folding tempearture
        pressure_list = [0.1]
        rg_list = [0.08, 0.2]
        memb_k_list=[1, 2, 4]
        force_list =[0.01, 0.02, 0.04, 0.08]
        temperature_list=[230, 300]
        start_from_list=["extended", "topology"]
        # start_from_list=["native"]
        simulation_model_list=["go"]
        force_ramp_rate_list=[2]
        variable_test(temperature_list=temperature_list,
                        memb_k_list=memb_k_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=4, force_ramp_rate_list=force_ramp_rate_list)
    if args.mode == 3:
        # folding tempearture
        pressure_list = [0.1]
        rg_list = [0.08]
        memb_k_list=[1]
        force_list =[0.2, 0.3, 0.35, 0.4]
        temperature_list=[230]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go"]
        force_ramp_rate_list=[2]
        variable_test(temperature_list=temperature_list,
                        memb_k_list=memb_k_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=20, force_ramp_rate_list=force_ramp_rate_list)
    if args.mode == 4:              # pending
        # folding tempearture
        pressure_list = [0.1]
        rg_list = [0.08]
        memb_k_list=[1, 2, 4]
        force_list =[0.2, 0.3, 0.35, 0.4]
        temperature_list=[230, 300]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go"]
        force_ramp_rate_list=[2]
        variable_test(temperature_list=temperature_list,
                        memb_k_list=memb_k_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=10, force_ramp_rate_list=force_ramp_rate_list)
if args.day == "sep09":
    if args.mode == 1:
        cd("simulation")
        do("pulling_prepare.py -m 7")
        cd("..")
        do("mkdir freeEnergy")
        cd("freeEnergy")
        do("make_metadata.py -m 3 -k 0.02 -t 230")
        do("pulling_analysis.py -m 6")
    if args.mode == 2:
        cmd = "gg.py -m 8"
        run_slurm = base_slurm.format(cmd)
        # folder_list = glob.glob("force_*")
        # folder_list = ['force_0.08', 'force_0.03', 'force_0.0']
        folder_list = ['force_0.055']
        # folder_list = ['force_0.07', 'force_0.02', 'force_0.045']
        # folder_list = ['force_0.06', 'force_0.04']
        print(folder_list)
        for folder in folder_list:
            cd(folder)
            cd("simulation")
            run_list = glob.glob("*")
            for run in run_list:
                cd(run)
                cd("0")
                with open("compute_angle.slurm", "w") as r:
                    r.write(run_slurm)
                do("sbatch compute_angle.slurm")
                cd("../..")
            cd("../..")
    if args.mode == 3:
        # folding tempearture
        pressure_list = [0.0, 0.1, 0.2]
        force_list =[0.0]
        rg_list = [0, 0.08, 0.1, 0.2]
        temperature_list=["ramp"]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go", "single"]
        variable_test(temperature_list=temperature_list,
                        simulation_model_list=simulation_model_list,
                        rg_list=rg_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=3,force_ramp_rate_list=[10, 1])
    if args.mode == 4:
        # folding tempearture
        pressure_list = [0, 0.1, 0.2]
        rg_list = [0, 0.08, 0.1, 0.2]
        force_list =[0.0]
        temperature_list=[230]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go", "single"]
        force_ramp_rate_list=[1, 10]
        variable_test(temperature_list=temperature_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=5,force_ramp_rate_list=force_ramp_rate_list)
    if(args.mode == 5):
        print("Extract qw and distance info.")
        for i in range(100):
            cd(str(i))
            cd("0")
            do("awk '{print $2}' wham.dat |  sed 's/,$//' | sed 1d > qw.dat")
            do("awk '{print $2}' addforce.dat |  sed 's/,$//' | sed 1d > distance.dat")
            cd("../..")
        print("create directory_list")
        with open("directory_list", "w") as f:
            for i in range(40):
                # print(os.getcwd())
                location = os.getcwd() + "/../"
                f.write(location+str(i)+"/0\n")
        do("cp ../../2xov/2xov.pdb .")
        do("python2 ~/opt/small_script/CalcLocalDistanceStats.py 2xov directory_list out")
    if args.mode == 6:
        # folding tempearture
        pressure_list = [0.1]
        rg_list = [0.1]
        memb_k_list=[1, 4]
        force_list =[0.1, 0.15, 0.05]
        temperature_list=[230, 300]
        start_from_list=["extended", "topology"]
        # start_from_list=["native"]
        simulation_model_list=["go", "single"]
        force_ramp_rate_list=[10]
        variable_test(temperature_list=temperature_list,
                        memb_k_list=memb_k_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=7,force_ramp_rate_list=force_ramp_rate_list)
    if args.mode == 7:
        # folding tempearture
        pressure_list = [0.4]
        rg_list = [0.4]
        memb_k_list=[4]
        force_list =[0.0]
        temperature_list=[230]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go", "single"]
        force_ramp_rate_list=[1]
        variable_test(temperature_list=temperature_list,
                        memb_k_list=memb_k_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=20,force_ramp_rate_list=force_ramp_rate_list)
    if args.mode == 8:
        # folding tempearture
        pressure_list = [0.4]
        rg_list = [0.4]
        memb_k_list=[4]
        force_list =[0.1, 0.15, 0.05]
        temperature_list=[230, 300]
        start_from_list=["extended", "topology"]
        # start_from_list=["native"]
        simulation_model_list=["go", "single"]
        force_ramp_rate_list=[2]
        variable_test(temperature_list=temperature_list,
                        memb_k_list=memb_k_list,
                        rg_list=rg_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=6,force_ramp_rate_list=force_ramp_rate_list)
if args.day == "sep08":
    if args.mode == 1:
        pressure_list = [0.1, 0.2, 0.4]
        force_list =[0.01, 0.02, 0.04, 0.08]
        temperature_list=[230]
        start_from_list=["extended", "topology"]
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=20,force_ramp_rate_list=[1])
    if args.mode == 2:
        pressure_list = [0.1]
        force_list =[0.01, 0.02]
        temperature_list=[230, 240]
        # start_from_list=["extended", "topology"]
        start_from_list=["extended"]
        variable_test(temperature_list=temperature_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=50,force_ramp_rate_list=[0.5])
    if args.mode == 3:
        # folding tempearture
        pressure_list = [0.1, 0.2]
        force_list =[0.0]
        temperature_list=["ramp"]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go", "single"]
        variable_test(temperature_list=temperature_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=5,force_ramp_rate_list=[10])
    if args.mode == 4:
        # folding tempearture
        pressure_list = [0.1, 0.2]
        force_list =[0.0]
        temperature_list=[200]
        # start_from_list=["extended", "topology"]
        start_from_list=["native"]
        simulation_model_list=["go", "single"]
        force_ramp_rate_list=[1, 10]
        variable_test(temperature_list=temperature_list,
                        simulation_model_list=simulation_model_list,
                        start_from_list=start_from_list,
                        force_list=force_list,
                        pressure_list=pressure_list,
                        repeat=20,force_ramp_rate_list=force_ramp_rate_list)
if args.day == "sep07":
    if args.mode == 9:
        # simulation_model_list=["go", "single"]
        force_list =[0.55, 0.5, 0.56, 0.57, 0.58, 0.59]
        variable_test(force_list=force_list, repeat=20)
    if args.mode == 8:
        # simulation_model_list=["go", "single"]
        # force_list =[0.6, 0.5, 0.4]
        # memb_k_list = [1, 2, 4]
        pressure_list = [0.1, 0.2, 0.4]
        # pressure_list = [0.8, 1.6]
        rg_list = [0.08, 0.04, 0.02]
        k_list = [0.1, 0.2, 0.3, 0.4]
        variable_test(k_list=k_list, pressure_list=pressure_list,
                        force_ramp_rate_list=[10], rg_list=rg_list, repeat=3)
    if args.mode == 7:
        cd("simulation")
        do("pulling_prepare.py -m 7")
        cd("..")
        do("mkdir freeEnergy")
        cd("freeEnergy")
        do("make_metadata.py -m 3 -k 0.02 -t 230")
        do("pulling_analysis.py -m 6")
    if args.mode == 1:
        memb_k_list = [1, 2, 4, 8, 16]
        variable_test(memb_k_list=memb_k_list, repeat=5)
    if args.mode == 2:
        force_ramp_rate_list = [1, 2, 4, 8, 16, 32]
        variable_test(force_ramp_rate_list=force_ramp_rate_list, repeat=5)
    if args.mode == 3:
        force_ramp_rate_list = [1, 2, 4, 8, 16, 32]
        variable_test(temperature_list=[200], force_ramp_rate_list=force_ramp_rate_list, repeat=5)
    if args.mode == 4:
        # simulation_model_list=["go", "single"]
        # force_list =[0.6, 0.5, 0.4]
        force_list =[0.8, 0.6, 0.4]
        memb_k_list = [1, 2, 4]
        variable_test(memb_k_list=memb_k_list, force_list=force_list, repeat=5)
    if args.mode == 6:
        # simulation_model_list=["go", "single"]
        # force_list =[0.6, 0.5, 0.4]
        force_list =[0.3]
        memb_k_list = [1, 2, 4]
        # pressure_list = [0.1, 0.2, 0.4]
        pressure_list = [0.8, 1.6]
        k_list = [10, 11, 12, 13]
        variable_test(k_list=k_list, pressure_list=pressure_list,
                        force_ramp_rate_list=[10], memb_k_list=memb_k_list, force_list=force_list, repeat=3)
if args.day == "sep06":
    if args.mode == 1:
        start_from_list=["extended", "topology"]
        simulation_model_list=["go", "single"]
        temperature_list = [200, 250, 300]
        pressure_list = [0, 0.1, 0.5, 1]
        variable_test(pressure_list=pressure_list,
                        start_from_list=start_from_list,
                        simulation_model_list=simulation_model_list,
                        repeat=5,
                        temperature_list=temperature_list,
                        commons=0)
    if args.mode == 2:
        start_from_list=["extended", "topology"]
        simulation_model_list=["go", "single"]
        temperature_list = [250, 300]
        memb_k_list = [0, 1, 2, 4]
        rg_list = [0, 0.1]
        variable_test(rg_list=rg_list, memb_k_list=memb_k_list,
                        start_from_list=start_from_list,
                        simulation_model_list=simulation_model_list,
                        repeat=3,
                        temperature_list=temperature_list,
                        commons=0)
if args.mode == 20:
    rg_list = [0]
    temperature_list = [200]
    variable_test(rg_list=rg_list, repeat=40, temperature_list=temperature_list, commons=True)


if args.mode == 19:
    rg_list = [0]
    temperature_list = [175, 200, 225, 250]
    variable_test(rg_list=rg_list, repeat=20, temperature_list=temperature_list, commons=True)

if args.mode == 18:
    rg_list = [0, 0.1, 0.2, 1]
    memb_k_list = [0, 1, 2, 4]
    pressure_list = [0, 0.1, 0.2, 0.4, 0.8, 1, 2]
    # rg_list = [0.1]
    # memb_k_list = [1]
    # pressure_list = [0.1, 1]
    variable_test(rg_list=rg_list, memb_k_list=memb_k_list, pressure_list=pressure_list, repeat=2)

if args.mode == 17:
    # protocol_list = ["er", "awsemer", "frag", "raptor"]
    protocol_list = ["awsemer", "frag"]
    protein_list = ["1occ"]
    for protein in protein_list:
        for protocol in protocol_list:
            print("Work on protein: {}, protocol: {}".format(protein, protocol))
            if protocol == "raptor":
                do("cp ~/opt/gremlin/protein/1occ/raptor/go_rnativeC* {}/".format(protein))
            else:
                do("cp ~/opt/gremlin/protein/1occ/gremlin/go_rnativeC* {}/".format(protein))
            do("mkdir -p {}".format(protocol))
            do("cp -r {} {}/".format(protein, protocol))
            cd(protocol)
            cd(protein)
            fileName = "{}_multi.in".format(protein)
            if protocol == "raptor":
                backbone_file = "fix_backbone_coeff_er.data"
                do("cp ~/opt/gremlin/protein/{}/raptor/go_rnativeC* .".format(protein))
            else:
                backbone_file = "fix_backbone_coeff_{}.data".format(protocol)
                do("cp ~/opt/gremlin/protein/{}/gremlin/go_rnativeC* .".format(protein))
            with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
                for line in file:
                    tmp = line
                    tmp = tmp.replace("fix_backbone_coeff_er.data", backbone_file)
                    print(tmp, end='')
            cd("..")
            do("run.py -m 0 -n 20 {}".format(protein))
            cd("..")

if args.mode == 16:
    rg_list = [0, 0.1, 0.2, 0.4, 0.5, 1, 2, 4]
    variable_test(rg_list=rg_list, repeat=1, commons=True)

if(args.mode == 15):
    print("create directory_list")
    with open("directory_list", "w") as f:
        for i in range(40):
            # print(os.getcwd())
            location = os.getcwd() + "/../"
            f.write(location+str(i)+"/0\n")
    do("cp ../../2xov/2xov.pdb .")
    do("python2 ~/opt/small_script/CalcLocalDistanceStats.py 2xov directory_list out")
if(args.mode == 14):
    print("Extract qw and distance info.")
    for i in range(100):
        cd(str(i))
        cd("0")
        do("awk '{print $2}' wham.dat |  sed 's/,$//' | sed 1d > qw.dat")
        do("awk '{print $2}' addforce.dat |  sed 's/,$//' | sed 1d > distance.dat")
        cd("../..")

if args.mode == 13:
    rg_list = [0, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2]
    memb_k_list = [0, 1, 2, 4, 8]
    variable_test(rg_list=rg_list, memb_k_list=memb_k_list)
if args.mode == 12:
    rg_list = [0.1, 0.2, 0.4, 0.8, 1.6, 3.2]
    variable_test(rg_list=rg_list)
if args.mode == 11:
    zim_type_list = ["aug04", "aug26"]
    membrane_width_list = [30, 28.8]
    for zim in zim_type_list:
        for width in membrane_width_list:
            folder = "zim_{}_width_{}".format(zim, width)
            do("mkdir -p {}".format(folder))
            cd(folder)
            do("cp -r ../2xov .")
            cd("2xov")
            fixFile = "fix_backbone_coeff_single.data"
            with fileinput.FileInput(fixFile, inplace=True, backup='.bak') as file:
                for line in file:
                    print(line.replace("WIDTH", str(width)), end='')
            do("cp zim_{} zim".format(zim))
            cd("..")
            do("run.py -n 2 2xov")
            cd("..")
if args.mode == 10:
    distance_list = np.linspace(166, 180, 15)
    for distance in distance_list:
        folder = "dis_{}".format(distance)
        cd(folder)
        do("sbatch run_0.slurm")
        cd("..")
# if args.mode == 9:
#     cmd = "python3 ~/opt/small_script/find_distance.py"
#     run_slurm = base_slurm.format(cmd)
#     folder_list = ['force_0.045']
#     print(folder_list)
#     for folder in folder_list:
#         cd(folder)
#         cd("simulation")
#         run_list = glob.glob("*")
#         for run in run_list:
#             cd(run)
#             cd("0")
#             with open("find_distance.slurm", "w") as r:
#                 r.write(run_slurm)
#             do("sbatch find_distance.slurm")
#             cd("../..")
#         cd("../..")
#
# if args.mode == 8:
#     cmd = "gg.py -m 8"
#     run_slurm = base_slurm.format(cmd)
#
#     # folder_list = glob.glob("force_*")
#     # folder_list = ['force_0.08', 'force_0.03', 'force_0.0']
#     folder_list = ['force_0.055']
#     # folder_list = ['force_0.07', 'force_0.02', 'force_0.045']
#     # folder_list = ['force_0.06', 'force_0.04']
#     print(folder_list)
#     for folder in folder_list:
#         cd(folder)
#         cd("simulation")
#         run_list = glob.glob("*")
#         for run in run_list:
#             cd(run)
#             cd("0")
#             with open("compute_angle.slurm", "w") as r:
#                 r.write(run_slurm)
#             do("sbatch compute_angle.slurm")
#             cd("../..")
#         cd("../..")

# if args.mode == 7:
#     for i in range(80):
#         do("mv {0} ../../../new_force_ramp/memb_0_force_ramp_rg_0_new/simulation/{1}".format(i,i+90))
# if args.mode == 6:
#     force_list = [0.55, 0.6, 0.65]
#     # force_list = [0.25, 0.35, 0.4, 0.45]
#     # force_list = [0.15, 0.2]
#     for force in force_list:
#         do("mkdir force_{}".format(force))
#         do("cp -r 2xov force_{}/".format(force))
#         cd("force_{}".format(force))
#         with fileinput.FileInput("2xov/2xov_multi.in", inplace=True, backup='.bak') as file:
#             for line in file:
#                 print(line.replace("MY_FORCE", str(force)), end='')
#         do("run.py -n 10 2xov/")
#         cd("..")
#
#
# if args.mode == 5:
#     # cd("start_misfolded")
#     distance_list = np.linspace(0, 30, 16)
#     for dis in distance_list:
#         do("mkdir -p dis_{}".format(dis))
#         do("cp -r ../2xov/ dis_{}".format(dis))
#         do("cp ../../freeEnergy/go_model_start_unfolded/simulation/dis_{0}/restart.25000000 dis_{0}/2xov/".format(dis))
#         cd("dis_{}".format(dis))
#         do("run.py -n 10 2xov/")
#         cd("..")

# if args.mode == 4:
#     do("rm data")
#     for i in range(100):
#         do("cat dis_{}/0/data >> data.dat".format(i))
#     do("awk '{print $1}' data.dat  > e.dat")
#     do("awk '{print $2}' data.dat  > p.dat")
#     do("awk '{print $3}' data.dat  > qw.dat")
#
# if args.mode == 1:
#     cd("simulation")
#     do("pulling_prepare.py")
#     cd("..")
#     do("mkdir freeEnergy")
#     cd("freeEnergy")
#     do("make_metadata.py -k 0.05 -t 300")
#     do("pulling_analysis.py -m 3 -p 2")
#
# if args.mode == 2:
#     print("80 bins.")
#     # cd("simulation")
#     # do("pulling_prepare.py")
#     # cd("..")
#     do("mkdir more_bin")
#     cd("more_bin")
#     do("make_metadata.py -k 0.05 -t 600")
#     do("pulling_analysis.py -m 3 -p 1")
#
# if args.mode == 3:
#     # cd("simulation")
#     # do("pulling_prepare.py")
#     # cd("..")
#     do("mkdir -p only_less_than_100")
#     cd("only_less_than_100")
#     do("make_metadata.py -k 0.05 -t 300 -m 2")
#     do("pulling_analysis.py -m 3 -p 1")
#     # for i in range(52, 70):
#     #     do("mv {}/{} .".format(i, i-40))
#     #     do("mv {} {}".format(i, i-20))
#     # for i in range(50):
#     #     do("mv force_0.8_2/{} force_0.8/{}".format(i, i+50))
#         # do("mv half_contact_force_0.8_memb1_rg1_2/{} half_contact_force_0.8_memb1_rg1/{}".format(i, i+20))
#

if(args.run):
    print("Hello World")

    name = "T0833"
    n = 21
    do("mkdir "+name)
    cd(name)
    for i in range(1, n):
        do("mkdir -p job.{}".format(i))
        do("cp ../preparation_files/myjob_nots.slurm job.{}".format(i))
        do("cp ../preparation_files/loopsubmit.bash job.{}".format(i))
        do("cp -r ../preparation_files/{0}_runpackage job.{1}/runpackage".format(name, i))
        do("cp ../preparation_files/{1}_tpr/run.{0}.tpr job.{0}/runpackage/run.tpr".format(i, name))
    for i in range(1, n):
        cd("job.{}".format(i))
        fileName = "myjob_nots.slurm"
        with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace("T0766", name), end='')

        do("bash loopsubmit.bash")
        cd("..")

    # for i in range(1, 6):
    #     do("mkdir job.{}".format(i))
    #     do("cp myjob_nots.slurm job.{}".format(i))
    #     do("cp loopsubmit.bash job.{}".format(i))
    #     do("cp -r runpackage job.{}".format(i))
    #     do("cp run.{0}.tpr job.{0}/runpackage/run.tpr".format(i))
    # for i in range(1, 6):
    #     cd("job.{}".format(i))
    #     fileName = "myjob_nots.slurm"
    #     name = "T0833"
    #     with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
    #         for line in file:
    #             print(line.replace("T0766", name), end='')
    #     do("bash loopsubmit.bash")
    #     cd("..")

    # force_list = [0.3, 0.5, 0.7, 0.9, 1.1, 1.3]
    # for force in force_list:
    #     folder = "1d_force_" + str(force)
    #     do("mkdir -p " + folder)
    #     cd(folder)
    #     do("cp ../metadatafile .")
    #     do("~/bin/python3/bin/python3 ~/opt/pulling_analysis.py -f -m 5 --force {}".format(force))
    #     do("sbatch freeEnergy.slurm")
    #     cd("..")
    # do("rm -r test")
    # do("cp -r 2xov test")
    # cd("test")
    # do("test_run.py test.in")

if(args.see):
    do("head test/0/addforce.dat")
def fix_error_run():
    n = args.number
    for i in range(n):
        os.chdir(str(i))
        os.system("cp ~/opt/pulling/qnqc.slurm .")
        os.system("sbatch qnqc.slurm")
        os.chdir("..")
    # os.system("grep 'srun: error: Application launch failed: Socket timed out on send/recv operation' . -r | cut -d':' -f1 | rev | cut -d"/" -f2- | rev > list")
    # array = []
    # cwd = os.getcwd()
    # print(cwd)
    # with open('list', 'r') as ins:
    #     for line in ins:
    #         target = line.strip('\n')
    #         array.append(target)
    # for i in array:
    #     os.chdir(i)
    #     os.system("rm slurm-*")
    #     os.system("sbatch rerun.slurm")
    #     sleep(0.5)  # Time in seconds.
    #     os.system("sbatch qnqc.slurm")
    #     sleep(0.5)  # Time in seconds.
    #     # os.system("pwd")
    #     os.chdir(cwd)

    # os.system("cut -d'/' -f2 list >")
# if(args.fix):
#     fix_error_run()


def rerun():
    n = args.number
    for i in range(n):
        os.system("cp -r {0} rerun_{0}".format(str(i)))
        # source = "~/opt/gagb/gagb_constant200_rerun.in"
        # target = " 2lhc.in"
        source = "~/opt/pulling/2xov_force_load_dis.in"
        target = " 2xov.in"
        os.system("cp "+source + target)
        os.system("cp "+target+" rerun_{}/".format(str(i)))
        os.chdir("rerun_"+str(i))
        os.system("rm slurm*")
        os.system("sbatch run.slurm")
        os.chdir("..")
# if(args.rerun):
#     rerun()


def continue_run():
    n = 10
    os.chdir("simulation")
    for i in range(n):
        os.system("cp -r {0} continue_{0}".format(str(i)))
        # source = "~/opt/gagb/gagb_constant200_rerun.in"
        # target = " 2lhc.in"
        source = "~/opt/pulling/2xov_continue_run.in"
        target = " 2xov.in"
        os.system("cp "+source + target)
        os.system("cp "+target+" continue_{}/".format(str(i)))
        os.chdir("continue_"+str(i))
        os.system("rm slurm*")
        os.system(  # replace RANDOM with a radnom number
            "sed -i.bak 's/RANDOM/'" +
            str(randint(1, 10**6)) +
            "'/g' 2xov.in")
        os.system("sbatch run.slurm")
        os.chdir("..")
# if(args.go):
#     continue_run()
# parser = argparse.ArgumentParser(
#         description="This is my playground for current project")
# parser.add_argument("protein", help="the name of protein")
# # parser.add_argument("template", help="the name of template file")
# args = parser.parse_args()
# # protein_name = args.template.split('_', 1)[-1].strip('/')
# protein_name = args.protein.strip('/')
    # name = "ga_2m"
## -------------Pulling--------
# os.system("cp ~/opt/small_script/springForce.plt .")
# os.system("cp ~/opt/small_script/springForce_smooth.plt .")
# os.system("gnuplot springForce.plt")
# os.system("gnuplot springForce_smooth.plt")
# os.system("open springForce.pdf")
# os.system("open springForce_smooth.pdf")
# SpringConstant_list = [3e-05, 5e-05, 1e-06, 3e-06, 5e-06, 1e-05, 1e-07]
# # SpringConstant_list = [3e-06, 5e-06]
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
#     os.system("run.py 2xov/ -s 8 -n 2")
#     os.chdir("..")

# number_of_run_list = [2, 4, 8, 16]
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
