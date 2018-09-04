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
# from run_parameter import *
parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    do see the difference variable make \
    run simulation")

parser.add_argument("template", help="the name of template file")
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("--rerun",
                    type=int, default=0)
parser.add_argument("-m", "--mode", type=int, default=1)
parser.add_argument("--commons", type=int, default=0)
parser.add_argument("--nick", action="store_true", default=False)

args = parser.parse_args()
protein_name = args.template.strip('/')
template = protein_name
if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

# protein_name = args.template.split('_', 1)[-1].strip('/')
# os.system("cp ~/opt/variable_test_run.py .")

# run_slurm = '''\
# #!/bin/bash
# #SBATCH --job-name=CTBP_WL
# #SBATCH --account=ctbp-common
# #SBATCH --partition=ctbp-common
# #SBATCH --ntasks=1
# #SBATCH --mem-per-cpu=1G
# #SBATCH --time=1-00:00:00
# #SBATCH --mail-user=luwei0917@gmail.com
# #SBATCH --mail-type=FAIL
# echo "My job ran on:"
# echo $SLURM_NODELIST
# srun ~/build/lammps_awsemmd_20161127/src/lmp_serial -in 2xov_{}.in
# '''

# run_slurm = '''\
# #!/bin/bash
# #SBATCH --job-name=CTBP_WL
# #SBATCH --account=ctbp-common
# #SBATCH --partition=ctbp-common
# #SBATCH --ntasks=1
# #SBATCH --threads-per-core=1
# #SBATCH --mem-per-cpu=1G
# #SBATCH --time=1-00:00:00
# #SBATCH --mail-user=luwei0917@gmail.com
# #SBATCH --mail-type=FAIL
# echo "My job ran on:"
# echo $SLURM_NODELIST
# srun ~/build/lammps_awsemmd_20161127/src/lmp_serial -in 2xov_{}.in
# '''

# run_slurm = '''\
# #!/bin/bash
# #SBATCH --job-name=CTBP_WL
# #SBATCH --account=ctbp-common
# #SBATCH --partition=ctbp-common
# #SBATCH --ntasks=1
# #SBATCH --threads-per-core=1
# #SBATCH --mem-per-cpu=1G
# #SBATCH --time=1-00:00:00
# #SBATCH --mail-user=luwei0917@gmail.com
# #SBATCH --mail-type=FAIL
# echo "My job ran on:"
# echo $SLURM_NODELIST
# srun /home/wl45/build/awsem_new_membrane/src/lmp_serial -in 2xov_{}.in
# '''
if args.mode == 2 or args.mode == 4:
    run_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun /home/wl45/build/awsem_lipid_fluctuations/src/lmp_serial -in 2xov_{}.in
'''

if args.mode == 3:
    run_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=12
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun /home/wl45/build/awsem_lipid_fluctuations/src/lmp_mpi -p 12x1 -in 2xov_{}.in
'''
# if args.mode == 5:
#     run_slurm = '''\
# #!/bin/bash
# #SBATCH --job-name=CTBP_WL
# #SBATCH --account=ctbp-common
# #SBATCH --partition=ctbp-common
# #SBATCH --ntasks=12
# #SBATCH --threads-per-core=1
# #SBATCH --constraint=opath
# #SBATCH --ntasks-per-node=24
# #SBATCH --mem-per-cpu=1G
# #SBATCH --time=1-00:00:00
# #SBATCH --mail-user=luwei0917@gmail.com
# #SBATCH --mail-type=FAIL
# echo "My job ran on:"
# echo $SLURM_NODELIST
# srun /home/wl45/build/awsem_lipid_fluctuations/src/lmp_mpi -p 12x1 -in 2xov_{}.in
# '''



if args.mode == 5:
    run_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun /home/wl45/build/awsem_lipid_fluctuations/src/lmp_mpi -p 12x1 -in 2xov_{}.in
'''

if args.mode == 6:
    run_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=10
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun /home/wl45/build/awsem_lipid_fluctuations/src/lmp_mpi -p 10x1 -in 2xov_{}.in
'''

if args.commons == 1:
    run_slurm = run_slurm.replace("ctbp-common", "commons")
if args.nick:
    run_slurm = run_slurm.replace("/home/wl45/build/awsem_lipid_fluctuations/src/lmp_mpi", "/home/ns24/lmp_mpi")
def change(fileName, from_str, to_str):
    with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
        for line in file:
            tmp = line
            tmp = tmp.replace(from_str, str(to_str))
            print(tmp, end='')

fileName = "2xov_multi.in"
if args.rerun == 0:
    start_from = "read_data data.2xov"
if args.rerun == 1:
    start_from = "read_restart restart.2000000"

if args.mode ==6:
    qbias_list = np.linspace(0.1,1,46)
    i = args.rerun
    do("mkdir simulation")
    cd("simulation")
    # qbias_list = [i*0.02 for i in range(50)]
    for qbias in qbias_list:
        folder_name = "qbias_{:.2f}".format(qbias)
        do("cp -r ../2xov " + folder_name)
        cd(folder_name)
        # fixFile = "fix_backbone_coeff_go.data"
        fixFile = "fix_qbias.data"
        with fileinput.FileInput(fixFile, inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace("MY_QBIAS", str(qbias)), end='')

        do("cp 2xov_multi.in 2xov_{}.in".format(i))
        # with fileinput.FileInput("2xov_{}.in".format(i), inplace=True, backup='.bak') as file:
        #     for line in file:
        #         print(line.replace("START_FROM", start_from), end='')
        do(  # replace SIMULATION_STEPS with specific steps
            "sed -i.bak 's/NUMBER/'" +
            str(int(i)) +
            "'/g' 2xov_{}.in".format(i))
        do("mkdir -p {}".format(i))
        do(  # replace RANDOM with a radnom number
            "sed -i.bak 's/RANDOM/'" +
            str(randint(1, 10**6)) +
            "'/g' *.in")
        with open("run_{}.slurm".format(i), "w") as r:
            r.write(run_slurm.format(i))
        do("sbatch " + "run_{}.slurm".format(i))
        cd("..")
if args.mode == 5:
    # distance_list = np.linspace(30, 230, 101)
    # distance_list = np.linspace(30, 180, 51)
    # distance_list = np.linspace(30, 140, 56)
    # distance_list = np.linspace(30, 100, 36)
    # distance_list = np.linspace(40, 110, 36)
    distance_list = np.linspace(110, 350, 41)
    # distance_list = np.linspace(30, )
# rg_list = [0, 1, 5, 10]
# force_list = [2.0]
# memb_k_list = [0, 1, 5, 10]

# rg_list = [0, 1, 2, 5]
# force_list = [0.0, 1.0, 2.0, 3.0]
# memb_k_list = [0, 1, 2, 5]

# rg_list = [0, 1, 2, 5]
# force_list = [0.0, 3.0]
# memb_k_list = [0, 1, 2, 5]

# rg_list = [0, 0.1, 1, 5, 10]
# force_list = [0.0, 3.0]
# memb_k_list = [0, 0.1, 1, 5, 10]

# rg_list = [0, 0.1, 1]
# force_list = ["ramp"]
# memb_k_list = [0, 0.1, 1]

# rg_list = [0, 0.1, 0.5, 1, 2]
# rg_list = [3, 4]
# force_list = ["ramp"]
# memb_k_list = [0, 0.1, 1, 2, 5, 10]
if args.mode == 4:
    i = args.rerun
    do("mkdir simulation")
    cd("simulation")
    # qbias_list = [i*0.02 for i in range(50)]
    for ii in range(50):
        qbias = ii*0.02
        folder_name = "qbias_{}".format(ii)
        do("cp -r ../2xov " + folder_name)
        cd(folder_name)
        # fixFile = "fix_backbone_coeff_go.data"
        fixFile = "fix_qbias_equil.data"
        with fileinput.FileInput(fixFile, inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace("QBIAS", str(qbias)), end='')

        do("cp 2xov_multi.in 2xov_{}.in".format(i))
        with fileinput.FileInput("2xov_{}.in".format(i), inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace("START_FROM", start_from), end='')
        do(  # replace SIMULATION_STEPS with specific steps
            "sed -i.bak 's/NUMBER/'" +
            str(int(i)) +
            "'/g' 2xov_{}.in".format(i))
        do("mkdir -p {}".format(i))
        do(  # replace RANDOM with a radnom number
            "sed -i.bak 's/RANDOM/'" +
            str(randint(1, 10**6)) +
            "'/g' *.in")
        with open("run_{}.slurm".format(i), "w") as r:
            r.write(run_slurm.format(i))
        do("sbatch " + "run_{}.slurm".format(i))
        cd("..")

if args.mode == 1:
    distance_list = np.linspace(0, 100, 51)
if args.mode == 2:
    distance_list = np.linspace(30, 180, 151)
if args.mode == 3:
    # distance_list = np.linspace(30, 230, 101)
    distance_list = np.linspace(30, 100, 36)
    distance_list = np.linspace(100, 340, 41)
    # distance_list = np.linspace(60.5, 90.5, 16)
    # distance_list = np.linspace(132, 232, 51)

    # distance_list = np.linspace(10, 180, 171)
    # distance_list = np.linspace(0, 1, 1)

if args.mode <= 3 or args.mode ==5:
    i = args.rerun
    i = 0
    do("mkdir simulation")
    cd("simulation")
    for dis in distance_list:
        folder_name = "dis_{}".format(dis)
        do(f"cp -r ../{template} " + folder_name)
        cd(folder_name)
        # fixFile = "fix_backbone_coeff_go.data"
        # fixFile = "colvars.x"
        fixFile = f"{template}_multi.in"
        with fileinput.FileInput(fixFile, inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace("DISTANCE", str(dis)), end='')
        do(f"cp {template}_multi.in {template}_{i}.in")
        with fileinput.FileInput(f"{template}_{i}.in", inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace("START_FROM", start_from), end='')
        do(  # replace SIMULATION_STEPS with specific steps
            "sed -i.bak 's/NUMBER/'" +
            str(int(i)) +
            f"'/g' {template}_{i}.in")
        do(f"mkdir -p {i}")
        do(  # replace RANDOM with a radnom number
            "sed -i.bak 's/RANDOM/'" +
            str(randint(1, 10**6)) +
            "'/g' *.in")
        with open(f"run_{i}.slurm", "w") as r:
            r.write(run_slurm.format(i).replace("2xov", template))

        do("sbatch " + "run_{}.slurm".format(i))
        cd("..")
