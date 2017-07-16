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
args = parser.parse_args()
protein_name = args.template.strip('/')

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
srun ~/build/lammps_awsemmd_20161127/src/lmp_serial -in 2xov_{}.in
'''

fileName = "2xov_multi.in"
if args.rerun == 0:
    start_from = "read_data data.2xov"
if args.rerun == 1:
    start_from = "read_restart restart.5000000"
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
distance_list = np.linspace(0, 400, 40)

i = args.rerun

for dis in distance_list:
    folder_name = "dis_{}".format(dis)
    do("cp -r 2xov " + folder_name)
    cd(folder_name)
    # fixFile = "fix_backbone_coeff_go.data"
    fixFile = "colvars.x"
    with fileinput.FileInput(fixFile, inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace("DISTANCE", str(dis)), end='')

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
