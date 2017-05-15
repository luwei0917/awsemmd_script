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
import numpy
# Useful codes
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# tr " " "\n"
# sed 1d
# sort -u -k 3
# sed -e 's/+T//'
mypath = os.environ["PATH"]
os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
my_env = os.environ.copy()

parser = argparse.ArgumentParser(description="Glue is for code needs constant changes to meet various needs\
                                    of each task")

# parser.add_argument("protein", help="the name of protein")
# parser.add_argument("template", help="the name of template file")
parser.add_argument("-t", "--test", help="test ", action="store_true", default=False)
parser.add_argument("--plot", action="store_true", default=False)
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("--protein", default="2xov")
parser.add_argument("--dimension", type=int, default=1)
parser.add_argument("-f", "--freeEnergy", action="store_true", default=False)
parser.add_argument("--move", action="store_true", default=False)
parser.add_argument("-m", "--mode", type=int, default=1)
parser.add_argument("-n", "--number", type=int, default=40,
                    help="Number of simulation run")
args = parser.parse_args()


if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir


if(args.freeEnergy):
    # sim_list = "t250 t275 t300 t325 t350"
    # temp_list = "250 275 300 325 350"
    sim_list = "t250 t300 t350"
    temp_list = "250 300 350"
    # do("mult_calc_cv.sc . '{}' 40 '{}' 100 200 10 30 400 0.12 0.9 2xov q".format(sim_list, temp_list))
    do("mult_calc_cv.sc . '{}' 66 '{}' 200 400 10 30 0.02 25 345 2xov q".format(sim_list, temp_list))

if(args.plot):
    do("plotcontour.py pmf-200.dat -xmax 0.8 -xmin 0.2 -ymin 0.2 -ymax 0.8")

# ls */halfdata | sort -g | xargs cat > all_halfdata
# ls */tinydata | sort -g | xargs cat > all_tinydata
# awk '{print $3}' all_halfdata > p_total
# awk '{print $4}' all_halfdata > e_total
# ls [0-9]* |sort -g | xargs cat > data


def move_data_to_wham(temp_list):
    for temp in temp_list:
        do("cp ../data/{}/data data".format(temp))
        do("awk '{print $1}' data > qn_t%i" % (temp))
        do("awk '{print $2}' data > qc_t%i" % (temp))
        do("awk '{print $3}' data > q_t%i" % (temp))
        do("awk '{print $4}' data > energy_t%i" % (temp))


def write_simulation_list(temp_list):
    with open("T_list", "w") as f:
        for temp in temp_list:
            f.write(str(temp)+"\n")
    with open("sim_list", "w") as f:
        for temp in temp_list:
            f.write("t"+str(temp)+"\n")


def get_total_x(temp_list):
    x_list = ["q", "qn", "qc", "energy"]
    for x in x_list:
        for temp in temp_list:
            do("cat {0}_t{1} >> {0}_total".format(x, temp))


def replace_random(file_name):
    do(  # replace RANDOM with a radnom number
        "sed -i.bak 's/RANDOM/'" +
        str(randint(1, 10**6)) +
        "'/g' {}".format(file_name))



if(args.move):
    if(args.mode == 14):
        n = 40
        do("mkdir -p analysis/data")
        run_list = [0, 1, 2, 3, 4, 5]
        run_list = [6]
        run_list = [7, 8, 9]
        run_list = [0, 1, 2]
        run_list = [0]
        run_list = [0, 1]
        for j in run_list:
            for i in range(n):
                do("cp simulation/{0}/{1}/data analysis/data/{0}_{1}.dat".format(i, j))
    if(args.mode == 13):
        replace_random("*.in")
    if(args.mode == 12):
        seed(datetime.now())
        n = args.number     # default is 40
        cwd = os.getcwd()
        # do("cp -r simulation back_up")
        for i in range(n):
            cd("simulation/" + str(i))
            do("cp -r ~/continue_run_addon/* .")
            replace_random("2xov_1.in")
            # do("cp 2xov_a206g.seq 2xov.seq")
            do("cp 2xov_l155a.seq 2xov.seq")
            do("mv addforce.dat energy.dat dump.lammpstrj wham.dat 0/")
            do("sbatch run_1.slurm")
            cd(cwd)
    if(args.mode == 11):
        n = 40
        do("mkdir -p analysis/data")
        for i in range(n):
            do("cp simulation/{0}/addforce.dat analysis/data/{0}.dat".format(i))
    if(args.mode == 10):
        folder_name = "multi_temp_2"
        do("mkdir "+folder_name)
        cd(folder_name)
        temp_list = [135, 160, 185, 210]
        move_data_to_wham(temp_list)
        write_simulation_list(temp_list)
        get_total_x(temp_list)
        sim_list = 't135 t160 t185 t210'
        temp_list = '135 160 185 210'
        do("mult_calc_cv.sc . '{}' 20 '{}' 150 350 10 30 200 0 0.95 2xov q".format(sim_list, temp_list))
        cd("..")
    if(args.mode == 9):
        do("mult_calc_cv.sc . 't135 t160 t185 t210' 20 '135 160 185 210' 150 350 10 30 200 0 0.95 2xov q")
    if(args.mode == 8):
        x_list = ["q", "qn", "qc", "energy"]
        temp_list = [135, 160, 185, 210]
        for x in x_list:
            for temp in temp_list:
                do("cat {0}_t{1} >> {0}_total".format(x, temp))
    if(args.mode == 7):
        temp_list = [135, 160, 185, 210]
        write_simulation_list(temp_list)
    if(args.mode == 6):
        temp_list = [135, 160, 185, 210]
        move_data_to_wham(temp_list)
    if(args.mode == 5):
        temp_list = [135, 160, 185, 210]
        for temp in temp_list:
            cd(str(temp))
            do("ls [0-9]* |sort -g | xargs cat > data")
            cd("..")
    if(args.mode == 4):
        n = 20
        temp_list = [135, 160, 185, 210]
        # temp_list = ['300', '200', '250']
        cwd = os.getcwd()
        do("mkdir -p analysis")
        for temp in temp_list:
            do("mkdir -p analysis/{}".format(temp))
            for i in range(n):
                print(str(i))
                do("cp simulation/{0}/{1}/data analysis/{0}/{1}".format(temp, i))
                # do("cp simulation/{0}/{1}/small_data analysis/{0}/first_2000_{1}".format(temp, i))
        cd("analysis")
        do("mkdir data")
        do("mv * data/")
    if(args.mode == 1):
        n = 40
        # temp_list = [250, 275, 300, 325, 350]
        temp_list = [200]
        cwd = os.getcwd()
        for temp in temp_list:
            for i in range(n):
                print(str(i))
                do("mkdir -p analysis/data/{}".format(i))
                do("cp simulation/{0}/{1}/halfdata analysis/data/{1}/".format(temp, i))
                cd("analysis/data/{}".format(i))
                do("tail -n 1000 halfdata > tinydata")
                cd(cwd)
    if(args.mode == 2):
        array = []
        cwd = os.getcwd()
        print(cwd)
        with open('folder_list', 'r') as ins:
            for line in ins:
                target = line.strip('\n')
                temp = target.split("_")[1]
                x = target.split("_")[3]
                t1 = "simulation/" + target + "/"
                cd(t1)
                do("pwd")
                do("cat halfdata >> ../t{}".format(temp))
                cd(cwd)
    if(args.mode == 3):
        cwd = os.getcwd()
        temp_list = [250, 300, 350]
        for temp in temp_list:
            do("cp ../data/t{} data".format(temp))
            do("awk '{print $1}' data > qn_t%i" % (temp))
            do("awk '{print $2}' data > qc_t%i" % (temp))
            do("awk '{print $3}' data > q_t%i" % (temp))
            do("awk '{print $4}' data > energy_t%i" % (temp))
# if(args.test):
    # force_list = [1.0, 1.2, 1.4, 1.6, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5]
    # force_list = [round(i*0.1,2) for i in range(10)]
    # force_list = [round(i*0.1,2) for i in range(20)]
    # # force_list = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # for force in force_list:
    #     do("mkdir {}".format(force))
    #     cd("{}".format(force))
    #     do("cp ../freeEnergy.slurm .")
    #     do("cp ../metadatafile .")
    #     do(
    #         "sed -i.bak 's/FORCE/" +
    #         str(force) +
    #         "/g' freeEnergy.slurm")
    #     do("sbatch freeEnergy.slurm")
    #     cd("..")
