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

parser = argparse.ArgumentParser(description="This is my playground for current project")

# parser.add_argument("protein", help="the name of protein")
# parser.add_argument("template", help="the name of template file")
parser.add_argument("--gagb", help="Project gagb", action="store_true", default=False)
parser.add_argument("-f", "--freeEnergy", help="free energy calculation", action="store_true", default=False)
parser.add_argument("-t", "--test", help="test ", action="store_true", default=False)
parser.add_argument("--ga", action="store_true", default=False)
parser.add_argument("--gb", action="store_true", default=False)
parser.add_argument("-d", "--debug", action="store_true", default=False)
args = parser.parse_args()

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

if(args.test):
    print("Welcome to test world")


if(args.gb):
    name_list = ["gb77", "gb88b", "gb91", "gb95", "gb", "ga", "ga95", "ga91", "ga88", "ga77"]
    for name in name_list:
        os.chdir(name)
        os.system("mkdir one_d_gb")
        os.chdir("one_d_gb")
        os.system("make_metadata.py --gagb")
        os.system("python2 ~/opt/gagb-compute-pmf.py")
        os.system("cp pmf-350.dat ~/Research/results/gb/{}-350.dat".format(name))
        os.chdir("..")
        os.chdir("..")


def ga():
    name_list = ["gb77", "gb", "ga", "ga95", "ga77"]
    for name in name_list:
        os.chdir(name)
        os.system("mkdir one_d_ga")
        os.chdir("one_d_ga")
        os.system("make_metadata.py --gagb")
        os.system("python2 ~/opt/gagb-compute-pmf.py -v1 2")
        os.system("cp pmf-350.dat ~/Research/results/ga/{}-350.dat".format(name))
        os.chdir("..")
        os.chdir("..")
if(args.ga):
    ga()

gagb_free_energy_config = "\n\
ndim = {} \n\
biasing_variable_column = 3\n\
energy_column = 6                  \n\
pmf_variable_column_1 = {}          \n\
pmf_variable_column_2 = {}          \n\
nbins1 = 20                        \n\
nbins2 = 20                        \n\
start_temperature = 200            \n\
end_temperature = 400              \n\
temperature_increment = 10         \n\
N_samples = 4000                   \n\
"
ndim = 1
pmf_variable_column_1 = 3
pmf_variable_column_2 = 2


def gagb_freeEnergy_calculation():
    # os.chdir("simulation")
    # os.system("gg.py -f")
    # os.chdir("..")
    print("gagb_freeEnergy_calculation")
    # wham one d on ga, gb
    name_list = ["ga", "gb", "two_d"]
    # name_list = ["test"]
    for name in name_list:
        if name == "two_d":
            folder = "two_d"
        else:
            folder = "one_d_"+name
        os.system("mkdir -p "+folder)
        os.chdir(folder)
        os.system("make_metadata.py --gagb")
        config = open('config.py', 'w')
        if name == "ga":
            ndim = 1
            pmf_variable_column_1 = 2
            pmf_variable_column_2 = 3
        if name == "gb":
            ndim = 1
            pmf_variable_column_1 = 3
            pmf_variable_column_2 = 2
        if name == "two_d":
            ndim = 2
            pmf_variable_column_1 = 3           # gb
            pmf_variable_column_2 = 2
        this_config = gagb_free_energy_config.format(ndim, pmf_variable_column_1, pmf_variable_column_2)
        config.write(this_config)
        # config.write("protein_name = '%s'\nnumber_of_run = %d\nsimulation_steps = %d\n\
        # warm_up_steps = %d\n" % (protein_name, n, simulation_steps, warm_up_steps))
        config.close()
        os.system("python2 ~/opt/compute-pmf.py")
        os.system("cp pmf-330.dat ../"+name+"_pmf-330.dat")
        script = "tail -n+2 cv-200-400-10.dat | sort -r -k 2 | head -n1"
        cv_peak_temp = subprocess.check_output(script, shell=True).decode("utf-8").split()[0]
        os.system("cp pmf-"+str(cv_peak_temp)+".dat ../"+name+"_cv_peak_"+str(cv_peak_temp)+".dat")
        os.system("cp cv-200-400-10.dat ../"+name+"_cv-200-400-10.dat")
        # os.system("myplot.py --gagb ../" + name + ".pdf")
        os.chdir("..")
    # wham two d


def gagb_freeEnergy_analysis():
    name_list = ["ga", "gb", "two_d"]
    # name_list = ["test"]
    for name in name_list:
        os.system("cp pmf-330.dat ../"+name+"_pmf-330.dat")
        script = "tail -n+2 cv-200-400-10.dat | sort -r -k 2 | head -n1"
        cv_peak_temp = subprocess.check_output(script, shell=True).decode("utf-8").split()[0]
        os.system("cp pmf-"+str(cv_peak_temp)+".dat ../"+name+"_cv_peak_"+str(cv_peak_temp)+".dat")
        os.system("cp cv-200-400-10.dat ../"+name+"_cv-200-400-10.dat")
        # os.system("myplot.py --gagb ../" + name + ".pdf")
        os.chdir("..")


def gagb_calall():
    # name_list = ["ga77", "gb77"]
    name_list = ["ga", "gb"]
    os.system("mkdir -p results")
    for name in name_list:
        os.chdir(name)
        skip = False
        if not skip:
            gagb_freeEnergy_calculation()
        # os.system("cp ga.pdf ../results/" + name+"_ga.pdf")
        # os.system("cp gb.pdf ../results/" + name+"_gb.pdf")
        # os.system("cp two_d.pdf ../results/"+name+"_two_d.pdf")
        # pdf_list = glob.glob("*.pdf")
        # for pdf in pdf_list:
        #     os.system("mv "+pdf+" ../results/"+name+pdf)
        data_list = glob.glob("*.dat")
        for data in data_list:
            os.system("cp "+data+" ../results/"+name+data)
        # os.system("cp pmf-330.dat ../results/"+name+"_pmf-330.dat")
        os.chdir("..")

if(args.gagb and args.freeEnergy):
    gagb_calall()
