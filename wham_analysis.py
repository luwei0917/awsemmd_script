#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import imp

mypath = os.environ["PATH"]
os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
my_env = os.environ.copy()

parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        automatically prepare for wham analysis")

parser.add_argument("protein", help="the name of protein")
args = parser.parse_args()

exec (open("config.py").read())
n = number_of_run
steps = simulation_steps

protein_name = args.protein.strip('/')

temp = 400
folder_name = "{}_t{}_q100".format(protein_name, str(temp))
print("all going to "+folder_name)
os.system("mkdir -p "+folder_name)
os.system("rm -f "+folder_name + "/*")
command = 'sed 1d simulation/{}/%d/wham.dat \
>> {}/all_wham.dat'.format(temp, folder_name)
# cal rmsd
# os.system("paste q_ga.dat q_ga_part.dat q_gb.dat p.dat > data.dat")
os.chdir("simulation/"+str(temp))
for i in range(n):
    os.chdir(str(i))
    os.system("cat q_ga.dat >> ../../../"+folder_name+"/p_total")
    os.system("cat q_ga_part.dat >> ../../../"+folder_name+"/q_ga_part_total")
    os.system("cat q_gb.dat >> ../../../"+folder_name+"/q_gb_total")
    os.system("cat p.dat >> ../../../"+folder_name+"/e_total")
    os.chdir("..")
# os.chdir("simulation/"+str(temp))
# for i in range(n):
#     os.chdir(str(i))
#     os.system("python2 ~/opt/script/CalcRMSD.py top7 dump.lammpstrj rmsd_temp")
#     os.system('cat rmsd_temp | tr " " "\n" |sed 1d > rmsd')
#     os.system("cat rmsd >> ../../../"+folder_name+"/rmsd_total")
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

# os.system("~/opt/script/wham/fused_calc_cv.sc top_he_t400_q100/ top7 50 400 350 450 5 50 100 0 0.98")
