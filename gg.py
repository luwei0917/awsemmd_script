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
        automatically analysis the simulation")

parser.add_argument("template", help="the name of template file")
args = parser.parse_args()

exec (open("config.py").read())
n = number_of_run
steps = simulation_steps

protein_name = args.template.strip('/')
n= 20
temp = 400
folder_name = "top7_t400_q100"
os.system("mkdir -p "+folder_name)
os.system("rm -f "+folder_name + "/*")
command = 'sed 1d simulation/{}/%d/wham.dat \
>> {}/all_wham.dat'.format(temp, folder_name)
#cal rmsd
os.chdir("simulation/"+str(temp))
for i in range(n):
    os.chdir(str(i))
    os.system("python2 ~/opt/script/CalcRMSD.py top7 dump.lammpstrj rmsd_temp")
    os.system('cat rmsd_temp | tr " " "\n" |sed 1d > rmsd')
    os.system("cat rmsd >> ../../../"+folder_name+"/rmsd_total")
    os.chdir("..")
os.chdir("../..")
for i in range(n):
    cmd = command % i
    os.system(cmd)
os.chdir(folder_name)
os.system("awk '{print $2}' all_wham.dat > Qw_total")
os.system("awk '{print $3}' all_wham.dat > rg_total")
os.system("awk '{print $4}' all_wham.dat > p_total")
os.system("awk '{print $5}' all_wham.dat > tc_total")
os.system("awk '{print $NF}' all_wham.dat > e_total")

os.chdir("..")
