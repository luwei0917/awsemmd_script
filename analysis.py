#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess

mypath = os.environ["PATH"]
os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
my_env = os.environ.copy()

parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        automatically analysis the simulation")

parser.add_argument("template", help="the name of template file")
parser.add_argument("-n", "--number", type=int, default=20,
                    help="Number of simulation run")
parser.add_argument("-m", "--movie", help="generate the movie",
                    action="store_true")
parser.add_argument("-p", "--plotOnly", help="only generate the plot",
                    action="store_true")
args = parser.parse_args()

list_of_max_q = []

n = args.number
protein_name = args.template.strip('/')

os.system("mkdir -p results")
for i in range(6, n):
    # analysis
    os.system("mkdir -p analysis/"+str(i))
    os.chdir("analysis/"+str(i))
    if not args.plotOnly:
        # move necessary file into analysis folder
        sys.stdout = open("chosen.txt", "w")
        os.system("mv ../../simulation/"+str(i)+"/dump.lammpstrj .")
        os.system("mv ../../simulation/"+str(i)+"/wham.dat .")
        os.system("mv ../../simulation/"+str(i)+"/energy.dat .")
        if(args.movie):
            os.system(
                "python2 ~/opt/script/BuildAllAtomsFromLammps.py \
                dump.lammpstrj movie")
        os.system(
            "python2 ~/opt/script/CalcRMSD.py "+protein_name+" \
            dump.lammpstrj rmsd")
        record_time = 0
        with open('wham.dat') as input_data:
            # Skips text before the beginning of the interesting block:
            record_time = 0
            max_q = 0
            next(input_data)
            for line in input_data:
                # time, q = line.strip().split()

                time = int(line.strip().split()[0])
                q = float(line.strip().split()[1])
                if(q > max_q):
                    record_time = time
                    max_q = q
            list_of_max_q += [max_q]
        time_step = record_time

        print('ITEM: TIMESTEP')
        with open('dump.lammpstrj') as input_data:
            # Skips text before the beginning of the interesting block:
            for line in input_data:
                if line.strip() == str(time_step):
                    print(line.strip())  # Or whatever test is needed
                    break
            # Reads text until the end of the block:
            for line in input_data:  # This keeps reading the file
                if line.strip() == 'ITEM: TIMESTEP':
                    break
                print(line.strip())
        sys.stdout.close()
        os.system(
            "python2 ~/opt/script/BuildAllAtomsFromLammps.py \
            chosen.txt chosen")
    # os.system("cp ~/opt/plot_scripts/energy.plt .")
    os.system("cp ~/opt/plot_scripts/*.plt .")

    os.system(  # replace PROTEIN with pdb name
            "sed -i.bak 's/NUMBER/'" +
            str(i) +
            "'/g' energy.plt")
    os.system("gnuplot energy.plt")

    os.system(  # replace PROTEIN with pdb name
            "sed -i.bak 's/NUMBER/'" +
            str(i) +
            "'/g' q_value.plt")
    os.system("gnuplot q_value.plt")

    os.system("gnuplot detail_energy.plt")
    # subprocess.Popen("gnuplot q_value.plt", env=my_env)
    os.system("cp ~/opt/plot_scripts/show.tcl .")
    os.system(  # replace PROTEIN with pdb name
            "sed -i.bak 's/PROTEIN/'" +
            protein_name +
            "'/g' show.tcl")
    os.chdir("../..")
    if not args.plotOnly:
        os.system("cp "+protein_name+"/*.pdb analysis/"+str(i))

if not args.plotOnly:
    sys.stdout = open("analysis/list_of_max_q", "w")
    for q in list_of_max_q:
        print(q)
    sys.stdout.close()
