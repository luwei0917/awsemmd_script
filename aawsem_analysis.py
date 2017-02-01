import os
import argparse
import sys
from time import sleep
import subprocess

parser = argparse.ArgumentParser(description="This is my aawsem project headquarter")
parser.add_argument("--dec25", help="Run code on Dec 25", action="store_true", default=False)
parser.add_argument("--jan03", help="Run code on Jan 03 ", action="store_true", default=False)
parser.add_argument("-t", "--test", help="Test run", action="store_true", default=False)
parser.add_argument("--jan16", help="Run code on Dec 25", action="store_true", default=False)
parser.add_argument("--fix", action="store_true", default=False)
parser.add_argument("--move", action="store_true", default=False)
parser.add_argument("-d", "--debug", action="store_true", default=False)
args = parser.parse_args()

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

if(args.test):
    name = "frag.mem"
    name = "test"
    with open(name) as f:
        for i in range(4):
            next(f)
        for line in f:
            # print(line)
            pos, start, target_start, length, weight = line.split()
            print(pos, start, target_start, length)


if(args.move):
    do("mkdir -p frag_Hybrid")
    folder_list = ["T0792", "T0815", "T0778", "T0766", "T0782", "T0833", "T0844", "T0842", "T0846", "T0803"]
    # folder_list = ["T0766", "T0792", "T0778", "T0782", "T0833", "T0844"]
    for protein_name in folder_list:
        do("mkdir -p frag_Hybrid/"+protein_name)
        do("cp aawsemJan16/{0}/{0}/{0}.seq frag_Hybrid/".format(protein_name))

if(args.fix):
    os.chdir("aawsemJan15")
    folder_list = ["T0766", "T0792", "T0778", "T0782", "T0833", "T0844"]
    for protein_name in folder_list:
        os.chdir(protein_name)
        n = 20
        os.chdir("simulation")
        for i in range(n):
            os.chdir(str(i))
            os.system("cp ~/opt/AASEM/aawsem.slurm .")
            os.system(
                "sed -i.bak 's/PROTEIN/'" +
                protein_name +
                "'/g' aawsem.slurm")
            os.system("sbatch aawsem.slurm")
            os.chdir("..")
        # os.system("run.py -n 20 -o "+protein_name+"\/")
        os.chdir("../..")

if(args.jan16):
    os.chdir("aawsemJan16")
    folder_list = ["T0815", "T0778", "T0766", "T0782", "T0833", "T0844", "T0842", "T0846", "T0803"]
    folder_list = ["T0792"]
    for protein_name in folder_list:
        os.chdir(protein_name)
        n = 20
        os.chdir("simulation_iteration_1")
        for i in range(n):
            os.chdir(str(i))
            os.system("cp ~/opt/AASEM/aawsem.slurm .")
            os.system(
                "sed -i.bak 's/PROTEIN/'" +
                protein_name +
                "'/g' aawsem.slurm")
            os.system("sbatch aawsem.slurm")
            os.chdir("..")
        # os.system("run.py -n 20 -o "+protein_name+"\/")
        os.chdir("../..")
    # for protein_name in folder_list:
    #     os.chdir(protein_name)
    #     n = 20
    #     os.system("mkdir -p analysis")
    #     os.chdir('analysis')
    #     for i in range(n):
    #         os.system("mkdir -p {}".format(str(i)))
    #         os.chdir(str(i))
    #         os.ssytem("cp ../../simulation/{}/dump.lammpstrj .".format(str(i)))
    #         os.system("cp ~/opt/AASEM/aawsem.slurm .")
    #         os.system(
    #             "sed -i.bak 's/PROTEIN/'" +
    #             protein_name +
    #             "'/g' aawsem.slurm")
    #         os.system("sbatch aawsem.slurm")
    #         os.chdir("..")
    #     os.chdir("../..")

if(args.jan03):
    folder_list = ["T0792", "T0778", "T0782", "T0833", "T0844"]
    os.chdir("aawsemJan10")
    for protein_name in folder_list:
        os.chdir(protein_name)
        n = 20
        os.chdir("simulation")
        for i in range(n):
            os.chdir(str(i))
            os.system("cp ~/opt/AASEM/aawsem.slurm .")
            os.system(
                "sed -i.bak 's/PROTEIN/'" +
                protein_name +
                "'/g' aawsem.slurm")
            os.system("sbatch aawsem.slurm")
            os.chdir("..")
        # os.system("run.py -n 20 -o "+protein_name+"\/")
        os.chdir("../..")
if(args.dec25):
    os.chdir("aawsemDec25")
    folder_list = ["T0815", "T0778", "T0766", "T0782", "T0833", "T0844", "T0842", "T0846", "T0803"]
    # folder_list = ["T0792", "T0778", "T0782", "T0833", "T0844"]
    # folder_list = ["T0792"]
    for protein_name in folder_list:
        os.chdir(protein_name)
        n = 20
        os.chdir("simulation")
        for i in range(n):
            os.chdir(str(i))
            os.system("cp ~/opt/AASEM/aawsem.slurm .")
            os.system(
                "sed -i.bak 's/PROTEIN/'" +
                protein_name +
                "'/g' aawsem.slurm")
            os.system("sbatch aawsem.slurm")
            os.chdir("..")
        # os.system("run.py -n 20 -o "+protein_name+"\/")
        os.chdir("../..")
