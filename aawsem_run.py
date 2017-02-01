import os
import argparse
import sys
from time import sleep
import subprocess

parser = argparse.ArgumentParser(description="This is my aawsem project headquarter")
parser.add_argument("--dec25", help="Run code on Dec 25", action="store_true", default=False)
parser.add_argument("--jan07", help="Run code on Jan 07 ", action="store_true", default=False)
parser.add_argument("--jan08", help="Run code on Jan 08 ", action="store_true", default=False)
parser.add_argument("--jan12", action="store_true", default=False)
parser.add_argument("--jan15", action="store_true", default=False)
parser.add_argument("--jan16", action="store_true", default=False)
parser.add_argument("-d", "--debug", action="store_true", default=False)
args = parser.parse_args()

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir


def move(protein_name):
    do("mv simulation simulation_iteration_1")
    do("mv list_of_lowest_potential_energy list_of_lowest_potential_energy_v1")
    do("cp -r {0} {0}_v1".format(protein_name))
    do("cp data.v1 {0}/data.{0}".format(protein_name))


def lowestEnergy(protein_name):
    do("lowest_energy.py {}".format(protein_name))
    do("python2 ~/opt/LammpsPDBToCoordinates.py global_lowest_energy v1")
    do("python2 ~/opt/script/CoordinatesToWorkLammpsDataFile.py v1.coord data.v1 -b")

if(args.jan16):
    # folder_list = ["T0792", "T0815", "T0778", "T0766", "T0782", "T0833", "T0844", "T0842", "T0846", "T0803"]
    # speical notes, T0782 is not run second run yet.
    # folder_list = ["T0792"]
    folder_list = ["T0766"]
    # cd("aawsemJan16")
    for protein_name in folder_list:
        cd(protein_name)
        lowestEnergy(protein_name=protein_name)
        move(protein_name=protein_name)
        do("run.py -n 20 -o "+protein_name+"\/")
        cd("..")
    cd("..")

if(args.jan15):
    os.chdir("aawsemJan15")
    folder_list = ["T0766", "T0792", "T0778", "T0782", "T0833", "T0844"]
    # folder_list = ["T0792"]
    for protein_name in folder_list:
        os.chdir(protein_name)
        os.system("cp ../../fix_backbone_coeff.data {}/".format(protein_name))
        os.system("run.py -n 20 -o "+protein_name+"\/")
        os.chdir("..")

if(args.jan12):
    os.chdir("aawsemJan12")
    folder_list1 = ["T0792", "T0815", "T0778", "T0766", "T0782", "T0833", "T0844", "T0842", "T0846", "T0803"]
    folder_list2 = ["T0792", "T0778", "T0782", "T0833", "T0844"]
    # folder_list = ["T0792"]
    folder_list = list(set(folder_list1) - set(folder_list2))
    # folder_list = ["T0792"]
    for protein_name in folder_list:
        os.chdir(protein_name)
        os.system("cp ../../fix_backbone_coeff.data {}/".format(protein_name))
        os.system("run.py -n 20 -o "+protein_name+"\/")
        os.chdir("..")

if(args.jan08):
    os.chdir("aawsemJan08")
    #os.chdir("test")
    folder_list = ["T0792", "T0778", "T0782", "T0833", "T0844"]
    #folder_list = ["T0792"]
    for protein_name in folder_list:
        os.system("mkdir -p "+protein_name)
        os.chdir(protein_name)
        os.system("cp -r ../../aawsemJan07/{0}/{0} .".format(protein_name))
        os.system("cp ../../fix_backbone_coeff.data {}/".format(protein_name))
        os.system("run.py -n 20 -o "+protein_name+"\/")
        os.chdir("..")

if(args.jan07):
    os.chdir("aawsemJan07")
    folder_list = ["T0792", "T0778", "T0782", "T0833", "T0844"]
    # folder_list = ["T0792"]
    for protein_name in folder_list:
        os.chdir(protein_name)
        os.system("run.py -n 20 -o "+protein_name+"\/")
        os.chdir("..")

if(args.dec25):
    os.chdir("aawsemJan16")
    # os.chdir("aawsemDec25")
    folder_list = ["T0815", "T0778", "T0766", "T0782", "T0833", "T0844", "T0842", "T0846", "T0803"]
    # folder_list = ["T0792"]
    for protein_name in folder_list:
        os.chdir(protein_name)
        os.system("run.py -n 20 -o "+protein_name+"\/")
        os.chdir("..")
