import os
import argparse
import sys
from time import sleep
import subprocess

parser = argparse.ArgumentParser(description="This is my aawsem project headquarter")
parser.add_argument("--jan05", help="Run code on Jan 03 ", action="store_true", default=False)

args = parser.parse_args()
if(args.jan05):
    os.chdir("gagbJan05")
    folder_list = ["T0815", "T0778", "T0766", "T0782", "T0833", "T0844", "T0842", "T0846", "T0803"]
    # folder_list = ["T0792"]
    for protein_name in folder_list:
        os.chdir(protein_name)
        os.system("run.py -n 20 -o "+protein_name+"\/")
        os.chdir("..")
