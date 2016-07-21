#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import imp
parser = argparse.ArgumentParser(
        description="written by Wei Lu."
)
parser.add_argument("protein", help="The name of the protein")

args = parser.parse_args()
protein_name = args.protein.split('.')[0]
os.system(
    "python2 ~/opt/script/BuildAllAtomsFromLammps.py \
    dump.lammpstrj movie")
os.system("cp ~/opt/plot_scripts/*.tcl .")
os.system(  # replace PROTEIN with pdb name
        "sed -i.bak 's/PROTEIN/'" +
        protein_name +
        "'/g' membrane_show.tcl")
