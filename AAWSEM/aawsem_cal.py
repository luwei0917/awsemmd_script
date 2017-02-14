#!/usr/bin/env python3
import os
import argparse

parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    automatic copy the template file, \
    run simulation and analysis")
parser.add_argument("template", help="the name of template file")
args = parser.parse_args()
protein_name = args.template.strip('/')
n = 20
cal = "~/opt/myCalcQValue_multi.py"

os.system("cp ~/opt/AASEM/"+protein_name+".pdb .")
os.system("python2 "+cal+" "+protein_name+".pdb dump.lammpstrj qw 0")
os.system("python2 "+cal+" "+protein_name+".pdb dump.lammpstrj qo 1")
os.system("python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie "+protein_name+".seq")
