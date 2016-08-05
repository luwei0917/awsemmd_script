#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import imp
import glob


parser = argparse.ArgumentParser(
        description="This moves all plot scripts here")

parser.add_argument("template", help="the name of template file")
args = parser.parse_args()
protein_name = args.template.split('.')[0]
os.system("cp -r ~/opt/plot_scripts .")
os.chdir("plot_scripts")
# os.system("/usr/local/bin/pymol -qc -r print_final.pml")
os.system(  # replace PROTEIN with pdb name
        "sed -i.bak 's/PROTEIN/'" +
        protein_name +
        "'/g' show_origin.pml")

# os.system(  # replace PROTEIN with pdb name
#         "sed -i.bak 's/NUMBER/'" +
#         str(i) +
#         "'/g' energy.plt")
# os.system("gnuplot energy.plt")

# os.system(  # replace PROTEIN with pdb name
#         "sed -i.bak 's/NUMBER/'" +
#         str(i) +
#         "'/g' q_value.plt")
# os.system("gnuplot q_value.plt")

# os.system("gnuplot detail_energy.plt")
# subprocess.Popen("gnuplot q_value.plt", env=my_env)

os.system(  # replace PROTEIN with pdb name
        "sed -i.bak 's/PROTEIN/'" +
        protein_name +
        "'/g' membraneProtein.tcl")

# os.system("cp ../../"+protein_name+"/*.pdb .")

# # os.chdir('')
# name_list = glob.glob('*.dat')
# # os.chdir('..')
# print(name_list)
# # os.system("gnuplot -e 'out_file_name="test"; in_file_name="T120_400"'")
# for name in name_list:
#     print("gnuplot -e \"out_file_name='%s'; in_file_name='%s'\" \
#         ~/opt/plot_scripts/list_q.gp" % (name[0:-4]+'.pdf', name))
#     os.system("gnuplot -e \"out_file_name='%s'; in_file_name='%s'\" \
#         ~/opt/plot_scripts/list_q.gp" % (name[0:-4]+'.pdf', name))
os.chdir("..")
